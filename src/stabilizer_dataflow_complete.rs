use std::{
    cmp::max,
    collections::{HashMap, HashSet},
    hash::Hash,
};

use bimap::BiHashMap;
use hugr::{
    extension::{prelude::qb_t, simple_op::MakeExtensionOp},
    ops::OpType,
    std_extensions::logic::LogicOp,
    types::TypeRV,
    HugrView, IncomingPort, OutgoingPort, PortIndex,
};
use hugr_core::hugr::internal::PortgraphNodeMap;
use itertools::{chain, Either, Itertools};
use petgraph::visit::{self as pv};
use tket::{
    extension::bool::{bool_type, BoolOp},
    TketOp,
};

use crate::{
    bit_vector::BitVector,
    symplectic_tableau::{PauliXZ, SymplecticTableau},
};

/// The level of detail included for a given node that is useful for phase folding in the enclosing scope
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum DataflowFlowDetail {
    /// No stabilizers used - treat the node as completely opaque
    None,
    /// Includes stabilizers that exist for any choice of control flow path within the node itself.
    /// Breaking this down in terms of types of nodes:
    /// - Gate: any known stabilizers across that type of gate.
    /// - Conditional: any stabilizers over the inputs and outputs that are common to all branches.
    /// - Loop: any stabilizers over the inputs and outputs that are loop invariants (holding regardless of the number of iterations performed).
    /// - Function call: any stabilizers over the inputs and outputs of the function body.
    Invariants,
    /// Same as Invariants, but Conditional nodes can now include stabilizers that hold on individual branches, using control qubits for branch selection.
    ConditionalInvs,
    /// Same as ConditionalInvs, with internal nodes included whenever it is possible to express them in terms of values in the outer scope. Specifically:
    /// - Conditional: the entire tableau for each branch is included, using control qubits for branch selection.
    /// - Loop: for any internal node which can be folded with itself between iterations of the loop, include stabilizers that express its action in terms of the values in the outer scope (as if we hoisted the node out of the loop)
    /// - Gate and function call are same as before. We do not hoist out of function calls, instead leaving it to the programmer or another compiler pass to inline functions.
    All,
}

/// A collection of settings that dictate, for each kind of node, what information we calculate and include in its dataflow analysis.
/// Changing the value of the settings can affect:
/// - How many opportunities for phase folding can be successfully identified.
/// - Whether there is enough information retained in the tableau to be able to resynthesise the circuit.
/// - Number of additional qubits used to represent boundary qubits or internal data.
/// - Number of control qubits used to handle branch selection or classical data.
/// - Computational cost of performing the analysis.
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct DataflowSettings {
    /// The level of detail useful for the enclosing scope.
    /// Needs to be at least DataflowFlowDetail::Invariants to not block phase folding across the node.
    /// Needs to be DataflowFlowDetail::All in order to perform gate hoisting.
    /// Having at least DataflowFlowDetail::ConditionalInvs allows branch-dependent optimisation across Conditional nodes.
    /// DataflowFlowDetail::All is needed for Conditional nodes to be resynthesised from the tableau for the enclosing scope (i.e. without maintaining a separate tableau for the nested blocks).
    pub(crate) flow_level: DataflowFlowDetail,
    /// Whether we track the external interface separately from any flow information.
    /// This is required in order to resynthesise from the tableau.
    /// Relationships between the interface qubits of different nodes can identify when phase folding (or similar gate merging) can be applied.
    /// Turning this on adds an additional qubit for each qubit input or output port on the given node.
    /// If this is on and flow_level is not DataflowFlowDetail::None, additional control qubits will be required to toggle between them.
    pub(crate) include_external_interface: bool,
    /// When qubits exist within a Sum type, do we attempt to track information about the qubits or leave it as a point of inaccuracy.
    /// Turning this on can enable more branch-dependent phase folding opportunities at the cost of more qubits for both the qubits themselves and control qubits for the branches of the sums.
    /// Examples of where this is useful include TryQAlloc gates and internalising the classical outcome of a measurement operation.
    pub(crate) include_sum_types: bool,
    /// For loops, whether or not we include a full copy of the loop body in the tableau.
    /// This is required in order to resynthesise loops from the tableau for the enclosing scope (i.e. without maintaining a separate tableau for the nested block).
    /// When flow_level is DataflowFlowDetail::All, the qubits for the internal nodes within the loop body are reused, but additional qubits will still be created for the input and output nodes within the loop body.
    pub(crate) include_loop_body: bool,
    /// For each gate, we typically want both the external interface to tell us where the gate occurs in the circuit and some flow information across the gate. For an Rz gate, this would take 2 qubits for the external interface and 4 control qubits to turn off the external interface when we want to focus on the flow information; giving a total of 8 qubits (2 ins, 2 outs, 4 controls) and 5 stabilizers.
    /// However, for Rz (and relatives like Rx, Ry, T), we can use 3 qubits (1 in, 1 out, 1 special) and 3 stabilizers to capture all the information we need for synthesis and flow - we essentially model it as a ZX spider with an open wire onto which we can post-select the angle of rotation, similar to its MBQC implementation.
    /// This exploits the fact that we don't need to load specific Pauli strings into both Z and X during resynthesis; an Rz gate only introduces a phase depending on the Pauli string loaded into Z.
    /// A stabilizer involving Z on the rotation's qubit tells us the Pauli string it acts around, and one involving the X says the rest of the stabilizer is dependent on the angle of rotation chosen, i.e. indicating the flow of causal influence through the circuit.
    /// This setting overrides all others for rotation gates to use this shorthand trick.
    /// The same trick of treating it as a ZX spider works for both Measure and MeasureFree, in which case the special qubit can just be treated as the classical output bit (DataflowPoint::SumOutPhControl) when include_sum_types == true.
    pub(crate) rotation_override: bool,
}

impl DataflowSettings {
    pub fn requires_role_control(self, optype: &OpType) -> bool {
        // Under most circumstances, we need role control if we have both external interface data and any flow data
        let default =
            self.flow_level != DataflowFlowDetail::None && self.include_external_interface;
        match optype {
            OpType::ExtensionOp(op) => {
                match TketOp::from_extension_op(&op) {
                    Ok(tkop) => {
                        match tkop {
                            TketOp::Rx
                            | TketOp::Ry
                            | TketOp::Rz
                            | TketOp::T
                            | TketOp::Tdg
                            | TketOp::Measure
                            | TketOp::MeasureFree => {
                                // For these gate types, rotation_override can allow us to avoid the need for role control
                                !self.rotation_override && default
                            }
                            _ => default,
                        }
                    }
                    Err(_) => default,
                }
            }
            OpType::Conditional(_) => {
                // With DataflowFlowDetail::All, we can completely resynthesise the conditional so there is never need to explicitly store the external interface
                self.flow_level != DataflowFlowDetail::All && default
            }
            OpType::TailLoop(_) => {
                // For TailLoop, we use the same option in role control for both the external interface and the loop body (since they use disjoint qubits), so we need role control if either of them is tracked alongside flow information
                self.flow_level != DataflowFlowDetail::None
                    && (self.include_external_interface || self.include_loop_body)
            }
            _ => default,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DataflowPoint<N: Copy + Eq + Hash> {
    /// Qubit is an IncomingPort to the given node.
    /// This includes being used as the output of a particular circuit when the node is the unique Output node.
    /// The third value is used when the type of wire is a Sum type, in which case the vector describes the position in the Sum type's AST of the qubit. We can unify qubits between different variants, so we only care about tracking the positional indices in the row vectors.
    /// For example:
    /// - [] says this qubit is not in a Sum.
    /// - [0, 1] says we have Sum(...| [Sum(...| [_, Qb, ...] |...), ...] |...).
    NodeIn(N, IncomingPort, Vec<usize>),
    /// Same as NodeIn but for OutgoingPorts of a node.
    NodeOut(N, OutgoingPort, Vec<usize>),
    /// As the analysis is primarily capturing the set of wires and connections through the program with the nodes cut out, each NodeIn is really an open output wire out of the system, onto which we would compose the node.
    /// When we give a description of a node, we really give a description for the collection of identity wires incident on it.
    /// A TempIn qubit is the other end of the identity wire to a NodeIn, representing a true input into the modelled system.
    /// These should be treated as temporary qubit indices we can use to express the semantics around a node or subgraph without having to know what they connect to.
    /// We then perform Bell post-selections to join these with the corresponding TempOut qubits or to give a canonical set of input qubits to refer to when taking controlled combinations of tableaus, as used in role combinations and conditionals.
    TempIn(N, IncomingPort, Vec<usize>),
    /// Analogous to TempIn.
    /// As we build the dataflow analysis by starting at the inputs and growing it as we consume the hugr in order, TempOut qubits are additionally used to capture the frontier of the portion of the program analysed so far.
    TempOut(N, OutgoingPort, Vec<usize>),
    /// A specialist qubit used for a rotation/measurement gate when DataflowSettings::rotation_override == true.
    /// Z is used to indicate the Pauli string on which the rotation acts or on which the rotation post-selects as its 0 (+1) outcome (with the 1 (-1) outcome given by the string's negation). This is the only detail needed for synthesising this gate.
    /// X is used to represent the causal influence of this rotation/measurement. Without this, we may find other gates that depend on this rotation/measurement outcome have could no longer be tracked.
    /// For measurements, this is tracked in addition to the SumOutPhControl for the classical value which only features Z on the stabilizer for the measurement effect and nothing for the causal influence.
    Rotation(N),
    /// A control qubit used for a sum type expression coming into a node (e.g. the branch condition on a Conditional).
    /// (node, port, row-tree-index, sum-index-bit, multiplicity-index)
    /// When taking a controlled combination of two tableaus, we will likely need multiple control qubits which will be distinguished using multiplicity-index (there is not particular semantics assigned to each multiplicity value).
    /// When taking controlled combinations of >2 tableaus, we do so pairwise, building a balanced binary tree of controls.
    /// If each branch (in this instance, a case in the sum) is labelled numerically in binary, the bits of that labels describes the path through the binary tree to reach the node.
    /// sum-index-bit indicates the bit in the label or depth in the tree where this control qubit is used for branch selection; 0 is used for the least significant label bit to make it easier to construct these multi-variant controlled conditionals by repeatedly chunking the list into pairs.
    /// When Sum types are nested, we assume other control qubits pick the right variants to get to the Sum we want, but we still need to care about the positions in the type rows; row-tree-index describes these.
    /// For example:
    /// - [] says this control is referring to the top scope Sum.
    /// - [0, 1] says we have Sum(...| [Sum(...| [_, Sum, ...] |...), ...] |...) with this control qubit relating to the inner Sum in this picture.
    /// We assume post-selecting Z on this qubit corresponds to branch selection for the 0 case, and post-selecting X takes the 1 case.
    SumInControl(N, IncomingPort, Vec<usize>, usize, usize),
    /// In controlled combinations, we can ensure at most one generator where each tableau includes the string but with different phases.
    /// This refers to the particular phase control qubit used at some level in the branch selection tree for a sum.
    /// (node, port, row-tree-index, sum-index-bit)
    /// We assume post-selecting +Z on this qubit corresponds to branch selection for the 0 case, and post-selecting -Z takes the 1 case.
    /// Stabilizers should never feature X components on this qubit, it should only contain I or Z.
    /// If there exists any SumInControl for a particular (node, port, row-tree-index, sum-index-bit), then there must also exist a SumInPhControl.
    SumInPhControl(N, IncomingPort, Vec<usize>, usize),
    /// Same as SumInControl but for outgoing ports, e.g. TryQAlloc gates.
    SumOutControl(N, OutgoingPort, Vec<usize>, usize, usize),
    /// Same as SumInPControl but for outgoing ports, e.g. the result bit of a Measure gate.
    SumOutPhControl(N, OutgoingPort, Vec<usize>, usize),
    /// Given a node where we are tracking both the external interface and flow information, these controls allow us to toggle between them.
    /// (node, multiplicity-index)
    /// We assume post-selecting Z on this qubit corresponds to selecting the external interface (also selects the loop body for loops when DataflowFlowDetail::include_loop_body == true), and post-selecting X selects the flow information.
    /// No stabilizers can be shared across these cases, so there is no need for a phase control.
    RoleControl(N, usize),
}

/// When making a controlled combination of two tableaus, we may choose the controls to be encoded either as inputs, outputs, or role controls.
#[derive(Clone, PartialEq, Eq, Hash)]
enum ControlSpec {
    SumInControl(IncomingPort, Vec<usize>, usize),
    SumOutControl(OutgoingPort, Vec<usize>, usize),
    Role,
}

pub struct StabilizerDataflow<H: HugrView> {
    /// Relational dataflow values are captured as stabilizer relations on the Choi-state of the circuit skeleton.
    /// The full tableau represents a summary of a region of the program.
    pub tab: SymplecticTableau,
    /// A bimap relating wires of the program (and supplementary control qubits) to the qubit indices in tab.
    /// For complexity analysis, we use the average complexity O(1) for lookups, insertions, and deletions.
    pub q_index_map: BiHashMap<DataflowPoint<H::Node>, usize>,
}

impl<H: HugrView> Default for StabilizerDataflow<H> {
    fn default() -> Self {
        StabilizerDataflow {
            tab: SymplecticTableau::new(0),
            q_index_map: BiHashMap::default(),
        }
    }
}

impl<H: HugrView> Clone for StabilizerDataflow<H> {
    // O(Q * S)
    fn clone(&self) -> Self {
        StabilizerDataflow {
            tab: self.tab.clone(),
            q_index_map: self.q_index_map.clone(),
        }
    }
}

impl<H: HugrView> StabilizerDataflow<H> {
    // O(self.S * other.Q + other.S * (self.Q + other.Q))
    // Q' = self.Q + other.Q
    // S' = self.S + other.S
    fn tensor_product_in_place(&mut self, other: &Self) {
        // Assumes self and other use disjoint qubits
        let offset = self.tab.add_qubits(other.tab.nb_qubits);
        for (dfp, index) in other.q_index_map.iter() {
            self.q_index_map.insert(dfp.clone(), offset + index);
        }
        for i in 0..other.tab.nb_stabs {
            let mut new_z = BitVector::new(offset);
            new_z.extend_vec(
                other.tab.z[i].get_sized_boolean_vec(other.tab.nb_qubits),
                offset,
            );
            let mut new_x = BitVector::new(offset);
            new_x.extend_vec(
                other.tab.x[i].get_sized_boolean_vec(other.tab.nb_qubits),
                offset,
            );
            self.tab.add_stab(new_z, new_x, other.tab.signs.get(i));
        }
    }

    // O((left.Q + right.Q) * (left.S + right.S))
    // Q' = left.Q + right.Q
    // S' = left.S + right.S
    fn tensor_product(left: &Self, right: &Self) -> Self {
        // Assumes left and right use disjoint qubits
        let mut res = left.clone();
        res.tensor_product_in_place(right);
        res
    }

    // When building identity wires for output interfaces, we may wish to use SumOutControls instead of SumInControls, and combining flow and interface data uses the same logic as for conditionals.
    // O(union_Q * (left.S + right.S)^2)
    //      Unifying qubit indexing O(right.Q)
    //      Cloning left O(left.Q * left.S)
    //      Extending left to new indexing O(right.Q * left.S)
    //      Reordering right to new indexing O(union_Q * right.S)
    //      Joint decomposition O(union_Q * (left.S + right.S)^2)
    //      Adding control qubits takes O(left.S + right.S) each; at most left.S + right.S controls, so O((left.S + right.S)^2) total
    //  Q' = union_Q + n_ac_pairs + just_left + just_right <= union_Q + left.S + right.S
    //  S' = left.S + right.S - join_ignoring_phases <= left.S + right.S
    fn controlled_conditional(
        left: &Self,
        right: &Self,
        control_node: H::Node,
        control_spec: &ControlSpec,
        always_add_phase_ctrl: bool,
    ) -> StabilizerDataflow<H> {
        // Unify qubit indexing then extend and reorder tableaus to match
        let mut q_ind_map: BiHashMap<DataflowPoint<H::Node>, usize> = left.q_index_map.clone();
        let mut l = left.tab.clone();
        for q in 0..right.tab.nb_qubits {
            let dfp = right.q_index_map.get_by_right(&q).unwrap();
            if !q_ind_map.contains_left(dfp) {
                q_ind_map.insert(dfp.clone(), l.add_qubit());
            }
        }
        let mut r = SymplecticTableau::new(l.nb_qubits);
        for i in 0..right.tab.nb_stabs {
            let mut new_z = BitVector::new(l.nb_qubits);
            let mut new_x = BitVector::new(l.nb_qubits);
            for r_index in right.tab.z[i].get_all_ones(right.tab.nb_qubits) {
                let j_index = q_ind_map
                    .get_by_left(right.q_index_map.get_by_right(&r_index).unwrap())
                    .unwrap();
                new_z.xor_bit(*j_index);
            }
            for r_index in right.tab.x[i].get_all_ones(right.tab.nb_qubits) {
                let j_index = q_ind_map
                    .get_by_left(right.q_index_map.get_by_right(&r_index).unwrap())
                    .unwrap();
                new_x.xor_bit(*j_index);
            }
            r.add_stab(new_z, new_x, right.tab.signs.get(i));
        }
        // Take the joint decomposition of the tableaus and add control qubits appropriately
        let mut jd = SymplecticTableau::joint_decomposition(&l, &r);
        if always_add_phase_ctrl || (jd.join_range().0 != 0 && *control_spec != ControlSpec::Role) {
            // If any controls are needed, always include a phase control as well.
            // Add this first so it is in a predictable qubit index.
            let control_index = jd.tab.add_qubit();
            let control_point = match control_spec {
                ControlSpec::SumInControl(in_port, ref row_index, sum_index) => {
                    DataflowPoint::SumInPhControl(
                        control_node,
                        *in_port,
                        row_index.clone(),
                        *sum_index,
                    )
                }
                ControlSpec::SumOutControl(out_port, ref row_index, sum_index) => {
                    DataflowPoint::SumOutPhControl(
                        control_node,
                        *out_port,
                        row_index.clone(),
                        *sum_index,
                    )
                }
                ControlSpec::Role => unreachable!(),
            };
            q_ind_map.insert(control_point, control_index);
            if let Some(i) = jd.phase_flip_term_index() {
                // Role combinations always have empty joins and no phase flip term
                jd.tab.z[i].xor_bit(control_index);
            }
        }
        let mut control_number = 0;
        for i in 0..jd.n_ac_pairs {
            let control_index = jd.tab.add_qubit();
            let control_point = match control_spec {
                ControlSpec::SumInControl(in_port, ref row_index, sum_index) => {
                    DataflowPoint::SumInControl(
                        control_node,
                        *in_port,
                        row_index.clone(),
                        *sum_index,
                        control_number,
                    )
                }
                ControlSpec::SumOutControl(out_port, ref row_index, sum_index) => {
                    DataflowPoint::SumOutControl(
                        control_node,
                        *out_port,
                        row_index.clone(),
                        *sum_index,
                        control_number,
                    )
                }
                ControlSpec::Role => DataflowPoint::RoleControl(control_node, control_number),
            };
            q_ind_map.insert(control_point, control_index);
            jd.tab.z[2 * i].xor_bit(control_index);
            jd.tab.x[2 * i + 1].xor_bit(control_index);
            control_number += 1;
        }
        let (l_begin, l_end) = jd.just_left_range();
        for i in l_begin..l_end {
            let control_index = jd.tab.add_qubit();
            let control_point = match control_spec {
                ControlSpec::SumInControl(in_port, ref row_index, sum_index) => {
                    DataflowPoint::SumInControl(
                        control_node,
                        *in_port,
                        row_index.clone(),
                        *sum_index,
                        control_number,
                    )
                }
                ControlSpec::SumOutControl(out_port, ref row_index, sum_index) => {
                    DataflowPoint::SumOutControl(
                        control_node,
                        *out_port,
                        row_index.clone(),
                        *sum_index,
                        control_number,
                    )
                }
                ControlSpec::Role => DataflowPoint::RoleControl(control_node, control_number),
            };
            q_ind_map.insert(control_point, control_index);
            jd.tab.z[i].xor_bit(control_index);
            control_number += 1;
        }
        let (r_begin, r_end) = jd.just_right_range();
        for i in r_begin..r_end {
            let control_index = jd.tab.add_qubit();
            let control_point = match control_spec {
                ControlSpec::SumInControl(in_port, ref row_index, sum_index) => {
                    DataflowPoint::SumInControl(
                        control_node,
                        *in_port,
                        row_index.clone(),
                        *sum_index,
                        control_number,
                    )
                }
                ControlSpec::SumOutControl(out_port, ref row_index, sum_index) => {
                    DataflowPoint::SumOutControl(
                        control_node,
                        *out_port,
                        row_index.clone(),
                        *sum_index,
                        control_number,
                    )
                }
                ControlSpec::Role => DataflowPoint::RoleControl(control_node, control_number),
            };
            q_ind_map.insert(control_point, control_index);
            jd.tab.x[i].xor_bit(control_index);
            control_number += 1;
        }
        StabilizerDataflow {
            tab: jd.tab,
            q_index_map: q_ind_map,
        }
    }

    // Unlike controlled_conditional, injections always use SumOutPhControls
    // O(Q + S)
    //      Adding new qubit is O(S)
    //      Adding new stabilizer is O(Q)
    // Q' = Q + 1
    // S' = S + 1
    fn sum_injection(
        &mut self,
        inj_variant: Either<(), ()>,
        control_node: H::Node,
        control_port: OutgoingPort,
        control_row_index: Vec<usize>,
        control_sum_index: usize,
    ) {
        let control_qb = self.tab.add_qubit();
        let dfp = DataflowPoint::SumOutPhControl(
            control_node,
            control_port,
            control_row_index,
            control_sum_index,
        );
        self.q_index_map.insert(dfp, control_qb);
        let mut z = BitVector::new(self.tab.nb_qubits);
        z.xor_bit(control_qb);
        self.tab.add_stab(
            z,
            BitVector::new(self.tab.nb_qubits),
            inj_variant.is_right(),
        );
    }

    // With nesting_depth and max_width, there are at most max_width^nesting_depth many possible values for the row_index of a DataflowPoint, which gives us at most 2*max_width^nesting_depth qubits (excluding controls); but the stabilizers from them might have different controls on each branch of a variant and thus necessitate adding extra controls at each point in the type's syntax tree
    // Complexity analysis reveals the actual bounds are pretty big...
    // Qubit base case
    //      O(1) time
    //      Q' = 2
    //      S' = 2
    // Bool base case
    //      O(1) time
    //      Q' = 1
    //      S' = 0
    // Sum case in terms of num_variants, max_width of variants, max_Q and max_S of recursive calls
    //      Time complexity dominated by adding control qubits
    //          Folding the tensor product within a single variant is O(max_width^3 * max_Q * max_S) time
    //              Each tensor product takes time equal to its resulting size
    //              Stacking rectangles for the tableaus after each fold gives a pyramid of base (max_width * max_Q) qubits by (max_width * max_S) stabilizers and height max_width; volume of a pyramid is 1/3 B h
    //              This could obviously be improved by computing the tensor product of all tableaus directly, saving a factor of max_width in the time complexity
    //          Adding control qubits is hard to analyse
    //              Assuming no overlapping qubit names (and hence no joint stabilizers) between any variants, each controlled_conditional gives Q' = 2*Q + 2*S, S' = 2*S
    //              Hence, in total,
    //                  Q' = 2^{log2 num_variants}*(max_width * max_Q) + 2^{log2 num_variants}*(log2 num_variants)*(max_width * max_S) < 2*num_variants*max_width*(max_Q + log2(num_variants) * max_S)
    //                  S' = 2^{log2 num_variants}*(max_width * max_S) < 2*num_variants * max_width * max_S
    //                  Round 1: O(2^3 QS^2)
    //                  Round 2: O(2^3 (2Q+2S)(2S)^2) = O(2^6 (QS^2 + S^3))
    //                  Round 3: O(2^3 (4Q+4S+4S)(4S)^2) = O(2^9 (QS^2 + 2S^3))
    //                  Round N: O(2^{3N} (QS^2 + (N-1)S^3))
    //                  Cumulative to Round N: 1*RN + 2*R(N-1) + ... + 2^(N-1)*R1
    //                                          = O(N2^{3N}QS^2 + N^2 2^{4N}S^3)
    //                  Total time O(num_variants^3 * (log2 num_variants) * max_width^3 * max_S^2 * (max_Q + max_S * num_variants * (log2 num_variants)))
    //              Obviously, this is a very conservative estimate and in practice we can get a lot of qubits overlapping and the sum types themselves are typically not very large
    // Iterating this for nesting_depth times...
    //          S' < (2*num_variants * max_width)^{nesting_depth} * max_S
    //          Q' < (2*num_variants * max_width)^{nesting_depth} * (max_Q + nesting_depth * log2(num_variants)^{nesting_depth} * max_S)
    //          Using N = log2 num_variants, W = max_width, Q = max_Q, S = max_S
    //          Round 1: O(N2^{3N}W^3 QS^2 + N^2 2^{4N}W^3 S^3)
    //          Round 2: O(N2^{3N}W^3 2^{N+1} W (Q + NS) (2^{N+1} W S)^2 + N^2 2^{4N}W^3 (2^{N+1} W S)^3)
    //          Round 2: O(N2^{6N+3}W^6 (QS^2 + NS^3) + N^2 2^{7N+3}W^6 S^3)
    //          Round 3: O(N2^{3N}W^3 2^{2N+2} W^2 (Q + 2NS) (2^{2N+2} W^2 S)^2 + N^2 2^{4N}W^3 (2^{2N+2} W^2 S)^3)
    //          Round 3: O(N2^{9N+6}W^9 (QS^2 + 2NS^3) + N^2 2^{10N+6}W^9 S^3)
    //          Round D: O(N2^{3DN+3D}W^{3D} QS^2 + DN^2 2^{3DN+3D}W^{3D} S^3 + N^2 2^{3DN + 3D + N}W^{3D} S^3)
    //          Round D: O(2^{3DN+3D}W^{3D} (NQS^2 + (D + 2^N) N^2 S^3))
    //          Cumulative to Round D: (W2^N)^{D-d} many calls to round d
    //                  \sum_d (W2^N)^{D-d} 2^{3dN+3d}W^{3d} (NQS^2 + (d + 2^N) N^2 S^3)
    //                  \sum_d 2^{DN + 2dN + 3d} W^{D + 2d} (NQS^2 + (d + 2^N) N^2 S^3)
    //                  O(2^{3DN + 3D}W^{3D} (NQS^2 + (D + 2^N) N^2 S^3))
    // Fixing initial Q and S as O(1) from base cases...
    // O((2 * num_variants * max_width)^{3*nesting_depth} * (nesting_depth + num_variants) * (log2 num_variants)^2)
    // With the same assumptions, the input syntax tree for the sum type itself can be O((num_variants * max_width)^{nesting_depth}), so the dimensions of the final tableau still grow reasonably wrt the actual input size, and the time complexity is a little over cubic in the size of the syntax tree
    // For complexity analysis of anything using this method, we will assume that the qubit count, stabilizer count, and time cost of all calls to this method are bound by some (potentially large) variables IdQ, IdS, IdT respectively, which are indicators of how complex the Sum types are in the target program
    fn identity_wire_recursive(
        hugr: &H,
        source: H::Node,
        source_port: Either<IncomingPort, OutgoingPort>,
        target: H::Node,
        target_port: Either<IncomingPort, OutgoingPort>,
        wire_type: &TypeRV,
        row_index: &Vec<usize>,
    ) -> Self {
        if *wire_type == qb_t() {
            // If type is a qubit, return an identity tableau
            let source_dfp = match source_port {
                Either::Left(in_port) => DataflowPoint::TempIn(source, in_port, row_index.clone()),
                Either::Right(out_port) => {
                    DataflowPoint::NodeOut(source, out_port, row_index.clone())
                }
            };
            let target_dfp = match target_port {
                Either::Left(in_port) => DataflowPoint::NodeIn(target, in_port, row_index.clone()),
                Either::Right(out_port) => {
                    DataflowPoint::TempOut(target, out_port, row_index.clone())
                }
            };
            let q_index_map: BiHashMap<DataflowPoint<H::Node>, usize> =
                BiHashMap::from_iter([(source_dfp, 0), (target_dfp, 1)]);
            let mut tab = SymplecticTableau::new(2);
            let mut both_qubits = BitVector::new(2);
            both_qubits.xor_bit(0);
            both_qubits.xor_bit(1);
            tab.add_stab(both_qubits.clone(), BitVector::new(2), false);
            tab.add_stab(BitVector::new(2), both_qubits, false);
            StabilizerDataflow {
                tab: tab,
                q_index_map: q_index_map,
            }
        } else if *wire_type == bool_type() {
            // Whilst prelude::bool_t is just Sum[unit, unit], tket::bool_type is an opaque boolean with its own set of logical operations.
            // Since TketOp::Measure produces a prelude::bool_t and TketOp::MeasureFree produces a tket::bool_type, we can expect them to be mixed in practice.
            // Whenever we encounter an opaque boolean, treat it the same as another boolean.
            // This routine for Sum[unit, unit] would produce a tableau with a single qubit and no stabilizers since the value of the boolean is unknown

            let control = match source_port {
                Either::Left(in_port) => {
                    DataflowPoint::SumInPhControl(source, in_port, row_index.clone(), 0)
                }
                Either::Right(out_port) => {
                    DataflowPoint::SumOutPhControl(source, out_port, row_index.clone(), 0)
                }
            };
            StabilizerDataflow {
                tab: SymplecticTableau::new(1),
                q_index_map: BiHashMap::from_iter([(control, 0)]),
            }
        } else if let Some(sum) = wire_type.as_sum() {
            // If type is a Sum type, iterate through the variants, and recursively call for each element of the row and tensor product, then combine the variants
            let mut variant_analyses = sum
                .variants()
                .map(|type_row| {
                    let row_analyses =
                        type_row
                            .iter()
                            .enumerate()
                            .map(|(element_i, element_type)| {
                                let mut new_row_index = row_index.clone();
                                new_row_index.push(element_i);
                                Self::identity_wire_recursive(
                                    hugr,
                                    source,
                                    source_port,
                                    target,
                                    target_port,
                                    element_type,
                                    &new_row_index,
                                )
                            });
                    row_analyses.fold(StabilizerDataflow::default(), |acc, val| {
                        StabilizerDataflow::tensor_product(&acc, &val)
                    })
                })
                .collect_vec();
            match sum.num_variants() {
                0 => StabilizerDataflow {
                    tab: SymplecticTableau::new(0),
                    q_index_map: BiHashMap::default(),
                },
                1 => variant_analyses[0].clone(),
                num => {
                    // ilog2 computes the floor of the log, we want ceiling of log so that e.g. 2 variants uses 1 bit, 3 or 4 variants uses 2 bits
                    let sum_index_bits = (num - 1).ilog2() as usize + 1;
                    for sum_index in 0..sum_index_bits {
                        let control_spec = match source_port {
                            Either::Left(in_port) => {
                                ControlSpec::SumInControl(in_port, row_index.clone(), sum_index)
                            }
                            Either::Right(out_port) => {
                                ControlSpec::SumOutControl(out_port, row_index.clone(), sum_index)
                            }
                        };
                        variant_analyses = variant_analyses
                            .chunks(2)
                            .into_iter()
                            .map(|chunk_it| {
                                if chunk_it.len() == 1 {
                                    chunk_it[0].clone()
                                } else {
                                    // ID(A+B) = L(ID(A)) + R(ID(B))
                                    // We don't need to care about the injections since we are using the same controls for the classical variables at the source and target, so we should only need to use Sum(In/Out)Controls to pick between different sets of qubit identities
                                    StabilizerDataflow::controlled_conditional(
                                        &chunk_it[0],
                                        &chunk_it[1],
                                        source,
                                        &control_spec,
                                        true,
                                    )
                                }
                            })
                            .collect_vec();
                    }
                    variant_analyses[0].clone()
                }
            }
        } else {
            // Don't know this type; return an empty tableau
            StabilizerDataflow::default()
        }
    }

    // O(IdT * inp.out_ports + IdQ * IdS * inp.out_ports^3)
    // Q' < IdQ * inp.out_ports
    // S' < IdS * inp.out_ports
    // Similar to identity_wire_recursive, we could reduce the factor of inp.out_ports by building the tensor product directly instead of with a fold
    fn initialise_from_input(hugr: &H, inp: H::Node, use_sums: bool) -> Self {
        hugr.out_value_types(inp)
            .fold(StabilizerDataflow::default(), |acc, (out, out_type)| {
                let new_wire = if use_sums || out_type == qb_t() {
                    Self::identity_wire_recursive(
                        hugr,
                        inp,
                        Either::Right(out),
                        inp,
                        Either::Right(out),
                        &out_type.into(),
                        &vec![],
                    )
                } else {
                    StabilizerDataflow::default()
                };
                Self::tensor_product(&acc, &new_wire)
            })
    }

    /// Projecting immediately to boundary information will lose flow information that relies on nodes with both flow and interface roles. This method post-selects all role controls to the flow information so we can then safely project to flow information over the boundaries with project_non_io
    /// O(Q * S * min(n_roles, S))
    /// Q' = Q
    fn project_to_flow(&mut self) {
        let post_selects: Vec<(usize, PauliXZ, bool)> = (0..self.tab.nb_qubits).filter_map(|q| {
            if let Some(DataflowPoint::RoleControl(_, _)) = self.q_index_map.get_by_right(&q) {
                Some((q, PauliXZ::X, false))
            } else { None }
        }).collect_vec();
        self.tab.post_select_1qs(&post_selects);
    }

    // O(Q * S^2)
    fn project_non_io(&mut self, inp: H::Node, out: H::Node) {
        let non_ios: Vec<usize> = (0..self.tab.nb_qubits).filter(|q| {
            match self.q_index_map.get_by_right(q).unwrap() {
                DataflowPoint::NodeIn(n, _, _)
                | DataflowPoint::NodeOut(n, _, _)
                | DataflowPoint::SumInControl(n, _, _, _, _)
                | DataflowPoint::SumInPhControl(n, _, _, _)
                | DataflowPoint::SumOutControl(n, _, _, _, _)
                | DataflowPoint::SumOutPhControl(n, _, _, _) => *n != inp && *n != out,
                _ => true,
            }
        }).collect_vec();
        let project_cols = chain!(
            non_ios.iter().map(|i| (*i, PauliXZ::X)),
            non_ios.iter().map(|i| (*i, PauliXZ::Z)),
        )
        .collect_vec();
        self.tab.project_cols_to_zero(&project_cols);
        for i in non_ios {
            let moved_qb = self.tab.delete_qubit(i);
            self.q_index_map.remove_by_right(&i);
            if let Some(mq) = moved_qb {
                // If the removal caused a qubit to have changed index, update it
                let (removed_dfp, _) = self.q_index_map.remove_by_right(&mq).unwrap();
                self.q_index_map.insert(removed_dfp, i);
            }
        }
    }

    /// Helper method for updating the tableau for a node operation.
    ///
    /// Looks up the tableau columns corresponding to the node inputs and updates the
    /// frontier with the node outputs provided by the `go` closure.
    fn apply_op_with<const IN: usize, const OUT: usize>(
        &mut self,
        hugr: &H,
        node: H::Node,
        go: impl FnOnce(
            &mut SymplecticTableau,
            &mut BiHashMap<DataflowPoint<H::Node>, usize>,
            [usize; IN],
        ) -> [usize; OUT],
    ) {
        // Collect tableau columns for node inputs
        let in_cols = (0..IN)
            .map(|i| {
                let (pred, pred_port) = hugr
                    .single_linked_output(node, IncomingPort::from(i))
                    .unwrap();
                let found = self.q_index_map.remove_by_left(&DataflowPoint::TempOut(
                    pred,
                    pred_port,
                    vec![],
                ));
                match found {
                    Some((_, col)) => col,
                    // If a qubit does not exist for the input, make one that doesn't feature in any of the stabilizers, i.e. a maximally mixed state
                    None => self.tab.add_qubit(),
                }
            })
            .collect_array()
            .unwrap();
        // Run closure and update frontier with the returned tableau columns
        let out_cols = go(&mut self.tab, &mut self.q_index_map, in_cols);
        for (out_port, out_col) in hugr.node_outputs(node).zip(out_cols) {
            self.q_index_map
                .insert(DataflowPoint::TempOut(node, out_port, vec![]), out_col);
        }
    }

    /// Helper method for updating the tableau for a Clifford operation.
    fn apply_clifford_with<const N: usize>(
        &mut self,
        hugr: &H,
        node: H::Node,
        go: impl FnOnce(&mut SymplecticTableau, [usize; N]),
    ) {
        self.apply_op_with(hugr, node, |tab, _, cols| {
            go(tab, cols);
            // Output columns for Cliffords are the same as the input columns
            cols
        });
    }

    /// Used when DataflowSettings::include_external_interface == true and DataflowSettings::flow_level == DataflowFlowDetail::None
    /// O(InQ * Q * S + node.out_ports * IdQ * S + node.out_ports * Q * IdS + node.out_ports^2 * IdQ * IdS) time where InQ is the number of qubits currently in the tableau associated to predecessors of node
    ///     Identifying qubits to be renamed/removed/copied is done in O(Q)
    ///     Each qubit rename is O(1)
    ///     Each qubit remove is O(Q * S)
    ///     Each copied qubit costs O(Q + S) to add both the qubit and the new stabilizer
    ///     Doing the tensor products in place takes O(S * node.out_ports * IdQ + node.out_ports * IdS * (Q + node.out_ports * IdQ))
    /// Q' < Q + node.out_ports * IdQ
    /// S' < S + InQ + node.out_ports * IdS
    fn apply_opaque(&mut self, hugr: &H, node: H::Node, use_sums: bool) {
        // Rename frontier qubits to the correct qubit labels
        // If use_sums == false, then we also get rid of any inputs that are held within sum types
        // If use_sums == true, then for any controls that exist from predecessors (which only exist if they have use_sums == true), copy the classical value by adding a corresponding Ph control and a ZZ stabilizer to connect them
        // The following map returns an IncomingPort if the qubit needs to be relabelled and None if it should be removed
        let mut pred_map: HashMap<(H::Node, OutgoingPort), Option<IncomingPort>> = HashMap::new();
        for (in_port, in_type) in hugr.in_value_types(node) {
            pred_map.insert(
                hugr.single_linked_output(node, in_port).unwrap(),
                if use_sums || in_type == qb_t() {
                    Some(in_port)
                } else {
                    None
                },
            );
        }
        let mut qbs_to_rename: Vec<(usize, DataflowPoint<H::Node>)> = vec![];
        let mut qbs_to_remove: Vec<(DataflowPoint<H::Node>, usize)> = vec![];
        let mut ctrls_to_copy: Vec<(usize, DataflowPoint<H::Node>)> = vec![];
        // Iterate through self.q_index_map in a deterministic order
        for qb in 0..self.tab.nb_qubits {
            let dfp = self.q_index_map.get_by_right(&qb).unwrap();
            match dfp {
                DataflowPoint::TempOut(source, out_port, row_index) => {
                    if let Some(found) = pred_map.get(&(*source, *out_port)) {
                        if let Some(in_port) = found {
                            qbs_to_rename.push((
                                qb,
                                DataflowPoint::NodeIn(node, *in_port, row_index.clone()),
                            ));
                        } else {
                            qbs_to_remove.push((dfp.clone(), qb));
                        }
                    }
                }
                DataflowPoint::SumOutPhControl(source, out_port, row_index, sum_index) => {
                    if let Some(found) = pred_map.get(&(*source, *out_port)) {
                        if let Some(in_port) = found {
                            ctrls_to_copy.push((
                                qb,
                                DataflowPoint::SumInPhControl(
                                    node,
                                    *in_port,
                                    row_index.clone(),
                                    *sum_index,
                                ),
                            ));
                        }
                    }
                }
                _ => {}
            }
        }
        for (i, dfp) in qbs_to_rename {
            self.q_index_map.remove_by_right(&i);
            self.q_index_map.insert(dfp, i);
        }
        self.tab.project_cols_to_zero(
            &chain!(
                qbs_to_remove.iter().map(|(_, q)| (*q, PauliXZ::Z)),
                qbs_to_remove.iter().map(|(_, q)| (*q, PauliXZ::X)),
            )
            .collect_vec(),
        );
        let new_ctrls_offset = self.tab.add_qubits(ctrls_to_copy.len());
        for (i, (ctrl, new_ctrl_dfp)) in ctrls_to_copy.iter().enumerate() {
            let new_ctrl = new_ctrls_offset + i;
            self.q_index_map.insert(new_ctrl_dfp.clone(), new_ctrl);
            let mut zz = BitVector::new(self.tab.nb_qubits);
            zz.xor_bit(*ctrl);
            zz.xor_bit(new_ctrl);
            self.tab
                .add_stab(zz, BitVector::new(self.tab.nb_qubits), false);
        }
        for (dfp, _) in qbs_to_remove {
            // qubit index might have been moved by an earlier removal, so grab it from the index map
            let (_, q) = self.q_index_map.remove_by_left(&dfp).unwrap();
            if let Some(moved_qb) = self.tab.delete_qubit(q) {
                let (moved_dfp, _) = self.q_index_map.remove_by_right(&moved_qb).unwrap();
                self.q_index_map.insert(moved_dfp, q);
            }
        }
        // Add identity wires for each output port
        for (out, out_type) in hugr.out_value_types(node) {
            let new_wire = if use_sums || out_type == qb_t() {
                Self::identity_wire_recursive(
                    hugr,
                    node,
                    Either::Right(out),
                    node,
                    Either::Right(out),
                    &out_type.into(),
                    &vec![],
                )
            } else {
                StabilizerDataflow::default()
            };
            self.tensor_product_in_place(&new_wire);
        }
    }

    /// Apply a quantum gate to the end of the summary.
    /// Used for nodes with TketOps when DataflowSettings::include_external_interface == false and DataflowSettings::flow_level != DataflowFlowDetail::None (so we assume flow_level passed in here is not None), i.e. we can apply some amount of flow information in-place
    /// O(Q * S) time (unless we add more TketOps that are handled as opaque)
    ///     Unitary Clifford gates are O(S)
    ///     QFree/CRz/Toffoli require projection of 2/2/3 columns, so O(Q * S)
    ///     Basic rotations and measurements require O(S) to add the new qubit(s) and extend existing stabilizers, then O(Q) to add the new one; if DataflowSettings::rotation_override==false, then O(Q * S) to project a column
    ///     (Try)QAlloc requies O(S) to add the qubit(s) and O(Q) to add the stabilizer
    ///     Reset combines an O(Q * S) projection and an O(Q) new stabilizer
    fn apply_quantum_gate_flow_only(
        &mut self,
        hugr: &H,
        node: H::Node,
        op: TketOp,
        settings: DataflowSettings,
    ) {
        match op {
            TketOp::H => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_h(col);
                });
            }
            TketOp::CX => {
                self.apply_clifford_with(hugr, node, |tab, [col0, col1]| {
                    tab.append_cx(col0, col1);
                });
            }
            TketOp::CY => {
                self.apply_clifford_with(hugr, node, |tab, [col0, col1]| {
                    tab.append_cy(col0, col1);
                });
            }
            TketOp::CZ => {
                self.apply_clifford_with(hugr, node, |tab, [col0, col1]| {
                    tab.append_cz(col0, col1);
                });
            }
            TketOp::CRz => {
                self.apply_op_with(hugr, node, |tab, _, [col_in0, col_in1]| {
                    // Z on each qubit commute through, so project away any stabilizer with X on either qubit
                    tab.project_cols_to_zero(&vec![(col_in0, PauliXZ::X), (col_in1, PauliXZ::X)]);
                    [col_in0, col_in1]
                });
            }
            TketOp::T | TketOp::Tdg | TketOp::Rz => {
                if settings.rotation_override {
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        // Add a new qubit for the rotation
                        let col_rot = tab.add_qubit();
                        q_index_map.insert(DataflowPoint::Rotation(node), col_rot);
                        // Z on input propagates through to output, so leave it unchanged
                        // X on input copies to X on rotation and output
                        // Applying a CX(col_in, col_rot) combines these two effects
                        tab.append_cx(col_in, col_rot);
                        // Add the remaining stabilizer between rotation and output
                        let mut zz = BitVector::new(tab.nb_qubits);
                        zz.xor_bit(col_in);
                        zz.xor_bit(col_rot);
                        tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                        [col_in]
                    });
                } else {
                    self.apply_op_with(hugr, node, |tab, _, [col_in]| {
                        // The only consistent stabilizer is that Z propagates, so project away Xs
                        tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X)]);
                        [col_in]
                    });
                }
            }
            TketOp::Measure => {
                if settings.rotation_override {
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        if settings.include_sum_types {
                            // IORSRole
                            // ZZ       Handled by qubit renaming
                            // XXX Z    Qubit renaming and gates to copy X
                            //  ZZ      Added as new stab
                            //  Z ZX    Added as new stab
                            let [col_rot, col_ctrl, col_role] = tab.add_n_qubits();
                            q_index_map.insert(DataflowPoint::Rotation(node), col_rot);
                            q_index_map.insert(
                                DataflowPoint::SumOutPhControl(
                                    node,
                                    OutgoingPort::from(1),
                                    vec![],
                                    0,
                                ),
                                col_ctrl,
                            );
                            q_index_map.insert(DataflowPoint::RoleControl(node, 0), col_role);
                            tab.append_cx(col_in, col_rot);
                            tab.append_cz(col_in, col_role);
                            let mut zz = BitVector::new(tab.nb_qubits);
                            zz.xor_bit(col_in);
                            zz.xor_bit(col_rot);
                            let mut zzx_z = BitVector::new(tab.nb_qubits);
                            zzx_z.xor_bit(col_in);
                            zzx_z.xor_bit(col_ctrl);
                            let mut zzx_x = BitVector::new(tab.nb_qubits);
                            zzx_x.xor_bit(col_role);
                            tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                            tab.add_stab(zzx_z, zzx_x, false);
                            [col_in]
                        } else {
                            // IOR
                            // ZZ   Handled by qubit renaming
                            // XXX  Qubit renaming and CX gate to copy X
                            //  ZZ  Added as new stab
                            let col_rot = tab.add_qubit();
                            q_index_map.insert(DataflowPoint::Rotation(node), col_rot);
                            tab.append_cx(col_in, col_rot);
                            let mut zz = BitVector::new(tab.nb_qubits);
                            zz.xor_bit(col_in);
                            zz.xor_bit(col_rot);
                            tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                            [col_in]
                        }
                    });
                } else {
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X)]);
                        if settings.include_sum_types {
                            // In addition to projecting away Xs, we make a record of how the measurement result affects the projection
                            let ctrl = tab.add_qubit();
                            q_index_map.insert(
                                DataflowPoint::SumOutPhControl(
                                    node,
                                    OutgoingPort::from(1),
                                    vec![],
                                    0,
                                ),
                                ctrl,
                            );
                            let mut zz = BitVector::new(tab.nb_qubits);
                            zz.xor_bit(col_in);
                            zz.xor_bit(ctrl);
                            tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                        }
                        [col_in]
                    });
                }
            }
            TketOp::S => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_s(col);
                });
            }
            TketOp::Sdg => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_s(col);
                    tab.append_z(col);
                });
            }
            TketOp::X => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_x(col);
                });
            }
            TketOp::Y => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_x(col);
                    tab.append_z(col);
                });
            }
            TketOp::Z => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_z(col);
                });
            }
            TketOp::Rx => {
                if settings.rotation_override {
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        // Add a new qubit for the rotation
                        let col_rot = tab.add_qubit();
                        q_index_map.insert(DataflowPoint::Rotation(node), col_rot);
                        // X on input propagates through to output, so leave it unchanged
                        // Z on input copies to X on rotation and output
                        // No single native gate combines these effects, so just conjugate the Rz case with Hs
                        tab.append_h(col_in);
                        tab.append_cx(col_in, col_rot);
                        // Add the remaining stabilizer between rotation and output
                        let mut zz = BitVector::new(tab.nb_qubits);
                        zz.xor_bit(col_in);
                        zz.xor_bit(col_rot);
                        tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                        tab.append_h(col_in);
                        [col_in]
                    });
                } else {
                    self.apply_op_with(hugr, node, |tab, _, [col_in]| {
                        // The only consistent stabilizer is that X propagates, so project away Zs
                        tab.project_cols_to_zero(&vec![(col_in, PauliXZ::Z)]);
                        [col_in]
                    });
                }
            }
            TketOp::Ry => {
                if settings.rotation_override {
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        // Same as for Rx, there is no native gate that combines the effects neatly, so just conjugate the Rz case
                        let col_rot = tab.add_qubit();
                        q_index_map.insert(DataflowPoint::Rotation(node), col_rot);
                        tab.append_v(col_in);
                        tab.append_cx(col_in, col_rot);
                        let mut zz = BitVector::new(tab.nb_qubits);
                        zz.xor_bit(col_in);
                        zz.xor_bit(col_rot);
                        tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                        tab.append_v(col_in);
                        tab.append_x(col_in);
                        [col_in]
                    });
                } else {
                    // Rather than projecting away anything anti-commuting with Y (which requires the general method for arbitrary Paulis), likely quicker to map in and out of Z
                    self.apply_op_with(hugr, node, |tab, _, [col_in]| {
                        tab.append_v(col_in);
                        tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X)]);
                        tab.append_v(col_in);
                        tab.append_x(col_in);
                        [col_in]
                    });
                }
            }
            TketOp::Toffoli => {
                self.apply_op_with(hugr, node, |tab, _, [col_in0, col_in1, col_in2]| {
                    // Z/Z/X on each qubit commute through, so project away any stabilizer with X/X/Z
                    tab.project_cols_to_zero(&vec![
                        (col_in0, PauliXZ::X),
                        (col_in1, PauliXZ::X),
                        (col_in2, PauliXZ::Z),
                    ]);
                    [col_in0, col_in1, col_in2]
                });
            }
            TketOp::MeasureFree => {
                if settings.rotation_override {
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        q_index_map.insert(DataflowPoint::Rotation(node), col_in);
                        if settings.include_sum_types {
                            // IRSRole
                            // ZZ       Handled by qubit renaming
                            // XX Z     Handled by qubit renaming and CZ to copy X to Z
                            //  ZZX     Added as a new stab
                            let [ctrl, role] = tab.add_n_qubits();
                            q_index_map.insert(
                                DataflowPoint::SumOutPhControl(
                                    node,
                                    OutgoingPort::from(0),
                                    vec![],
                                    0,
                                ),
                                ctrl,
                            );
                            q_index_map.insert(DataflowPoint::RoleControl(node, 0), role);
                            tab.append_cz(col_in, role);
                            let mut zzx_z = BitVector::new(tab.nb_qubits);
                            zzx_z.xor_bit(col_in);
                            zzx_z.xor_bit(ctrl);
                            let mut zzx_x = BitVector::new(tab.nb_qubits);
                            zzx_x.xor_bit(role);
                            tab.add_stab(zzx_z, zzx_x, false);
                        }
                        // Without Sum types, we just need to track the input and the rotation; we can treat this as an identity wire so it is already implemented by qubit renaming
                        []
                    });
                } else if settings.include_sum_types {
                    // Either post-selection would project away Xs and store the outcome in +-Z, so it is enough to project the Xs and rename the qubit to the measurement outcome
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X)]);
                        q_index_map.insert(
                            DataflowPoint::SumOutPhControl(node, OutgoingPort::from(0), vec![], 0),
                            col_in,
                        );
                        []
                    });
                } else {
                    // When we don't track the measurement outcome, no stabilizer information remains; essentially the same as QFree
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        // Project out non-commuting rows and remove column from tableau
                        tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X), (col_in, PauliXZ::Z)]);
                        let moved_qb = tab.delete_qubit(col_in);
                        if moved_qb.is_some() {
                            let mq = moved_qb.unwrap();
                            let (removed_dfp, _) = q_index_map.remove_by_right(&mq).unwrap();
                            q_index_map.insert(removed_dfp, col_in);
                        }
                        []
                    });
                }
            }
            TketOp::QAlloc => {
                self.apply_op_with(hugr, node, |tab, _, []| {
                    let col_front: usize = tab.add_qubit();
                    // Add row for Z over col_front
                    let mut front_bv = BitVector::new(tab.nb_qubits);
                    front_bv.xor_bit(col_front);
                    tab.add_stab(front_bv, BitVector::new(tab.nb_qubits), false);
                    [col_front]
                });
            }
            TketOp::QFree => {
                self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                    // Project out non-commuting rows and remove column from tableau
                    tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X), (col_in, PauliXZ::Z)]);
                    let moved_qb = tab.delete_qubit(col_in);
                    if moved_qb.is_some() {
                        let mq = moved_qb.unwrap();
                        let (removed_dfp, _) = q_index_map.remove_by_right(&mq).unwrap();
                        q_index_map.insert(removed_dfp, col_in);
                    }
                    []
                });
            }
            TketOp::Reset => {
                self.apply_op_with(hugr, node, |tab, _, [col_in]| {
                    // Project out non-commuting rows
                    tab.project_cols_to_zero(&vec![(col_in, PauliXZ::X), (col_in, PauliXZ::Z)]);
                    // Reuse col_in for the output qubit
                    // Add row for Z over col_in
                    let mut bv = BitVector::new(tab.nb_qubits);
                    bv.xor_bit(col_in);
                    tab.add_stab(bv, BitVector::new(tab.nb_qubits), false);
                    [col_in]
                });
            }
            TketOp::V => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_v(col);
                });
            }
            TketOp::Vdg => {
                self.apply_clifford_with(hugr, node, |tab, [col]| {
                    tab.append_v(col);
                    tab.append_x(col);
                });
            }
            TketOp::TryQAlloc => {
                if settings.include_sum_types {
                    // Add the row for Z conditional on the classical result
                    let [out_q, ctrl] = self.tab.add_n_qubits();
                    let mut zz = BitVector::new(self.tab.nb_qubits);
                    zz.xor_bit(out_q);
                    zz.xor_bit(ctrl);
                    self.tab
                        .add_stab(zz, BitVector::new(self.tab.nb_qubits), false);
                    self.q_index_map.insert(
                        DataflowPoint::NodeOut(node, OutgoingPort::from(0), vec![0]),
                        out_q,
                    );
                    self.q_index_map.insert(
                        DataflowPoint::SumOutControl(node, OutgoingPort::from(0), vec![], 0, 0),
                        ctrl,
                    );
                }
                // If we aren't tracking data within sum types, then there is nothing to do
            }
            _ => {
                // In case other options are added later on, handle them as opaque unless we explicitly add a custom handler for them
                self.apply_opaque(hugr, node, settings.include_sum_types);
            }
        }
    }

    /// Used when DataflowSettings::include_external_interface == false and DataflowSettings::flow_level == DataflowFlowDetail::None, i.e. we want to treat this node as a completely decoherent channel that destroys all information incident on it
    /// O(InQ * Q * S) time, where InQ is the number of existing qubits associated with predecessors of node
    fn apply_decoherent(&mut self, hugr: &H, node: H::Node) {
        let preds: HashSet<(H::Node, OutgoingPort)> =
            HashSet::from_iter(hugr.all_linked_outputs(node));
        // Collect qubit indices to be removed in reverse order to maintain consistency when deleting qubits
        let qbs_to_remove = (0..self.tab.nb_qubits)
            .rev()
            .filter(|q| {
                if let DataflowPoint::TempOut(n, out_port, _) =
                    self.q_index_map.get_by_right(q).unwrap()
                {
                    preds.contains(&(*n, *out_port))
                } else {
                    false
                }
            })
            .collect_vec();
        // Remove any stabilizers involving the input qubits
        self.tab.project_cols_to_zero(
            &chain!(
                qbs_to_remove.iter().map(|q| (*q, PauliXZ::X)),
                qbs_to_remove.iter().map(|q| (*q, PauliXZ::Z)),
            )
            .collect_vec(),
        );
        for qb in qbs_to_remove {
            self.q_index_map.remove_by_right(&qb);
            if let Some(moved_qb) = self.tab.delete_qubit(qb) {
                let (moved_dfp, _) = self.q_index_map.remove_by_right(&moved_qb).unwrap();
                self.q_index_map.insert(moved_dfp, qb);
            }
        }
        // Add no new stabilizers for outputs
    }

    /// Takes the summary of a given node and appends it to the current tableau.
    /// Naively takes the tensor product of the two tableaus and performs a Bell post-selection to join TempOut/TempIn qubits, as well as adding ZZ stabilizers to propagate SumInPhControl/SumOutPhControl qubits.
    /// We assume flow has already been constructed according to the DataflowSettings for the node.
    /// In future, it may be worth avoiding the Bell post-selection by instead finding a set of gate applications that extend the existing stabilizers appropriately and then only adding the stabilizers independent of TempIns, since this could allow us to only use operations that are fast in qubit-major tableau implementations during summary construction.
    /// O(InQ * (self.Q + summary.Q) * (self.S + summary.S)) time for InQ qubits used in the input interface of node
    ///     Tensor product takes O(self.S * summary.Q + summary.S * (self.Q + summary.Q))
    ///     O(self.Q + summary.Q) to identify the post-selections, discards, and propagations to perform
    ///     For each qubit in the type signature of an input, it must be either post-selected or projected, so O(InQ * (self.Q + summary.Q) * (self.S + summary.S))
    ///     Adding stabilizers to propagate controls is O(InQ * (self.Q + summary.Q))
    fn append_node_summary(&mut self, hugr: &H, node: H::Node, summary: &StabilizerDataflow<H>) {
        self.tensor_product_in_place(summary);
        let preds: HashMap<(H::Node, OutgoingPort), IncomingPort> = HashMap::from_iter(
            hugr.in_value_types(node)
                .map(|(in_port, _)| (hugr.single_linked_output(node, in_port).unwrap(), in_port)),
        );
        // Find all Bell post-selections to be performed
        let mut bells: Vec<(DataflowPoint<H::Node>, DataflowPoint<H::Node>)> = vec![];
        // Sometimes only one side of a Bell post-selection is available (the other is implicitly the maximally mixed state) so we just discard them
        let mut halfbells: HashSet<DataflowPoint<H::Node>> = HashSet::new();
        // We only propagate the classical values if both exist
        let mut classical_propagations: Vec<(usize, usize)> = vec![];
        // Iterate through self.q_index_map in a deterministic order so each of the above happen in a fixed order
        for col in 0..self.tab.nb_qubits {
            let dfp = self.q_index_map.get_by_right(&col).unwrap();
            match dfp {
                DataflowPoint::TempIn(n, in_port, row_index) => {
                    if *n == node {
                        let (pred_node, pred_port) =
                            hugr.single_linked_output(node, *in_port).unwrap();
                        let counterpart =
                            DataflowPoint::TempOut(pred_node, pred_port, row_index.clone());
                        if halfbells.remove(&counterpart) {
                            bells.push((counterpart, dfp.clone()));
                        } else {
                            halfbells.insert(dfp.clone());
                        }
                    }
                }
                DataflowPoint::TempOut(n, out_port, row_index) => {
                    if let Some(in_port) = preds.get(&(*n, *out_port)) {
                        let counterpart = DataflowPoint::TempIn(node, *in_port, row_index.clone());
                        if halfbells.remove(&counterpart) {
                            bells.push((dfp.clone(), counterpart));
                        } else {
                            halfbells.insert(dfp.clone());
                        }
                    }
                }
                DataflowPoint::SumOutPhControl(n, out_port, row_index, sum_index) => {
                    if let Some(in_port) = preds.get(&(*n, *out_port)) {
                        if let Some(counterpart_col) =
                            self.q_index_map.get_by_left(&DataflowPoint::SumInPhControl(
                                node,
                                *in_port,
                                row_index.clone(),
                                *sum_index,
                            ))
                        {
                            classical_propagations.push((col, *counterpart_col));
                        }
                    }
                }
                _ => {}
            }
        }
        let mut post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        let mut delete_mask = BitVector::new(self.tab.nb_qubits);
        // Sort halfbells by qubit index so we get a deterministic order of the post-selections
        for q in halfbells {
            let col = self.q_index_map.get_by_left(&q).unwrap();
            delete_mask.xor_bit(*col);
        }
        let halfbell_indices: Vec<usize> = delete_mask.get_all_ones(self.tab.nb_qubits);
        // Transform Bell post-selections into post-selections on individual qubits
        for (l, r) in bells {
            let l_col = self.q_index_map.get_by_left(&l).unwrap();
            let r_col = self.q_index_map.get_by_left(&r).unwrap();
            self.tab.append_cx(*l_col, *r_col);
            post_selects.push((*l_col, PauliXZ::X, false));
            post_selects.push((*r_col, PauliXZ::Z, false));
            delete_mask.xor_bit(*l_col);
            delete_mask.xor_bit(*r_col);
        }
        for q in halfbell_indices {
            post_selects.push((q, PauliXZ::X, false));
            post_selects.push((q, PauliXZ::Z, false));
        }
        self.tab.post_select_1qs(&post_selects);
        // Propagate classical values
        for (l, r) in classical_propagations {
            let mut zz = BitVector::new(self.tab.nb_qubits);
            zz.xor_bit(l);
            zz.xor_bit(r);
            self.tab
                .add_stab(zz, BitVector::new(self.tab.nb_qubits), false);
        }
        // Remove post-selected qubits in order so to minimise issues from relabelling qubits
        for to_delete in delete_mask.get_all_ones(self.tab.nb_qubits).iter().rev() {
            self.q_index_map.remove_by_right(to_delete);
            if let Some(moved_qb) = self.tab.delete_qubit(*to_delete) {
                let (moved_dfp, _) = self.q_index_map.remove_by_right(&moved_qb).unwrap();
                self.q_index_map.insert(moved_dfp, *to_delete);
            }
        }
    }

    /// Calls apply_quantum_gate_flow_only to generate a standalone summary of the flow for a quantum gate.
    /// The resulting summary describes the flow between TempIns and TempOuts of the corresponding node.
    /// The settings configures usage of sum types and rotation overrides.
    /// O(node.in_ports^3 * IdQ * IdS)
    ///     Combining input wire descriptions gives node.in_ports * IdQ qubits and node.in_ports * IdS stabilizers
    ///     Building this using a fold operation takes O(node.in_ports^3 * IdQ * IdS) time (again, we can cut this down by building in-place)
    ///     Using the complexity of apply_quantum_gate_flow_only assuming a TketOp that is not handled as opaque, applying it takes O(node.in_ports^2 * IdQ * IdS) time
    fn generate_quantum_gate_flow(
        hugr: &H,
        node: H::Node,
        op: TketOp,
        settings: DataflowSettings,
    ) -> Self {
        let mut res = hugr.in_value_types(node).fold(
            StabilizerDataflow::default(),
            |acc, (in_port, in_type)| {
                let new_wire = if settings.include_sum_types || in_type == qb_t() {
                    let (pred, out_port) = hugr.single_linked_output(node, in_port).unwrap();
                    Self::identity_wire_recursive(
                        hugr,
                        node,
                        Either::Left(in_port),
                        pred,
                        Either::Right(out_port),
                        &in_type.into(),
                        &vec![],
                    )
                } else {
                    StabilizerDataflow::default()
                };
                Self::tensor_product(&acc, &new_wire)
            },
        );
        res.apply_quantum_gate_flow_only(hugr, node, op, settings);
        res
    }

    /// Generates a standalone summary for the interface to a given node, i.e. a set of identity wires for any qubit inputs and outputs (optionally including any of those within sum types).
    /// O((node.in_ports + node.out_ports)^3 * IdQ * IdS) time (could be reduced if we build the tensor product in-place)
    fn generate_interface(hugr: &H, node: H::Node, settings: DataflowSettings) -> Self {
        let in_interface = hugr.in_value_types(node).fold(
            StabilizerDataflow::default(),
            |acc, (in_port, in_type)| {
                let new_wire = if settings.include_sum_types || in_type == qb_t() {
                    Self::identity_wire_recursive(
                        hugr,
                        node,
                        Either::Left(in_port),
                        node,
                        Either::Left(in_port),
                        &in_type.into(),
                        &vec![],
                    )
                } else {
                    StabilizerDataflow::default()
                };
                Self::tensor_product(&acc, &new_wire)
            },
        );
        hugr.out_value_types(node)
            .fold(in_interface, |acc, (out_port, out_type)| {
                let new_wire = if settings.include_sum_types || out_type == qb_t() {
                    Self::identity_wire_recursive(
                        hugr,
                        node,
                        Either::Right(out_port),
                        node,
                        Either::Right(out_port),
                        &out_type.into(),
                        &vec![],
                    )
                } else {
                    StabilizerDataflow::default()
                };
                Self::tensor_product(&acc, &new_wire)
            })
    }

    // O(Q*S*mult + Q*S^2 + mult * log(mult)) time
    //      Post-selecting the controls takes O(Q*S*mult)
    //      Ordering the controls for deletion takes O(mult * log(mult))
    //      Deletion of qubits takes O(S*mult)
    //      When recombining, they have the same qubits, so O(Q*S^2)
    pub fn recalculate_control(
        &self,
        control_node: H::Node,
        control_port: Either<IncomingPort, OutgoingPort>,
        control_row_index: &Vec<usize>,
        control_sum_index: usize,
    ) -> Self {
        let mut positive_copy = self.tab.clone();
        let mut positive_post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        let mut negative_copy = self.tab.clone();
        let mut negative_post_selects: Vec<(usize, PauliXZ, bool)> = vec![];

        let phase_control_dfp = match control_port {
            Either::Left(in_port) => DataflowPoint::SumInPhControl(
                control_node,
                in_port,
                control_row_index.clone(),
                control_sum_index,
            ),
            Either::Right(out_port) => DataflowPoint::SumOutPhControl(
                control_node,
                out_port,
                control_row_index.clone(),
                control_sum_index,
            ),
        };
        let phase_control_index = self.q_index_map.get_by_left(&phase_control_dfp).unwrap();
        positive_post_selects.push((*phase_control_index, PauliXZ::Z, false));
        negative_post_selects.push((*phase_control_index, PauliXZ::Z, true));
        let mut to_delete: Vec<usize> = vec![*phase_control_index];
        let mut mult = 0;
        loop {
            let control_dfp = match control_port {
                Either::Left(in_port) => DataflowPoint::SumInControl(
                    control_node,
                    in_port,
                    control_row_index.clone(),
                    control_sum_index,
                    mult,
                ),
                Either::Right(out_port) => DataflowPoint::SumOutControl(
                    control_node,
                    out_port,
                    control_row_index.clone(),
                    control_sum_index,
                    mult,
                ),
            };
            match self.q_index_map.get_by_left(&control_dfp) {
                Some(control_index) => {
                    positive_post_selects.push((*control_index, PauliXZ::Z, false));
                    negative_post_selects.push((*control_index, PauliXZ::X, false));
                    to_delete.push(*control_index);
                }
                None => {
                    break;
                }
            }
            mult += 1;
        }

        positive_copy.post_select_1qs(&positive_post_selects);
        negative_copy.post_select_1qs(&negative_post_selects);

        // Remove post-selected qubits in order so to minimise issues from relabelling qubits
        to_delete.sort();
        let mut new_q_index_map = self.q_index_map.clone();
        for to_delete in to_delete.iter().rev() {
            new_q_index_map.remove_by_right(to_delete);
            let positive_moved = positive_copy.delete_qubit(*to_delete);
            let negative_moved = negative_copy.delete_qubit(*to_delete);
            if let Some(moved_qb) = positive_moved {
                let (moved_dfp, _) = new_q_index_map.remove_by_right(&moved_qb).unwrap();
                new_q_index_map.insert(moved_dfp, *to_delete);
            } else if let Some(moved_qb) = negative_moved {
                let (moved_dfp, _) = new_q_index_map.remove_by_right(&moved_qb).unwrap();
                new_q_index_map.insert(moved_dfp, *to_delete);
            }
        }

        let control_spec = match control_port {
            Either::Left(in_port) => {
                ControlSpec::SumInControl(in_port, control_row_index.clone(), control_sum_index)
            }
            Either::Right(out_port) => {
                ControlSpec::SumOutControl(out_port, control_row_index.clone(), control_sum_index)
            }
        };
        StabilizerDataflow::controlled_conditional(
            &StabilizerDataflow {
                tab: positive_copy,
                q_index_map: new_q_index_map.clone(),
            },
            &StabilizerDataflow {
                tab: negative_copy,
                q_index_map: new_q_index_map,
            },
            control_node,
            &control_spec,
            true,
        )
    }
}

/// Accumulates summaries of regions of a hugr
pub struct SDFAnalysis<H: HugrView> {
    targets: HashSet<H::Node>,
    summaries: HashMap<H::Node, StabilizerDataflow<H>>,
}

impl<H: HugrView> SDFAnalysis<H> {
    pub fn summarise_targets(
        hugr: &H,
        targets: &Vec<H::Node>,
        detail_callback: &mut impl FnMut(H::Node) -> DataflowSettings,
    ) -> HashMap<H::Node, StabilizerDataflow<H>> {
        let mut res = SDFAnalysis {
            targets: HashSet::from_iter(targets.iter().copied()),
            summaries: HashMap::new(),
        };
        for t in targets {
            // Some targets may have already been constructed when recursively building others
            if !res.summaries.contains_key(t) {
                match hugr.get_optype(*t) {
                    OpType::FuncDefn(_)
                    | OpType::DFG(_)
                    | OpType::DataflowBlock(_)
                    | OpType::Case(_)
                    | OpType::TailLoop(_) => {
                        // Treat all of these as just a DFG
                        let summary = res.run_dfg(hugr, *t, detail_callback);
                        res.summaries.insert(*t, summary);
                    }
                    _ => {
                        // Not a DFG node; just return an empty summary
                        res.summaries.insert(*t, StabilizerDataflow::default());
                    }
                }
            }
        }
        res.summaries
    }

    // Time complexity breakdown:
    //      Per node:
    //          If DataflowSettings::flow_level == DataflowFlowDetail::None, O(InQ * Q * S + node.out_ports * IdQ * S + node.out_ports * Q * IdS + node.out_ports^2 * IdQ * IdS)
    //          For TketOps or unknown ops, O(Q*S) in place or the breakdown below with interface
    //              O((node.in_ports + node.out_ports)^3 * IdQ * IdS) to generate the interface
    //              O(node.in_ports^3 * IdQ * IdS) to generate flow (trivial if unknown op)
    //              O((node.in_ports + node.out_ports)^3 * IdQ * IdS^2) to combine
    //                  union_Q = O((node.in_ports + node.out_ports) * IdQ)
    //                  left.S = O((node.in_ports + node.out_ports) * IdS)
    //                  right.S = O(node.in_ports * IdS)
    //              O(InQ * (Q + (node.in_ports + node.out_ports) * IdQ + (node.in_ports + node.out_ports) * IdS) * (S + (node.in_ports + node.out_ports) * IdS)) to apply the summary
    //          For LogicOps/BoolOps, O(Q + S) to add new qubits and stabilizers
    //          For conditionals
    //              Recursive calls to compute the branch summaries
    //              O(num_variants*branch.Q*branch.S^2) to project summaries if needed
    //              O(num_variants*branch.Q) to update row_index vectors of qubits
    //              O(num_variants*branch.Q*branch.S) to reorder qubits
    //              O(num_variants*branch.Q*branch.S^2) to compute invariants when DataflowSettings::flow_level == DataflowFlowDetail::Invariants
    //              When DataflowSettings::flow_level == DataflowFlowDetail::ConditionalInvs, branches have the same qubits but possibly completely different stabilizers
    //                  S' < num_variants * branch.S
    //                  Q' < branch.Q + 3 * num_variants * branch.S
    //                      The origin of the 3 is weird
    //                      When num_variants = 2^k, need up to (2^{k+1}-2)*branch.S control qubits
    //                      When num_variants = 2^k +1, need up to (2^{k+1}-2)*branch.S to merge into groups of 2^k, then the final merge has (2^k +1)*branch.S stabs to uniquely identify
    //                      (2^{k+1}-2)+(2^k +1) = 3*2^k - 1 < 3*(2^k +1)
    //                      Every new variant until the next power of 2 just adds up to branch.S new control qubits
    //                  Total time O(num_variants^3 * (log2 num_variants) * branch.Q * branch.S^2 + num_variants^4 * (log2 num_variants)^2 * branch.S^3)
    //                      Round 1: O(branch.Q * (2*branch.S)^2)
    //                      Round 2: O((branch.Q + 2*branch.S) * (4*branch.S)^2)
    //                      Round 3: O((branch.Q + 6*branch.S) * (8*branch.S)^2)
    //                      Round 4: O((branch.Q + 14*branch.S) * (16*branch.S)^2)
    //                      Round N: O(2^{2N} * branch.Q * branch.S^2 + 2^{3N} * branch.S^3)
    //                      Cumulative to Round N: O(N * 2^{3N} * branch.Q * branch.S^2 + N^2 * 2^{4N} * branch.S^3)
    //              When DataflowSettings::flow_level == DataflowFlowDetail::All, we assume the majority of qubits in each branch are completely disjoint, and hence so are the stabilizers
    //                  S' < num_variants * max_S
    //                  Q' < num_variants * max_Q + 3 * num_variants * max_S
    //                  Total time O(num_variants^3 * (log2 num_variants) * max_Q * max_S^2 + num_variants^4 * (log2 num_variants)^2 * max_S^3)
    //              O(InQ * (Q + num_variants * max_Q + 3 * num_variants * max_S) * (S + num_variants * max_S)) time to apply the summary
    //          For loops
    //              Recursive call to compute the body summary
    //              O(body.Q * body.S^2) to project to the continue/break branch, and project to flow over io
    //              O(InQ^3) to join the projected tableau (2*InQ qubits, up to same number of stabilizers) and the identity (same size)
    //              O(InQ * (body.Q + body.S + InQ)) to tensor product
    //              O(InQ * (body.Q + InQ) * (body.S + InQ)) to compose via Bell post-selection; same for eliminating the adjoining qubits
    //              Worst case includes adding role controls for both interface and body summary
    //                  O((node.in_ports + node.out_ports)^3 * IdQ * IdS) to generate interface tableau
    //                  O((body.Q + (node.in_ports + node.out_ports) * IdQ) * (body.S + (node.in_ports + node.out_ports) * IdS)) to tensor product
    //                  O((body.Q + (node.in_ports + node.out_ports) * IdQ) * (body.S + (node.in_ports + node.out_ports) * IdS)^2) to take conditional
    //              O(InQ * (Q + body.Q + body.S + (node.in_ports + node.out_ports) * (IdQ + IdS)) * (S + body.S + (node.in_ports + node.out_ports) * IdS)) to apply overall loop summary to the context summary
    //          For function calls
    //              Recursive call to compute the summary
    //              O(fun.Q*fun.S^2) to project the function summary to flow over IO
    //              Same cost as TketOps to apply a summary
    //          For tags O(InQ*Q + S)
    //              O(Q+S) to add the injections
    //              O(Q) to update all of the row_index vectors of qubits
    //              O(InQ*Q) to propagate classical values
    //
    // The worst node kinds for time complexity are conditional nodes; where size of the hugr will be roughly O(num_variants * max_Q), building the summary is close to quartic in this size
    fn run_dfg(
        &mut self,
        hugr: &H,
        parent: H::Node,
        detail_callback: &mut impl FnMut(H::Node) -> DataflowSettings,
    ) -> StabilizerDataflow<H> {
        let inp = hugr
            .children(parent)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let mut summary = StabilizerDataflow::initialise_from_input(
            hugr,
            inp,
            detail_callback(inp).include_sum_types,
        );
        let (region, node_map) = hugr.region_portgraph(parent);
        let mut topo = pv::Topo::new(&region);
        while let Some(pgnode) = topo.next(&region) {
            let node = node_map.from_portgraph(pgnode);
            let node_settings = detail_callback(node);
            if node_settings.flow_level == DataflowFlowDetail::None {
                // With no flow information, the action does not depend on the semantics of the node
                if node_settings.include_external_interface {
                    summary.apply_opaque(hugr, node, node_settings.include_sum_types);
                } else {
                    summary.apply_decoherent(hugr, node);
                }
            } else {
                let optype: &OpType = hugr.get_optype(node);
                match optype {
                    OpType::ExtensionOp(op) => {
                        if let Ok(tkop) = TketOp::from_extension_op(op) {
                            if node_settings.requires_role_control(optype) {
                                let interface = StabilizerDataflow::generate_interface(
                                    hugr,
                                    node,
                                    node_settings,
                                );
                                let flow = StabilizerDataflow::generate_quantum_gate_flow(
                                    hugr,
                                    node,
                                    tkop,
                                    node_settings,
                                );
                                let controlled = StabilizerDataflow::controlled_conditional(
                                    &interface,
                                    &flow,
                                    node,
                                    &ControlSpec::Role,
                                    false,
                                );
                                summary.append_node_summary(hugr, node, &controlled);
                            } else {
                                summary.apply_quantum_gate_flow_only(
                                    hugr,
                                    node,
                                    tkop,
                                    node_settings,
                                );
                            }
                        } else if let Ok(logic_op) = LogicOp::from_extension_op(op) {
                            if node_settings.include_sum_types {
                                match logic_op {
                                    LogicOp::Not => {
                                        let [sumin, sumout] = summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        let mut zz = BitVector::new(summary.tab.nb_qubits);
                                        zz.xor_bit(sumin);
                                        zz.xor_bit(sumout);
                                        summary.tab.add_stab(
                                            zz,
                                            BitVector::new(summary.tab.nb_qubits),
                                            true,
                                        );
                                        let (pred_node, pred_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred_node,
                                                pred_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                    }
                                    LogicOp::Xor | LogicOp::Eq => {
                                        // Eq is just Xor followed by negation
                                        let [sumin0, sumin1, sumout] = summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin0,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(1),
                                                vec![],
                                                0,
                                            ),
                                            sumin1,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        let mut zzz = BitVector::new(summary.tab.nb_qubits);
                                        zzz.xor_bit(sumin0);
                                        zzz.xor_bit(sumin1);
                                        zzz.xor_bit(sumout);
                                        summary.tab.add_stab(
                                            zzz,
                                            BitVector::new(summary.tab.nb_qubits),
                                            logic_op == LogicOp::Eq,
                                        );
                                        let (pred0_node, pred0_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred0_node,
                                                pred0_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin0);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                        let (pred1_node, pred1_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(1))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred1_node,
                                                pred1_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin1);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                    }
                                    LogicOp::And | LogicOp::Or => {
                                        // And and Or require the same sets of control qubits based on short-circuit definitions
                                        let [sumin0, sumin1, sumout, sumctrl0, sumctrl1] =
                                            summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin0,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(1),
                                                vec![],
                                                0,
                                            ),
                                            sumin1,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                                0,
                                            ),
                                            sumctrl0,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                                1,
                                            ),
                                            sumctrl1,
                                        );
                                        // AND(false, _) = false
                                        // AND(true, x) = x
                                        // OR(false, x) = x
                                        // OR(true, _) = true
                                        let mut false_z = BitVector::new(summary.tab.nb_qubits);
                                        let mut true_z = BitVector::new(summary.tab.nb_qubits);
                                        let mut true_x = BitVector::new(summary.tab.nb_qubits);
                                        false_z.xor_bit(sumout);
                                        false_z.xor_bit(sumctrl0);
                                        true_z.xor_bit(sumout);
                                        true_x.xor_bit(sumctrl1);
                                        let true_phase = if logic_op == LogicOp::And {
                                            true_z.xor_bit(sumin1);
                                            false
                                        } else {
                                            false_z.xor_bit(sumin1);
                                            true
                                        };
                                        summary.tab.add_stab(
                                            false_z,
                                            BitVector::new(summary.tab.nb_qubits),
                                            false,
                                        );
                                        summary.tab.add_stab(true_z, true_x, true_phase);
                                        let (pred0_node, pred0_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred0_node,
                                                pred0_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin0);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                        let (pred1_node, pred1_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(1))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred1_node,
                                                pred1_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin1);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                    }
                                    _ => {
                                        // LogicOps will always only involve classical data, so we just omit the classical values if we don't know the semantics of the op
                                    }
                                }
                            }
                        } else if let Ok(bool_op) = BoolOp::from_extension_op(op) {
                            if node_settings.include_sum_types {
                                match bool_op {
                                    BoolOp::not => {
                                        let [sumin, sumout] = summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        let mut zz = BitVector::new(summary.tab.nb_qubits);
                                        zz.xor_bit(sumin);
                                        zz.xor_bit(sumout);
                                        summary.tab.add_stab(
                                            zz,
                                            BitVector::new(summary.tab.nb_qubits),
                                            true,
                                        );
                                        let (pred_node, pred_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred_node,
                                                pred_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                    }
                                    BoolOp::xor | BoolOp::eq => {
                                        // Eq is just Xor followed by negation
                                        let [sumin0, sumin1, sumout] = summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin0,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(1),
                                                vec![],
                                                0,
                                            ),
                                            sumin1,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        let mut zzz = BitVector::new(summary.tab.nb_qubits);
                                        zzz.xor_bit(sumin0);
                                        zzz.xor_bit(sumin1);
                                        zzz.xor_bit(sumout);
                                        summary.tab.add_stab(
                                            zzz,
                                            BitVector::new(summary.tab.nb_qubits),
                                            bool_op == BoolOp::eq,
                                        );
                                        let (pred0_node, pred0_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred0_node,
                                                pred0_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin0);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                        let (pred1_node, pred1_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(1))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred1_node,
                                                pred1_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin1);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                    }
                                    BoolOp::and | BoolOp::or => {
                                        // And and Or require the same sets of control qubits based on short-circuit definitions
                                        let [sumin0, sumin1, sumout, sumctrl0, sumctrl1] =
                                            summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin0,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(1),
                                                vec![],
                                                0,
                                            ),
                                            sumin1,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                                0,
                                            ),
                                            sumctrl0,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                                1,
                                            ),
                                            sumctrl1,
                                        );
                                        // AND(false, _) = false
                                        // AND(true, x) = x
                                        // OR(false, x) = x
                                        // OR(true, _) = true
                                        let mut false_z = BitVector::new(summary.tab.nb_qubits);
                                        let mut true_z = BitVector::new(summary.tab.nb_qubits);
                                        let mut true_x = BitVector::new(summary.tab.nb_qubits);
                                        false_z.xor_bit(sumout);
                                        false_z.xor_bit(sumctrl0);
                                        true_z.xor_bit(sumout);
                                        true_x.xor_bit(sumctrl1);
                                        let true_phase = if bool_op == BoolOp::and {
                                            true_z.xor_bit(sumin1);
                                            false
                                        } else {
                                            false_z.xor_bit(sumin1);
                                            true
                                        };
                                        summary.tab.add_stab(
                                            false_z,
                                            BitVector::new(summary.tab.nb_qubits),
                                            false,
                                        );
                                        summary.tab.add_stab(true_z, true_x, true_phase);
                                        let (pred0_node, pred0_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred0_node,
                                                pred0_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin0);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                        let (pred1_node, pred1_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(1))
                                            .unwrap();
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred1_node,
                                                pred1_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin1);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                    }
                                    BoolOp::read | BoolOp::make_opaque => {
                                        // Both of these act as an identity, just converting between the opaque boolean and Sum[unit, unit].
                                        // Because of the way we handle classical values via propagations, we need to introduce two new sum controls (one for the input, one for the output) and propagate the values all the way along
                                        let (pred_node, pred_port) = hugr
                                            .single_linked_output(node, IncomingPort::from(0))
                                            .unwrap();
                                        let [sumin, sumout] = summary.tab.add_n_qubits();
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumInPhControl(
                                                node,
                                                IncomingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumin,
                                        );
                                        summary.q_index_map.insert(
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                vec![],
                                                0,
                                            ),
                                            sumout,
                                        );
                                        if let Some(q) = summary.q_index_map.get_by_left(
                                            &DataflowPoint::SumOutPhControl(
                                                pred_node,
                                                pred_port,
                                                vec![],
                                                0,
                                            ),
                                        ) {
                                            // Propagate value from predecessor into the input
                                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                                            zz.xor_bit(*q);
                                            zz.xor_bit(sumin);
                                            summary.tab.add_stab(
                                                zz,
                                                BitVector::new(summary.tab.nb_qubits),
                                                false,
                                            );
                                        }
                                        // Propagate the value from the input to the output
                                        let mut zz = BitVector::new(summary.tab.nb_qubits);
                                        zz.xor_bit(sumin);
                                        zz.xor_bit(sumout);
                                        summary.tab.add_stab(
                                            zz,
                                            BitVector::new(summary.tab.nb_qubits),
                                            false,
                                        );
                                    }
                                    _ => {
                                        // BoolOps will always only involve classical data, so we just omit the classical values if we don't know the semantics of the op
                                    }
                                }
                            }
                        } else {
                            if node_settings.include_external_interface {
                                let interface = StabilizerDataflow::generate_interface(
                                    hugr,
                                    node,
                                    node_settings,
                                );
                                let controlled = StabilizerDataflow::controlled_conditional(
                                    &interface,
                                    &StabilizerDataflow::default(),
                                    node,
                                    &ControlSpec::Role,
                                    false,
                                );
                                summary.append_node_summary(hugr, node, &controlled);
                            } else {
                                summary.apply_decoherent(hugr, node);
                            }
                        }
                    }
                    OpType::Conditional(_) => {
                        summary.append_node_summary(
                            hugr,
                            node,
                            &self.run_conditional(hugr, node, node_settings, detail_callback),
                        );
                    }
                    OpType::TailLoop(_) => {
                        summary.append_node_summary(
                            hugr,
                            node,
                            &self.run_tail_loop(hugr, node, node_settings, detail_callback),
                        );
                    }
                    OpType::Call(_) => {
                        // Even if node_settings.flow_level == DataflowFlowDetail::All, we do not hoist out of function calls, so we only ever look at invariants over the function body
                        let call_port = optype.static_input_port().unwrap();
                        let (fun_def_node, _) = hugr
                            .linked_outputs(node, call_port)
                            .exactly_one()
                            .ok()
                            .unwrap();
                        let mut fun_body = if self.targets.contains(&fun_def_node) {
                            if let Some(fun_body) = self.summaries.get(&fun_def_node) {
                                fun_body.clone()
                            } else {
                                let fun_body = self.run_dfg(hugr, fun_def_node, detail_callback);
                                self.summaries.insert(fun_def_node, fun_body.clone());
                                fun_body
                            }
                        } else {
                            self.run_dfg(hugr, fun_def_node, detail_callback)
                        };
                        let fun_inp = hugr
                            .children(parent)
                            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
                            .exactly_one()
                            .ok()
                            .unwrap();
                        let fun_outp = hugr
                            .children(parent)
                            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
                            .exactly_one()
                            .ok()
                            .unwrap();
                        fun_body.project_to_flow();
                        fun_body.project_non_io(fun_inp, fun_outp);
                        summary.append_node_summary(hugr, node, &fun_body);
                    }
                    OpType::Input(_) => {
                        // Already handled during setup
                    }
                    OpType::Output(_) => {
                        // Qubit relabelling is handled by treating this node as opaque
                        summary.apply_opaque(hugr, node, node_settings.include_sum_types);
                    }
                    OpType::Tag(tag) => {
                        // Inject
                        let mut sum_index = 0;
                        let mut sum_bit_mask = 1;
                        while sum_bit_mask < tag.variants.len() {
                            summary.sum_injection(
                                if tag.tag & sum_bit_mask == 0 {
                                    Either::Left(())
                                } else {
                                    Either::Right(())
                                },
                                node,
                                OutgoingPort::from(0),
                                vec![],
                                sum_index,
                            );
                            sum_index += 1;
                            sum_bit_mask *= 2;
                        }
                        // Rename the inputs to the single output port, pushing the input port to their row index
                        let pred_map: HashMap<(H::Node, OutgoingPort), IncomingPort> =
                            HashMap::from_iter(hugr.in_value_types(node).map(|(in_port, _)| {
                                (hugr.single_linked_output(node, in_port).unwrap(), in_port)
                            }));
                        let mut renames: Vec<(usize, DataflowPoint<H::Node>)> = vec![];
                        let mut classical_propagations: Vec<(usize, DataflowPoint<H::Node>)> =
                            vec![];
                        // Iterate through summary.q_index_map in a deterministic order so classical_propagations add stabilizers in a fixed order
                        for col in 0..summary.tab.nb_qubits {
                            let dfp = summary.q_index_map.get_by_right(&col).unwrap();
                            match dfp {
                                DataflowPoint::TempOut(n, out_port, ri) => {
                                    if let Some(in_port) = pred_map.get(&(*n, *out_port)) {
                                        let mut new_ri = vec![in_port.index()];
                                        new_ri.extend(ri);
                                        renames.push((
                                            col,
                                            DataflowPoint::TempOut(
                                                node,
                                                OutgoingPort::from(0),
                                                new_ri,
                                            ),
                                        ));
                                    }
                                }
                                DataflowPoint::SumOutPhControl(n, out_port, ri, si) => {
                                    if let Some(in_port) = pred_map.get(&(*n, *out_port)) {
                                        let mut new_ri = vec![in_port.index()];
                                        new_ri.extend(ri);
                                        classical_propagations.push((
                                            col,
                                            DataflowPoint::SumOutPhControl(
                                                node,
                                                OutgoingPort::from(0),
                                                new_ri,
                                                *si,
                                            ),
                                        ));
                                    }
                                }
                                _ => {}
                            }
                        }
                        for (col, dfp) in renames {
                            summary.q_index_map.remove_by_right(&col);
                            summary.q_index_map.insert(dfp, col);
                        }
                        // Propagate any classical values from the inputs to the outputs, pushing the input port to their row index
                        for (prev_col, new_dfp) in classical_propagations {
                            let new_q = summary.tab.add_qubit();
                            summary.q_index_map.insert(new_dfp, new_q);
                            let mut zz = BitVector::new(summary.tab.nb_qubits);
                            zz.xor_bit(prev_col);
                            zz.xor_bit(new_q);
                            summary
                                .tab
                                .add_stab(zz, BitVector::new(summary.tab.nb_qubits), false);
                        }
                    }
                    _ => {
                        if node_settings.include_external_interface {
                            // Even though no flow information is known, the caller has still asked for the ability to toggle the interface information off
                            let interface =
                                StabilizerDataflow::generate_interface(hugr, node, node_settings);
                            let controlled = StabilizerDataflow::controlled_conditional(
                                &interface,
                                &StabilizerDataflow::default(),
                                node,
                                &ControlSpec::Role,
                                false,
                            );
                            summary.append_node_summary(hugr, node, &controlled);
                        } else {
                            summary.apply_decoherent(hugr, node);
                        }
                    }
                }
            }
        }
        summary
    }

    /// Constructs the summary for a conditional node.
    /// Assumes cond_settings.flow_level != DataflowFlowDetail::None.
    fn run_conditional(
        &mut self,
        hugr: &H,
        node: H::Node,
        cond_settings: DataflowSettings,
        detail_callback: &mut impl FnMut(H::Node) -> DataflowSettings,
    ) -> StabilizerDataflow<H> {
        let cond_op = hugr.get_optype(node).as_conditional().unwrap();
        let mut bodies_to_join: Vec<StabilizerDataflow<H>> = hugr
            .children(node)
            .zip(cond_op.sum_rows.iter())
            .map(|(cond_node, case_row)| {
                let mut cond_body = if self.targets.contains(&cond_node) {
                    if let Some(body) = self.summaries.get(&cond_node) {
                        body.clone()
                    } else {
                        let body = self.run_dfg(hugr, cond_node, detail_callback);
                        self.summaries.insert(cond_node, body.clone());
                        body
                    }
                } else {
                    self.run_dfg(hugr, cond_node, detail_callback)
                };
                let body_in = hugr
                    .children(cond_node)
                    .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
                    .exactly_one()
                    .ok()
                    .unwrap();
                let body_out = hugr
                    .children(cond_node)
                    .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
                    .exactly_one()
                    .ok()
                    .unwrap();
                if cond_settings.flow_level != DataflowFlowDetail::All {
                    cond_body.project_to_flow();
                    cond_body.project_non_io(body_in, body_out);
                }
                let new_qubit_labels: Vec<(DataflowPoint<H::Node>, usize)> = cond_body
                    .q_index_map
                    .iter()
                    .filter_map(|(dfp, q)| match dfp {
                        DataflowPoint::NodeIn(n, i, ri) => {
                            if *n == body_out {
                                Some((
                                    DataflowPoint::TempOut(
                                        node,
                                        OutgoingPort::from(i.index()),
                                        ri.clone(),
                                    ),
                                    *q,
                                ))
                            } else {
                                None
                            }
                        }
                        DataflowPoint::SumInControl(n, i, ri, si, m) => {
                            if *n == body_out {
                                Some((
                                    DataflowPoint::SumOutControl(
                                        node,
                                        OutgoingPort::from(i.index()),
                                        ri.clone(),
                                        *si,
                                        *m,
                                    ),
                                    *q,
                                ))
                            } else {
                                None
                            }
                        }
                        DataflowPoint::SumInPhControl(n, i, ri, si) => {
                            if *n == body_out {
                                Some((
                                    DataflowPoint::SumOutPhControl(
                                        node,
                                        OutgoingPort::from(i.index()),
                                        ri.clone(),
                                        *si,
                                    ),
                                    *q,
                                ))
                            } else {
                                None
                            }
                        }
                        DataflowPoint::NodeOut(n, o, ri) => {
                            if *n == body_in {
                                if o.index() < case_row.len() {
                                    let mut new_ri: Vec<usize> = vec![o.index()];
                                    new_ri.extend(ri);
                                    Some((
                                        DataflowPoint::TempIn(node, IncomingPort::from(0), new_ri),
                                        *q,
                                    ))
                                } else {
                                    Some((
                                        DataflowPoint::TempIn(
                                            node,
                                            IncomingPort::from(o.index() + 1),
                                            ri.clone(),
                                        ),
                                        *q,
                                    ))
                                }
                            } else {
                                None
                            }
                        }
                        DataflowPoint::SumOutControl(n, o, ri, si, m) => {
                            if *n == body_in {
                                if o.index() < case_row.len() {
                                    let mut new_ri: Vec<usize> = vec![o.index()];
                                    new_ri.extend(ri);
                                    Some((
                                        DataflowPoint::SumInControl(
                                            node,
                                            IncomingPort::from(0),
                                            new_ri,
                                            *si,
                                            *m,
                                        ),
                                        *q,
                                    ))
                                } else {
                                    Some((
                                        DataflowPoint::SumInControl(
                                            node,
                                            IncomingPort::from(o.index() + 1),
                                            ri.clone(),
                                            *si,
                                            *m,
                                        ),
                                        *q,
                                    ))
                                }
                            } else {
                                None
                            }
                        }
                        DataflowPoint::SumOutPhControl(n, o, ri, si) => {
                            if *n == body_in {
                                if o.index() < case_row.len() {
                                    let mut new_ri: Vec<usize> = vec![o.index()];
                                    new_ri.extend(ri);
                                    Some((
                                        DataflowPoint::SumInPhControl(
                                            node,
                                            IncomingPort::from(0),
                                            new_ri,
                                            *si,
                                        ),
                                        *q,
                                    ))
                                } else {
                                    Some((
                                        DataflowPoint::SumInPhControl(
                                            node,
                                            IncomingPort::from(o.index() + 1),
                                            ri.clone(),
                                            *si,
                                        ),
                                        *q,
                                    ))
                                }
                            } else {
                                None
                            }
                        }
                        _ => None,
                    })
                    .collect_vec();
                for (new_dfp, q) in new_qubit_labels {
                    cond_body.q_index_map.remove_by_right(&q);
                    cond_body.q_index_map.insert(new_dfp, q);
                }
                cond_body
            })
            .collect_vec();
        if bodies_to_join.is_empty() {
            // I'm not sure if hugr permits empty Sums to be plugged into conditional nodes, but just in case...
            return StabilizerDataflow::default();
        }
        // Reorder qubits in all of the bodies to a unified order before we begin making joins
        let mut unified_q_index: BiHashMap<DataflowPoint<H::Node>, usize> = BiHashMap::new();
        for body in bodies_to_join.iter() {
            for q in 0..(body.tab.nb_qubits) {
                let dfp = body.q_index_map.get_by_right(&q).unwrap();
                if !unified_q_index.contains_left(dfp) {
                    unified_q_index.insert(dfp.clone(), unified_q_index.len());
                }
            }
        }
        let mut qubit_count = unified_q_index.len();
        bodies_to_join = bodies_to_join
            .iter()
            .map(|body| {
                let mut reordered_body = StabilizerDataflow {
                    tab: SymplecticTableau::new(qubit_count),
                    q_index_map: unified_q_index.clone(),
                };
                for i in 0..body.tab.nb_stabs {
                    let mut new_z = BitVector::new(qubit_count);
                    let mut new_x = BitVector::new(qubit_count);
                    for z_index in body.tab.z[i].get_all_ones(body.tab.nb_qubits) {
                        let new_index = unified_q_index
                            .get_by_left(body.q_index_map.get_by_right(&z_index).unwrap())
                            .unwrap();
                        new_z.xor_bit(*new_index);
                    }
                    for x_index in body.tab.x[i].get_all_ones(body.tab.nb_qubits) {
                        let new_index = unified_q_index
                            .get_by_left(body.q_index_map.get_by_right(&x_index).unwrap())
                            .unwrap();
                        new_x.xor_bit(*new_index);
                    }
                    reordered_body
                        .tab
                        .add_stab(new_z, new_x, body.tab.signs.get(i));
                }
                reordered_body
            })
            .collect_vec();
        // Now take intersections or controlled conditionals of each pair until a single summary remains
        let final_flow = if cond_settings.flow_level == DataflowFlowDetail::Invariants {
            while bodies_to_join.len() > 1 {
                bodies_to_join = bodies_to_join
                    .chunks(2)
                    .map(|body_array| {
                        if body_array.len() == 1 {
                            body_array[0].clone()
                        } else {
                            // NOTE:: Suppose we have some classical variable X flowing into our branches, and all branches contain P when X is true (similarly if all false). P may not exist in the intersection here since it may be associated with a different subset of the control qubits for X between the different branches. In this sense, even though the information is practically invariant on our choice of branches, we accept under DataflowFlowDetail::Invariants that we may lose information dependent on any classical variable, so this is an acceptable loss.
                            let mut jd = SymplecticTableau::joint_decomposition(
                                &body_array[0].tab,
                                &body_array[1].tab,
                            );
                            let (join_begin, _) = jd.join_range();
                            for i in (0..join_begin).rev() {
                                jd.tab.delete_stab(i);
                            }
                            StabilizerDataflow {
                                tab: jd.tab,
                                q_index_map: unified_q_index.clone(),
                            }
                        }
                    })
                    .collect_vec();
            }
            bodies_to_join[0].clone()
        } else {
            // At each sum-index value, the various controlled_conditional's may introduce different numbers of control qubits.
            // These controls are always added in the same order (following the multiplicity value) with the same control spec, so to unify the qubit indexing again we will only need to pad them out to the same number of control qubits.
            let mut sum_index = 0;
            while bodies_to_join.len() > 1 {
                let control_spec =
                    ControlSpec::SumInControl(IncomingPort::from(0), vec![], sum_index);
                let mut max_qbs = 0;
                bodies_to_join = bodies_to_join
                    .chunks(2)
                    .map(|body_array| {
                        if body_array.len() == 1 {
                            max_qbs = max(max_qbs, body_array[0].tab.nb_qubits);
                            body_array[0].clone()
                        } else {
                            let combined = StabilizerDataflow::controlled_conditional(
                                &body_array[0],
                                &body_array[1],
                                node,
                                &control_spec,
                                true,
                            );
                            max_qbs = max(max_qbs, combined.tab.nb_qubits);
                            combined
                        }
                    })
                    .collect_vec();
                if max_qbs != qubit_count {
                    // At least one control was added somewhere
                    for body in bodies_to_join.iter_mut() {
                        if body.tab.nb_qubits == qubit_count {
                            // This body did not have any controls, and in particular does not have a phase control, so add that first
                            body.tab.add_qubit();
                            body.q_index_map.insert(
                                DataflowPoint::SumInPhControl(
                                    node,
                                    IncomingPort::from(0),
                                    vec![],
                                    sum_index,
                                ),
                                qubit_count,
                            );
                        }
                        if body.tab.nb_qubits < max_qbs {
                            for i in body.tab.nb_qubits..max_qbs {
                                body.q_index_map.insert(
                                    DataflowPoint::SumInControl(
                                        node,
                                        IncomingPort::from(0),
                                        vec![],
                                        sum_index,
                                        i - (qubit_count + 1),
                                    ),
                                    i,
                                );
                            }
                            body.tab.add_qubits(max_qbs - body.tab.nb_qubits);
                        }
                    }
                    qubit_count = max_qbs;
                }
                sum_index += 1;
            }
            bodies_to_join[0].clone()
        };
        if cond_settings.include_external_interface {
            let interface = StabilizerDataflow::generate_interface(hugr, node, cond_settings);
            StabilizerDataflow::controlled_conditional(
                &interface,
                &final_flow,
                node,
                &ControlSpec::Role,
                false,
            )
        } else {
            final_flow
        }
    }

    /// Constructs the summary for a loop node.
    /// Assumes cond_settings.flow_level != DataflowFlowDetail::None.
    /// Current implementation does not gather any information for hoisting, just loop invariants
    fn run_tail_loop(
        &mut self,
        hugr: &H,
        node: H::Node,
        loop_settings: DataflowSettings,
        detail_callback: &mut impl FnMut(H::Node) -> DataflowSettings,
    ) -> StabilizerDataflow<H> {
        let loop_body = if self.targets.contains(&node) {
            if let Some(body) = self.summaries.get(&node) {
                body.clone()
            } else {
                let body = self.run_dfg(hugr, node, detail_callback);
                self.summaries.insert(node, body.clone());
                body
            }
        } else {
            self.run_dfg(hugr, node, detail_callback)
        };
        let loop_op = hugr.get_optype(node).as_tail_loop().unwrap();
        println!("{}", loop_body.tab);
        for i in 0..loop_body.tab.nb_qubits {
            println!("{}\t{:?}", i, loop_body.q_index_map.get_by_right(&i).unwrap());
        }
        // Gates can be hoisted out of a loop in multiple ways:
        // - Identifying that it is only run on the first/final loop iteration, when it has a canonical effect wrt the inputs/outputs of the loop (because we use TailLoops, at least one iteration is always performed, so this gate is performed precisely once)
        // - Folding all copies of a rotation gate into one another (this always implies it has a canonical effect wrt the inputs to the loop, and therefore can only be hoisted out the end of the loop when there is a loop invariant)
        // - More complicated variants of this where the phase of the rotation gate might vary between iterations (so the angle accumulator may need to add or subtract the angle on each iteration), or where the gate alternates between effects on a small number of terms
        // - Folding invertible gates (e.g. CCX) to hoist it as a conditional gate applied at most one time
        // - Probably more to come
        // Each of these might require more information to be included in the tableau generated here. For the time being, we will limit ourselves to just folding basic rotation gates (any gate type that uses the rotation overload).

        // We always want to find the loop invariants, so find those first.
        // The loop execution is always zero or more copies of the loop body post-selected to the "continue" branch, followed by exactly one copy of the loop body post-selected to the "break" branch.
        // We will calculate invariants here by finding the invariants of the "continue" branch as a region-based summary (i.e. the stabilizers it shares with the identity), then compose this with the flow information over the "break" branch.
        // Post-selections in the "break" branch might propagate backwards through the "continue"s in an alternating fashion which are proper invariants but would not be detected through this method. They may be detectable if we replace this with an iterative procedure in the future.
        let body_in = hugr
            .children(node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let body_out = hugr
            .children(node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let mut continue_flow = loop_body.clone();
        let mut continue_post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        let mut break_post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        if let Some(q) = continue_flow
            .q_index_map
            .get_by_left(&DataflowPoint::SumInPhControl(
                body_out,
                IncomingPort::from(0),
                vec![],
                0,
            ))
        {
            continue_post_selects.push((*q, PauliXZ::Z, false));
            break_post_selects.push((*q, PauliXZ::Z, true));
        }
        let mut multiplicity = 0;
        while let Some(q) = continue_flow
            .q_index_map
            .get_by_left(&DataflowPoint::SumInControl(
                body_out,
                IncomingPort::from(0),
                vec![],
                0,
                multiplicity,
            ))
        {
            continue_post_selects.push((*q, PauliXZ::Z, false));
            break_post_selects.push((*q, PauliXZ::X, false));
            multiplicity += 1;
        }
        continue_flow.tab.post_select_1qs(&continue_post_selects);
        continue_flow.project_to_flow();
        // This projection might remove information on how classical values propagate, preventing us from spotting some of the classical invariants of the loop. We leave the expensive propagation of classical information for future development
        continue_flow.project_non_io(body_in, body_out);
        let mut simple_identity_tab = SymplecticTableau::new(continue_flow.tab.nb_qubits);
        // These sets will be useful for us to compose onto the break case
        let mut bells: Vec<(usize, DataflowPoint<H::Node>)> = vec![];
        let mut classical_propagations: Vec<(usize, DataflowPoint<H::Node>)> = vec![];
        for q in 0..continue_flow.tab.nb_qubits {
            let input_dfp: Option<DataflowPoint<H::Node>> =
                match continue_flow.q_index_map.get_by_right(&q).unwrap() {
                    DataflowPoint::NodeIn(_, in_port, ri) => {
                        // Since we have already called project_non_io, all NodeIns must be from body_out
                        let counterpart = if in_port.index() == 0 {
                            DataflowPoint::NodeOut(
                                body_in,
                                OutgoingPort::from(ri[0]),
                                ri.as_slice()[1..].to_vec(),
                            )
                        } else {
                            DataflowPoint::NodeOut(
                                body_in,
                                OutgoingPort::from(in_port.index() + loop_op.just_inputs.len() - 1),
                                ri.clone(),
                            )
                        };
                        bells.push((q, counterpart.clone()));
                        Some(counterpart)
                    }
                    DataflowPoint::SumInPhControl(_, in_port, ri, si) => {
                        // Similarly, the SumInPhControl must be from body_out
                        if in_port.index() == 0 {
                            if ri.is_empty() {
                                // Controls for the top-level sum have no correspondant
                                None
                            } else {
                                let counterpart = DataflowPoint::SumOutPhControl(
                                    body_in,
                                    OutgoingPort::from(ri[0]),
                                    ri.as_slice()[1..].to_vec(),
                                    *si,
                                );
                                classical_propagations.push((q, counterpart.clone()));
                                Some(counterpart)
                            }
                        } else {
                            let counterpart = DataflowPoint::SumOutPhControl(
                                body_in,
                                OutgoingPort::from(in_port.index() + loop_op.just_inputs.len() - 1),
                                ri.clone(),
                                *si,
                            );
                            classical_propagations.push((q, counterpart.clone()));
                            Some(counterpart)
                        }
                    }
                    _ => None,
                };
            if let Some(in_dfp) = input_dfp {
                let mut zz = BitVector::new(simple_identity_tab.nb_qubits);
                zz.xor_bit(q);
                let input_q = continue_flow.q_index_map.get_by_left(&in_dfp).unwrap();
                zz.xor_bit(*input_q);
                simple_identity_tab.add_stab(
                    zz.clone(),
                    BitVector::new(simple_identity_tab.nb_qubits),
                    false,
                );
                simple_identity_tab.add_stab(
                    BitVector::new(simple_identity_tab.nb_qubits),
                    zz,
                    false,
                );
            }
        }
        let mut invariants_tab = SymplecticTableau::join(&continue_flow.tab, &simple_identity_tab);

        // Now to compose this with the break branch to give the invariants across the full loop
        let mut break_flow = loop_body.clone();
        break_flow.tab.post_select_1qs(&break_post_selects);
        break_flow.project_to_flow();
        break_flow.project_non_io(body_in, body_out);
        // Tensor product
        let offset = invariants_tab.add_qubits(break_flow.tab.nb_qubits);
        for i in 0..break_flow.tab.nb_stabs {
            let mut new_z = BitVector::new(offset);
            new_z.extend_vec(
                break_flow.tab.z[i].get_sized_boolean_vec(break_flow.tab.nb_qubits),
                offset,
            );
            let mut new_x = BitVector::new(offset);
            new_x.extend_vec(
                break_flow.tab.x[i].get_sized_boolean_vec(break_flow.tab.nb_qubits),
                offset,
            );
            invariants_tab.add_stab(new_z, new_x, break_flow.tab.signs.get(i));
        }
        let mut post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        // Transform Bell post-selections into post-selections on individual qubits
        for (l, r_dfp) in bells {
            let r = offset + break_flow.q_index_map.get_by_left(&r_dfp).unwrap();
            invariants_tab.append_cx(l, r);
            post_selects.push((l, PauliXZ::X, false));
            post_selects.push((r, PauliXZ::Z, false));
        }
        invariants_tab.post_select_1qs(&post_selects);
        // Propagate classical values
        for (l, r_dfp) in classical_propagations {
            let r = offset + break_flow.q_index_map.get_by_left(&r_dfp).unwrap();
            let mut zz = BitVector::new(invariants_tab.nb_qubits);
            zz.xor_bit(l);
            zz.xor_bit(r);
            invariants_tab.add_stab(zz, BitVector::new(invariants_tab.nb_qubits), false);
        }
        let mut invariants_q_map: BiHashMap<DataflowPoint<H::Node>, usize> = BiHashMap::new();
        let mut cols_to_project: Vec<(usize, PauliXZ)> = vec![];
        let mut cols_to_delete: Vec<usize> = vec![];
        for i in (offset..invariants_tab.nb_qubits).rev() {
            match break_flow.q_index_map.get_by_right(&(i - offset)).unwrap() {
                DataflowPoint::NodeIn(_, in_port, ri) => {
                    let new_dfp = if in_port.index() == 0 {
                        DataflowPoint::TempOut(
                            node,
                            OutgoingPort::from(ri[0]),
                            ri.as_slice()[1..].to_vec(),
                        )
                    } else {
                        DataflowPoint::TempOut(
                            node,
                            OutgoingPort::from(in_port.index() + loop_op.just_outputs.len() - 1),
                            ri.clone(),
                        )
                    };
                    invariants_q_map.insert(new_dfp, i);
                }
                DataflowPoint::SumInControl(_, in_port, ri, si, m) => {
                    if in_port.index() == 0 {
                        if ri.is_empty() {
                            // Controls for the top-level sum are removed
                            cols_to_project.push((i, PauliXZ::X));
                            cols_to_project.push((i, PauliXZ::Z));
                            cols_to_delete.push(i);
                        } else {
                            invariants_q_map.insert(
                                DataflowPoint::SumOutControl(
                                    node,
                                    OutgoingPort::from(ri[0]),
                                    ri.as_slice()[1..].to_vec(),
                                    *si,
                                    *m,
                                ),
                                i,
                            );
                        }
                    } else {
                        invariants_q_map.insert(
                            DataflowPoint::SumOutControl(
                                node,
                                OutgoingPort::from(
                                    in_port.index() + loop_op.just_outputs.len() - 1,
                                ),
                                ri.clone(),
                                *si,
                                *m,
                            ),
                            i,
                        );
                    }
                }
                DataflowPoint::SumInPhControl(_, in_port, ri, si) => {
                    if in_port.index() == 0 {
                        if ri.is_empty() {
                            // Controls for the top-level sum are removed
                            cols_to_project.push((i, PauliXZ::X));
                            cols_to_project.push((i, PauliXZ::Z));
                            cols_to_delete.push(i);
                        } else {
                            invariants_q_map.insert(
                                DataflowPoint::SumOutPhControl(
                                    node,
                                    OutgoingPort::from(ri[0]),
                                    ri.as_slice()[1..].to_vec(),
                                    *si,
                                ),
                                i,
                            );
                        }
                    } else {
                        invariants_q_map.insert(
                            DataflowPoint::SumOutPhControl(
                                node,
                                OutgoingPort::from(
                                    in_port.index() + loop_op.just_outputs.len() - 1,
                                ),
                                ri.clone(),
                                *si,
                            ),
                            i,
                        );
                    }
                }
                _ => {
                    cols_to_project.push((i, PauliXZ::X));
                    cols_to_project.push((i, PauliXZ::Z));
                    cols_to_delete.push(i);
                }
            }
        }
        for i in (0..offset).rev() {
            match continue_flow.q_index_map.get_by_right(&i).unwrap() {
                DataflowPoint::NodeOut(_, out_port, ri) => {
                    invariants_q_map.insert(
                        DataflowPoint::TempIn(
                            node,
                            IncomingPort::from(out_port.index()),
                            ri.clone(),
                        ),
                        i,
                    );
                }
                DataflowPoint::SumOutControl(_, out_port, ri, si, m) => {
                    invariants_q_map.insert(
                        DataflowPoint::SumInControl(
                            node,
                            IncomingPort::from(out_port.index()),
                            ri.clone(),
                            *si,
                            *m,
                        ),
                        i,
                    );
                }
                DataflowPoint::SumOutPhControl(_, out_port, ri, si) => {
                    invariants_q_map.insert(
                        DataflowPoint::SumInPhControl(
                            node,
                            IncomingPort::from(out_port.index()),
                            ri.clone(),
                            *si,
                        ),
                        i,
                    );
                }
                _ => {
                    cols_to_project.push((i, PauliXZ::X));
                    cols_to_project.push((i, PauliXZ::Z));
                    cols_to_delete.push(i);
                }
            }
        }
        invariants_tab.project_cols_to_zero(&cols_to_project);
        for q in cols_to_delete {
            if let Some(moved_qb) = invariants_tab.delete_qubit(q) {
                let (moved_dfp, _) = invariants_q_map.remove_by_right(&moved_qb).unwrap();
                invariants_q_map.insert(moved_dfp, q);
            }
        }
        let invariants = StabilizerDataflow {
            tab: invariants_tab,
            q_index_map: invariants_q_map,
        };
        match (
            loop_settings.include_external_interface,
            loop_settings.include_loop_body,
        ) {
            (false, false) => invariants,
            (false, true) => StabilizerDataflow::controlled_conditional(
                &loop_body,
                &invariants,
                node,
                &ControlSpec::Role,
                false,
            ),
            (true, false) => {
                let interface = StabilizerDataflow::generate_interface(hugr, node, loop_settings);
                StabilizerDataflow::controlled_conditional(
                    &interface,
                    &invariants,
                    node,
                    &ControlSpec::Role,
                    false,
                )
            }
            (true, true) => {
                let interface = StabilizerDataflow::generate_interface(hugr, node, loop_settings);
                let product = StabilizerDataflow::tensor_product(&interface, &loop_body);
                StabilizerDataflow::controlled_conditional(
                    &product,
                    &invariants,
                    node,
                    &ControlSpec::Role,
                    false,
                )
            }
        }
    }
}

#[cfg(test)]
mod test {
    use hugr::{
        Extension, Hugr, HugrView, IncomingPort, OutgoingPort, Wire, builder::{
            BuildHandle, Container, DFGWrapper, Dataflow, DataflowHugr, DataflowSubContainer, FunctionBuilder, HugrBuilder, SubContainer, endo_sig
        }, extension::{
            Version, prelude::{ConstError, UnwrapBuilder, bool_t, qb_t, usize_t}, simple_op::MakeExtensionOp
        }, ops::{OpType, Tag, Value, handle::{FuncID, NodeHandle}}, std_extensions::logic::LogicOp, type_row, types::{Signature, Type, TypeBase, TypeRow}
    };
    use itertools::{Itertools, chain, interleave};
    use tket::{TketOp, extension::{bool::BoolOp, rotation::{ConstRotation, rotation_type}}};

    use crate::{
        stabilizer_dataflow_complete::{
            DataflowFlowDetail, DataflowPoint, DataflowSettings, SDFAnalysis,
        },
        symplectic_tableau::PauliXZ,
    };

    fn default_settings(_: hugr::Node) -> DataflowSettings {
        DataflowSettings {
            flow_level: DataflowFlowDetail::All,
            include_external_interface: false,
            include_sum_types: true,
            include_loop_body: false,
            rotation_override: true,
        }
    }

    fn max_settings(_: hugr::Node) -> DataflowSettings {
        DataflowSettings {
            flow_level: DataflowFlowDetail::All,
            include_external_interface: true,
            include_sum_types: true,
            include_loop_body: true,
            rotation_override: true,
        }
    }

    #[test]
    fn test_print_circuit() {
        let mut builder = FunctionBuilder::new("cond", Signature::new(vec![qb_t(), qb_t(), bool_t()], vec![qb_t(), qb_t()])).unwrap();
        let [qb0, qb1, b] = builder.input_wires_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b), [(qb_t(), qb0), (qb_t(), qb1)], vec![qb_t(); 2].into()).unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let [c0q0, c0q1] = cond0_builder.input_wires_arr();
        let cx = cond0_builder.add_dataflow_op(TketOp::CX, [c0q0, c0q1]).unwrap();
        cond0_builder.finish_with_outputs([cx.out_wire(0), cx.out_wire(1)]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [c1q0, c1q1] = cond1_builder.input_wires_arr();
        let t = cond1_builder.add_dataflow_op(TketOp::T, [c1q0]).unwrap();
        let z = cond1_builder.add_dataflow_op(TketOp::Z, [c1q1]).unwrap();
        cond1_builder.finish_with_outputs([t.out_wire(0), z.out_wire(0)]).ok();
        let cond = cond_builder.finish_sub_container().unwrap();
        let hugr = builder.finish_hugr_with_outputs([cond.out_wire(0), cond.out_wire(1)]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_print_constant_propagation() {
        let mut builder = FunctionBuilder::new("bool_prop", endo_sig(vec![qb_t()])).unwrap();
        let [qb] = builder.input_wires_arr();
        let tag = builder.add_dataflow_op(Tag::new(1, vec![type_row![]; 2]), []).unwrap();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], tag.out_wire(0)), [(qb_t(), qb)], vec![qb_t()].into()).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [c0qb] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([c0qb]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [c1qb] = cond1_builder.input_wires_arr();
        let z = cond1_builder.add_dataflow_op(TketOp::Z, [c1qb]).unwrap();
        cond1_builder.finish_with_outputs([z.out_wire(0)]).ok();
        let cond = cond_builder.finish_sub_container().unwrap();
        let hugr = builder.finish_hugr_with_outputs([cond.out_wire(0)]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_multi_conditional() {
        let mut builder = FunctionBuilder::new("multi_conditional", Signature::new(vec![TypeBase::new_unit_sum(5)], vec![qb_t()])).unwrap();
        let [five] = builder.input_wires_arr();
        let mut cond_builder = builder.conditional_builder((vec![type_row![]; 5], five), [], vec![qb_t()].into()).unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let alloc0 = cond0_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        cond0_builder.finish_with_outputs([alloc0.out_wire(0)]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let alloc1 = cond1_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let x = cond1_builder.add_dataflow_op(TketOp::X, [alloc1.out_wire(0)]).unwrap();
        cond1_builder.finish_with_outputs([x.out_wire(0)]).ok();
        let mut cond2_builder = cond_builder.case_builder(2).unwrap();
        let alloc2 = cond2_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let h2 = cond2_builder.add_dataflow_op(TketOp::H, [alloc2.out_wire(0)]).unwrap();
        cond2_builder.finish_with_outputs([h2.out_wire(0)]).ok();
        let mut cond3_builder = cond_builder.case_builder(3).unwrap();
        let alloc3 = cond3_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let h3 = cond3_builder.add_dataflow_op(TketOp::H, [alloc3.out_wire(0)]).unwrap();
        let z = cond3_builder.add_dataflow_op(TketOp::Z, [h3.out_wire(0)]).unwrap();
        cond3_builder.finish_with_outputs([z.out_wire(0)]).ok();
        let mut cond4_builder = cond_builder.case_builder(4).unwrap();
        let alloc4 = cond4_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let vdg = cond4_builder.add_dataflow_op(TketOp::Vdg, [alloc4.out_wire(0)]).unwrap();
        cond4_builder.finish_with_outputs([vdg.out_wire(0)]).ok();
        let cond = cond_builder.finish_sub_container().unwrap();
        let hugr = builder.finish_hugr_with_outputs([cond.out_wire(0)]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_print_nested_sums() {
        let qb_row = TypeRow::from(qb_t());
        let qb2_row = TypeRow::from(vec![qb_t(); 2]);
        let inner_sum = Type::new_sum([type_row![], qb_row.clone()]);
        let sum_qb_row = TypeRow::from(vec![inner_sum, qb_t()]);
        let outer_sum = Type::new_sum([qb2_row.clone(), sum_qb_row.clone()]);
        let mut builder = FunctionBuilder::new("nested_conds", Signature::new(vec![outer_sum], vec![qb_t(); 2])).unwrap();
        let [s] = builder.input_wires_arr();
        let mut cond_builder = builder.conditional_builder((vec![qb2_row.clone(), sum_qb_row], s), [], qb2_row).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [qb0, qb1] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([qb0, qb1]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [s, qb1] = cond1_builder.input_wires_arr();
        let mut inner_cond_builder = cond1_builder.conditional_builder(([type_row![], qb_row.clone()], s), [], qb_row).unwrap();
        let mut inner_cond0_builder = inner_cond_builder.case_builder(0).unwrap();
        let [qb0] = inner_cond0_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        inner_cond0_builder.finish_with_outputs([qb0]).ok();
        let inner_cond1_builder = inner_cond_builder.case_builder(1).unwrap();
        let [qb0] = inner_cond1_builder.input_wires_arr();
        inner_cond1_builder.finish_with_outputs([qb0]).ok();
        let [qb0] = inner_cond_builder.finish_sub_container().unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([qb0, qb1]).ok();
        let [qb0, qb1] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_recombine_ghz() {
        let mut builder = FunctionBuilder::new("ghz_two_ways", Signature::new(vec![bool_t()], vec![qb_t(); 3])).unwrap();
        let [b] = builder.input_wires_arr();
        let [q0] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q1] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q0] = builder.add_dataflow_op(TketOp::H, [q0]).unwrap().outputs_arr();
        let [q0, q1] = builder.add_dataflow_op(TketOp::CX, [q0, q1]).unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder((vec![type_row![], type_row![]], b), [(qb_t(), q0), (qb_t(), q1)], vec![qb_t(); 3].into()).unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q0, q1] = cond0_builder.input_wires_arr();
        let [q2] = cond0_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q0, q2] = cond0_builder.add_dataflow_op(TketOp::CX, [q0, q2]).unwrap().outputs_arr();
        cond0_builder.finish_with_outputs([q0, q1, q2]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [q0, q1] = cond1_builder.input_wires_arr();
        let [q2] = cond1_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q1, q2] = cond1_builder.add_dataflow_op(TketOp::CX, [q1, q2]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([q0, q1, q2]).ok();
        let [q0, q1, q2] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([q0, q1, q2]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_loop_invariant() {
        let mut builder = FunctionBuilder::new("loop_example", endo_sig(vec![qb_t(); 2])).unwrap();
        let [qb0, qb1] = builder.input_wires_arr();
        let mut loop_builder = builder.tail_loop_builder([], [(qb_t(), qb0), (qb_t(), qb1)], type_row![]).unwrap();
        let [qb0, qb1] = loop_builder.input_wires_arr();
        let [qb0, qb1] = loop_builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
        let [qb1, qb0] = loop_builder.add_dataflow_op(TketOp::CX, [qb1, qb0]).unwrap().outputs_arr();
        let [qb1, b] = loop_builder.add_dataflow_op(TketOp::Measure, [qb1]).unwrap().outputs_arr();
        let [qb0, qb1] = loop_builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
        let [qb0, qb1] = loop_builder.finish_with_outputs(b, [qb0, qb1]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        assert_eq!(summary.tab.nb_qubits, 4);
        assert_eq!(summary.tab.nb_stabs, 4);
        assert_eq!(summary.q_index_map.len(), 4);
        let body_in = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let body_out = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(body_in, OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::NodeIn(body_out, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(body_in, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::NodeIn(body_out, IncomingPort::from(0), vec![])
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // XX commutes through
        assert_eq!(summary.tab.stab_as_string(0), "+XXXX");
        // ZI -> -ZI
        assert_eq!(summary.tab.stab_as_string(1), "-Z  Z");
        // -ZZ output stabilizer
        assert_eq!(summary.tab.stab_as_string(2), "- Z Z");
        // IZ -> ZI
        assert_eq!(summary.tab.stab_as_string(3), "+  ZZ");
    }

    #[test]
    fn test_conditional_hoisting() {
        let mut builder = FunctionBuilder::new("cond_hoist", Signature::new(vec![qb_t(), qb_t(), bool_t(), rotation_type()], vec![qb_t(), qb_t()])).unwrap();
        let [qb0, qb1, b, r] = builder.input_wires_arr();
        let [qb0] = builder.add_dataflow_op(TketOp::Rx, [qb0, r]).unwrap().outputs_arr();
        let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b), [(qb_t(), qb0), (qb_t(), qb1), (rotation_type(), r)], vec![qb_t(), qb_t()].into()).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [qb0, qb1, _ri] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([qb0, qb1]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [qb0, qb1, ri] = cond1_builder.input_wires_arr();
        let [qb1, qb0] = cond1_builder.add_dataflow_op(TketOp::CX, [qb1, qb0]).unwrap().outputs_arr();
        let [qb0] = cond1_builder.add_dataflow_op(TketOp::Rz, [qb0, ri]).unwrap().outputs_arr();
        let [qb1] = cond1_builder.add_dataflow_op(TketOp::Rx, [qb1, ri]).unwrap().outputs_arr();
        let [qb1, qb0] = cond1_builder.add_dataflow_op(TketOp::CX, [qb1, qb0]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([qb0, qb1]).ok();
        let [qb0, qb1] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
        let [qb1] = builder.add_dataflow_op(TketOp::Rz, [qb1, r]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", hugr.mermaid_string());
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_repeat_until_success() {
        let mut builder = FunctionBuilder::new("rus", endo_sig(vec![qb_t()])).unwrap();
        let [qb0] = builder.input_wires_arr();
        let mut loop_builder = builder.tail_loop_builder([], [(qb_t(), qb0)], type_row![]).unwrap();
        let [qb0] = loop_builder.input_wires_arr();
        let mut inner_loop_builder = loop_builder.tail_loop_builder([], [], vec![qb_t()].into()).unwrap();
        let [qb1] = inner_loop_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [qb2] = inner_loop_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [qb1] = inner_loop_builder.add_dataflow_op(TketOp::H, [qb1]).unwrap().outputs_arr();
        let [qb2] = inner_loop_builder.add_dataflow_op(TketOp::H, [qb2]).unwrap().outputs_arr();
        let [qb2] = inner_loop_builder.add_dataflow_op(TketOp::Tdg, [qb2]).unwrap().outputs_arr();
        let [qb1, qb2] = inner_loop_builder.add_dataflow_op(TketOp::CX, [qb1, qb2]).unwrap().outputs_arr();
        let [qb2] = inner_loop_builder.add_dataflow_op(TketOp::T, [qb2]).unwrap().outputs_arr();
        let [qb2] = inner_loop_builder.add_dataflow_op(TketOp::H, [qb2]).unwrap().outputs_arr();
        let [b2] = inner_loop_builder.add_dataflow_op(TketOp::MeasureFree, [qb2]).unwrap().outputs_arr();
        let [b2] = inner_loop_builder.add_dataflow_op(BoolOp::read, [b2]).unwrap().outputs_arr();
        let [b2] = inner_loop_builder.add_dataflow_op(LogicOp::Not, [b2]).unwrap().outputs_arr();
        let mut cond_builder = inner_loop_builder.conditional_builder(([type_row![], type_row![]], b2), [(qb_t(), qb1)], TypeBase::new_sum([type_row![], vec![qb_t()].into()]).into()).unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let [qb1] = cond0_builder.input_wires_arr();
        cond0_builder.add_dataflow_op(TketOp::QFree, [qb1]).ok();
        let [v] = cond0_builder.add_dataflow_op(Tag::new(0, [type_row![], vec![qb_t()].into()].into()), []).unwrap().outputs_arr();
        cond0_builder.finish_with_outputs([v]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [qb1] = cond1_builder.input_wires_arr();
        let [v] = cond1_builder.add_dataflow_op(Tag::new(1, [type_row![], vec![qb_t()].into()].into()), [qb1]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([v]).ok();
        let [v] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let [qb1] = inner_loop_builder.finish_with_outputs(v, []).unwrap().outputs_arr();
        let [qb0] = loop_builder.add_dataflow_op(TketOp::T, [qb0]).unwrap().outputs_arr();
        let [qb0] = loop_builder.add_dataflow_op(TketOp::Z, [qb0]).unwrap().outputs_arr();
        let [qb0, qb1] = loop_builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
        let [qb1] = loop_builder.add_dataflow_op(TketOp::T, [qb1]).unwrap().outputs_arr();
        let [qb1] = loop_builder.add_dataflow_op(TketOp::H, [qb1]).unwrap().outputs_arr();
        let [b1] = loop_builder.add_dataflow_op(TketOp::MeasureFree, [qb1]).unwrap().outputs_arr();
        let [b1] = loop_builder.add_dataflow_op(BoolOp::read, [b1]).unwrap().outputs_arr();
        let [b1] = loop_builder.add_dataflow_op(LogicOp::Not, [b1]).unwrap().outputs_arr();
        let [qb0] = loop_builder.finish_with_outputs(b1, [qb0]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0]).unwrap();
        println!("{}", hugr.mermaid_string());
        let func_node = hugr.first_child(hugr.module_root()).unwrap();

        let mut custom_settings = |n: hugr::Node| {
            let ot = hugr.get_optype(n);
            let include_interface = if let Some(extop) = ot.as_extension_op() {
                if let Ok(tkop) = TketOp::from_extension_op(extop) {
                    match tkop {
                        TketOp::QAlloc | TketOp::H | TketOp::CX | TketOp::Z => false,
                        _ => true
                    }
                }
                else { true }
            } else { !(ot.is_tail_loop() || ot.is_conditional()) };
            DataflowSettings {
                flow_level: DataflowFlowDetail::All,
                include_external_interface: include_interface,
                include_sum_types: true,
                include_loop_body: false,
                rotation_override: true,
            }
        };
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut custom_settings);
        let summary = summary_map.get(&func_node).unwrap();
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    fn add_pauli_error(builder: &mut DFGWrapper<Hugr, BuildHandle<FuncID<true>>>, qb : Wire, b: Wire, gate: TketOp) -> BuildHandle<hugr::ops::handle::ConditionalID> {
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b), [(qb_t(), qb)], vec![qb_t()].into()).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([q]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [q] = cond1_builder.input_wires_arr();
        let [q] = cond1_builder.add_dataflow_op(gate, [q]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([q]).ok();
        cond_builder.finish_sub_container().unwrap()
    }

    #[test]
    fn test_fault_tolerant() {
        // A circuit for ZZ phase between the first qubits from two iceberg blocks
        // Expected to not be fault-tolerant
        let mut builder = FunctionBuilder::new("iceberg_zz", Signature::new(chain![vec![qb_t(); 4], vec![rotation_type()], vec![bool_t(); 36]].collect_vec(), vec![qb_t(); 8])).unwrap();
        let [a0, a1, b0, b1, r, xat, zat, xa0, za0, xa1, za1, xab, zab, xbt, zbt, xb0, zb0, xb1, zb1, xbb, zbb, xab_cx, zab_cx, xbb_cx, zbb_cx, xa0_cx, za0_cx, xb0_cx, zb0_cx, xbb_zz, zbb_zz, xb0_zz, zb0_zz, xab_out, zab_out, xbb_out, zbb_out, xa0_out, za0_out, xb0_out, zb0_out] = builder.input_wires_arr();
        // Encoder (assume no noise)
        // Stabilizers ZZZZ, XXXX
        // Logicals XXII, IZIZ, XIXI, IIZZ
        let [at] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [ab] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        // IIIZ -> IIIX
        let [ab] = builder.add_dataflow_op(TketOp::H, [ab]).unwrap().outputs_arr();
        // IZII -> IZIZ
        // IIIX -> IXIX
        let [ab, a0] = builder.add_dataflow_op(TketOp::CX, [ab, a0]).unwrap().outputs_arr();
        // IIZI -> IIZZ
        // IXIX -> IXXX
        let [ab, a1] = builder.add_dataflow_op(TketOp::CX, [ab, a1]).unwrap().outputs_arr();
        // IXII -> XXII
        // ZIII -> ZZII
        let [a0, at] = builder.add_dataflow_op(TketOp::CX, [a0, at]).unwrap().outputs_arr();
        // IIXI -> XIXI
        // ZZII -> ZZZI
        let [a1, at] = builder.add_dataflow_op(TketOp::CX, [a1, at]).unwrap().outputs_arr();
        // IXXX -> XXXX
        // ZZZI -> ZZZZ
        let [ab, at] = builder.add_dataflow_op(TketOp::CX, [ab, at]).unwrap().outputs_arr();
        // Same for the b block
        let [bt] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [bb] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [bb] = builder.add_dataflow_op(TketOp::H, [bb]).unwrap().outputs_arr();
        let [bb, b0] = builder.add_dataflow_op(TketOp::CX, [bb, b0]).unwrap().outputs_arr();
        let [bb, b1] = builder.add_dataflow_op(TketOp::CX, [bb, b1]).unwrap().outputs_arr();
        let [b0, bt] = builder.add_dataflow_op(TketOp::CX, [b0, bt]).unwrap().outputs_arr();
        let [b1, bt] = builder.add_dataflow_op(TketOp::CX, [b1, bt]).unwrap().outputs_arr();
        let [bb, bt] = builder.add_dataflow_op(TketOp::CX, [bb, bt]).unwrap().outputs_arr();
        // Errors at start of circuit
        let [at] = add_pauli_error(&mut builder, at, xat, TketOp::X).outputs_arr();
        let [at] = add_pauli_error(&mut builder, at, zat, TketOp::Z).outputs_arr();
        let [a0] = add_pauli_error(&mut builder, a0, xa0, TketOp::X).outputs_arr();
        let [a0] = add_pauli_error(&mut builder, a0, za0, TketOp::Z).outputs_arr();
        let [a1] = add_pauli_error(&mut builder, a1, xa1, TketOp::X).outputs_arr();
        let [a1] = add_pauli_error(&mut builder, a1, za1, TketOp::Z).outputs_arr();
        let [ab] = add_pauli_error(&mut builder, ab, xab, TketOp::X).outputs_arr();
        let [ab] = add_pauli_error(&mut builder, ab, zab, TketOp::Z).outputs_arr();
        let [bt] = add_pauli_error(&mut builder, bt, xbt, TketOp::X).outputs_arr();
        let [bt] = add_pauli_error(&mut builder, bt, zbt, TketOp::Z).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, xb0, TketOp::X).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, zb0, TketOp::Z).outputs_arr();
        let [b1] = add_pauli_error(&mut builder, b1, xb1, TketOp::X).outputs_arr();
        let [b1] = add_pauli_error(&mut builder, b1, zb1, TketOp::Z).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, xbb, TketOp::X).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, zbb, TketOp::Z).outputs_arr();
        // CX[a.bottom, b.bottom]
        let [ab, bb] = builder.add_dataflow_op(TketOp::CX, [ab, bb]).unwrap().outputs_arr();
        // Errors from CX
        let [ab] = add_pauli_error(&mut builder, ab, xab_cx, TketOp::X).outputs_arr();
        let [ab] = add_pauli_error(&mut builder, ab, zab_cx, TketOp::Z).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, xbb_cx, TketOp::X).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, zbb_cx, TketOp::Z).outputs_arr();
        // CX[a.q[a_target], b.q[b_target]]
        let [a0, b0] = builder.add_dataflow_op(TketOp::CX, [a0, b0]).unwrap().outputs_arr();
        // Errors from CX
        let [a0] = add_pauli_error(&mut builder, a0, xa0_cx, TketOp::X).outputs_arr();
        let [a0] = add_pauli_error(&mut builder, a0, za0_cx, TketOp::Z).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, xb0_cx, TketOp::X).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, zb0_cx, TketOp::Z).outputs_arr();
        // ZZPhase[b.bottom, b.q[b_target]]
        let [bb, b0] = builder.add_dataflow_op(TketOp::CX, [bb, b0]).unwrap().outputs_arr();
        let [b0] = builder.add_dataflow_op(TketOp::Rz, [b0, r]).unwrap().outputs_arr();
        let [bb, b0] = builder.add_dataflow_op(TketOp::CX, [bb, b0]).unwrap().outputs_arr();
        // Errors from ZZPhase
        let [bb] = add_pauli_error(&mut builder, bb, xbb_zz, TketOp::X).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, zbb_zz, TketOp::Z).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, xb0_zz, TketOp::X).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, zb0_zz, TketOp::Z).outputs_arr();
        // CX[a.bottom, b.bottom]
        let [ab, bb] = builder.add_dataflow_op(TketOp::CX, [ab, bb]).unwrap().outputs_arr();
        // Errors from CX
        let [ab] = add_pauli_error(&mut builder, ab, xab_out, TketOp::X).outputs_arr();
        let [ab] = add_pauli_error(&mut builder, ab, zab_out, TketOp::Z).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, xbb_out, TketOp::X).outputs_arr();
        let [bb] = add_pauli_error(&mut builder, bb, zbb_out, TketOp::Z).outputs_arr();
        // CX[a.q[a_target], b.q[b_target]]
        let [a0, b0] = builder.add_dataflow_op(TketOp::CX, [a0, b0]).unwrap().outputs_arr();
        // Errors from CX
        let [a0] = add_pauli_error(&mut builder, a0, xa0_out, TketOp::X).outputs_arr();
        let [a0] = add_pauli_error(&mut builder, a0, za0_out, TketOp::Z).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, xb0_out, TketOp::X).outputs_arr();
        let [b0] = add_pauli_error(&mut builder, b0, zb0_out, TketOp::Z).outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([at, a0, a1, ab, bt, b0, b1, bb]).unwrap();
        println!("{}", hugr.mermaid_string());
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        
        let mut custom_settings = |n: hugr::Node| {
            let ot = hugr.get_optype(n);
            DataflowSettings {
                flow_level: DataflowFlowDetail::All,
                include_external_interface: false,
                include_sum_types: !ot.is_input(),
                include_loop_body: false,
                rotation_override: true,
            }
        };
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut custom_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        // let echelon_order: Vec<usize> = vec![11, 1, 3, 10, 9, 5, 7, 8, 32, 0, 2, 4, 6];
        let echelon_order: Vec<usize> = vec![0, 2, 4, 6, 32, 11, 1, 3, 10, 9, 5, 7, 8];
        // let echelon_order: Vec<usize> = vec![26, 27, 20, 21, 39, 40, 28, 29, 14, 15, 12, 13, 18, 19, 16, 17, 37, 38, 30, 31, 24, 25, 22, 23, 35, 36, 33, 34];
        let echelon_columns = interleave(echelon_order.iter().map(|c| (*c, PauliXZ::X)), echelon_order.iter().map(|c| (*c, PauliXZ::Z))).collect_vec();
        summary.tab.echelon(&echelon_columns);
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_five_q_code_syndromes() {
        let mut builder = FunctionBuilder::new("five_q_code", Signature::new(chain!(vec![qb_t()], vec![bool_t(); 10]).collect_vec(), vec![qb_t(); 5])).unwrap();
        let [q0, x0, z0, x1, z1, x2, z2, x3, z3, x4, z4] = builder.input_wires_arr();
        // Encoder (not the most efficient, don't care, it just needs to prepare the right state)
        // Logicals XXXXX, ZZZZZ
        // Stabilizers XZZXI, IXZZX, XIXZZ, ZXIXZ
        let [q1] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q2] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q3] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q4] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q1] = builder.add_dataflow_op(TketOp::H, [q1]).unwrap().outputs_arr();
        let [q2] = builder.add_dataflow_op(TketOp::H, [q2]).unwrap().outputs_arr();
        let [q3] = builder.add_dataflow_op(TketOp::H, [q3]).unwrap().outputs_arr();
        let [q4] = builder.add_dataflow_op(TketOp::H, [q4]).unwrap().outputs_arr();
        let [q1, q2] = builder.add_dataflow_op(TketOp::CZ, [q1, q2]).unwrap().outputs_arr();
        let [q1, q3] = builder.add_dataflow_op(TketOp::CZ, [q1, q3]).unwrap().outputs_arr();
        let [q1, q4] = builder.add_dataflow_op(TketOp::CZ, [q1, q4]).unwrap().outputs_arr();
        let [q2, q3] = builder.add_dataflow_op(TketOp::CZ, [q2, q3]).unwrap().outputs_arr();
        let [q2, q4] = builder.add_dataflow_op(TketOp::CZ, [q2, q4]).unwrap().outputs_arr();
        let [q3, q4] = builder.add_dataflow_op(TketOp::CZ, [q3, q4]).unwrap().outputs_arr();
        let [q3, q1] = builder.add_dataflow_op(TketOp::CX, [q3, q1]).unwrap().outputs_arr();
        let [q1, q4] = builder.add_dataflow_op(TketOp::CX, [q1, q4]).unwrap().outputs_arr();
        let [q4, q2] = builder.add_dataflow_op(TketOp::CX, [q4, q2]).unwrap().outputs_arr();
        let [q0, q1] = builder.add_dataflow_op(TketOp::CX, [q0, q1]).unwrap().outputs_arr();
        let [q0, q3] = builder.add_dataflow_op(TketOp::CX, [q0, q3]).unwrap().outputs_arr();
        let [q2, q0] = builder.add_dataflow_op(TketOp::CX, [q2, q0]).unwrap().outputs_arr();
        let [q4, q0] = builder.add_dataflow_op(TketOp::CX, [q4, q0]).unwrap().outputs_arr();
        let [q1, q2] = builder.add_dataflow_op(TketOp::CX, [q1, q2]).unwrap().outputs_arr();
        let [q3, q4] = builder.add_dataflow_op(TketOp::CX, [q3, q4]).unwrap().outputs_arr();
        // Error channels
        let [q0] = add_pauli_error(&mut builder, q0, x0, TketOp::X).outputs_arr();
        let [q0] = add_pauli_error(&mut builder, q0, z0, TketOp::Z).outputs_arr();
        let [q1] = add_pauli_error(&mut builder, q1, x1, TketOp::X).outputs_arr();
        let [q1] = add_pauli_error(&mut builder, q1, z1, TketOp::Z).outputs_arr();
        let [q2] = add_pauli_error(&mut builder, q2, x2, TketOp::X).outputs_arr();
        let [q2] = add_pauli_error(&mut builder, q2, z2, TketOp::Z).outputs_arr();
        let [q3] = add_pauli_error(&mut builder, q3, x3, TketOp::X).outputs_arr();
        let [q3] = add_pauli_error(&mut builder, q3, z3, TketOp::Z).outputs_arr();
        let [q4] = add_pauli_error(&mut builder, q4, x4, TketOp::X).outputs_arr();
        let [q4] = add_pauli_error(&mut builder, q4, z4, TketOp::Z).outputs_arr();
        // Syndrome measurements
        let [s0] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [s1] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [s2] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [s3] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [s0] = builder.add_dataflow_op(TketOp::H, [s0]).unwrap().outputs_arr();
        let [s1] = builder.add_dataflow_op(TketOp::H, [s1]).unwrap().outputs_arr();
        let [s2] = builder.add_dataflow_op(TketOp::H, [s2]).unwrap().outputs_arr();
        let [s3] = builder.add_dataflow_op(TketOp::H, [s3]).unwrap().outputs_arr();
        let [s0, q0] = builder.add_dataflow_op(TketOp::CX, [s0, q0]).unwrap().outputs_arr();
        let [s0, q1] = builder.add_dataflow_op(TketOp::CZ, [s0, q1]).unwrap().outputs_arr();
        let [s0, q2] = builder.add_dataflow_op(TketOp::CZ, [s0, q2]).unwrap().outputs_arr();
        let [s0, q3] = builder.add_dataflow_op(TketOp::CX, [s0, q3]).unwrap().outputs_arr();
        let [s1, q1] = builder.add_dataflow_op(TketOp::CX, [s1, q1]).unwrap().outputs_arr();
        let [s1, q2] = builder.add_dataflow_op(TketOp::CZ, [s1, q2]).unwrap().outputs_arr();
        let [s1, q3] = builder.add_dataflow_op(TketOp::CZ, [s1, q3]).unwrap().outputs_arr();
        let [s1, q4] = builder.add_dataflow_op(TketOp::CX, [s1, q4]).unwrap().outputs_arr();
        let [s2, q2] = builder.add_dataflow_op(TketOp::CX, [s2, q2]).unwrap().outputs_arr();
        let [s2, q3] = builder.add_dataflow_op(TketOp::CZ, [s2, q3]).unwrap().outputs_arr();
        let [s2, q4] = builder.add_dataflow_op(TketOp::CZ, [s2, q4]).unwrap().outputs_arr();
        let [s2, q0] = builder.add_dataflow_op(TketOp::CX, [s2, q0]).unwrap().outputs_arr();
        let [s3, q3] = builder.add_dataflow_op(TketOp::CX, [s3, q3]).unwrap().outputs_arr();
        let [s3, q4] = builder.add_dataflow_op(TketOp::CZ, [s3, q4]).unwrap().outputs_arr();
        let [s3, q0] = builder.add_dataflow_op(TketOp::CZ, [s3, q0]).unwrap().outputs_arr();
        let [s3, q1] = builder.add_dataflow_op(TketOp::CX, [s3, q1]).unwrap().outputs_arr();
        let [s0] = builder.add_dataflow_op(TketOp::H, [s0]).unwrap().outputs_arr();
        let [s1] = builder.add_dataflow_op(TketOp::H, [s1]).unwrap().outputs_arr();
        let [s2] = builder.add_dataflow_op(TketOp::H, [s2]).unwrap().outputs_arr();
        let [s3] = builder.add_dataflow_op(TketOp::H, [s3]).unwrap().outputs_arr();
        let [_m0] = builder.add_dataflow_op(TketOp::MeasureFree, [s0]).unwrap().outputs_arr();
        let [_m1] = builder.add_dataflow_op(TketOp::MeasureFree, [s1]).unwrap().outputs_arr();
        let [_m2] = builder.add_dataflow_op(TketOp::MeasureFree, [s2]).unwrap().outputs_arr();
        let [_m3] = builder.add_dataflow_op(TketOp::MeasureFree, [s3]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([q0, q1, q2, q3, q4]).unwrap();
        println!("{}", hugr.mermaid_string());
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        
        let mut custom_settings = |n: hugr::Node| {
            let ot = hugr.get_optype(n);
            DataflowSettings {
                flow_level: DataflowFlowDetail::All,
                include_external_interface: false,
                include_sum_types: !ot.is_input(),
                include_loop_body: false,
                rotation_override: true,
            }
        };
        let summary_map = SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut custom_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        // let echelon_order: Vec<usize> = vec![11, 1, 3, 10, 9, 5, 7, 8, 32, 0, 2, 4, 6];
        let echelon_order: Vec<usize> = vec![0, 1, 9, 8, 7, 6, 5, 26, 27, 4, 24, 25, 3, 22, 23, 2, 20, 21];
        // let echelon_order: Vec<usize> = vec![26, 27, 20, 21, 39, 40, 28, 29, 14, 15, 12, 13, 18, 19, 16, 17, 37, 38, 30, 31, 24, 25, 22, 23, 35, 36, 33, 34];
        let echelon_columns = interleave(echelon_order.iter().map(|c| (*c, PauliXZ::X)), echelon_order.iter().map(|c| (*c, PauliXZ::Z))).collect_vec();
        summary.tab.echelon(&echelon_columns);
        println!("{}", summary.tab);
        for i in 0..summary.tab.nb_qubits {
            println!("{}\t{:?}", i, summary.q_index_map.get_by_right(&i).unwrap());
        }
        // assert!(false);
    }

    #[test]
    fn test_print_preliminary_example() {
        let mut builder = FunctionBuilder::new("my_func", Signature::new(vec![qb_t(), qb_t()], vec![qb_t(), qb_t()])).unwrap();
        let [q0, q1] = builder.input_wires_arr();
        let [anc] = builder.add_dataflow_op(TketOp::TryQAlloc, []).unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], vec![qb_t()].into()], anc), [(qb_t(), q0)], vec![qb_t(), bool_t()].into()).unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q0] = cond0_builder.input_wires_arr();
        let [q0] = cond0_builder.add_dataflow_op(TketOp::S, [q0]).unwrap().outputs_arr();
        let [tag] = cond0_builder.add_dataflow_op(Tag::new(0, vec![type_row![]; 2]), []).unwrap().outputs_arr();
        cond0_builder.finish_with_outputs([q0, tag]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [anc, q0] = cond1_builder.input_wires_arr();
        let [anc] = cond1_builder.add_dataflow_op(TketOp::H, [anc]).unwrap().outputs_arr();
        let [anc] = cond1_builder.add_dataflow_op(TketOp::S, [anc]).unwrap().outputs_arr();
        let [anc, q0] = cond1_builder.add_dataflow_op(TketOp::CX, [anc, q0]).unwrap().outputs_arr();
        let [q0, res] = cond1_builder.add_dataflow_op(TketOp::Measure, [q0]).unwrap().outputs_arr();
        cond1_builder.add_dataflow_op(TketOp::QFree, [q0]).ok();
        cond1_builder.finish_with_outputs([anc, res]).ok();
        let [q0, correct] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let [q0, q1] = builder.add_dataflow_op(TketOp::CX, [q0, q1]).unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], correct), [(qb_t(), q0)], vec![qb_t()].into()).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q0] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([q0]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [q0] = cond1_builder.input_wires_arr();
        let [q0] = cond1_builder.add_dataflow_op(TketOp::Y, [q0]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([q0]).ok();
        let [q0] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], correct), [(qb_t(), q1)], vec![qb_t()].into()).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q1] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([q1]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [q1] = cond1_builder.input_wires_arr();
        let [q1] = cond1_builder.add_dataflow_op(TketOp::X, [q1]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([q1]).ok();
        let [q1] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([q0, q1]).unwrap();
        println!("{}", hugr.mermaid_string());
        // assert!(false);
    }

    #[test]
    fn test_print_assertion() {
        let mut builder = FunctionBuilder::new("assertion", endo_sig(vec![qb_t(); 2])).unwrap();
        let [q0, q1] = builder.input_wires_arr();
        let [anc] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [anc] = builder.add_dataflow_op(TketOp::H, [anc]).unwrap().outputs_arr();
        let [anc, q0] = builder.add_dataflow_op(TketOp::CZ, [anc, q0]).unwrap().outputs_arr();
        let [anc, q1] = builder.add_dataflow_op(TketOp::CY, [anc, q1]).unwrap().outputs_arr();
        let [anc] = builder.add_dataflow_op(TketOp::H, [anc]).unwrap().outputs_arr();
        let [anc, b] = builder.add_dataflow_op(TketOp::Measure, [anc]).unwrap().outputs_arr();
        builder.add_dataflow_op(TketOp::QFree, [anc]).ok();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b), [], type_row![]).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        cond0_builder.finish_with_outputs([]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let err = ConstError::new(0, "Error detected");
        cond1_builder.add_panic(err, [], []).ok();
        cond1_builder.finish_with_outputs([]).ok();
        cond_builder.finish_sub_container().ok();
        let hugr = builder.finish_hugr_with_outputs([q0, q1]).unwrap();
        println!("{}", hugr.mermaid_string());
        // assert!(false);
    }

    #[test]
    fn test_loop_internal() {
        let mut builder = FunctionBuilder::new("loop", Signature::new(vec![qb_t(); 2], vec![qb_t()])).unwrap();
        let [q0, q1] = builder.input_wires_arr();
        let mut tl_builder = builder.tail_loop_builder([(qb_t(), q0)], [(qb_t(), q1)], type_row![]).unwrap();
        let [q0, q1] = tl_builder.input_wires_arr();
        let [q1, q0] = tl_builder.add_dataflow_op(TketOp::CX, [q1, q0]).unwrap().outputs_arr();
        let [q0] = tl_builder.add_dataflow_op(TketOp::T, [q0]).unwrap().outputs_arr();
        let [q1] = tl_builder.add_dataflow_op(TketOp::Tdg, [q1]).unwrap().outputs_arr();
        let [q0] = tl_builder.add_dataflow_op(TketOp::H, [q0]).unwrap().outputs_arr();
        let [q0, b] = tl_builder.add_dataflow_op(TketOp::Measure, [q0]).unwrap().outputs_arr();
        tl_builder.add_dataflow_op(TketOp::QFree, [q0]).ok();
        let mut cond_builder = tl_builder.conditional_builder(([type_row![], type_row![]], b), [(qb_t(), q1)], vec![Type::new_sum([vec![qb_t()].into(), type_row![]]), qb_t()].into()).unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q1] = cond0_builder.input_wires_arr();
        let [q0] = cond0_builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q0] = cond0_builder.add_dataflow_op(Tag::new(0, vec![vec![qb_t()].into(), type_row![]]), [q0]).unwrap().outputs_arr();
        cond0_builder.finish_with_outputs([q0, q1]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [q1] = cond1_builder.input_wires_arr();
        let [q1] = cond1_builder.add_dataflow_op(TketOp::H, [q1]).unwrap().outputs_arr();
        let [q1] = cond1_builder.add_dataflow_op(TketOp::T, [q1]).unwrap().outputs_arr();
        let [q0] = cond1_builder.add_dataflow_op(Tag::new(1, vec![vec![qb_t()].into(), type_row![]]), []).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([q0, q1]).ok();
        let [q0, q1] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let [q1] = tl_builder.finish_with_outputs(q0, [q1]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([q1]).unwrap();
        println!("{}", hugr.mermaid_string());
        // assert!(false);
    }

    #[test]
    fn test_back_and_forth() {
        let mut builder = FunctionBuilder::new("bidirectional", Signature::new(vec![], vec![qb_t()])).unwrap();
        let [q0] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q1] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q2] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
        let [q0] = builder.add_dataflow_op(TketOp::X, [q0]).unwrap().outputs_arr();
        let [q1] = builder.add_dataflow_op(TketOp::H, [q1]).unwrap().outputs_arr();
        let [q1, q2] = builder.add_dataflow_op(TketOp::CX, [q1, q2]).unwrap().outputs_arr();
        let [q0, b0] = builder.add_dataflow_op(TketOp::Measure, [q0]).unwrap().outputs_arr();
        let [] = builder.add_dataflow_op(TketOp::QFree, [q0]).unwrap().outputs_arr();
        let [q1, b1] = builder.add_dataflow_op(TketOp::Measure, [q1]).unwrap().outputs_arr();
        let [] = builder.add_dataflow_op(TketOp::QFree, [q1]).unwrap().outputs_arr();
        let [b2] = builder.add_dataflow_op(LogicOp::And, [b0, b1]).unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b2), [], type_row![]).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        cond0_builder.finish_with_outputs([]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let err = ConstError::new(0, "Error detected");
        cond1_builder.add_panic(err, [], []).ok();
        cond1_builder.finish_with_outputs([]).ok();
        cond_builder.finish_sub_container().ok();
        let hugr = builder.finish_hugr_with_outputs([q2]).unwrap();
        println!("{}", hugr.mermaid_string());
        // assert!(false);
    }

    #[test]
    fn test_dynamic_merge() {
        let mut builder = FunctionBuilder::new("dynamic_merge", Signature::new(vec![qb_t(); 2], vec![bool_t(); 2])).unwrap();
        let [q0, q1] = builder.input_wires_arr();
        let [q0] = builder.add_dataflow_op(TketOp::T, [q0]).unwrap().outputs_arr();
        let [q1] = builder.add_dataflow_op(TketOp::T, [q1]).unwrap().outputs_arr();
        let [q0, q1] = builder.add_dataflow_op(TketOp::CX, [q0, q1]).unwrap().outputs_arr();
        let [q0] = builder.add_dataflow_op(TketOp::H, [q0]).unwrap().outputs_arr();
        let [b0] = builder.add_dataflow_op(TketOp::MeasureFree, [q0]).unwrap().outputs_arr();
        let [b0] = builder.add_dataflow_op(BoolOp::read, [b0]).unwrap().outputs_arr();
        let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b0), [(qb_t(), q1)], vec![qb_t()].into()).unwrap();
        let cond0_builder = cond_builder.case_builder(0).unwrap();
        let [q1] = cond0_builder.input_wires_arr();
        cond0_builder.finish_with_outputs([q1]).ok();
        let mut cond1_builder = cond_builder.case_builder(1).unwrap();
        let [q1] = cond1_builder.input_wires_arr();
        let [q1] = cond1_builder.add_dataflow_op(TketOp::X, [q1]).unwrap().outputs_arr();
        cond1_builder.finish_with_outputs([q1]).ok();
        let [q1] = cond_builder.finish_sub_container().unwrap().outputs_arr();
        let [b1] = builder.add_dataflow_op(TketOp::MeasureFree, [q1]).unwrap().outputs_arr();
        let [b1] = builder.add_dataflow_op(BoolOp::read, [b1]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([b0, b1]).unwrap();
        println!("{}", hugr.mermaid_string());
        // assert!(false);
    }

    #[test]
    fn test_dynamic_merge2() {
        let mut builder = FunctionBuilder::new("dynamic_merge", Signature::new(vec![qb_t(); 2], vec![bool_t(); 2])).unwrap();
        let [q0, q1] = builder.input_wires_arr();
        let [q1] = builder.add_dataflow_op(TketOp::T, [q1]).unwrap().outputs_arr();
        let [q0, q1] = builder.add_dataflow_op(TketOp::CX, [q0, q1]).unwrap().outputs_arr();
        let [q0] = builder.add_dataflow_op(TketOp::T, [q0]).unwrap().outputs_arr();
        let [q0] = builder.add_dataflow_op(TketOp::H, [q0]).unwrap().outputs_arr();
        let [b1] = builder.add_dataflow_op(TketOp::MeasureFree, [q1]).unwrap().outputs_arr();
        let [b1] = builder.add_dataflow_op(BoolOp::read, [b1]).unwrap().outputs_arr();
        let [b0] = builder.add_dataflow_op(TketOp::MeasureFree, [q0]).unwrap().outputs_arr();
        let [b0] = builder.add_dataflow_op(BoolOp::read, [b0]).unwrap().outputs_arr();
        let [b1] = builder.add_dataflow_op(LogicOp::Xor, [b0, b1]).unwrap().outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([b0, b1]).unwrap();
        println!("{}", hugr.mermaid_string());
        assert!(false);
    }

    #[test]
    fn test_empty_analysis() {
        let builder = FunctionBuilder::new("empty", endo_sig(vec![])).unwrap();
        let hugr = builder.finish_hugr().unwrap();
        let empty_map = SDFAnalysis::summarise_targets(&hugr, &vec![], &mut default_settings);
        assert!(empty_map.is_empty());
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        assert_eq!(summary_map.len(), 1);
        let summary = summary_map.get(&func_node).unwrap();
        assert_eq!(summary.tab.nb_qubits, 0);
        assert_eq!(summary.tab.nb_stabs, 0);
    }

    #[test]
    fn test_identity_analysis() {
        // Add an extra integer input to make sure we only track the qubits and sum types
        let builder = FunctionBuilder::new(
            "identity",
            endo_sig(vec![usize_t(), qb_t(), qb_t(), bool_t()]),
        )
        .unwrap();
        let [i, qb0, qb1, b] = builder.input_wires_arr();
        let hugr = builder.finish_hugr_with_outputs([i, qb0, qb1, b]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        assert_eq!(summary.tab.nb_qubits, 6);
        assert_eq!(summary.tab.nb_stabs, 5);
        // Check the right ports are stored for tracking the qubits
        assert_eq!(summary.q_index_map.len(), 6);
        let body_in = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let body_out = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(body_in, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::NodeIn(body_out, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(body_in, OutgoingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::NodeIn(body_out, IncomingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::SumOutPhControl(body_in, OutgoingPort::from(3), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::SumInPhControl(body_out, IncomingPort::from(3), vec![], 0)
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Check that the rows correspond to the identity operations
        // Note that BitVector assigns index 0 to the least significant bit and always adds an extra block than needed
        assert_eq!(summary.tab.stab_as_string(0), "+XX    ");
        assert_eq!(summary.tab.stab_as_string(1), "+ZZ    ");
        assert_eq!(summary.tab.stab_as_string(2), "+  XX  ");
        assert_eq!(summary.tab.stab_as_string(3), "+  ZZ  ");
        assert_eq!(summary.tab.stab_as_string(4), "+    ZZ");
    }

    #[test]
    fn test_bell_state() {
        let mut builder =
            FunctionBuilder::new("bell", Signature::new(vec![], vec![qb_t(), qb_t()])).unwrap();
        let [qb0] = builder
            .add_dataflow_op(TketOp::QAlloc, [])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::QAlloc, [])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::H, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        assert_eq!(summary.tab.nb_qubits, 2);
        assert_eq!(summary.tab.nb_stabs, 2);
        // Check that the rows correspond to the Bell state stabilizers
        summary.tab.echelon(&summary.tab.all_columns());
        assert_eq!(summary.tab.stab_as_string(0), "+XX");
        assert_eq!(summary.tab.stab_as_string(1), "+ZZ");
    }

    #[test]
    fn test_opaque() {
        let ext = Extension::new_arc(
            "ext".try_into().unwrap(),
            Version::new(0, 0, 0),
            |ext, extension_ref| {
                ext.add_op(
                    "op".into(),
                    String::new(),
                    Signature::new_endo(vec![qb_t()]),
                    extension_ref,
                )
                .unwrap();
            },
        );
        let mut builder =
            FunctionBuilder::new("opaque_test", Signature::new(vec![], vec![qb_t(), qb_t()]))
                .unwrap();
        let [qb0] = builder
            .add_dataflow_op(TketOp::QAlloc, [])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::QAlloc, [])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::H, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let op = ext.instantiate_extension_op("op", []).unwrap();
        let opaque_op = builder.add_dataflow_op(op, [qb1]).unwrap();
        let [qb1] = opaque_op.outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::H, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let mut detail_settings = |node| {
            if node == opaque_op.node() {
                // Ask for the external interface of the opaque node to be tracked and no flow
                DataflowSettings {
                    flow_level: DataflowFlowDetail::None,
                    include_external_interface: true,
                    include_sum_types: true,
                    include_loop_body: true,
                    rotation_override: true,
                }
            } else {
                default_settings(node)
            }
        };
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut detail_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        assert_eq!(summary.tab.nb_qubits, 4);
        assert_eq!(summary.tab.nb_stabs, 4);
        // Reduce summary.tab to row echelon form with qubit ordering [op_in, out0, op_out, out1] (because the second QAlloc occurs first in the topological sort)
        let out_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeIn(opaque_op.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(opaque_op.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(1), vec![])
        );
        // Check the rows
        summary.tab.echelon(&summary.tab.all_columns());
        assert_eq!(summary.tab.stab_as_string(0), "+XX X");
        assert_eq!(summary.tab.stab_as_string(1), "+Z XZ");
        assert_eq!(summary.tab.stab_as_string(2), "+ ZXZ");
        assert_eq!(summary.tab.stab_as_string(3), "+  ZX");
    }

    #[test]
    fn test_clifford_gates() {
        // Need to cover H, CX, CY, CZ, S, Sdg, X, Y, Z, V, Vdg
        let mut builder =
            FunctionBuilder::new("cliff", endo_sig(vec![qb_t(), qb_t(), qb_t()])).unwrap();
        let [qb0, qb1, qb2] = builder.input_wires_arr();
        // H;S;V;S = I
        let [qb0] = builder
            .add_dataflow_op(TketOp::H, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::S, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::V, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::S, [qb0])
            .unwrap()
            .outputs_arr();
        // CX;IH;CZ;IH = II
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::H, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CZ, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::H, [qb1])
            .unwrap()
            .outputs_arr();
        // CX;IS;CY;ISdg = II
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::S, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CY, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::Sdg, [qb1])
            .unwrap()
            .outputs_arr();
        // S;Sdg = II
        let [qb0] = builder
            .add_dataflow_op(TketOp::S, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::Sdg, [qb0])
            .unwrap()
            .outputs_arr();
        // V;Vdh = II
        let [qb1] = builder
            .add_dataflow_op(TketOp::V, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::Vdg, [qb1])
            .unwrap()
            .outputs_arr();
        // Test Paulis explicitly
        let [qb0] = builder
            .add_dataflow_op(TketOp::Z, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::X, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb2] = builder
            .add_dataflow_op(TketOp::Y, [qb2])
            .unwrap()
            .outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1, qb2]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        assert_eq!(summary.tab.nb_qubits, 6);
        assert_eq!(summary.tab.nb_stabs, 6);
        // Reduce summary.tab to row echelon form with qubit ordering [in0, out0, in1, out1, in2, out2]
        let in_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let out_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(2), vec![])
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Check the rows
        assert_eq!(summary.tab.stab_as_string(0), "-XX    ");
        assert_eq!(summary.tab.stab_as_string(1), "+ZZ    ");
        assert_eq!(summary.tab.stab_as_string(2), "+  XX  ");
        assert_eq!(summary.tab.stab_as_string(3), "-  ZZ  ");
        assert_eq!(summary.tab.stab_as_string(4), "-    XX");
        assert_eq!(summary.tab.stab_as_string(5), "-    ZZ");
    }

    #[test]
    fn test_nonclifford() {
        // Need to cover the separate logic for CRz, T/Tdg/Rz/Measure, Rx, Ry, Toffoli
        let mut builder =
            FunctionBuilder::new("non_cliff", endo_sig(vec![qb_t(), qb_t(), qb_t()])).unwrap();
        let [qb0, qb1, qb2] = builder.input_wires_arr();
        let constant = builder.add_constant(Value::extension(ConstRotation::new(0.137).unwrap()));
        let loaded_const = builder.load_const(&constant);
        let t = builder.add_dataflow_op(TketOp::T, [qb0]).unwrap();
        let tdg = builder
            .add_dataflow_op(TketOp::Tdg, [t.out_wire(0)])
            .unwrap();
        let rz = builder
            .add_dataflow_op(TketOp::Rz, [tdg.out_wire(0), loaded_const])
            .unwrap();
        let meas = builder
            .add_dataflow_op(TketOp::Measure, [rz.out_wire(0)])
            .unwrap();
        let ry = builder
            .add_dataflow_op(TketOp::Ry, [qb1, loaded_const])
            .unwrap();
        let rx = builder
            .add_dataflow_op(TketOp::Rx, [qb2, loaded_const])
            .unwrap();
        let crz = builder
            .add_dataflow_op(
                TketOp::CRz,
                [meas.out_wire(0), ry.out_wire(0), loaded_const],
            )
            .unwrap();
        let toffoli = builder
            .add_dataflow_op(
                TketOp::Toffoli,
                [crz.out_wire(0), crz.out_wire(1), rx.out_wire(0)],
            )
            .unwrap();
        let hugr = builder
            .finish_hugr_with_outputs(toffoli.outputs_arr::<3>())
            .unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut max_settings);
        let summary = summary_map.get(&func_node).unwrap().clone();
        // 3 inputs, 3 outputs, 6 basic rotations, 1 classical value + 1 role for measure, 4 interface + 8 role for crz, 6 interface + 12 role for toffoli
        // 3+3+6+2+4+8+6+12=44
        assert_eq!(summary.tab.nb_qubits, 44);
        // One stabilizer per qubit in the choi state (3+3+6+4+6=22), plus one per flow stabilizer over measure, crz and toffoli (1+2+3=6)
        assert_eq!(summary.tab.nb_stabs, 28);
        let in_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let out_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        // The weird ordering of the qubits here can be explained by (a) topological sort of the hugr, and (b) some qubits being reordered when we append the crz or toffoli by Bell post-selection which removes the intermediate qubits
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::Rotation(rx.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::Rotation(ry.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::Rotation(t.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::Rotation(tdg.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::Rotation(rz.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::Rotation(meas.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::SumOutPhControl(meas.node(), OutgoingPort::from(1), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::RoleControl(meas.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::NodeIn(crz.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&16).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&17).unwrap(),
            DataflowPoint::NodeIn(crz.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&18).unwrap(),
            DataflowPoint::NodeOut(crz.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&19).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&20).unwrap(),
            DataflowPoint::NodeOut(crz.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&21).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 8)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&22).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&23).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&24).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&25).unwrap(),
            DataflowPoint::RoleControl(crz.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&26).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 9)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&27).unwrap(),
            DataflowPoint::NodeIn(toffoli.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&28).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 10)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&29).unwrap(),
            DataflowPoint::NodeIn(toffoli.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&30).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 11)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&31).unwrap(),
            DataflowPoint::NodeIn(toffoli.node(), IncomingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&32).unwrap(),
            DataflowPoint::NodeOut(toffoli.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&33).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&34).unwrap(),
            DataflowPoint::NodeOut(toffoli.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&35).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&36).unwrap(),
            DataflowPoint::NodeOut(toffoli.node(), OutgoingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&37).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&38).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&39).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&40).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&41).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&42).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&43).unwrap(),
            DataflowPoint::RoleControl(toffoli.node(), 5)
        );

        // To make checking this against our intentions easier, we will post-select onto flow information and onto interface information
        let mut flow_tab = summary.tab.clone();
        // List the measurement rc, followed by crz rcs, then toffoli rcs
        let rcs = vec![
            13, 22, 23, 24, 25, 1, 3, 14, 16, 38, 39, 40, 41, 42, 43, 5, 19, 21, 26, 28, 30,
        ];
        let flow_post_selects: Vec<(usize, PauliXZ, bool)> =
            rcs.iter().map(|q| (*q, PauliXZ::X, false)).collect_vec();
        flow_tab.post_select_1qs(&flow_post_selects);
        // Pick the echelon order to make interpreting the stabilizers easier
        // outs, rx, ry, meas, meas.out, rz, tdg, t, ins
        let flow_order = vec![33, 35, 37, 6, 7, 11, 12, 10, 9, 8, 0, 2, 4]
            .iter()
            .map(|q| [(*q, PauliXZ::Z), (*q, PauliXZ::X)])
            .flatten()
            .collect_vec();
        flow_tab.echelon(&flow_order);
        // The remaining choi state for flow analysis has 3 ins, 3 outs, one classical value, and 6 basic rotations, so a complete set of stabilizers would require 13, but we lose 3 from the Toffoli (the 2 lost from the CRz coincide with some already lost from the Toffoli)
        assert_eq!(flow_tab.nb_stabs, 10);
        // Zin0 Zout0
        assert_eq!(
            flow_tab.stab_as_string(0),
            "+Z                                Z          "
        );
        // Zin1 Xry Zout1
        assert_eq!(
            flow_tab.stab_as_string(1),
            "+  Z    X                           Z        "
        );
        // Xin2 Xout2
        assert_eq!(
            flow_tab.stab_as_string(2),
            "+    X                                X      "
        );
        // Xin2 Zrx
        assert_eq!(
            flow_tab.stab_as_string(3),
            "+    X Z                                     "
        );
        // -Yin1 Zry (the negation is because of the transposition for inputs)
        assert_eq!(
            flow_tab.stab_as_string(4),
            "-  Y    Z                                    "
        );
        // Zin0 Zmeas
        assert_eq!(
            flow_tab.stab_as_string(5),
            "+Z          Z                                "
        );
        // Zin0 Zmeas.out
        assert_eq!(
            flow_tab.stab_as_string(6),
            "+Z           Z                               "
        );
        // Zin0 Zrz
        assert_eq!(
            flow_tab.stab_as_string(7),
            "+Z         Z                                 "
        );
        // Zin0 Ztdg
        assert_eq!(
            flow_tab.stab_as_string(8),
            "+Z        Z                                  "
        );
        // Zin0 Zt
        assert_eq!(
            flow_tab.stab_as_string(9),
            "+Z       Z                                   "
        );

        let mut interface_tab = summary.tab.clone();
        let interface_post_selects: Vec<(usize, PauliXZ, bool)> =
            rcs.iter().map(|q| (*q, PauliXZ::Z, false)).collect_vec();
        interface_tab.post_select_1qs(&interface_post_selects);
        // Pick the echelon order to make interpreting the stabilizers easier
        // toffoli.outs, toffoli.ins, crz.outs, crz.ins, [flow_order]
        let interface_order = vec![
            32, 34, 36, 27, 29, 31, 18, 20, 15, 17, 33, 35, 37, 6, 7, 11, 12, 10, 9, 8, 0, 2, 4,
        ]
        .iter()
        .map(|q| [(*q, PauliXZ::Z), (*q, PauliXZ::X)])
        .flatten()
        .collect_vec();
        interface_tab.echelon(&interface_order);
        // The choi state for the interface has 22 qubits and we have a complete set of stabilizers
        assert_eq!(interface_tab.nb_stabs, 22);
        // Ztof.out0 Zout0
        assert_eq!(
            interface_tab.stab_as_string(0),
            "+                                ZZ          "
        );
        // Xtof.out0 Xout0
        assert_eq!(
            interface_tab.stab_as_string(1),
            "+                                XX          "
        );
        // Ztof.out1 Zout1
        assert_eq!(
            interface_tab.stab_as_string(2),
            "+                                  ZZ        "
        );
        // Xtof.out1 Xout1
        assert_eq!(
            interface_tab.stab_as_string(3),
            "+                                  XX        "
        );
        // Ztof.out2 Zout2
        assert_eq!(
            interface_tab.stab_as_string(4),
            "+                                    ZZ      "
        );
        // Xtof.out2 Xout2
        assert_eq!(
            interface_tab.stab_as_string(5),
            "+                                    XX      "
        );
        // Zcrz.out0 Ztof.in0
        assert_eq!(
            interface_tab.stab_as_string(6),
            "+                  Z        Z                "
        );
        // Xcrz.out0 Xtof.in0
        assert_eq!(
            interface_tab.stab_as_string(7),
            "+                  X        X                "
        );
        // Zcrz.out1 Ztof.in1
        assert_eq!(
            interface_tab.stab_as_string(8),
            "+                    Z        Z              "
        );
        // Xcrz.out1 Xtof.in1
        assert_eq!(
            interface_tab.stab_as_string(9),
            "+                    X        X              "
        );
        // Zin2 Xrx Ztof.in2
        assert_eq!(
            interface_tab.stab_as_string(10),
            "+    Z X                        Z            "
        );
        // Xin2 Xtof.in2
        assert_eq!(
            interface_tab.stab_as_string(11),
            "+    X                          X            "
        );
        // Zin0 Zcrz.in0
        assert_eq!(
            interface_tab.stab_as_string(12),
            "+Z              Z                            "
        );
        // Xin0 Xt Xtdg Xrz Xmeas Xcrz.in0
        assert_eq!(
            interface_tab.stab_as_string(13),
            "+X       XXXX   X                            "
        );
        // Zin1 Xry Zcrz.in1
        assert_eq!(
            interface_tab.stab_as_string(14),
            "+  Z    X         Z                          "
        );
        // Xin1 Xry Xcrz.in1
        assert_eq!(
            interface_tab.stab_as_string(15),
            "+  X    X         X                          "
        );
        // Xin2 Zrx
        assert_eq!(
            interface_tab.stab_as_string(16),
            "+    X Z                                     "
        );
        // -Yin1 Zry (the negation is because of the transposition for inputs)
        assert_eq!(
            interface_tab.stab_as_string(17),
            "-  Y    Z                                    "
        );
        // Zin0 Zmeas
        assert_eq!(
            interface_tab.stab_as_string(18),
            "+Z          Z                                "
        );
        // Zin Zrz
        assert_eq!(
            interface_tab.stab_as_string(19),
            "+Z         Z                                 "
        );
        // Zin Ztdag
        assert_eq!(
            interface_tab.stab_as_string(20),
            "+Z        Z                                  "
        );
        // Zin Zt
        assert_eq!(
            interface_tab.stab_as_string(21),
            "+Z       Z                                   "
        );
    }

    #[test]
    fn test_mixed_clifford() {
        // Need to cover MeasureFree, QFree and Reset
        let mut builder = FunctionBuilder::new(
            "mixed_cliff",
            Signature::new(vec![qb_t(), qb_t()], vec![qb_t()]),
        )
        .unwrap();
        let [qb0, qb1] = builder.input_wires_arr();
        let qalloc = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let cx0 = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap();
        let reset = builder
            .add_dataflow_op(TketOp::Reset, [cx0.out_wire(0)])
            .unwrap();
        let meas = builder
            .add_dataflow_op(TketOp::MeasureFree, [cx0.out_wire(1)])
            .unwrap();
        let h = builder
            .add_dataflow_op(TketOp::H, [reset.out_wire(0)])
            .unwrap();
        let cx1 = builder
            .add_dataflow_op(TketOp::CX, [h.out_wire(0), qalloc.out_wire(0)])
            .unwrap();
        let qfree = builder
            .add_dataflow_op(TketOp::QFree, [cx1.out_wire(1)])
            .unwrap();
        let hugr = builder.finish_hugr_with_outputs([cx1.out_wire(0)]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut max_settings);
        let summary = summary_map.get(&func_node).unwrap().clone();
        let in_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let out_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();

        // 2 inputs, 1 output, 1 rotation for MeasureFree, 1 classical for MeasureFree, 14 interface (1+4+2+0+2+4+1, 0 needed for MeasureFree), 1+14*2=29 role controls
        // 2+1+1+1+14+28 = 48
        assert_eq!(summary.tab.nb_qubits, 48);
        // One stabilizer per qubit in the choi state (2+1+1+14=18), plus one per flow stabilizer added (1+4+1+1+2+4=13)
        assert_eq!(summary.tab.nb_stabs, 31);
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::NodeOut(qalloc.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::RoleControl(qalloc.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::RoleControl(qalloc.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::NodeIn(cx0.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::NodeIn(cx0.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::NodeOut(cx0.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::NodeOut(cx0.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::Rotation(meas.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&16).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&17).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&18).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&19).unwrap(),
            DataflowPoint::RoleControl(cx0.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&20).unwrap(),
            DataflowPoint::SumOutPhControl(meas.node(), OutgoingPort::from(0), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&21).unwrap(),
            DataflowPoint::RoleControl(meas.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&22).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&23).unwrap(),
            DataflowPoint::NodeIn(reset.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&24).unwrap(),
            DataflowPoint::NodeOut(reset.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&25).unwrap(),
            DataflowPoint::RoleControl(h.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&26).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&27).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&28).unwrap(),
            DataflowPoint::RoleControl(h.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&29).unwrap(),
            DataflowPoint::NodeIn(h.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&30).unwrap(),
            DataflowPoint::NodeOut(h.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&31).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&32).unwrap(),
            DataflowPoint::RoleControl(h.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&33).unwrap(),
            DataflowPoint::RoleControl(h.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&34).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&35).unwrap(),
            DataflowPoint::NodeIn(cx1.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&36).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&37).unwrap(),
            DataflowPoint::NodeIn(cx1.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&38).unwrap(),
            DataflowPoint::NodeOut(cx1.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&39).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&40).unwrap(),
            DataflowPoint::NodeOut(cx1.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&41).unwrap(),
            DataflowPoint::RoleControl(qfree.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&42).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&43).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&44).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&45).unwrap(),
            DataflowPoint::RoleControl(cx1.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&46).unwrap(),
            DataflowPoint::RoleControl(qfree.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&47).unwrap(),
            DataflowPoint::NodeIn(qfree.node(), IncomingPort::from(0), vec![])
        );

        // Post-select onto flow information
        let mut flow_tab = summary.tab.clone();
        // List the rcs to post-select (grouped by node)
        let rcs = vec![
            6, 7, 16, 17, 18, 19, 1, 3, 8, 10, 26, 27, 13, 22, 21, 32, 33, 25, 28, 42, 43, 44, 45,
            5, 31, 34, 36, 41, 46,
        ];
        let flow_post_selects: Vec<(usize, PauliXZ, bool)> =
            rcs.iter().map(|q| (*q, PauliXZ::X, false)).collect_vec();
        flow_tab.post_select_1qs(&flow_post_selects);
        // Pick the echelon order to make interpreting the stabilizers easier
        // outs, meas, meas.out, ins
        let flow_order = vec![39, 15, 20, 0, 2]
            .iter()
            .map(|q| [(*q, PauliXZ::Z), (*q, PauliXZ::X)])
            .flatten()
            .collect_vec();
        flow_tab.echelon(&flow_order);
        // The remaining state can be purified onto 2 inputs, 1 output, 2 discards (QFree and Reset), 1 basic rotation and 1 classical value (MeasureFree), so a complete set of stabilizers would have 7, but the MeasureFree is short one and the discards remove 2 each, leaving 2 left
        assert_eq!(flow_tab.nb_stabs, 2);
        // Zin0 Zin1 Zmeas
        assert_eq!(
            flow_tab.stab_as_string(0),
            "+Z Z            Z                                "
        );
        // Xin1 Xmeas
        assert_eq!(
            flow_tab.stab_as_string(1),
            "+Z Z                 Z                           "
        );

        // Post-select onto interface information
        let mut interface_tab = summary.tab.clone();
        let interface_post_selects: Vec<(usize, PauliXZ, bool)> =
            rcs.iter().map(|q| (*q, PauliXZ::Z, false)).collect_vec();
        interface_tab.post_select_1qs(&interface_post_selects);
        // Pick the echelon order to make interpreting the stabilizers easier
        // qfree.ins, cx1.outs, cx1.ins, h.outs, h.ins, reset.outs, reset.ins, cx0.outs, cx0.ins, qalloc.outs, [flow_order]
        let interface_order = vec![
            47, 38, 40, 35, 37, 29, 30, 23, 24, 12, 14, 9, 11, 4, 6, 7, 16, 17, 18, 19, 1, 3, 8,
            10, 26, 27, 13, 22, 32, 33, 25, 28, 42, 43, 44, 45, 5, 31, 34, 36, 41, 46,
        ]
        .iter()
        .map(|q| [(*q, PauliXZ::Z), (*q, PauliXZ::X)])
        .flatten()
        .collect_vec();
        interface_tab.echelon(&interface_order);
        // 18 qubits left in the Choi state, so we get a complete set of 18 stabilizers
        // Note that 2+18 != 30, i.e. in projecting into these two cases, we have lost 10 stabilizers. Those 10 lost are those that depend on taking the flow data of some gates and the interface data of others. For example, the fact that X flows through the target of cx0 requires us to post-select onto flow for cx0 and onto interface for reset, since the main inference from it (Xin1 Xreset.in0) is removed when we factor in the discarding effect of the reset. Testing the full tableau might be preferable, but it is definitely much harder to a human to decipher enough to write the damn test
        assert_eq!(interface_tab.nb_stabs, 18);
        // Zcx1.out1 Zqfree.in0
        assert_eq!(
            interface_tab.stab_as_string(0),
            "+                                        Z      Z"
        );
        // Xcx1.out1 Xqfree.in0
        assert_eq!(
            interface_tab.stab_as_string(1),
            "+                                        X      X"
        );
        // Zcx1.out0 Zout0
        assert_eq!(
            interface_tab.stab_as_string(2),
            "+                                      ZZ        "
        );
        // Xcx1.out0 Xout0
        assert_eq!(
            interface_tab.stab_as_string(3),
            "+                                      XX        "
        );
        // Zh.out0 Zcx1.in0
        assert_eq!(
            interface_tab.stab_as_string(4),
            "+                              Z    Z            "
        );
        // Xh.out0 Xcx1.in0
        assert_eq!(
            interface_tab.stab_as_string(5),
            "+                              X    X            "
        );
        // Zqalloc.out0 Zcx1.in1
        assert_eq!(
            interface_tab.stab_as_string(6),
            "+    Z                                Z          "
        );
        // Xqalloc.out0 Xcx1.in1
        assert_eq!(
            interface_tab.stab_as_string(7),
            "+    X                                X          "
        );
        // Zreset.out0 Zh.in0
        assert_eq!(
            interface_tab.stab_as_string(8),
            "+                        Z    Z                  "
        );
        // Xreset.out0 Xh.in0
        assert_eq!(
            interface_tab.stab_as_string(9),
            "+                        X    X                  "
        );
        // Zcx0.out0 Zreset.in0
        assert_eq!(
            interface_tab.stab_as_string(10),
            "+            Z          Z                        "
        );
        // Xcx0.out0 Xreset.in0
        assert_eq!(
            interface_tab.stab_as_string(11),
            "+            X          X                        "
        );
        // Zcx0.out1 Zmeas
        assert_eq!(
            interface_tab.stab_as_string(12),
            "+              ZZ                                "
        );
        // Xcx0.out1 Xmeas
        assert_eq!(
            interface_tab.stab_as_string(13),
            "+              XX                                "
        );
        // Zin0 Zcx0.in0
        assert_eq!(
            interface_tab.stab_as_string(14),
            "+Z        Z                                      "
        );
        // Xin0 Xcx0.in0
        assert_eq!(
            interface_tab.stab_as_string(15),
            "+X        X                                      "
        );
        // Zin1 Zcx0.in1
        assert_eq!(
            interface_tab.stab_as_string(16),
            "+  Z        Z                                    "
        );
        // Xin1 Xcx0.in1
        assert_eq!(
            interface_tab.stab_as_string(17),
            "+  X        X                                    "
        );
    }

    #[test]
    fn test_if_simple() {
        let mut builder =
            FunctionBuilder::new("if_simple", endo_sig(vec![qb_t(), qb_t(), bool_t()])).unwrap();
        let [qb0, qb1, b] = builder.input_wires_arr();
        let t = builder.add_dataflow_op(TketOp::T, [qb0]).unwrap();
        let [qb0] = t.outputs_arr();
        let mut cond_builder = builder
            .conditional_builder(
                ([type_row![], type_row![]], b),
                [(qb_t(), qb0), (qb_t(), qb1)],
                vec![qb_t(); 2].into(),
            )
            .unwrap();
        let mut cond0_builder = cond_builder.case_builder(0).unwrap();
        let [c0q0, c0q1] = cond0_builder.input_wires_arr();
        let cx = cond0_builder
            .add_dataflow_op(TketOp::CX, [c0q0, c0q1])
            .unwrap();
        cond0_builder
            .finish_with_outputs([cx.out_wire(0), cx.out_wire(1)])
            .ok();
        let cond1_builder = cond_builder.case_builder(1).unwrap();
        let [c1c0, c1q1] = cond1_builder.input_wires_arr();
        cond1_builder.finish_with_outputs([c1c0, c1q1]).ok();
        let cond = cond_builder.finish_sub_container().unwrap();
        let [qb0, qb1] = cond.outputs_arr();
        let tdg = builder.add_dataflow_op(TketOp::Tdg, [qb0]).unwrap();
        let [qb0] = tdg.outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut max_settings);
        let summary = summary_map.get(&func_node).unwrap().clone();
        // 3 inputs, 3 outputs, 2 basic rotations, 4 interface + 8 role for CX, 5 interface + 12 sum controls (case tableaus have 12 and 4 stabs which can be paired with anti-commuting ones) + 16 role (8 interface stabs and 12 flow stabs, can find 4 anti-commuting pairs) for cond
        assert_eq!(summary.tab.nb_qubits, 53);
        // CX summary has 12 stabs (8 interface, 4 flow)
        // Case summaries have 12 and 4 stabs, so 16 total in flow for conditional
        // Overall choi state is 2 input qubits, 2 output qubits, 2 basic rotations, 4 interface qubits for 10 stabs + the 16 from the conditional flow + 2 from propagating classical bits
        assert_eq!(summary.tab.nb_stabs, 28);
        let in_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let out_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();

        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 12)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 13)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::SumOutPhControl(in_node, OutgoingPort::from(2), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::Rotation(t.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::SumInPhControl(cond.node(), IncomingPort::from(0), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 14)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::NodeIn(cond.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 15)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::NodeIn(cond.node(), IncomingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::NodeOut(cond.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::NodeOut(cond.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&16).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&17).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&18).unwrap(),
            DataflowPoint::NodeIn(cx.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&19).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&20).unwrap(),
            DataflowPoint::NodeIn(cx.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&21).unwrap(),
            DataflowPoint::NodeOut(cx.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&22).unwrap(),
            DataflowPoint::NodeOut(cx.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&23).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&24).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&25).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&26).unwrap(),
            DataflowPoint::RoleControl(cx.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&27).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&28).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&29).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&30).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&31).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&32).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&33).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&34).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&35).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 8)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&36).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 9)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&37).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 10)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&38).unwrap(),
            DataflowPoint::SumInControl(cond.node(), IncomingPort::from(0), vec![], 0, 11)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&39).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&40).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&41).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&42).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&43).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&44).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&45).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&46).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&47).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 8)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&48).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 9)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&49).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 10)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&50).unwrap(),
            DataflowPoint::RoleControl(cond.node(), 11)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&51).unwrap(),
            DataflowPoint::Rotation(tdg.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&52).unwrap(),
            DataflowPoint::SumInPhControl(out_node, IncomingPort::from(2), vec![], 0)
        );

        let mut top_level_tab = summary.tab.clone();
        let cond_rcs = vec![39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 1, 3, 7, 9];
        let cond_interface_post_selects: Vec<(usize, PauliXZ, bool)> = cond_rcs
            .iter()
            .map(|q| (*q, PauliXZ::Z, false))
            .collect_vec();
        top_level_tab.post_select_1qs(&cond_interface_post_selects);
        // 2in qubits, 2out qubits, 2 basic rotations, 4 interface qubits on conditional, 2 classical propagations
        assert_eq!(top_level_tab.nb_stabs, 12);
        top_level_tab.echelon(&top_level_tab.all_columns());
        // Xin0 Xt Xcond.in1
        assert_eq!(
            top_level_tab.stab_as_string(0),
            "+X    X  X                                            "
        );
        // Zin0 Zcond.in1
        assert_eq!(
            top_level_tab.stab_as_string(1),
            "+Z       Z                                            "
        );
        // Xin1 Xcond.in2
        assert_eq!(
            top_level_tab.stab_as_string(2),
            "+  X       X                                          "
        );
        // Zin1 Zcond.in2
        assert_eq!(
            top_level_tab.stab_as_string(3),
            "+  Z       Z                                          "
        );
        // Zin.sumoutph2 Zout.suminph2
        assert_eq!(
            top_level_tab.stab_as_string(4),
            "+    Z                                               Z"
        );
        // Zt Zcond.in1
        assert_eq!(
            top_level_tab.stab_as_string(5),
            "+     Z  Z                                            "
        );
        // Zcond.suminph0 Zout.suminph2
        assert_eq!(
            top_level_tab.stab_as_string(6),
            "+      Z                                             Z"
        );
        // Xcond.out0 Xout0 Xtdg
        assert_eq!(
            top_level_tab.stab_as_string(7),
            "+           XX                                      X "
        );
        // Zcond.out0 Ztdg
        assert_eq!(
            top_level_tab.stab_as_string(8),
            "+           Z                                       Z "
        );
        // Zout0 Ztdg
        assert_eq!(
            top_level_tab.stab_as_string(9),
            "+            Z                                      Z "
        );
        // Xcond.out1 Xout1
        assert_eq!(
            top_level_tab.stab_as_string(10),
            "+             XX                                      "
        );
        // Zcond.out1 Zout1
        assert_eq!(
            top_level_tab.stab_as_string(11),
            "+             ZZ                                      "
        );

        let mut if_cx_interface_tab = summary.tab.clone();
        let cond_controls = vec![27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38];
        let cx_rcs = vec![23, 24, 25, 26, 15, 16, 17, 19];
        let if_cx_interface_post_selects: Vec<(usize, PauliXZ, bool)> = chain!(
            cond_rcs.iter().map(|q| (*q, PauliXZ::X, false)),
            [(6, PauliXZ::Z, false)], // Phase control
            cond_controls.iter().map(|q| (*q, PauliXZ::Z, false)),
            cx_rcs.iter().map(|q| (*q, PauliXZ::Z, false)),
        )
        .collect_vec();
        if_cx_interface_tab.post_select_1qs(&if_cx_interface_post_selects);
        // 2in qubits, 2out qubits, 2 basic rotations, 4 interface qubits on cx, 2 deduced classical values
        assert_eq!(if_cx_interface_tab.nb_stabs, 12);
        if_cx_interface_tab.echelon(&if_cx_interface_tab.all_columns());
        // Xin0 Xt Xcx.in0
        assert_eq!(
            if_cx_interface_tab.stab_as_string(0),
            "+X    X            X                                  "
        );
        // Zin0 Zcx.in0
        assert_eq!(
            if_cx_interface_tab.stab_as_string(1),
            "+Z                 Z                                  "
        );
        // Xin1 Xcx.in1
        assert_eq!(
            if_cx_interface_tab.stab_as_string(2),
            "+  X                 X                                "
        );
        // Zin1 Zcx.in1
        assert_eq!(
            if_cx_interface_tab.stab_as_string(3),
            "+  Z                 Z                                "
        );
        // Zin.sumoutph2
        assert_eq!(
            if_cx_interface_tab.stab_as_string(4),
            "+    Z                                                "
        );
        // Zt Zcx.in0
        assert_eq!(
            if_cx_interface_tab.stab_as_string(5),
            "+     Z            Z                                  "
        );
        // Xout0 Xcx.out0 Xtdg
        assert_eq!(
            if_cx_interface_tab.stab_as_string(6),
            "+            X        X                             X "
        );
        // Zout0 Ztdg
        assert_eq!(
            if_cx_interface_tab.stab_as_string(7),
            "+            Z                                      Z "
        );
        // Xout1 Xcx.out1
        assert_eq!(
            if_cx_interface_tab.stab_as_string(8),
            "+              X       X                              "
        );
        // Zout1 Zcx.out1
        assert_eq!(
            if_cx_interface_tab.stab_as_string(9),
            "+              Z       Z                              "
        );
        // Zcx.out0 Ztdg
        assert_eq!(
            if_cx_interface_tab.stab_as_string(10),
            "+                     Z                             Z "
        );
        // Zout.suminph2
        assert_eq!(
            if_cx_interface_tab.stab_as_string(11),
            "+                                                    Z"
        );

        let mut if_cx_flow_tab = summary.tab.clone();
        let if_cx_flow_post_selects: Vec<(usize, PauliXZ, bool)> = chain!(
            cond_rcs.iter().map(|q| (*q, PauliXZ::X, false)),
            [(6, PauliXZ::Z, false)], // Phase control
            cond_controls.iter().map(|q| (*q, PauliXZ::Z, false)),
            cx_rcs.iter().map(|q| (*q, PauliXZ::X, false)),
        )
        .collect_vec();
        if_cx_flow_tab.post_select_1qs(&if_cx_flow_post_selects);
        // 2in qubits, 2out qubits, 2 basic rotations, 2 deduced classical values
        assert_eq!(if_cx_flow_tab.nb_stabs, 8);
        if_cx_flow_tab.echelon(&if_cx_flow_tab.all_columns());
        // Xin0 Xt Xout0 Xout1 Xtdg
        assert_eq!(
            if_cx_flow_tab.stab_as_string(0),
            "+X    X      X X                                    X "
        );
        // Zin0 Ztdg
        assert_eq!(
            if_cx_flow_tab.stab_as_string(1),
            "+Z                                                  Z "
        );
        // Xin1 Xout1
        assert_eq!(
            if_cx_flow_tab.stab_as_string(2),
            "+  X           X                                      "
        );
        // Zin1 Zout1 Ztdg
        assert_eq!(
            if_cx_flow_tab.stab_as_string(3),
            "+  Z           Z                                    Z "
        );
        // Zin.sumoutph2
        assert_eq!(
            if_cx_flow_tab.stab_as_string(4),
            "+    Z                                                "
        );
        // Zt Ztdg
        assert_eq!(
            if_cx_flow_tab.stab_as_string(5),
            "+     Z                                             Z "
        );
        // Zout0 Ztdg
        assert_eq!(
            if_cx_flow_tab.stab_as_string(6),
            "+            Z                                      Z "
        );
        // Zout.suminph2
        assert_eq!(
            if_cx_flow_tab.stab_as_string(7),
            "+                                                    Z"
        );

        let mut else_tab = summary.tab.clone();
        let else_post_selects: Vec<(usize, PauliXZ, bool)> = chain!(
            cond_rcs.iter().map(|q| (*q, PauliXZ::X, false)),
            [(6, PauliXZ::Z, true)], // Phase control
            cond_controls.iter().map(|q| (*q, PauliXZ::X, false)),
        )
        .collect_vec();
        else_tab.post_select_1qs(&else_post_selects);
        // 2in qubits, 2out qubits, 2 basic rotations, 2 deduced classical values
        assert_eq!(else_tab.nb_stabs, 8);
        else_tab.echelon(&else_tab.all_columns());
        // Xin0 Xt Xout0 Xtdg
        assert_eq!(
            else_tab.stab_as_string(0),
            "+X    X      X                                      X "
        );
        // Zin0 Ztdg
        assert_eq!(
            else_tab.stab_as_string(1),
            "+Z                                                  Z "
        );
        // Xin1 Xout1
        assert_eq!(
            else_tab.stab_as_string(2),
            "+  X           X                                      "
        );
        // Zin1 Zout1
        assert_eq!(
            else_tab.stab_as_string(3),
            "+  Z           Z                                      "
        );
        // -Zin.sumoutph2
        assert_eq!(
            else_tab.stab_as_string(4),
            "-    Z                                                "
        );
        // Zt Ztdg
        assert_eq!(
            else_tab.stab_as_string(5),
            "+     Z                                             Z "
        );
        // Zout0 Ztdg
        assert_eq!(
            else_tab.stab_as_string(6),
            "+            Z                                      Z "
        );
        // -Zout.suminph2
        assert_eq!(
            else_tab.stab_as_string(7),
            "-                                                    Z"
        );
    }

    #[test]
    fn test_loop_null() {
        let mut builder =
            FunctionBuilder::new("loop_null", endo_sig(vec![qb_t(), qb_t(), bool_t()])).unwrap();
        let [qb0, qb1, b] = builder.input_wires_arr();
        let reset = builder.add_dataflow_op(TketOp::Reset, [qb0]).unwrap();
        let [qb0] = reset.outputs_arr();
        let t = builder.add_dataflow_op(TketOp::T, [qb1]).unwrap();
        let [qb1] = t.outputs_arr();
        let mut loop_builder = builder
            .tail_loop_builder([(bool_t(), b)], [(qb_t(), qb0), (qb_t(), qb1)], type_row![])
            .unwrap();
        let [_, loop_qb0, loop_qb1] = loop_builder.input_wires_arr();
        let loop_t = loop_builder.add_dataflow_op(TketOp::T, [loop_qb0]).unwrap();
        let [loop_qb0] = loop_t.outputs_arr();
        let loop_tdg = loop_builder
            .add_dataflow_op(TketOp::Tdg, [loop_qb1])
            .unwrap();
        let [loop_qb1] = loop_tdg.outputs_arr();
        let loop_b = loop_builder
            .make_break(loop_builder.loop_signature().unwrap().clone(), [])
            .unwrap();
        let tl = loop_builder
            .finish_with_outputs(loop_b, [loop_qb0, loop_qb1])
            .unwrap();
        let [qb0, qb1] = tl.outputs_arr();
        let tdg = builder.add_dataflow_op(TketOp::Tdg, [qb1]).unwrap();
        let [qb1] = tdg.outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut max_settings);
        let summary = summary_map.get(&func_node).unwrap().clone();
        let in_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let out_node = hugr
            .children(func_node)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let loop_in_node = hugr
            .children(tl.node())
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        let loop_out_node = hugr
            .children(tl.node())
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_)))
            .exactly_one()
            .ok()
            .unwrap();
        // 3ins, 3outs, 4 basic rotations, 2 interface + 4 role for reset, 3loop-in, 3loop-out, 1tag, 5 interface + 16 role for loop
        assert_eq!(summary.tab.nb_qubits, 44);
        // For top-level, a complete set of stabilizers for 2in qubits, 2out qubits, 2 basic rotations, 6 interface (reset and loop), and 2 classical propagations
        // Plus 1 flow from reset
        // For loop body, a complete set for 2in qubits, 2out qubits, 2 basic rotations, 1 classical constant, 1 classical propagation
        // 2 proper loop invariants (Z commutes through on each qubit) and a redundant stabilizer from deterministic post-selection during the break case
        // (2+2+2+6+2)+1+(2+2+2+1+1)+(2+1) = 14+1+8+3 = 26
        assert_eq!(summary.tab.nb_stabs, 26);
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::NodeOut(in_node, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 12)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::SumOutPhControl(in_node, OutgoingPort::from(2), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::Rotation(t.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::NodeIn(reset.node(), IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::NodeOut(reset.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 13)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::RoleControl(reset.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::SumInPhControl(tl.node(), IncomingPort::from(0), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 14)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::NodeIn(tl.node(), IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 15)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&16).unwrap(),
            DataflowPoint::NodeIn(tl.node(), IncomingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&17).unwrap(),
            DataflowPoint::NodeOut(tl.node(), OutgoingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&18).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(0), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&19).unwrap(),
            DataflowPoint::NodeOut(tl.node(), OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&20).unwrap(),
            DataflowPoint::NodeIn(out_node, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&21).unwrap(),
            DataflowPoint::SumOutPhControl(loop_in_node, OutgoingPort::from(0), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&22).unwrap(),
            DataflowPoint::NodeOut(loop_in_node, OutgoingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&23).unwrap(),
            DataflowPoint::NodeIn(loop_out_node, IncomingPort::from(1), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&24).unwrap(),
            DataflowPoint::NodeOut(loop_in_node, OutgoingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&25).unwrap(),
            DataflowPoint::NodeIn(loop_out_node, IncomingPort::from(2), vec![])
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&26).unwrap(),
            DataflowPoint::SumOutPhControl(loop_b.node(), OutgoingPort::from(0), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&27).unwrap(),
            DataflowPoint::Rotation(loop_tdg.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&28).unwrap(),
            DataflowPoint::Rotation(loop_t.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&29).unwrap(),
            DataflowPoint::SumInPhControl(loop_out_node, IncomingPort::from(0), vec![], 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&30).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 0)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&31).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 1)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&32).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 2)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&33).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 3)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&34).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 4)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&35).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 5)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&36).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 6)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&37).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 7)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&38).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 8)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&39).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 9)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&40).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 10)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&41).unwrap(),
            DataflowPoint::RoleControl(tl.node(), 11)
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&42).unwrap(),
            DataflowPoint::Rotation(tdg.node())
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&43).unwrap(),
            DataflowPoint::SumInPhControl(out_node, IncomingPort::from(2), vec![], 0)
        );

        // Post-select onto loop interface
        let mut loop_interface_tab = summary.tab.clone();
        let loop_rcs = vec![30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 3, 9, 13, 15];
        let loop_interface_post_selects: Vec<(usize, PauliXZ, bool)> = loop_rcs
            .iter()
            .map(|q| (*q, PauliXZ::Z, false))
            .collect_vec();
        loop_interface_tab.post_select_1qs(&loop_interface_post_selects);
        // From previous calculations, we expect 14+1+8+1 stabilizers (include both roles from reset, loop body, and the trivial stabilizer)
        assert_eq!(loop_interface_tab.nb_stabs, 24);
        loop_interface_tab.echelon(&loop_interface_tab.all_columns());
        // Xin0 Zreset.rc2 Xreset.in Zreset.rc1
        assert_eq!(
            loop_interface_tab.stab_as_string(0),
            "+XZ     X   Z                                "
        );
        // Zin0 Zreset.rc2 Zreset.in
        assert_eq!(
            loop_interface_tab.stab_as_string(1),
            "+ZZ     Z                                    "
        );
        // Xin1 Xt Xloop.in2
        assert_eq!(
            loop_interface_tab.stab_as_string(2),
            "+  X  X          X                           "
        );
        // Zin1 Zloop.in2
        assert_eq!(
            loop_interface_tab.stab_as_string(3),
            "+  Z             Z                           "
        );
        // Zin.sumoutph2 Zout.suminph2
        assert_eq!(
            loop_interface_tab.stab_as_string(4),
            "+    Z                                      Z"
        );
        // Zt Zloop.in2
        assert_eq!(
            loop_interface_tab.stab_as_string(5),
            "+     Z          Z                           "
        );
        // Zreset.rc3 Zreset.out Zloop.in1
        assert_eq!(
            loop_interface_tab.stab_as_string(6),
            "+      Z Z     Z                             "
        );
        // Xreset.out Zreset.rc0 Xloop.in1
        assert_eq!(
            loop_interface_tab.stab_as_string(7),
            "+        X Z   X                             "
        );
        // Xreset.rc0 Zloop.in1
        assert_eq!(
            loop_interface_tab.stab_as_string(8),
            "+          X   Z                             "
        );
        // Zloop.suminph0 Zout2
        assert_eq!(
            loop_interface_tab.stab_as_string(9),
            "+            Z                              Z"
        );
        // Xloop.out0 Xout0
        assert_eq!(
            loop_interface_tab.stab_as_string(10),
            "+                 XX                         "
        );
        // Zloop.out0 Zout0
        assert_eq!(
            loop_interface_tab.stab_as_string(11),
            "+                 ZZ                         "
        );
        // Xloop.out1 Xout1 Xtdg
        assert_eq!(
            loop_interface_tab.stab_as_string(12),
            "+                   XX                     X "
        );
        // Zloop.out1 Ztdg
        assert_eq!(
            loop_interface_tab.stab_as_string(13),
            "+                   Z                      Z "
        );
        // Zout1 Ztdg
        assert_eq!(
            loop_interface_tab.stab_as_string(14),
            "+                    Z                     Z "
        );
        // Xloopin1 Xloopout1 Xloopt
        assert_eq!(
            loop_interface_tab.stab_as_string(15),
            "+                      XX    X               "
        );
        // Zloopin1 Zloopt
        assert_eq!(
            loop_interface_tab.stab_as_string(16),
            "+                      Z     Z               "
        );
        // Zloopout1 Zloopt
        assert_eq!(
            loop_interface_tab.stab_as_string(17),
            "+                       Z    Z               "
        );
        // Xloopin2 Xloopout2 Xlooptdg
        assert_eq!(
            loop_interface_tab.stab_as_string(18),
            "+                        XX X                "
        );
        // Zloopin2 Zlooptdg
        assert_eq!(
            loop_interface_tab.stab_as_string(19),
            "+                        Z  Z                "
        );
        // Zloopout2 Zlooptdp
        assert_eq!(
            loop_interface_tab.stab_as_string(20),
            "+                         Z Z                "
        );
        // -Zloopbreak
        assert_eq!(
            loop_interface_tab.stab_as_string(21),
            "-                          Z                 "
        );
        // -Zloopout0
        assert_eq!(
            loop_interface_tab.stab_as_string(22),
            "-                             Z              "
        );
        // Trivial stabilizer (expected because post-selection of the break branch during invariant finding trivially succeeds but we don't remove any stabilizers; then this propagates through joins and controlled conditionals because it is trivially generated by any other tableau)
        assert_eq!(
            loop_interface_tab.stab_as_string(23),
            "+                                            "
        );

        // Post-select onto loop flow (and reset flow)
        let mut loop_flow_tab = summary.tab.clone();
        let reset_rcs = vec![10, 11, 1, 6];
        let loop_flow_post_selects: Vec<(usize, PauliXZ, bool)> = chain!(
            loop_rcs.iter().map(|q| (*q, PauliXZ::X, false)),
            reset_rcs.iter().map(|q| (*q, PauliXZ::X, false)),
        )
        .collect_vec();
        loop_flow_tab.post_select_1qs(&loop_flow_post_selects);
        // Remaining qubits are 2inputs, 2outputs, 2 basic rotations; we lose 1 stabilizer from the reset, 2 from the loop invariants (though one of these was already lost from the reset); plus 2 classical propagations and the trivial stabilizer
        // 2+2+2-(1+1)+2+1 = 7
        assert_eq!(loop_flow_tab.nb_stabs, 7);
        loop_flow_tab.echelon(&loop_flow_tab.all_columns());
        // Zin1 Ztdg
        assert_eq!(
            loop_flow_tab.stab_as_string(0),
            "+  Z                                       Z "
        );
        // Zin2 Zout2
        assert_eq!(
            loop_flow_tab.stab_as_string(1),
            "+    Z                                      Z"
        );
        // Zt Ztdg
        assert_eq!(
            loop_flow_tab.stab_as_string(2),
            "+     Z                                    Z "
        );
        // Zloop.in0 Zout2
        assert_eq!(
            loop_flow_tab.stab_as_string(3),
            "+            Z                              Z"
        );
        // Zout0
        assert_eq!(
            loop_flow_tab.stab_as_string(4),
            "+                  Z                         "
        );
        // Zout1 Ztdg
        assert_eq!(
            loop_flow_tab.stab_as_string(5),
            "+                    Z                     Z "
        );
        // Trivial stabilizer
        assert_eq!(
            loop_flow_tab.stab_as_string(6),
            "+                                            "
        );
    }
}
