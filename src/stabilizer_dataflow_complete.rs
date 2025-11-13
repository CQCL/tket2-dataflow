use std::{
    cmp::max,
    collections::{HashMap, HashSet},
    convert::identity,
    hash::Hash,
};

use bimap::BiHashMap;
use hugr::{
    extension::{prelude::qb_t, simple_op::MakeExtensionOp},
    hugr::patch::outline_cfg,
    ops::{sum, OpTag, OpTrait, OpType},
    std_extensions::logic::LogicOp,
    types::{TypeBase, TypeEnum, TypeRV},
    HugrView, IncomingPort, OutgoingPort, Port, PortIndex,
};
use hugr_core::hugr::internal::PortgraphNodeMap;
use itertools::{chain, Either, Itertools};
use petgraph::visit as pv;
use tket::{extension::rotation, TketOp};

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
    flow_level: DataflowFlowDetail,
    /// Whether we track the external interface separately from any flow information.
    /// This is required in order to resynthesise from the tableau.
    /// Relationships between the interface qubits of different nodes can identify when phase folding (or similar gate merging) can be applied.
    /// Turning this on adds an additional qubit for each qubit input or output port on the given node.
    /// If this is on and flow_level is not DataflowFlowDetail::None, additional control qubits will be required to toggle between them.
    include_external_interface: bool,
    /// When qubits exist within a Sum type, do we attempt to track information about the qubits or leave it as a point of inaccuracy.
    /// Turning this on can enable more branch-dependent phase folding opportunities at the cost of more qubits for both the qubits themselves and control qubits for the branches of the sums.
    /// Examples of where this is useful include TryQAlloc gates and internalising the classical outcome of a measurement operation.
    include_sum_types: bool,
    /// For loops, whether or not we include a full copy of the loop body in the tableau.
    /// This is required in order to resynthesise loops from the tableau for the enclosing scope (i.e. without maintaining a separate tableau for the nested block).
    /// When flow_level is DataflowFlowDetail::All, the qubits for the internal nodes within the loop body are reused, but additional qubits will still be created for the input and output nodes within the loop body.
    include_loop_body: bool,
    /// For each gate, we typically want both the external interface to tell us where the gate occurs in the circuit and some flow information across the gate. For an Rz gate, this would take 2 qubits for the external interface and 4 control qubits to turn off the external interface when we want to focus on the flow information; giving a total of 8 qubits (2 ins, 2 outs, 4 controls) and 5 stabilizers.
    /// However, for Rz (and relatives like Rx, Ry, T), we can use 3 qubits (1 in, 1 out, 1 special) and 3 stabilizers to capture all the information we need for synthesis and flow - we essentially model it as a ZX spider with an open wire onto which we can post-select the angle of rotation, similar to its MBQC implementation.
    /// This exploits the fact that we don't need to load specific Pauli strings into both Z and X during resynthesis; an Rz gate only introduces a phase depending on the Pauli string loaded into Z.
    /// A stabilizer involving Z on the rotation's qubit tells us the Pauli string it acts around, and one involving the X says the rest of the stabilizer is dependent on the angle of rotation chosen, i.e. indicating the flow of causal influence through the circuit.
    /// This setting overrides all others for rotation gates to use this shorthand trick.
    /// The same trick of treating it as a ZX spider works for both Measure and MeasureFree, in which case the special qubit can just be treated as the classical output bit (DataflowPoint::SumOutPhControl) when include_sum_types == true.
    rotation_override: bool,
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

#[derive(Clone, PartialEq, Eq, Hash)]
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
    // /// Similar to the distinction between NodeIn and TempIn, we use SumInControl for the classical values incident on a node compared to ExtInControl for its connection to the environment (the classical value we have to match with ExtOutControls from the predecessor).
    // ExtInControl(N, IncomingPort, Vec<usize>, usize, usize),
    // /// Analogous to ExtInControl but for SumInPhControl
    // ExtInPhControl(N, IncomingPort, Vec<usize>, usize),
    // /// Analogous to ExtInControl but for SumOutControl
    // ExtOutControl(N, OutgoingPort, Vec<usize>, usize, usize),
    // /// Analogous to ExtInControl but for SumOutPhControl
    // ExtOutPhControl(N, OutgoingPort, Vec<usize>, usize),
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
    fn clone(&self) -> Self {
        StabilizerDataflow {
            tab: self.tab.clone(),
            q_index_map: self.q_index_map.clone(),
        }
    }
}

impl<H: HugrView> StabilizerDataflow<H> {
    fn tensor_product_in_place(&mut self, other: &Self) {
        // Assumes self and other use disjoint qubits
        let offset = self.tab.add_qubits(other.tab.nb_qubits);
        for (dfp, index) in other.q_index_map.iter() {
            self.q_index_map.insert(dfp.clone(), offset + index);
        }
        for i in 0..other.tab.nb_stabs {
            let mut new_z = BitVector::new(offset);
            new_z.extend_vec(other.tab.z[i].get_boolean_vec(), offset);
            let mut new_x = BitVector::new(offset);
            new_x.extend_vec(other.tab.x[i].get_boolean_vec(), offset);
            self.tab.add_stab(new_z, new_x, other.tab.signs.get(i));
        }
    }

    fn tensor_product(left: &Self, right: &Self) -> Self {
        // Assumes left and right use disjoint qubits
        let mut res = left.clone();
        res.tensor_product_in_place(right);
        res
    }

    // When building identity wires for output interfaces, we may wish to use SumOutControls instead of SumInControls, and combining flow and interface data uses the same logic as for conditionals.
    fn controlled_conditional(
        left: &Self,
        right: &Self,
        control_node: H::Node,
        control_spec: &ControlSpec,
    ) -> StabilizerDataflow<H> {
        // Unify qubit indexing then extend and reorder tableaus to match
        let mut q_ind_map: BiHashMap<DataflowPoint<H::Node>, usize> = left.q_index_map.clone();
        let mut l = left.tab.clone();
        for (point, _) in right.q_index_map.iter() {
            if !q_ind_map.contains_left(point) {
                q_ind_map.insert(point.clone(), l.add_qubit());
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
        if jd.join_range().0 != 0 && *control_spec != ControlSpec::Role {
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
    fn sum_injection(
        &mut self,
        inj_variant: Either<(), ()>,
        control_node: H::Node,
        // control_port: Either<IncomingPort, OutgoingPort>,
        control_port: OutgoingPort,
        control_row_index: Vec<usize>,
        control_sum_index: usize,
    ) {
        let control_qb = self.tab.add_qubit();
        // let dfp = match control_port {
        //     Either::Left(in_port) => DataflowPoint::ExtInPhControl(
        //         control_node,
        //         in_port,
        //         control_row_index,
        //         control_sum_index,
        //     ),
        //     Either::Right(out_port) => DataflowPoint::ExtOutPhControl(
        //         control_node,
        //         out_port,
        //         control_row_index,
        //         control_sum_index,
        //     ),
        // };
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

    // fn sum_injection_old(
    //     &mut self,
    //     data: Vec<DataflowPoint<H::Node>>,
    //     control_node: H::Node,
    //     control_port: IncomingPort,
    //     control_row_index: Vec<usize>,
    //     control_sum_index: usize,
    //     control_value: PauliXZ,
    // ) {
    //     // Echelon over data to obtain the minimal set of generators with non-identity values over them
    //     let cols: Vec<(usize, PauliXZ)> = data
    //         .iter()
    //         .map(|dfp| {
    //             let index = *self.q_index_map.get_by_left(dfp).unwrap();
    //             [(index, PauliXZ::Z), (index, PauliXZ::X)]
    //         })
    //         .flatten()
    //         .collect_vec();
    //     self.tab.echelon(&cols);
    //     // For each generator with non-identity over data, add a control qubit
    //     let mut col_iter = cols.iter();
    //     for s in 0..self.tab.nb_stabs {
    //         // Find the next pivot qubit
    //         let mut empty = true;
    //         while let Some((ind, xz)) = col_iter.next() {
    //             match xz {
    //                 PauliXZ::X => {
    //                     if self.tab.x[s].get(*ind) {
    //                         empty = false;
    //                         break;
    //                     }
    //                 }
    //                 PauliXZ::Z => {
    //                     if self.tab.z[s].get(*ind) {
    //                         empty = false;
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //         if empty {
    //             break;
    //         } else {
    //             let ctrl_dfp = DataflowPoint::SumInControl(
    //                 control_node,
    //                 control_port,
    //                 control_row_index.clone(),
    //                 control_sum_index,
    //                 s,
    //             );
    //             let new_index = self.tab.add_qubit();
    //             match control_value {
    //                 PauliXZ::X => {
    //                     self.tab.x[s].xor_bit(new_index);
    //                 }
    //                 PauliXZ::Z => {
    //                     self.tab.z[s].xor_bit(new_index);
    //                 }
    //             }
    //             self.q_index_map.insert(ctrl_dfp, new_index);
    //         }
    //     }
    // }

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
                                    // // ID(A+B) = L(ID(A)) + R(ID(B))
                                    // // Inject all data from both variants
                                    // let mut left = chunk_it[0].clone();
                                    // StabilizerDataflow::sum_injection(
                                    //     &mut left,
                                    //     Either::Left(()),
                                    //     target,
                                    //     target_port,
                                    //     row_index.clone(),
                                    //     sum_index,
                                    // );
                                    // let mut right = chunk_it[1].clone();
                                    // StabilizerDataflow::sum_injection(
                                    //     &mut right,
                                    //     Either::Right(()),
                                    //     target,
                                    //     target_port,
                                    //     row_index.clone(),
                                    //     sum_index,
                                    // );
                                    // // Combine injected copies
                                    // StabilizerDataflow::controlled_conditional(
                                    //     &left,
                                    //     &right,
                                    //     source,
                                    //     source_port,
                                    //     row_index.clone(),
                                    //     sum_index,
                                    // )

                                    // ID(A+B) = L(ID(A)) + R(ID(B))
                                    // We don't need to care about the injections since we are using the same controls for the classical variables at the source and target, so we should only need to use Sum(In/Out)Controls to pick between different sets of qubit identities
                                    StabilizerDataflow::controlled_conditional(
                                        &chunk_it[0],
                                        &chunk_it[1],
                                        source,
                                        &control_spec,
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

    // fn identity_wire(hugr: &H, source: H::Node, source_port: OutgoingPort) -> Self {
    //     let (target, target_port) = hugr.single_linked_input(source, source_port).unwrap();
    //     let (_, wire_type) = hugr
    //         .out_value_types(source)
    //         .find(|(p, _)| *p == source_port)
    //         .unwrap();
    //     Self::identity_wire_recursive(
    //         hugr,
    //         source,
    //         source_port,
    //         target,
    //         target_port,
    //         &wire_type.into(),
    //         &vec![],
    //     )
    // }

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
    fn project_to_flow(&mut self) {
        let post_selects: Vec<(usize, PauliXZ, bool)> = self
            .q_index_map
            .iter()
            .filter_map(|(dfp, q)| {
                if let DataflowPoint::RoleControl(_, _) = dfp {
                    Some((*q, PauliXZ::X, false))
                } else {
                    None
                }
            })
            .collect_vec();
        self.tab.post_select_1qs(&post_selects);
    }

    fn project_non_io(&mut self, inp: H::Node, out: H::Node) {
        let non_ios: Vec<usize> = self
            .q_index_map
            .iter()
            .filter(|(dfp, _)| match dfp {
                DataflowPoint::NodeIn(n, _, _)
                | DataflowPoint::NodeOut(n, _, _)
                | DataflowPoint::SumInControl(n, _, _, _, _)
                | DataflowPoint::SumInPhControl(n, _, _, _)
                | DataflowPoint::SumOutControl(n, _, _, _, _)
                | DataflowPoint::SumOutPhControl(n, _, _, _) => *n != inp && *n != out,
                _ => true,
            })
            .map(|(_, q)| *q)
            .collect_vec();
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
        for (dfp, qb) in self.q_index_map.iter() {
            if let DataflowPoint::TempOut(source, out_port, row_index) = dfp {
                if let Some(found) = pred_map.get(&(*source, *out_port)) {
                    if let Some(in_port) = found {
                        qbs_to_rename.push((
                            *qb,
                            DataflowPoint::NodeIn(node, *in_port, row_index.clone()),
                        ));
                    } else {
                        qbs_to_remove.push((dfp.clone(), *qb));
                    }
                }
            }
            // No need to rename control qubits as they are all strictly associated with the previous node and may be used for studying of the classical information flow

            match dfp {
                DataflowPoint::TempOut(source, out_port, row_index) => {
                    if let Some(found) = pred_map.get(&(*source, *out_port)) {
                        if let Some(in_port) = found {
                            qbs_to_rename.push((
                                *qb,
                                DataflowPoint::NodeIn(node, *in_port, row_index.clone()),
                            ));
                        } else {
                            qbs_to_remove.push((dfp.clone(), *qb));
                        }
                    }
                }
                DataflowPoint::SumOutPhControl(source, out_port, row_index, sum_index) => {
                    if let Some(found) = pred_map.get(&(*source, *out_port)) {
                        if let Some(in_port) = found {
                            ctrls_to_copy.push((
                                *qb,
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
                    tab.append_s(col1);
                    tab.append_z(col1);
                    tab.append_cx(col0, col1);
                    tab.append_s(col1);
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
                        if settings.include_sum_types {
                            // Include the classical output which copies the value in the rotation Pauli
                            let col_ctrl = tab.add_qubit();
                            q_index_map.insert(
                                DataflowPoint::SumOutPhControl(
                                    node,
                                    OutgoingPort::from(1),
                                    vec![],
                                    0,
                                ),
                                col_ctrl,
                            );
                            zz.xor_bit(col_ctrl);
                        }
                        tab.add_stab(zz, BitVector::new(tab.nb_qubits), false);
                        [col_in]
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
                    // Rather than needing a 3-ary Z spider, for MeasureFree we can use a 2-ary Z spider, i.e. an identity wire.
                    // So just rename the qubit
                    self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                        q_index_map.insert(DataflowPoint::Rotation(node), col_in);
                        if settings.include_sum_types {
                            // Copy the result to the phase control; a CX will copy the Z component to Z on another qubit
                            let ctrl = tab.add_qubit();
                            q_index_map.insert(
                                DataflowPoint::SumOutPhControl(
                                    node,
                                    OutgoingPort::from(0),
                                    vec![],
                                    0,
                                ),
                                ctrl,
                            );
                            tab.append_cx(ctrl, col_in);
                        }
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
        for (dfp, col) in self.q_index_map.iter() {
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
                            classical_propagations.push((*col, *counterpart_col));
                        }
                    }
                }
                _ => {}
            }
        }
        let mut post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        let mut delete_mask = BitVector::new(self.tab.nb_qubits);
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
        for q in halfbells {
            let col = self.q_index_map.get_by_left(&q).unwrap();
            post_selects.push((*col, PauliXZ::X, false));
            post_selects.push((*col, PauliXZ::Z, false));
            delete_mask.xor_bit(*col);
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

    // /// Used for nodes with TketOps when DataflowSettings::include_external_interface == false and DataflowSettings::flow_level != DataflowFlowDetail::None (so we assume flow_level passed in here is not None), i.e. we can apply some amount of flow information in-place
    // fn apply_flow(&mut self, hugr: &H, node: H::Node, flow: &mut StabilizerDataflow<H>, settings: DataflowSettings) {
    //     todo!("Tensor product; bell post-selection on qubit wires when the predecessor exists (removing the adjoining qubits; add ZZ to join classical values on inputs whenever predecessor exists... Given the similarity between this and apply_combined, consider whether or not to just create an \"append\" method which can either take the flow for a node or the combination of its flow and interface")
    // }

    // /// Used when DataflowSettings::include_external_interface == true and DataflowSettings::flow_level != DataflowFlowDetail::None.
    // /// This implementation handles this in a generic fashion, though for TketOps it may be faster in practice to define hand-written routines to do this for each gate type.
    // /// Assumes flow has been defined across TempIns and TempOuts to describe the node on its own.
    // fn apply_combined(&mut self, hugr: &H, node: H::Node, flow: &mut StabilizerDataflow<H>, settings: DataflowSettings) {
    //     // Start by analysing flow to find the maximal number of stabilizers independent of the interface; we will need one role control qubit for each of these
    //     let boundary_qbs: Vec<usize> = flow.q_index_map.iter().filter_map(|(dfp, q)| {
    //         match dfp {
    //             DataflowPoint::TempIn(_, _, _) | DataflowPoint::TempOut(_, _, _) => { Some(*q) }
    //             // Any content on Ext(In/Out)(Ph)Controls will commute with the interface since the identity wires only use the phase controls and everything on phase controls is guaranteed to be in the Z basis
    //             _ => { None }
    //         }
    //     }).collect_vec();
    //     let cols_to_reduce = chain!(
    //         boundary_qbs.iter().map(|q| (*q, PauliXZ::X)),
    //         boundary_qbs.iter().map(|q| (*q, PauliXZ::Z)),
    //     ).collect_vec();
    //     flow.tab.echelon(&cols_to_reduce);
    //     let mut first_non_interface_stab = 0;
    //     for (q, pauli) in cols_to_reduce {
    //         match pauli {
    //             PauliXZ::X => {
    //                 if flow.tab.x[first_non_interface_stab].get(q) {
    //                     first_non_interface_stab += 1;
    //                 }
    //             }
    //             PauliXZ::Z => {
    //                 if flow.tab.z[first_non_interface_stab].get(q) {
    //                     first_non_interface_stab += 1;
    //                 }
    //             }
    //         }
    //     }
    //     // As we now know how many role qubits we need, we can allocate them
    //     let num_commuting_flow_stabs = flow.tab.nb_stabs - first_non_interface_stab;
    //     todo!("Continue checking whether or not classical data messes up with the procedure here; may want to backtrack through everything else to incorporate identities over phase controls and the ability for appends to transfer ExtOut(Ph)Controls as they would for qubits and TempOuts.")
    // }
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
                        for (dfp, col) in summary.q_index_map.iter() {
                            match dfp {
                                DataflowPoint::TempOut(n, out_port, ri) => {
                                    if let Some(in_port) = pred_map.get(&(*n, *out_port)) {
                                        let mut new_ri = vec![in_port.index()];
                                        new_ri.extend(ri);
                                        renames.push((
                                            *col,
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
                                            *col,
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
        let n_common_ins = cond_op.other_inputs.len();
        let mut bodies_to_join = hugr
            .children(node)
            .map(|cond_node| {
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
                    .enumerate()
                    .filter_map(|(cond_i, (dfp, q))| match dfp {
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
                                if o.index() < n_common_ins {
                                    Some((
                                        DataflowPoint::TempIn(
                                            node,
                                            IncomingPort::from(o.index()),
                                            ri.clone(),
                                        ),
                                        *q,
                                    ))
                                } else {
                                    let mut new_ri: Vec<usize> = vec![cond_i];
                                    new_ri.extend(ri);
                                    Some((
                                        DataflowPoint::TempIn(
                                            node,
                                            IncomingPort::from(n_common_ins),
                                            new_ri,
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
                                if o.index() < n_common_ins {
                                    Some((
                                        DataflowPoint::SumInControl(
                                            node,
                                            IncomingPort::from(o.index()),
                                            ri.clone(),
                                            *si,
                                            *m,
                                        ),
                                        *q,
                                    ))
                                } else {
                                    let mut new_ri: Vec<usize> = vec![cond_i];
                                    new_ri.extend(ri);
                                    Some((
                                        DataflowPoint::SumInControl(
                                            node,
                                            IncomingPort::from(n_common_ins),
                                            new_ri,
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
                                if o.index() < n_common_ins {
                                    Some((
                                        DataflowPoint::SumInPhControl(
                                            node,
                                            IncomingPort::from(o.index()),
                                            ri.clone(),
                                            *si,
                                        ),
                                        *q,
                                    ))
                                } else {
                                    let mut new_ri: Vec<usize> = vec![cond_i];
                                    new_ri.extend(ri);
                                    Some((
                                        DataflowPoint::SumInPhControl(
                                            node,
                                            IncomingPort::from(n_common_ins),
                                            new_ri,
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
            for (dfp, _) in &body.q_index_map {
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
        if cond_settings.flow_level == DataflowFlowDetail::Invariants {
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
                                jd.tab.delete_qubit(i);
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
                    ControlSpec::SumInControl(IncomingPort::from(n_common_ins), vec![], sum_index);
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
                                    IncomingPort::from(n_common_ins),
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
                                        IncomingPort::from(n_common_ins),
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
        }
    }

    /// Constructs the summary for a loop node.
    /// Assumes cond_settings.flow_level != DataflowFlowDetail::None.
    /// Current implementation does not gather any information for hoisting, just loop invariants
    fn run_tail_loop(
        &mut self,
        hugr: &H,
        node: H::Node,
        _loop_settings: DataflowSettings,
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
        let loop_condition_port = loop_op.just_outputs.len();
        let mut continue_flow = loop_body.clone();
        let mut continue_post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        let mut break_post_selects: Vec<(usize, PauliXZ, bool)> = vec![];
        if let Some(q) = continue_flow
            .q_index_map
            .get_by_left(&DataflowPoint::SumInPhControl(
                body_out,
                IncomingPort::from(loop_condition_port),
                vec![],
                0,
            ))
        {
            break_post_selects.push((*q, PauliXZ::Z, true));
        }
        let mut multiplicity = 0;
        while let Some(q) = continue_flow
            .q_index_map
            .get_by_left(&DataflowPoint::SumInControl(
                body_out,
                IncomingPort::from(loop_condition_port),
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
        for (dfp, q) in continue_flow.q_index_map.clone() {
            let input_dfp: Option<DataflowPoint<H::Node>> = match dfp {
                DataflowPoint::NodeIn(_, in_port, ri) => {
                    // Since we have already called project_non_io, all NodeIns must be from body_out
                    let counterpart = if in_port.index() == loop_condition_port {
                        DataflowPoint::NodeOut(
                            body_in,
                            OutgoingPort::from(loop_condition_port + ri[0]),
                            ri.as_slice()[1..].to_vec(),
                        )
                    } else {
                        DataflowPoint::NodeOut(
                            body_in,
                            OutgoingPort::from(in_port.index()),
                            ri.clone(),
                        )
                    };
                    bells.push((q, counterpart.clone()));
                    Some(counterpart)
                }
                DataflowPoint::SumInPhControl(_, in_port, ri, si) => {
                    // Similarly, the SumInPhControl must be from body_out
                    let counterpart = if in_port.index() == loop_condition_port {
                        DataflowPoint::SumOutPhControl(
                            body_in,
                            OutgoingPort::from(loop_condition_port + ri[0]),
                            ri.as_slice()[1..].to_vec(),
                            si,
                        )
                    } else {
                        DataflowPoint::SumOutPhControl(
                            body_in,
                            OutgoingPort::from(in_port.index()),
                            ri.clone(),
                            si,
                        )
                    };
                    classical_propagations.push((q, counterpart.clone()));
                    Some(counterpart)
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
            new_z.extend_vec(break_flow.tab.z[i].get_boolean_vec(), offset);
            let mut new_x = BitVector::new(offset);
            new_x.extend_vec(break_flow.tab.x[i].get_boolean_vec(), offset);
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
                    let new_dfp = if in_port.index() == loop_condition_port {
                        DataflowPoint::TempOut(
                            node,
                            OutgoingPort::from(loop_condition_port + ri[0]),
                            ri.as_slice()[1..].to_vec(),
                        )
                    } else {
                        DataflowPoint::TempOut(
                            node,
                            OutgoingPort::from(in_port.index()),
                            ri.clone(),
                        )
                    };
                    invariants_q_map.insert(new_dfp, i);
                }
                DataflowPoint::SumInControl(_, in_port, ri, si, m) => {
                    let new_dfp = if in_port.index() == loop_condition_port {
                        DataflowPoint::SumOutControl(
                            node,
                            OutgoingPort::from(loop_condition_port + ri[0]),
                            ri.as_slice()[1..].to_vec(),
                            *si,
                            *m,
                        )
                    } else {
                        DataflowPoint::SumOutControl(
                            node,
                            OutgoingPort::from(in_port.index()),
                            ri.clone(),
                            *si,
                            *m,
                        )
                    };
                    invariants_q_map.insert(new_dfp, i);
                }
                DataflowPoint::SumInPhControl(_, in_port, ri, si) => {
                    let new_dfp = if in_port.index() == loop_condition_port {
                        DataflowPoint::SumOutPhControl(
                            node,
                            OutgoingPort::from(loop_condition_port + ri[0]),
                            ri.as_slice()[1..].to_vec(),
                            *si,
                        )
                    } else {
                        DataflowPoint::SumOutPhControl(
                            node,
                            OutgoingPort::from(in_port.index()),
                            ri.clone(),
                            *si,
                        )
                    };
                    invariants_q_map.insert(new_dfp, i);
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
        StabilizerDataflow {
            tab: invariants_tab,
            q_index_map: invariants_q_map,
        }
    }
}
