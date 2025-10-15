// Going to base this on Mark's python phase folding implementation instead of the hugr dataflow framework which I struggle to see how to adapt to relational values since we can't easily attribute them to individual wires

use crate::bit_vector::BitVector;
use crate::symplectic_tableau::{PauliXZ, SymplecticTableau};
use bimap::BiHashMap;
use hugr::extension::prelude::qb_t;
use hugr::ops::{DataflowOpTrait, OpTag, OpTrait};
use hugr::PortIndex;
use hugr_core::hugr::internal::PortgraphNodeMap;
use hugr_core::ops::OpType;
use hugr_core::{HugrView, IncomingPort, OutgoingPort};
use itertools::{chain, Itertools};
use petgraph::visit as pv;
use std::collections::HashMap;
use std::hash::Hash;
use tket::hugr::extension::simple_op::MakeExtensionOp;
use tket::TketOp;

/// Sets behaviour for function calls in dataflow analysis
#[derive(Clone, PartialEq, Eq)]
pub enum FunctionOpacity {
    /// Function calls are completely opaque and admit no information across them
    Opaque,
    /// Function bodies are analysed but we only pass on the stabilizers over the boundary qubits
    Boundary,
    /// Function bodies are analysed and fully inserted into the parent graph, allowing e.g. phase folding between gates inside the body and around the call site by inlining the function body
    Inline,
}

/// For each qubit in a StabilizerDataflow analysis, decides its role and attaches it to a node if relevant
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DataflowPoint<N: Copy + Eq + Hash> {
    /// Qubit is an OutgoingPort of the unique Input node
    Input(OutgoingPort),
    /// Qubit is an IncomingPort of the unique Output node
    Output(IncomingPort),
    /// Qubit is part of the frontier still being evaluated; recorded as the IncomingPort to the next node
    Frontier(N, IncomingPort),
    /// Qubit is an input to an internal non-Clifford node
    InternalIn(N, IncomingPort),
    /// Qubit is an output of an internal non-Clifford node
    InternalOut(N, OutgoingPort),
    /// Qubit is an input wire within a hierarchical node
    NestedIn(N, OutgoingPort),
    /// Qubit is an output wire within a hierarchical node
    NestedOut(N, IncomingPort),
}

pub struct StabilizerDataflow<H: HugrView> {
    /// Relational dataflow value captured as a set of stabilizer relations on the Choi-state of the circuit skeleton
    pub(crate) tab: SymplecticTableau,
    /// Maps from wires of the program to columns of the tableau. We separately need to track columns for:
    /// - Each input qubit (indexed by OutgoingPorts of the unique Input node)
    /// - Each output qubit (indexed by IncomingPorts of the unique Output node)
    /// - A frontier that moves forward through the program (eventually becoming the output qubits and being removed from here)
    /// - For any internal non-Clifford (or opaque) node, we use columns for each input and output qubit separately; for nodes with stabilizers across them (e.g. Rz has Z_i Z_o), we impose these via projections on the tableau rather than reducing the number of qubits used as this allows every node kind to be handled identically and preventing more tableau management from column elimination
    /// - For any hierarchical node, we use additional columns for each input and output port within their internal representation that we compose to "internal" columns here by projections on the tableau, again so we don't fuss with column elimination
    pub(crate) q_index_map: BiHashMap<DataflowPoint<H::Node>, usize>,
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
    fn new(hugr: &H, parent: H::Node) -> Self {
        let mut q_ind_map: BiHashMap<DataflowPoint<H::Node>, usize> = BiHashMap::default();
        let mut n_in_qubits = 0;
        let inp = hugr
            .children(parent)
            .filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_)))
            .exactly_one()
            .ok()
            .unwrap();
        for (out, out_type) in hugr.out_value_types(inp) {
            if out_type == qb_t() {
                q_ind_map.insert(DataflowPoint::Input(out), 2 * n_in_qubits);
                let (next, next_p) = hugr.single_linked_input(inp, out).unwrap();
                q_ind_map.insert(DataflowPoint::Frontier(next, next_p), 2 * n_in_qubits + 1);
                n_in_qubits += 1;
            }
        }
        let mut tab = SymplecticTableau::new(2 * n_in_qubits);
        // For each input wire, add an identity operation to the tableau
        for i in 0..n_in_qubits {
            let mut ii = BitVector::new(2 * n_in_qubits);
            ii.xor_bit(2 * i);
            ii.xor_bit(2 * i + 1);
            // XX
            tab.add_stab(BitVector::new(2 * n_in_qubits), ii.clone(), false);
            // ZZ
            tab.add_stab(ii, BitVector::new(2 * n_in_qubits), false);
        }
        Self {
            tab,
            q_index_map: q_ind_map,
        }
    }

    fn remove_non_io_qubits(&mut self) {
        // Project out non-IO columns
        let non_ios: Vec<usize> = self
            .q_index_map
            .iter()
            .filter(|(dfp, _)| {
                matches!(dfp, DataflowPoint::Input(_)) || matches!(dfp, DataflowPoint::Output(_))
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
            match moved_qb {
                Some(mq) => {
                    // If the removal caused an IO qubit to have changed index, update it
                    let (removed_dfp, _) = self.q_index_map.remove_by_right(&mq).unwrap();
                    if matches!(removed_dfp, DataflowPoint::Input(_))
                        || matches!(removed_dfp, DataflowPoint::Output(_))
                    {
                        self.q_index_map.insert(removed_dfp, i);
                    }
                }
                None => {
                    self.q_index_map.remove_by_right(&i);
                }
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
                self.q_index_map
                    .remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(i)))
                    .unwrap()
                    .1
            })
            .collect_array()
            .unwrap();
        // Run closure and update frontier with the returned tableau columns
        let out_cols = go(&mut self.tab, &mut self.q_index_map, in_cols);
        for (out_port, out_col) in hugr.node_outputs(node).zip(out_cols) {
            let (next_node, next_port) = hugr.single_linked_input(node, out_port).unwrap();
            self.q_index_map
                .insert(DataflowPoint::Frontier(next_node, next_port), out_col);
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

    fn apply_quantum_gate(&mut self, hugr: &H, node: H::Node, op: TketOp) {
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
                self.apply_op_with(hugr, node, |tab, q_index_map, [col_in0, col_in1]| {
                    let [col_out0, col_out1, col_front0, col_front1] = tab.add_n_qubits();
                    // Add rows for identities col_out0/1--col_front0/1
                    let mut out_front_0 = BitVector::new(tab.nb_qubits);
                    out_front_0.xor_bit(col_out0);
                    out_front_0.xor_bit(col_front0);
                    tab.add_stab(BitVector::new(tab.nb_qubits), out_front_0.clone(), false);
                    tab.add_stab(out_front_0, BitVector::new(tab.nb_qubits), false);
                    let mut out_front_1 = BitVector::new(tab.nb_qubits);
                    out_front_1.xor_bit(col_out1);
                    out_front_1.xor_bit(col_front1);
                    tab.add_stab(BitVector::new(tab.nb_qubits), out_front_1.clone(), false);
                    tab.add_stab(out_front_1, BitVector::new(tab.nb_qubits), false);
                    // Add rows for ZZ over col_in0/1--col_out0/1 and project to commuting
                    let mut in_out_0 = BitVector::new(tab.nb_qubits);
                    in_out_0.xor_bit(col_in0);
                    in_out_0.xor_bit(col_out0);
                    let mut in_out_1 = BitVector::new(tab.nb_qubits);
                    in_out_1.xor_bit(col_in1);
                    in_out_1.xor_bit(col_out1);
                    tab.project_commuting_with(&in_out_0, &BitVector::new(tab.nb_qubits));
                    tab.project_commuting_with(&in_out_1, &BitVector::new(tab.nb_qubits));
                    tab.add_stab(in_out_0, BitVector::new(tab.nb_qubits), false);
                    tab.add_stab(in_out_1, BitVector::new(tab.nb_qubits), false);
                    q_index_map.insert(
                        DataflowPoint::InternalIn(node, IncomingPort::from(0)),
                        col_in0,
                    );
                    q_index_map.insert(
                        DataflowPoint::InternalIn(node, IncomingPort::from(1)),
                        col_in1,
                    );
                    q_index_map.insert(
                        DataflowPoint::InternalOut(node, OutgoingPort::from(0)),
                        col_out0,
                    );
                    q_index_map.insert(
                        DataflowPoint::InternalOut(node, OutgoingPort::from(1)),
                        col_out1,
                    );
                    [col_front0, col_front1]
                });
            }
            TketOp::T | TketOp::Tdg | TketOp::Rz | TketOp::Measure => {
                self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                    let [col_out, col_front] = tab.add_n_qubits();
                    // Add rows for identity col_out--col_front
                    let mut out_front = BitVector::new(tab.nb_qubits);
                    out_front.xor_bit(col_out);
                    out_front.xor_bit(col_front);
                    tab.add_stab(BitVector::new(tab.nb_qubits), out_front.clone(), false);
                    tab.add_stab(out_front, BitVector::new(tab.nb_qubits), false);
                    // Add row for ZZ over col_in--col_out and project to commuting
                    let mut in_out = BitVector::new(tab.nb_qubits);
                    in_out.xor_bit(col_in);
                    in_out.xor_bit(col_out);
                    tab.project_commuting_with(&in_out, &BitVector::new(tab.nb_qubits));
                    tab.add_stab(in_out, BitVector::new(tab.nb_qubits), false);
                    q_index_map.insert(
                        DataflowPoint::InternalIn(node, IncomingPort::from(0)),
                        col_in,
                    );
                    q_index_map.insert(
                        DataflowPoint::InternalOut(node, OutgoingPort::from(0)),
                        col_out,
                    );
                    [col_front]
                });
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
                self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                    let [col_out, col_front] = tab.add_n_qubits();
                    // Add rows for identity col_out--col_front
                    let mut out_front = BitVector::new(tab.nb_qubits);
                    out_front.xor_bit(col_out);
                    out_front.xor_bit(col_front);
                    tab.add_stab(BitVector::new(tab.nb_qubits), out_front.clone(), false);
                    tab.add_stab(out_front, BitVector::new(tab.nb_qubits), false);
                    // Add row for XX over col_in--col_out and project to commuting
                    let mut in_out = BitVector::new(tab.nb_qubits);
                    in_out.xor_bit(col_in);
                    in_out.xor_bit(col_out);
                    tab.project_commuting_with(&BitVector::new(tab.nb_qubits), &in_out);
                    tab.add_stab(BitVector::new(tab.nb_qubits), in_out, false);
                    q_index_map.insert(
                        DataflowPoint::InternalIn(node, IncomingPort::from(0)),
                        col_in,
                    );
                    q_index_map.insert(
                        DataflowPoint::InternalOut(node, OutgoingPort::from(0)),
                        col_out,
                    );
                    [col_front]
                });
            }
            TketOp::Ry => {
                self.apply_op_with(hugr, node, |tab, q_index_map, [col_in]| {
                    let [col_out, col_front] = tab.add_n_qubits();
                    // Add rows for identity col_out--col_front
                    let mut out_front = BitVector::new(tab.nb_qubits);
                    out_front.xor_bit(col_out);
                    out_front.xor_bit(col_front);
                    tab.add_stab(BitVector::new(tab.nb_qubits), out_front.clone(), false);
                    tab.add_stab(out_front, BitVector::new(tab.nb_qubits), false);
                    // Add row for -YY (negative because of the partial transpose under the Choi isomorphism) over col_in--col_out and project to commuting
                    let mut in_out = BitVector::new(tab.nb_qubits);
                    in_out.xor_bit(col_in);
                    in_out.xor_bit(col_out);
                    tab.project_commuting_with(&in_out, &in_out);
                    tab.add_stab(in_out.clone(), in_out, true);
                    q_index_map.insert(
                        DataflowPoint::InternalIn(node, IncomingPort::from(0)),
                        col_in,
                    );
                    q_index_map.insert(
                        DataflowPoint::InternalOut(node, OutgoingPort::from(0)),
                        col_out,
                    );
                    [col_front]
                });
            }
            TketOp::Toffoli => {
                self.apply_op_with(
                    hugr,
                    node,
                    |tab, q_index_map, [col_in0, col_in1, col_in2]| {
                        let [col_out0, col_out1, col_out2, col_front0, col_front1, col_front2] =
                            tab.add_n_qubits();
                        // Add rows for identities col_out0/1/2--col_front0/1/2
                        let mut out_front_0 = BitVector::new(tab.nb_qubits);
                        out_front_0.xor_bit(col_out0);
                        out_front_0.xor_bit(col_front0);
                        tab.add_stab(BitVector::new(tab.nb_qubits), out_front_0.clone(), false);
                        tab.add_stab(out_front_0, BitVector::new(tab.nb_qubits), false);
                        let mut out_front_1 = BitVector::new(tab.nb_qubits);
                        out_front_1.xor_bit(col_out1);
                        out_front_1.xor_bit(col_front1);
                        tab.add_stab(BitVector::new(tab.nb_qubits), out_front_1.clone(), false);
                        tab.add_stab(out_front_1, BitVector::new(tab.nb_qubits), false);
                        let mut out_front_2 = BitVector::new(tab.nb_qubits);
                        out_front_2.xor_bit(col_out2);
                        out_front_2.xor_bit(col_front2);
                        tab.add_stab(BitVector::new(tab.nb_qubits), out_front_2.clone(), false);
                        tab.add_stab(out_front_2, BitVector::new(tab.nb_qubits), false);
                        // Add rows for ZZ/ZZ/XX over col_in0/1/2--col_out0/1/2 and project to commuting
                        let mut in_out_0 = BitVector::new(tab.nb_qubits);
                        in_out_0.xor_bit(col_in0);
                        in_out_0.xor_bit(col_out0);
                        let mut in_out_1 = BitVector::new(tab.nb_qubits);
                        in_out_1.xor_bit(col_in1);
                        in_out_1.xor_bit(col_out1);
                        let mut in_out_2 = BitVector::new(tab.nb_qubits);
                        in_out_2.xor_bit(col_in2);
                        in_out_2.xor_bit(col_out2);
                        tab.project_commuting_with(&in_out_0, &BitVector::new(tab.nb_qubits));
                        tab.project_commuting_with(&in_out_1, &BitVector::new(tab.nb_qubits));
                        tab.project_commuting_with(&BitVector::new(tab.nb_qubits), &in_out_2);
                        tab.add_stab(in_out_0, BitVector::new(tab.nb_qubits), false);
                        tab.add_stab(in_out_1, BitVector::new(tab.nb_qubits), false);
                        tab.add_stab(BitVector::new(tab.nb_qubits), in_out_2, false);
                        q_index_map.insert(
                            DataflowPoint::InternalIn(node, IncomingPort::from(0)),
                            col_in0,
                        );
                        q_index_map.insert(
                            DataflowPoint::InternalIn(node, IncomingPort::from(1)),
                            col_in1,
                        );
                        q_index_map.insert(
                            DataflowPoint::InternalIn(node, IncomingPort::from(2)),
                            col_in2,
                        );
                        q_index_map.insert(
                            DataflowPoint::InternalOut(node, OutgoingPort::from(0)),
                            col_out0,
                        );
                        q_index_map.insert(
                            DataflowPoint::InternalOut(node, OutgoingPort::from(1)),
                            col_out1,
                        );
                        q_index_map.insert(
                            DataflowPoint::InternalOut(node, OutgoingPort::from(2)),
                            col_out2,
                        );
                        [col_front0, col_front1, col_front2]
                    },
                );
            }
            TketOp::MeasureFree => {
                self.apply_op_with(hugr, node, |_, q_index_map, [col_in]| {
                    q_index_map.insert(
                        DataflowPoint::InternalIn(node, IncomingPort::from(0)),
                        col_in,
                    );
                    []
                });
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
            _ => {
                // Only other remaining TketOp option at time of writing is TryQAlloc which has no qubits in its signature (the output is a Sum and therefore we currently don't track any relations involving it)
                // In case other options are added later on, handle them as opaque unless we explicitly add a custom handler for them
                self.apply_opaque(hugr, node)
            }
        }
    }

    fn apply_opaque(&mut self, hugr: &H, node: H::Node) {
        // For each Qubit input, move the column from frontier_cols to internal_in_cols
        for (p, t) in hugr.in_value_types(node) {
            if t == qb_t() {
                let (_, col) = self
                    .q_index_map
                    .remove_by_left(&DataflowPoint::Frontier(node, p))
                    .unwrap();
                self.q_index_map
                    .insert(DataflowPoint::InternalIn(node, p), col);
            }
        }
        // For each Qubit output, create a pair of columns with the identity for internal_out_cols and frontier_cols
        for (p, t) in hugr.out_value_types(node) {
            if t == qb_t() {
                let [col_out, col_front] = self.tab.add_n_qubits();
                // Add rows for identity col_out--col_front
                let mut out_front = BitVector::new(self.tab.nb_qubits);
                out_front.xor_bit(col_out);
                out_front.xor_bit(col_front);
                self.tab
                    .add_stab(BitVector::new(self.tab.nb_qubits), out_front.clone(), false);
                self.tab
                    .add_stab(out_front, BitVector::new(self.tab.nb_qubits), false);
                self.q_index_map
                    .insert(DataflowPoint::InternalOut(node, p), col_out);
                let (next_node, next_port) = hugr.single_linked_input(node, p).unwrap();
                self.q_index_map
                    .insert(DataflowPoint::Frontier(next_node, next_port), col_front);
            }
        }
    }

    /// Suppose we have already recursively calculated a StabilizerDataflow for node, factoring in any computation from the node itself (e.g. Kleene closure for TailLoop or projections for function calls); performs sequential composition to append it to the appropriate qubits here
    fn apply_summary(&mut self, hugr: &H, node: H::Node, node_summary: &StabilizerDataflow<H>) {
        let old_n_qbs = self.tab.nb_qubits;
        let n_added_qbs = node_summary.tab.nb_qubits;
        self.tab.add_qubits(n_added_qbs);
        for (dfp, col) in node_summary.q_index_map.iter() {
            match dfp {
                DataflowPoint::Input(p) => {
                    self.q_index_map
                        .insert(DataflowPoint::NestedIn(node, *p), col + old_n_qbs);
                }
                DataflowPoint::Output(p) => {
                    self.q_index_map
                        .insert(DataflowPoint::NestedOut(node, *p), col + old_n_qbs);
                }
                _ => {
                    // It is decided at analysis construction whether or not we project to IO or keep internals, so always copy any internals that have remained (e.g. for inlining function calls)
                    self.q_index_map.insert(dfp.clone(), col + old_n_qbs);
                }
            }
        }
        for i in 0..node_summary.tab.nb_stabs {
            let mut new_z = BitVector::new(old_n_qbs);
            new_z.extend_vec(
                node_summary.tab.z[i].get_sized_boolean_vec(n_added_qbs),
                old_n_qbs,
            );
            let mut new_x = BitVector::new(old_n_qbs);
            new_x.extend_vec(
                node_summary.tab.x[i].get_sized_boolean_vec(n_added_qbs),
                old_n_qbs,
            );
            self.tab
                .add_stab(new_z, new_x, node_summary.tab.signs.get(i));
        }
        for (port, in_type) in hugr.in_value_types(node) {
            if in_type == qb_t() {
                let out_port = OutgoingPort::from(port.index());
                let (_, internal_col) = self
                    .q_index_map
                    .remove_by_left(&DataflowPoint::Frontier(node, port))
                    .unwrap();
                self.q_index_map
                    .insert(DataflowPoint::InternalIn(node, port), internal_col);
                let nested_col = self
                    .q_index_map
                    .get_by_left(&DataflowPoint::NestedIn(node, out_port))
                    .unwrap();
                // Project ZZ and XX to compose nested_col and internal_col
                let mut nested_internal = BitVector::new(self.tab.nb_qubits);
                nested_internal.xor_bit(*nested_col);
                nested_internal.xor_bit(internal_col);
                self.tab
                    .project_commuting_with(&BitVector::new(self.tab.nb_qubits), &nested_internal);
                self.tab
                    .project_commuting_with(&nested_internal, &BitVector::new(self.tab.nb_qubits));
                self.tab.add_stab(
                    BitVector::new(self.tab.nb_qubits),
                    nested_internal.clone(),
                    false,
                );
                self.tab
                    .add_stab(nested_internal, BitVector::new(self.tab.nb_qubits), false);
            }
        }
        for (port, out_type) in hugr.out_value_types(node) {
            if out_type == qb_t() {
                let in_port = IncomingPort::from(port.index());
                let nested_col = *self
                    .q_index_map
                    .get_by_left(&DataflowPoint::NestedOut(node, in_port))
                    .unwrap();
                let [internal_col, front_col] = self.tab.add_n_qubits();
                self.q_index_map
                    .insert(DataflowPoint::InternalOut(node, port), internal_col);
                let (next_node, next_port) = hugr.single_linked_input(node, port).unwrap();
                self.q_index_map
                    .insert(DataflowPoint::Frontier(next_node, next_port), front_col);
                // Add rows for identity internal_col--front_col
                let mut internal_front = BitVector::new(self.tab.nb_qubits);
                internal_front.xor_bit(internal_col);
                internal_front.xor_bit(front_col);
                self.tab.add_stab(
                    BitVector::new(self.tab.nb_qubits),
                    internal_front.clone(),
                    false,
                );
                self.tab
                    .add_stab(internal_front, BitVector::new(self.tab.nb_qubits), false);
                // Project ZZ and XX to compose nested_col and internal_col
                let mut nested_internal = BitVector::new(self.tab.nb_qubits);
                nested_internal.xor_bit(nested_col);
                nested_internal.xor_bit(internal_col);
                self.tab
                    .project_commuting_with(&BitVector::new(self.tab.nb_qubits), &nested_internal);
                self.tab
                    .project_commuting_with(&nested_internal, &BitVector::new(self.tab.nb_qubits));
                self.tab.add_stab(
                    BitVector::new(self.tab.nb_qubits),
                    nested_internal.clone(),
                    false,
                );
                self.tab
                    .add_stab(nested_internal, BitVector::new(self.tab.nb_qubits), false);
            }
        }
    }
}

// For any control-flow region or hierarchical node, store the analysis for its internal calculations
// For TailLoop and function calls, this is the summary of the body and not of the invariants or projecting away interior information
pub struct SDFAnalysis<H: HugrView>(pub(crate) HashMap<H::Node, StabilizerDataflow<H>>);

impl<H: HugrView> SDFAnalysis<H> {
    pub fn run_hugr(hugr: &H, fun_op: &FunctionOpacity) -> SDFAnalysis<H> {
        let mut res = SDFAnalysis(HashMap::default());
        for n in hugr.nodes() {
            if OpTag::DataflowParent.is_superset(hugr.get_optype(n).tag())
                && !res.0.contains_key(&n)
            {
                let summary = res.run_dfg(hugr, n, fun_op);
                res.0.insert(n, summary);
            }
        }
        res
    }

    fn run_dfg(
        &mut self,
        hugr: &H,
        parent: H::Node,
        fun_op: &FunctionOpacity,
    ) -> StabilizerDataflow<H> {
        let mut summary = StabilizerDataflow::new(hugr, parent);
        let (region, node_map) = hugr.region_portgraph(parent);
        let mut topo = pv::Topo::new(&region);
        while let Some(pgnode) = topo.next(&region) {
            let node = node_map.from_portgraph(pgnode);
            let optype: &OpType = hugr.get_optype(node);
            match optype {
                OpType::ExtensionOp(op) => match TketOp::from_extension_op(op) {
                    Ok(tkop) => summary.apply_quantum_gate(hugr, node, tkop),
                    Err(_) => summary.apply_opaque(hugr, node),
                },
                OpType::Conditional(_) => {
                    if self.0.contains_key(&node) {
                        summary.apply_summary(hugr, node, self.0.get(&node).unwrap());
                    } else {
                        let cond_summary = self.run_conditional(hugr, node, fun_op);
                        summary.apply_summary(hugr, node, &cond_summary);
                        self.0.insert(node, cond_summary);
                    }
                }
                OpType::TailLoop(_) => {
                    if !self.0.contains_key(&node) {
                        let body_summary = self.run_dfg(hugr, node, fun_op);
                        self.0.insert(node, body_summary);
                    }
                    let loop_summary = self.run_tail_loop(hugr, node, fun_op);
                    summary.apply_summary(hugr, node, &loop_summary);
                }
                OpType::Call(_) => match *fun_op {
                    FunctionOpacity::Opaque => {
                        summary.apply_opaque(hugr, node);
                    }
                    FunctionOpacity::Boundary => {
                        let call_port = optype.static_input_port().unwrap();
                        let (fun_def_node, _) = hugr
                            .linked_outputs(node, call_port)
                            .exactly_one()
                            .ok()
                            .unwrap();
                        if !self.0.contains_key(&fun_def_node) {
                            let body_summary = self.run_dfg(hugr, fun_def_node, fun_op);
                            self.0.insert(fun_def_node, body_summary);
                        }
                        let mut fun_summary = (*self.0.get(&fun_def_node).unwrap()).clone();
                        fun_summary.remove_non_io_qubits();
                        summary.apply_summary(hugr, node, &fun_summary);
                    }
                    FunctionOpacity::Inline => {
                        let call_port = optype.static_input_port().unwrap();
                        let (fun_def_node, _) = hugr
                            .linked_outputs(node, call_port)
                            .exactly_one()
                            .ok()
                            .unwrap();
                        if self.0.contains_key(&fun_def_node) {
                            let body_summary = self.0.get(&fun_def_node).unwrap();
                            summary.apply_summary(hugr, fun_def_node, body_summary);
                        } else {
                            let body_summary = self.run_dfg(hugr, fun_def_node, fun_op);
                            summary.apply_summary(hugr, node, &body_summary);
                            self.0.insert(fun_def_node, body_summary);
                        }
                    }
                },
                OpType::Input(_) => {
                    // Already handled during setup
                }
                OpType::Output(_) => {
                    // Frontier finished, move it to out_cols
                    let to_move: Vec<(DataflowPoint<H::Node>, usize)> = summary
                        .q_index_map
                        .iter()
                        .filter(|(dfp, _)| matches!(dfp, DataflowPoint::Frontier(_, _)))
                        .map(|(dfp, q)| (dfp.clone(), *q))
                        .collect_vec();
                    for (dfp, q) in to_move {
                        if let DataflowPoint::Frontier(_node, port) = dfp {
                            summary.q_index_map.remove_by_right(&q);
                            summary.q_index_map.insert(DataflowPoint::Output(port), q);
                        }
                    }
                }
                _ => summary.apply_opaque(hugr, node),
            }
        }
        summary
    }

    fn run_conditional(
        &mut self,
        hugr: &H,
        node: H::Node,
        fun_op: &FunctionOpacity,
    ) -> StabilizerDataflow<H> {
        // Assume no information is passed about Qubits within the Sum types, so our summary only incorporates the Qubits in the other args
        let cond = hugr.get_optype(node).as_conditional().unwrap();
        let sig = cond.signature();
        // Determins consistent column indexing for inputs and outputs
        let mut unified_q_index: BiHashMap<DataflowPoint<H::Node>, usize> = BiHashMap::default();
        let mut n_unified_qbs = 0;
        for in_port in sig.input_ports() {
            if *sig.in_port_type(in_port).unwrap() == qb_t() {
                unified_q_index.insert(
                    DataflowPoint::Input(OutgoingPort::from(in_port.index())),
                    n_unified_qbs,
                );
                n_unified_qbs += 1;
            }
        }
        for out_port in sig.output_ports() {
            if *sig.out_port_type(out_port).unwrap() == qb_t() {
                unified_q_index.insert(
                    DataflowPoint::Output(IncomingPort::from(out_port.index())),
                    n_unified_qbs,
                );
                n_unified_qbs += 1;
            }
        }
        let mut summary: Option<StabilizerDataflow<H>> = None;
        for (cond_i, cond_node) in hugr.children(node).enumerate() {
            if !self.0.contains_key(&cond_node) {
                let cond_summary = self.run_dfg(hugr, cond_node, fun_op);
                self.0.insert(cond_node, cond_summary);
            }
            let cond_summary = self.0.get(&cond_node).unwrap();
            let mut projected_tab = cond_summary.tab.clone();
            // Number of ports from the condition row; given port p on input, corresponds to IncomingPort::from(p + 1 - cond_len) to the Conditional
            let cond_len = cond.sum_rows.get(cond_i).unwrap().len();
            // Project out non-IO columns (including any qubits from the condition row)
            let non_ios: Vec<usize> = cond_summary
                .q_index_map
                .iter()
                .filter(|(dfp, _)| match dfp {
                    DataflowPoint::Input(port) => port.index() < cond_len,
                    DataflowPoint::Output(_) => false,
                    _ => true,
                })
                .map(|(_, q)| *q)
                .collect_vec();
            let project_cols: Vec<(usize, PauliXZ)> = chain!(
                non_ios.iter().map(|i| (*i, PauliXZ::X)),
                non_ios.iter().map(|i| (*i, PauliXZ::Z)),
            )
            .collect_vec();
            projected_tab.project_cols_to_zero(&project_cols);
            // Rebuild projected_tab with the column order given by unified_X_cols
            let mut unified_order_tab = SymplecticTableau::new(n_unified_qbs);
            for i in 0..projected_tab.nb_stabs {
                let mut z = BitVector::new(n_unified_qbs);
                let mut x = BitVector::new(n_unified_qbs);
                for (dfp, col) in &unified_q_index {
                    match dfp {
                        DataflowPoint::Input(port) => {
                            let old_col = cond_summary
                                .q_index_map
                                .get_by_left(&DataflowPoint::Input(OutgoingPort::from(
                                    port.index() + cond_len - 1,
                                )))
                                .unwrap();
                            if projected_tab.z[i].get(*old_col) {
                                z.xor_bit(*col);
                            }
                            if projected_tab.x[i].get(*old_col) {
                                x.xor_bit(*col);
                            }
                        }
                        DataflowPoint::Output(port) => {
                            let old_col = cond_summary
                                .q_index_map
                                .get_by_left(&DataflowPoint::Output(*port))
                                .unwrap();
                            if projected_tab.z[i].get(*old_col) {
                                z.xor_bit(*col);
                            }
                            if projected_tab.x[i].get(*old_col) {
                                x.xor_bit(*col);
                            }
                        }
                        _ => {}
                    }
                }
                unified_order_tab.add_stab(z, x, projected_tab.signs.get(i));
            }
            // Build summary
            match summary {
                Some(ref mut summ) => {
                    // Compute join of unified_order_tab and summ.tab
                    summ.tab = SymplecticTableau::join(&unified_order_tab, &summ.tab);
                }
                None => {
                    summary = Some(StabilizerDataflow {
                        tab: unified_order_tab,
                        q_index_map: unified_q_index.clone(),
                    });
                }
            }
        }
        summary.unwrap()
    }

    fn run_tail_loop(
        &mut self,
        hugr: &H,
        node: H::Node,
        fun_op: &FunctionOpacity,
    ) -> StabilizerDataflow<H> {
        // The output of the loop body includes a Sum[just_inputs, just_outputs] to dictate whether to loop; region-based analysis would scale exponentially in the number of Sum types, so we assume any qubits included in it have been projected out and therefore we have no information about qubits in just_outputs, or about qubits in just_inputs passed to the next iteration
        // Tail loops are run at least once; despite this we still reach a fixpoint with a single join:
        // Suppose PCQ and PCCQ; the latter suggests there is some R s.t. PCR and RCQ; combining, we get CQR, PRC, RCR; then for any number of iterations we can go PCR,RCR,RCR,...,RCQ
        // Even if we interpose them with some black-box initialisation B over [just_inputs], PCBCQ still implies a split by R which commutes with B (i.e. is on disjoint qubits) so we can still identify that R is an invariant of BC
        if !self.0.contains_key(&node) {
            let body_summary = self.run_dfg(hugr, node, fun_op);
            self.0.insert(node, body_summary);
        }
        let mut body_summary = self.0.get(&node).unwrap().clone();
        let tl = hugr.get_optype(node).as_tail_loop().unwrap();
        let mut summary = StabilizerDataflow {
            tab: SymplecticTableau::new(0),
            q_index_map: BiHashMap::default(),
        };
        // Build q_index_map for target qubit structure
        for (in_port, in_type) in tl.just_inputs.iter().enumerate() {
            if *in_type == qb_t() {
                let new_col = summary.tab.add_qubit();
                summary
                    .q_index_map
                    .insert(DataflowPoint::Input(OutgoingPort::from(in_port)), new_col);
            }
        }
        for (out_port, out_type) in tl.just_outputs.iter().enumerate() {
            if *out_type == qb_t() {
                let new_col = summary.tab.add_qubit();
                summary
                    .q_index_map
                    .insert(DataflowPoint::Output(IncomingPort::from(out_port)), new_col);
            }
        }
        for (io_port, io_type) in tl.rest.iter().enumerate() {
            if *io_type == qb_t() {
                let [first_col, second_col] = summary.tab.add_n_qubits();
                summary.q_index_map.insert(
                    DataflowPoint::Input(OutgoingPort::from(io_port + tl.just_inputs.len())),
                    first_col,
                );
                summary.q_index_map.insert(
                    DataflowPoint::Output(IncomingPort::from(io_port + tl.just_outputs.len())),
                    second_col,
                );
            }
        }
        // Project body_summary.tab and reorder to match target q_index_map
        let non_ios: Vec<usize> = body_summary
            .q_index_map
            .iter()
            .filter(|(dfp, _)| {
                !(matches!(dfp, DataflowPoint::Input(_)) || matches!(dfp, DataflowPoint::Output(_)))
            })
            .map(|(_, q)| *q)
            .collect_vec();
        let project_cols: Vec<(usize, PauliXZ)> = chain!(
            non_ios.iter().map(|i| (*i, PauliXZ::X)),
            non_ios.iter().map(|i| (*i, PauliXZ::Z))
        )
        .collect_vec();
        body_summary.tab.project_cols_to_zero(&project_cols);
        for i in 0..body_summary.tab.nb_stabs {
            let mut z = BitVector::new(summary.q_index_map.len());
            let mut x = BitVector::new(summary.q_index_map.len());
            for (dfp, new_col) in summary.q_index_map.iter() {
                match dfp {
                    DataflowPoint::Input(_) => {
                        // Port indexing inside and outside the loop match at the inputs
                        let old_col = body_summary.q_index_map.get_by_left(dfp).unwrap();
                        if body_summary.tab.z[i].get(*old_col) {
                            z.xor_bit(*new_col)
                        }
                        if body_summary.tab.x[i].get(*old_col) {
                            x.xor_bit(*new_col)
                        }
                    }
                    DataflowPoint::Output(port) => {
                        if port.index() >= tl.just_outputs.len() {
                            let old_col = body_summary
                                .q_index_map
                                .get_by_left(&DataflowPoint::Output(IncomingPort::from(
                                    port.index() + 1 - tl.just_outputs.len(),
                                )))
                                .unwrap();
                            if body_summary.tab.z[i].get(*old_col) {
                                z.xor_bit(*new_col)
                            }
                            if body_summary.tab.x[i].get(*old_col) {
                                x.xor_bit(*new_col)
                            }
                        }
                    }
                    _ => {}
                }
            }
            summary.tab.add_stab(z, x, body_summary.tab.signs.get(i));
        }
        // Compose body_summary.tab with itself and reorder to match target q_index_map
        let mut iter_2_tab = summary.tab.clone();
        iter_2_tab.add_qubits(iter_2_tab.nb_qubits);
        for i in 0..summary.tab.nb_stabs {
            let mut z = BitVector::new(summary.tab.nb_qubits);
            let mut x = BitVector::new(summary.tab.nb_qubits);
            z.extend_vec(
                summary.tab.z[i].get_sized_boolean_vec(summary.tab.nb_qubits),
                summary.tab.nb_qubits,
            );
            x.extend_vec(
                summary.tab.x[i].get_sized_boolean_vec(summary.tab.nb_qubits),
                summary.tab.nb_qubits,
            );
            iter_2_tab.add_stab(z, x, summary.tab.signs.get(i));
        }
        for (rest_index, rest_type) in tl.rest.iter().enumerate() {
            if *rest_type == qb_t() {
                let iter1_out = summary
                    .q_index_map
                    .get_by_left(&DataflowPoint::Output(IncomingPort::from(
                        tl.just_outputs.len() + rest_index,
                    )))
                    .unwrap();
                let iter2_in = summary.tab.nb_qubits
                    + summary
                        .q_index_map
                        .get_by_left(&DataflowPoint::Input(OutgoingPort::from(
                            tl.just_inputs.len() + rest_index,
                        )))
                        .unwrap();
                // Project ZZ and XX to compose
                let mut both_qubits = BitVector::new(iter_2_tab.nb_qubits);
                both_qubits.xor_bit(*iter1_out);
                both_qubits.xor_bit(iter2_in);
                iter_2_tab
                    .project_commuting_with(&BitVector::new(iter_2_tab.nb_qubits), &both_qubits);
                iter_2_tab
                    .project_commuting_with(&both_qubits, &BitVector::new(iter_2_tab.nb_qubits));
                iter_2_tab.add_stab(
                    BitVector::new(iter_2_tab.nb_qubits),
                    both_qubits.clone(),
                    false,
                );
                iter_2_tab.add_stab(both_qubits, BitVector::new(iter_2_tab.nb_qubits), false);
                // Swap qubits to get all final qubits in their intended position
                let iter2_out = iter2_in + 1;
                iter_2_tab.append_swap(*iter1_out, iter2_out);
            }
        }
        let iter_2_project_cols: Vec<(usize, PauliXZ)> = chain!(
            (summary.tab.nb_qubits..iter_2_tab.nb_qubits).map(|i| (i, PauliXZ::X)),
            (summary.tab.nb_qubits..iter_2_tab.nb_qubits).map(|i| (i, PauliXZ::Z))
        )
        .collect_vec();
        iter_2_tab.project_cols_to_zero(&iter_2_project_cols);
        for q in (summary.tab.nb_qubits..iter_2_tab.nb_qubits).rev() {
            // This removes the joined qubits and any of just_inputs from iteration 2
            // Since no information is obtained for just_outputs, it doesn't matter which copy of the qubits we remove
            iter_2_tab.delete_qubit(q);
        }
        // Take join
        summary.tab = SymplecticTableau::join(&summary.tab, &iter_2_tab);
        summary
    }
}

#[cfg(test)]
mod test {
    use hugr::{
        builder::{
            endo_sig, Container, Dataflow, DataflowHugr, DataflowSubContainer, FunctionBuilder,
            HugrBuilder, SubContainer,
        },
        extension::{
            prelude::{bool_t, qb_t, usize_t},
            Version,
        },
        ops::{handle::NodeHandle, Value},
        type_row,
        types::Signature,
        Extension, HugrView, IncomingPort, OutgoingPort,
    };
    use tket::{extension::rotation::ConstRotation, TketOp};

    use crate::stabilizer_dataflow::{DataflowPoint, FunctionOpacity, SDFAnalysis};

    #[test]
    fn test_empty_analysis() {
        let builder = FunctionBuilder::new("empty", endo_sig(vec![])).unwrap();
        let hugr = builder.finish_hugr().unwrap();
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap();
        assert_eq!(summary.tab.nb_qubits, 0);
        assert_eq!(summary.tab.nb_stabs, 0);
    }

    #[test]
    fn test_identity_analysis() {
        // Add an extra integer input to make sure we only track the qubits
        let builder =
            FunctionBuilder::new("identity", endo_sig(vec![usize_t(), qb_t(), qb_t()])).unwrap();
        let [i, qb0, qb1] = builder.input_wires_arr();
        let hugr = builder.finish_hugr_with_outputs([i, qb0, qb1]).unwrap();
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 4);
        assert_eq!(summary.tab.nb_stabs, 4);
        // Check the right ports are stored for tracking the qubits
        assert_eq!(summary.q_index_map.len(), 4);
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::Output(IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::Output(IncomingPort::from(2))
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Check that the rows correspond to the identity operations
        // Note that BitVector assigns index 0 to the least significant bit and always adds an extra block than needed
        assert_eq!(summary.tab.x[0].get_integer_vec()[0], 0b0011i128);
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(0), false);
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.z[1].get_integer_vec()[0], 0b0011i128);
        assert_eq!(summary.tab.signs.get(1), false);
        assert_eq!(summary.tab.x[2].get_integer_vec()[0], 0b1100i128);
        assert_eq!(summary.tab.z[2].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(2), false);
        assert_eq!(summary.tab.x[3].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.z[3].get_integer_vec()[0], 0b1100i128);
        assert_eq!(summary.tab.signs.get(3), false);
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
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 2);
        assert_eq!(summary.tab.nb_stabs, 2);
        // Check that the rows correspond to the Bell state stabilizers
        summary.tab.echelon(&summary.tab.all_columns());
        assert_eq!(summary.tab.x[0].get_integer_vec()[0], 0b11i128);
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0b00i128);
        assert_eq!(summary.tab.signs.get(0), false);
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0b00i128);
        assert_eq!(summary.tab.z[1].get_integer_vec()[0], 0b11i128);
        assert_eq!(summary.tab.signs.get(1), false);
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
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 4);
        assert_eq!(summary.tab.nb_stabs, 4);
        // Reduce summary.tab to row echelon form with qubit ordering [op_in, out0, op_out, out1] (because the second QAlloc occurs first in the topological sort)
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::InternalIn(opaque_op.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::Output(IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::InternalOut(opaque_op.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::Output(IncomingPort::from(1))
        );
        // Check the rows
        summary.tab.echelon(&summary.tab.all_columns());
        // Xop_in Xout0 Xout1
        // Zop_in Xop_out Zout1
        // Zout0 Xop_out Zout1
        // Zop_out Xout1
        assert_eq!(summary.tab.x[0].get_integer_vec()[0], 0b1011i128);
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(0), false);
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0b0100i128);
        assert_eq!(summary.tab.z[1].get_integer_vec()[0], 0b1001i128);
        assert_eq!(summary.tab.signs.get(1), false);
        assert_eq!(summary.tab.x[2].get_integer_vec()[0], 0b0100i128);
        assert_eq!(summary.tab.z[2].get_integer_vec()[0], 0b1010i128);
        assert_eq!(summary.tab.signs.get(2), false);
        assert_eq!(summary.tab.x[3].get_integer_vec()[0], 0b1000i128);
        assert_eq!(summary.tab.z[3].get_integer_vec()[0], 0b0100i128);
        assert_eq!(summary.tab.signs.get(3), false);
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
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 6);
        assert_eq!(summary.tab.nb_stabs, 6);
        // Reduce summary.tab to row echelon form with qubit ordering [in0, out0, in1, out1, in2, out2]
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::Output(IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::Output(IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::Output(IncomingPort::from(2))
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Check the rows
        assert_eq!(summary.tab.x[0].get_integer_vec()[0], 0b000011i128);
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(0), true);
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.z[1].get_integer_vec()[0], 0b000011i128);
        assert_eq!(summary.tab.signs.get(1), false);
        assert_eq!(summary.tab.x[2].get_integer_vec()[0], 0b001100i128);
        assert_eq!(summary.tab.z[2].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(2), false);
        assert_eq!(summary.tab.x[3].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.z[3].get_integer_vec()[0], 0b001100i128);
        assert_eq!(summary.tab.signs.get(3), true);
        assert_eq!(summary.tab.x[4].get_integer_vec()[0], 0b110000i128);
        assert_eq!(summary.tab.z[4].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(4), true);
        assert_eq!(summary.tab.x[5].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.z[5].get_integer_vec()[0], 0b110000i128);
        assert_eq!(summary.tab.signs.get(5), true);
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
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 28);
        assert_eq!(summary.tab.nb_stabs, 28);
        // Reduce summary.tab to row echelon form with qubit ordering:
        // [in0, t.in, in1, ry.in, in2, rx.in, rx.out, toffoli.in2, ry.out, crz.in1, t.out, tdg.in, tdg.out, rz.in,
        // rz.out, meas.in, meas.out, crz.in0, crz.out0, crz.out1, toffoli.in0, toffoli.in1, toffoli.out0,
        // toffoli.out1, toffoli.out2, out0, out1, out2]
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::InternalIn(t.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::InternalIn(ry.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::InternalIn(rx.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::InternalOut(rx.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::InternalIn(toffoli.node(), IncomingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::InternalOut(ry.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::InternalIn(crz.node(), IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::InternalOut(t.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::InternalIn(tdg.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::InternalOut(tdg.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::InternalIn(rz.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::InternalOut(rz.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::InternalIn(meas.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&16).unwrap(),
            DataflowPoint::InternalOut(meas.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&17).unwrap(),
            DataflowPoint::InternalIn(crz.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&18).unwrap(),
            DataflowPoint::InternalOut(crz.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&19).unwrap(),
            DataflowPoint::InternalOut(crz.node(), OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&20).unwrap(),
            DataflowPoint::InternalIn(toffoli.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&21).unwrap(),
            DataflowPoint::InternalIn(toffoli.node(), IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&22).unwrap(),
            DataflowPoint::InternalOut(toffoli.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&23).unwrap(),
            DataflowPoint::InternalOut(toffoli.node(), OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&24).unwrap(),
            DataflowPoint::InternalOut(toffoli.node(), OutgoingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&25).unwrap(),
            DataflowPoint::Output(IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&26).unwrap(),
            DataflowPoint::Output(IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&27).unwrap(),
            DataflowPoint::Output(IncomingPort::from(2))
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Check the rows
        // Xin0
        assert_eq!(
            summary.tab.x[0].get_integer_vec()[0],
            0b0010010101111111110000000011i128
        );
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(0), false);
        // Zin0
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[1].get_integer_vec()[0],
            0b0010000000000000000000000001i128
        );
        assert_eq!(summary.tab.signs.get(1), false);
        // Zt.in
        assert_eq!(summary.tab.x[2].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[2].get_integer_vec()[0],
            0b0010000000000000000000000010i128
        );
        assert_eq!(summary.tab.signs.get(2), false);
        // Xin1
        assert_eq!(
            summary.tab.x[3].get_integer_vec()[0],
            0b0100101010000000001000000100i128
        );
        assert_eq!(
            summary.tab.z[3].get_integer_vec()[0],
            0b0000000000000000000100001000i128
        );
        assert_eq!(summary.tab.signs.get(3), false);
        // Zin1
        assert_eq!(summary.tab.x[4].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[4].get_integer_vec()[0],
            0b0100000000000000000100001100i128
        );
        assert_eq!(summary.tab.signs.get(4), false);
        // Yry.in
        assert_eq!(
            summary.tab.x[5].get_integer_vec()[0],
            0b0100101010000000001000001000i128
        );
        assert_eq!(
            summary.tab.z[5].get_integer_vec()[0],
            0b0100000000000000000000001000i128
        );
        assert_eq!(summary.tab.signs.get(5), false);
        // Xin2
        assert_eq!(
            summary.tab.x[6].get_integer_vec()[0],
            0b1000000000000000000000010000i128
        );
        assert_eq!(summary.tab.z[6].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(6), false);
        // Zin2
        assert_eq!(summary.tab.x[7].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[7].get_integer_vec()[0],
            0b1001000000000000000011110000i128
        );
        assert_eq!(summary.tab.signs.get(7), false);
        // Xrx.in
        assert_eq!(
            summary.tab.x[8].get_integer_vec()[0],
            0b1000000000000000000000100000i128
        );
        assert_eq!(summary.tab.z[8].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(8), false);
        // Xrx.out
        assert_eq!(
            summary.tab.x[9].get_integer_vec()[0],
            0b1000000000000000000001000000i128
        );
        assert_eq!(summary.tab.z[9].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(9), false);
        // Xtoffoli.in2
        assert_eq!(
            summary.tab.x[10].get_integer_vec()[0],
            0b1000000000000000000010000000i128
        );
        assert_eq!(summary.tab.z[10].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(10), false);
        // Yry.out
        assert_eq!(
            summary.tab.x[11].get_integer_vec()[0],
            0b0100101010000000001100000000i128
        );
        assert_eq!(
            summary.tab.z[11].get_integer_vec()[0],
            0b0100000000000000000100000000i128
        );
        assert_eq!(summary.tab.signs.get(11), true);
        // Zcrz.in1
        assert_eq!(summary.tab.x[12].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[12].get_integer_vec()[0],
            0b0100000000000000001000000000i128
        );
        assert_eq!(summary.tab.signs.get(12), false);
        // Zt.out
        assert_eq!(summary.tab.x[13].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[13].get_integer_vec()[0],
            0b0010000000000000010000000000i128
        );
        assert_eq!(summary.tab.signs.get(13), false);
        // Ztdg.in
        assert_eq!(summary.tab.x[14].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[14].get_integer_vec()[0],
            0b0010000000000000100000000000i128
        );
        assert_eq!(summary.tab.signs.get(14), false);
        // Ztdg.out
        assert_eq!(summary.tab.x[15].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[15].get_integer_vec()[0],
            0b0010000000000001000000000000i128
        );
        assert_eq!(summary.tab.signs.get(15), false);
        // Zrz.in
        assert_eq!(summary.tab.x[16].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[16].get_integer_vec()[0],
            0b0010000000000010000000000000i128
        );
        assert_eq!(summary.tab.signs.get(16), false);
        // Zrz.out
        assert_eq!(summary.tab.x[17].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[17].get_integer_vec()[0],
            0b0010000000000100000000000000i128
        );
        assert_eq!(summary.tab.signs.get(17), false);
        // Zmeas.in
        assert_eq!(summary.tab.x[18].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[18].get_integer_vec()[0],
            0b0010000000001000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(18), false);
        // Zmeas.out
        assert_eq!(summary.tab.x[19].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[19].get_integer_vec()[0],
            0b0010000000010000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(19), false);
        // Zcrz.in0
        assert_eq!(summary.tab.x[20].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[20].get_integer_vec()[0],
            0b0010000000100000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(20), false);
        // Zcrz.out0
        assert_eq!(summary.tab.x[21].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[21].get_integer_vec()[0],
            0b0010000001000000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(21), false);
        // Zcrz.out1
        assert_eq!(summary.tab.x[22].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[22].get_integer_vec()[0],
            0b0100000010000000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(22), false);
        // Ztoffoli.in0
        assert_eq!(summary.tab.x[23].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[23].get_integer_vec()[0],
            0b0010000100000000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(23), false);
        // Ztoffoli.in1
        assert_eq!(summary.tab.x[24].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[24].get_integer_vec()[0],
            0b0100001000000000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(24), false);
        // Ztoffoli.out0
        assert_eq!(summary.tab.x[25].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[25].get_integer_vec()[0],
            0b0010010000000000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(25), false);
        // Ztoffoli.out1
        assert_eq!(summary.tab.x[26].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[26].get_integer_vec()[0],
            0b0100100000000000000000000000i128
        );
        assert_eq!(summary.tab.signs.get(26), false);
        // Xtoffoli.out2
        assert_eq!(
            summary.tab.x[27].get_integer_vec()[0],
            0b1001000000000000000000000000i128
        );
        assert_eq!(summary.tab.z[27].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(27), false);
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
        let [qb2] = builder
            .add_dataflow_op(TketOp::QAlloc, [])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb0] = builder
            .add_dataflow_op(TketOp::Reset, [qb0])
            .unwrap()
            .outputs_arr();
        let meas = builder
            .add_dataflow_op(TketOp::MeasureFree, [qb1])
            .unwrap()
            .node();
        let [qb0] = builder
            .add_dataflow_op(TketOp::H, [qb0])
            .unwrap()
            .outputs_arr();
        let [qb0, qb2] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb2])
            .unwrap()
            .outputs_arr();
        builder.add_dataflow_op(TketOp::QFree, [qb2]).ok();
        let hugr = builder.finish_hugr_with_outputs([qb0]).unwrap();
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 4);
        // Input wires, alloc, and reset-alloc give 6 qubits/stabs
        // Reset-free and QFree remove 2 each
        // MeasureFree just acts as an opaque gate so doesn't remove any
        assert_eq!(summary.tab.nb_stabs, 2);
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::Output(IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::InternalIn(meas, IncomingPort::from(0))
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Zin0 Zin1 Zmeas
        assert_eq!(summary.tab.x[0].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0b1101i128);
        assert_eq!(summary.tab.signs.get(0), false);
        // Xin1 Xmeas
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0b1100i128);
        assert_eq!(summary.tab.z[1].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(1), false);
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
        let [c0q0, c0q1] = cond0_builder
            .add_dataflow_op(TketOp::CX, [c0q0, c0q1])
            .unwrap()
            .outputs_arr();
        cond0_builder.finish_with_outputs([c0q0, c0q1]).ok();
        let cond1_builder = cond_builder.case_builder(1).unwrap();
        let [c1c0, c1q1] = cond1_builder.input_wires_arr();
        cond1_builder.finish_with_outputs([c1c0, c1q1]).ok();
        let cond = cond_builder.finish_sub_container().unwrap();
        let [qb0, qb1] = cond.outputs_arr();
        let tdg = builder.add_dataflow_op(TketOp::Tdg, [qb0]).unwrap();
        let [qb0] = tdg.outputs_arr();
        let hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 16);
        assert_eq!(summary.tab.nb_stabs, 14);
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::InternalIn(t.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::InternalIn(cond.node(), IncomingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::InternalOut(t.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::InternalIn(cond.node(), IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::NestedIn(cond.node(), OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::NestedIn(cond.node(), OutgoingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::NestedOut(cond.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::NestedOut(cond.node(), IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::InternalOut(cond.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::InternalIn(tdg.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::InternalOut(cond.node(), OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::Output(IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::InternalOut(tdg.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::Output(IncomingPort::from(0))
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Zin0
        assert_eq!(summary.tab.x[0].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[0].get_integer_vec()[0],
            0b1000000000000001i128
        );
        assert_eq!(summary.tab.signs.get(0), false);
        // Zt.in
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[1].get_integer_vec()[0],
            0b1000000000000010i128
        );
        assert_eq!(summary.tab.signs.get(1), false);
        // Xin1
        assert_eq!(
            summary.tab.x[2].get_integer_vec()[0],
            0b0010000000000100i128
        );
        assert_eq!(summary.tab.z[2].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(2), false);
        // Note that after projecting ZZ and XX, the internals and nested boundaries of the conditional block have been disconnected from the rest of the circuit, now sharing Bell states
        // Xcond.in1
        assert_eq!(
            summary.tab.x[3].get_integer_vec()[0],
            0b0000000010001000i128
        );
        assert_eq!(summary.tab.z[3].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(3), false);
        // Zcond.in1
        assert_eq!(summary.tab.x[4].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[4].get_integer_vec()[0],
            0b0000000010001000i128
        );
        assert_eq!(summary.tab.signs.get(4), false);
        // Zt.out
        assert_eq!(summary.tab.x[5].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[5].get_integer_vec()[0],
            0b1000000000010000i128
        );
        assert_eq!(summary.tab.signs.get(5), false);
        // Xcond.in0
        assert_eq!(
            summary.tab.x[6].get_integer_vec()[0],
            0b0000000001100000i128
        );
        assert_eq!(summary.tab.z[6].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(6), false);
        // Zcond.in0
        assert_eq!(summary.tab.x[7].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[7].get_integer_vec()[0],
            0b0000000001100000i128
        );
        assert_eq!(summary.tab.signs.get(7), false);
        // Xcond.nout0
        assert_eq!(
            summary.tab.x[8].get_integer_vec()[0],
            0b0000010100000000i128
        );
        assert_eq!(summary.tab.z[8].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(8), false);
        // Zcond.nout0
        assert_eq!(summary.tab.x[9].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[9].get_integer_vec()[0],
            0b0000010100000000i128
        );
        assert_eq!(summary.tab.signs.get(9), false);
        // Xcond.nout1
        assert_eq!(
            summary.tab.x[10].get_integer_vec()[0],
            0b0001001000000000i128
        );
        assert_eq!(summary.tab.z[10].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(10), false);
        // Zcond.nout1
        assert_eq!(summary.tab.x[11].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[11].get_integer_vec()[0],
            0b0001001000000000i128
        );
        assert_eq!(summary.tab.signs.get(11), false);
        // Ztdg.in
        assert_eq!(summary.tab.x[12].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[12].get_integer_vec()[0],
            0b1000100000000000i128
        );
        assert_eq!(summary.tab.signs.get(12), false);
        // Ztdg.out
        assert_eq!(summary.tab.x[13].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[13].get_integer_vec()[0],
            0b1100000000000000i128
        );
        assert_eq!(summary.tab.signs.get(13), false);
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
        let analysis = SDFAnalysis::run_hugr(&hugr, &FunctionOpacity::Opaque);
        let mut summary = analysis
            .0
            .get(&hugr.first_child(hugr.module_root()).unwrap())
            .unwrap()
            .clone();
        assert_eq!(summary.tab.nb_qubits, 16);
        assert_eq!(summary.tab.nb_stabs, 14);
        assert_eq!(
            *summary.q_index_map.get_by_right(&0).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&1).unwrap(),
            DataflowPoint::InternalIn(tl.node(), IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&2).unwrap(),
            DataflowPoint::Input(OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&3).unwrap(),
            DataflowPoint::InternalIn(t.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&4).unwrap(),
            DataflowPoint::InternalOut(t.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&5).unwrap(),
            DataflowPoint::InternalIn(tl.node(), IncomingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&6).unwrap(),
            DataflowPoint::NestedIn(tl.node(), OutgoingPort::from(1))
        );
        // Ports from NestedIn and NestedOut are wrt the nested summary; in this case respecting the final signature of the TailLoop and not the signature of the body
        assert_eq!(
            *summary.q_index_map.get_by_right(&7).unwrap(),
            DataflowPoint::NestedOut(tl.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&8).unwrap(),
            DataflowPoint::NestedIn(tl.node(), OutgoingPort::from(2))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&9).unwrap(),
            DataflowPoint::NestedOut(tl.node(), IncomingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&10).unwrap(),
            DataflowPoint::InternalOut(tl.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&11).unwrap(),
            DataflowPoint::Output(IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&12).unwrap(),
            DataflowPoint::InternalOut(tl.node(), OutgoingPort::from(1))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&13).unwrap(),
            DataflowPoint::InternalIn(tdg.node(), IncomingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&14).unwrap(),
            DataflowPoint::InternalOut(tdg.node(), OutgoingPort::from(0))
        );
        assert_eq!(
            *summary.q_index_map.get_by_right(&15).unwrap(),
            DataflowPoint::Output(IncomingPort::from(1))
        );
        summary.tab.echelon(&summary.tab.all_columns());
        // Xtl.in1
        assert_eq!(
            summary.tab.x[0].get_integer_vec()[0],
            0b0000000001000010i128
        );
        assert_eq!(summary.tab.z[0].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(0), false);
        // Ztl.in1
        assert_eq!(summary.tab.x[1].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[1].get_integer_vec()[0],
            0b0000000001000010i128
        );
        assert_eq!(summary.tab.signs.get(1), false);
        // Zin1
        assert_eq!(summary.tab.x[2].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[2].get_integer_vec()[0],
            0b1000000000000100i128
        );
        assert_eq!(summary.tab.signs.get(2), false);
        // Zt.in
        assert_eq!(summary.tab.x[3].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[3].get_integer_vec()[0],
            0b1000000000001000i128
        );
        assert_eq!(summary.tab.signs.get(3), false);
        // Zt.out
        assert_eq!(summary.tab.x[4].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[4].get_integer_vec()[0],
            0b1000000000010000i128
        );
        assert_eq!(summary.tab.signs.get(4), false);
        // Xtl.in2
        assert_eq!(
            summary.tab.x[5].get_integer_vec()[0],
            0b0000000100100000i128
        );
        assert_eq!(summary.tab.z[5].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(5), false);
        // Ztl.in2
        assert_eq!(summary.tab.x[6].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[6].get_integer_vec()[0],
            0b0000000100100000i128
        );
        assert_eq!(summary.tab.signs.get(6), false);
        // Xtl.nout0
        assert_eq!(
            summary.tab.x[7].get_integer_vec()[0],
            0b0000010010000000i128
        );
        assert_eq!(summary.tab.z[7].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(7), false);
        // Ztl.nout0
        assert_eq!(summary.tab.x[8].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[8].get_integer_vec()[0],
            0b0000010010000000i128
        );
        assert_eq!(summary.tab.signs.get(8), false);
        // Xtl.nout1
        assert_eq!(
            summary.tab.x[9].get_integer_vec()[0],
            0b0001001000000000i128
        );
        assert_eq!(summary.tab.z[9].get_integer_vec()[0], 0i128);
        assert_eq!(summary.tab.signs.get(9), false);
        // Ztl.nout1
        assert_eq!(summary.tab.x[10].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[10].get_integer_vec()[0],
            0b0001001000000000i128
        );
        assert_eq!(summary.tab.signs.get(10), false);
        // Zout0
        assert_eq!(summary.tab.x[11].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[11].get_integer_vec()[0],
            0b0000100000000000i128
        );
        assert_eq!(summary.tab.signs.get(11), false);
        // Ztdg.in
        assert_eq!(summary.tab.x[12].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[12].get_integer_vec()[0],
            0b1010000000000000i128
        );
        assert_eq!(summary.tab.signs.get(12), false);
        // Ztdg.out
        assert_eq!(summary.tab.x[13].get_integer_vec()[0], 0i128);
        assert_eq!(
            summary.tab.z[13].get_integer_vec()[0],
            0b1100000000000000i128
        );
        assert_eq!(summary.tab.signs.get(13), false);
    }
}
