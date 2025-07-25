// Going to base this on Mark's python phase folding implementation instead of the hugr dataflow framework which I struggle to see how to adapt to relational values since we can't easily attribute them to individual wires

use std::collections::HashMap;
use hugr_core::{HugrView, IncomingPort, OutgoingPort};
use hugr_core::ops::OpType;
use hugr::extension::prelude::qb_t;
use petgraph::visit as pv;
use tket2::Tk2Op;
use crate::tableau::ChoiTableau;

pub struct StabilizerDataflow<H: HugrView> {
    /// Relational dataflow value captured as a set of stabilizer relations on the Choi-state of the circuit skeleton
    tab: ChoiTableau,
    /// Maps from wires of the program to columns of the tableau. We separately need to track columns for:
    /// - Each input qubit (indexed by OutgoingPorts of the unique Input node)
    /// - A frontier that moves forward through the program (eventually becoming the output qubits)
    /// - For any internal non-Clifford (or opaque) node, we use columns for each input and output qubit separately; for nodes with stabilizers across them (e.g. Rz has Z_i Z_o), we impose these via projections on the tableau rather than reducing the number of qubits used as this allows every node kind to be handled identically
    in_cols: HashMap<OutgoingPort, usize>,
    frontier_cols: HashMap<(H::Node, IncomingPort), usize>,
    internal_in_cols: HashMap<(H::Node, IncomingPort), usize>,
    internal_out_cols: HashMap<(H::Node, OutgoingPort), usize>,

    // For any control-flow region or hierarchical node, store the analysis for its internal calculations
    nested_analysis: HashMap<H::Node, StabilizerDataflow>,

    // Maintain a HugrView and parent node for convenience
    hugr: H,
    parent: H::Node,
}

impl<H: HugrView> StabilizerDataflow<H> {
    fn new(hugr: H, parent: H::Node) -> Self {
        let mut in_cols = HashMap<H::Node, HashMap<OutgoingPort, usize>>;
        let mut frontier_cols = HashMap<H::Node, HashMap<IncomingPort, usize>>;
        let mut n_in_qubits = 0;
        let inp = find_unique(hugr.children(parent), |n| matches!(hugr.get_optype(*n), OpType::Input(_)));
        for (out, out_type) in hugr.out_value_types(inp) {
            if out_type == qb_t() {
                in_cols.insert(out, 2*n_in_qubits);
                let next, next_p = hugr.single_linked_input(inp, out).unwrap();
                frontier_cols.insert((next, next_p), 2*n_in_qubits + 1)
                n_in_qubits = n_in_qubits + 1;
            }
        }
        let tab = ChoiTableau(2*n_in_qubits);
        //TODO:: Add rows to tableau
        Self(tab, in_cols, frontier_cols, HashMap<(H::Node, IncomingPort), usize>::default(), HashMap<(H::Node, OutgoingPort), usize>::default(), HashMap<H::Node, StabilizerDataflow>::default(), hugr, parent)
    }

    pub fn run_dfg(hugr: H, parent: H::Node) -> StabilizerDataflow {
        let analysis = StabilizerDataflow(hugr, parent);
        let region: SiblingGraph = SiblingGraph::try_new(hugr, parent).unwrap();
        let mut topo: pv::Topo = pv::Topo::new(&region);
        while let Some(node) = topo.next(&region) {
            let optype: OpType = hugr.get_optype(node);
            match optype {
                OpType::ExtensionOp(op) => {
                    match Tk2Op::from_extension_op(op) {
                        Ok(tkop) => analysis.apply_quantum_gate(node, tkop),
                        Err(_) => analysis.apply_opaque(node)
                    }
                }
                OpType::Conditional(_) => {
                    //TODO:: Handle conditionals
                }
                OpType::TailLoop(_) => {
                    //TODO:: Handle loops
                }
                OpType::Call(_) => {
                    //TODO:: Handle function calls
                }
                OpType::Output(_) => {
                    // No need to do anything
                }
                _ => {
                    analysis.apply_opaque(node)
                }
            }
        }
        analysis
    }

    fn run_conditional(hugr: H, node: H::Node) -> StabilizerDataflow {
        let mut summary: Option<StabilizerDataflow> = None;
        for cond_node in hugr.children(node) {
            let analysis = run_dfg(hugr, cond_node);
            let mut tab = analysis.tab.clone();
            //TODO:: Project out non-IO columns
            match summary {
                Some(pfa) => {
                    //TODO:: Reorder and remove columns of tab to match pfa
                    //TODO:: Compute join of tab and pfa.tab
                }
                None => {
                    // Determine consistent column indexing for inputs and outputs
                    let mut in_cols = HashMap<OutgoingPort, usize>::default();
                    let mut out_cols = HashMap<(H::Node, IncomingPort), usize>::default();
                    let mut n_qbs = 0;
                    let inp = find_unique(self.hugr.children(self.cond_node), |n| matches!(self.hugr.get_optype(*n), OpType::Input(_)));
                    for (out_port, out_type) in self.hugr.out_value_types(inp) {
                        if out_type == qb_t() {
                            in_cols.insert(out_port, n_qbs);
                            n_qbs = n_qbs + 1;
                        }
                    }
                    let outp = find_unique(self.hugr.children(self.cond_node), |n| matches!(self.hugr.get_optype(*n), OpType::Output(_)));
                    for (in_port, in_type) in self.hugr.in_value_types(outp) {
                        if in_type == qb_t() {
                            out_cols.insert((outp, in_port), n_qbs);
                            n_qbs = n_qbs + 1;
                        }
                    }
                    //TODO:: Reorder and remove columns of tab to match this column indexing
                    summary = Some(StabilizerDataflow(tab, in_cols, out_cols, HashMap<(H::Node, IncomingPort), usize>::default(), HashMap<(H::Node, OutgoingPort), usize>::default(), HashMap<H::Node, StabilizerDataflow>::default(), self.hugr, self.parent));
                    summary.unwrap().nested_analysis.insert(cond_node, analysis);
                }
            }
        }
        summary.unwrap()
    }

    fn apply_quantum_gate(&mut self, node: H::Node, op: Tk2Op) {
        match op {
            Tk2Op::H => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_h(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::CX => {
                let col0: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col1: usize = self.frontier_cols.remove((node, IncomingPort::from(1))).unwrap();
                self.tab.append_cx(col0, col1);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
            }
            Tk2Op::CY => {
                let col0: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col1: usize = self.frontier_cols.remove((node, IncomingPort::from(1))).unwrap();
                self.tab.append_cy(col0, col1);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
            }
            Tk2Op::CZ => {
                let col0: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col1: usize = self.frontier_cols.remove((node, IncomingPort::from(1))).unwrap();
                self.tab.append_cz(col0, col1);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
            }
            Tk2Op::CRz => {
                let col_in0: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col_in1: usize = self.frontier_cols.remove((node, IncomingPort::from(1))).unwrap();
                let col_out0: usize = self.tab.add_col();
                let col_out1: usize = self.tab.add_col();
                let col_front0: usize = self.tab.add_col();
                let col_front1: usize = self.tab.add_col();
                //TODO:: Add rows for identities col_out0/1--col_front0/1
                //TODO:: Add rows for ZZ over col_in0/1--col_out0/1 and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in0);
                self.internal_in_cols.insert((node, IncomingPort::from(1)), col_in1);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out0);
                self.internal_out_cols.insert((node, OutgoingPort::from(1)), col_out1);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front0);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col_front1);
            }
            Tk2Op::T, Tk2Op::Tdg, Tk2Op::Rz, Tk2Op::Measure => {
                let col_in: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_col();
                let col_front: usize = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                //TODO:: Add row for ZZ over col_in--col_out and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            Tk2Op::S => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_s(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::Sdg => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_s(col);
                self.tab.append_z(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::X => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_x(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::Y => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_x(col);
                self.tab.append_z(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::Z => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_z(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::Rx => {
                let col_in: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_col();
                let col_front: usize = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                //TODO:: Add row for XX over col_in--col_out and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            Tk2Op::Ry => {
                let col_in: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_col();
                let col_front: usize = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                //TODO:: Add row for YY over col_in--col_out and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            Tk2Op::MeasureFree => {
                let col_in: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
            }
            Tk2Op::QAlloc => {
                let col_front: usize = self.tab.add_col();
                //TODO:: Add row for Z over col_front
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            Tk2Op::QFree => {
                let col_in: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                //TODO:: Project out non-commuting rows and remove column from tableau
            }
            Tk2Op::Reset => {
                let col_in: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                //TODO:: Project out non-commuting rows
                // Reuse col_in for the output qubit
                //TODO:: Add row for Z over col_in
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_in);
            }
            Tk2Op::V => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_v(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            Tk2Op::Vdg => {
                let col: usize = self.frontier_cols.remove((node, IncomingPort::from(0))).unwrap();
                self.tab.append_v(col);
                self.tab.append_x(col);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            _ => {
              apply_opaque(node)
            }
        }
    }

    fn apply_opaque(&mut self, node: H::Node) {
        // For each Qubit input, move the column from frontier_cols to internal_in_cols
        for (p, t) in self.hugr.in_value_types(node) {
            if t == qb_t() {
                let col: usize = self.frontier_cols.remove((node, p)).unwrap();
                self.internal_in_cols.insert((node, p), col);
            }
        }
        // For each Qubit output, create a pair of columns with the identity for internal_out_cols and frontier_cols
        for (p, t) in self.hugr.out_value_types(node) {
            if t == qb_t() {
                let col_out = self.tab.add_col();
                let col_front = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                self.internal_out_cols.insert((node, p), col_out);
                self.frontier_cols.insert(self.hugr.single_linked_input(node, p).unwrap(), col_front);
            }
        }
    }

    /// Suppose we have already recursively calculated a StabilizerDataflow for node and stored it in nested_analysis; performs sequential composition to append it to the appropriate qubits here
    fn apply_analysis(&mut self, node: H::Node) {
        //TODO:: LEFT OFF HERE
    }

}