// Going to base this on Mark's python phase folding implementation instead of the hugr dataflow framework which I struggle to see how to adapt to relational values since we can't easily attribute them to individual wires

use std::collections::HashMap;
use hugr::PortIndex;
use hugr_core::hugr::internal::PortgraphNodeMap;
use hugr_core::{HugrView, IncomingPort, OutgoingPort};
use hugr_core::ops::OpType;
use hugr::extension::prelude::qb_t;
use itertools::Itertools;
use petgraph::visit as pv;
use tket::hugr::extension::simple_op::MakeExtensionOp;
use tket::TketOp;
use crate::bit_vector::BitVector;
use crate::pauli_product::PauliProduct;
use crate::tableau::Tableau;

/// Sets behaviour for function calls in dataflow analysis
#[derive(Clone)]
pub enum FunctionOpacity {
    /// Function calls are completely opaque and admit no information across them
    Opaque,
    /// Function bodies are analysed but we only pass on the stabilizers over the boundary qubits
    Boundary,
    /// Function bodies are analysed and fully inserted into the parent graph, allowing e.g. phase folding between gates inside the body and around the call site by inlining the function body
    Inline,
}

pub struct StabilizerDataflow<H: HugrView> {
    /// Relational dataflow value captured as a set of stabilizer relations on the Choi-state of the circuit skeleton
    tab: Tableau,
    /// Maps from wires of the program to columns of the tableau. We separately need to track columns for:
    /// - Each input qubit (indexed by OutgoingPorts of the unique Input node)
    /// - Each output qubit (indexed by IncomingPorts of the unique Output node)
    /// - A frontier that moves forward through the program (eventually becoming the output qubits and being removed from here)
    /// - For any internal non-Clifford (or opaque) node, we use columns for each input and output qubit separately; for nodes with stabilizers across them (e.g. Rz has Z_i Z_o), we impose these via projections on the tableau rather than reducing the number of qubits used as this allows every node kind to be handled identically and preventing more tableau management from column elimination
    /// - For any hierarchical node, we use additional columns for each input and output port within their internal representation that we compose to "internal" columns here by projections on the tableau, again so we don't fuss with column elimination
    in_cols: HashMap<OutgoingPort, usize>,
    out_cols: HashMap<IncomingPort, usize>,
    frontier_cols: HashMap<(H::Node, IncomingPort), usize>,
    internal_in_cols: HashMap<(H::Node, IncomingPort), usize>,
    internal_out_cols: HashMap<(H::Node, OutgoingPort), usize>,
    nested_in_cols: HashMap<(H::Node, OutgoingPort), usize>,
    nested_out_cols: HashMap<(H::Node, IncomingPort), usize>,

    // For any control-flow region or hierarchical node, store the analysis for its internal calculations
    nested_analysis: HashMap<H::Node, StabilizerDataflow<H>>,
}

impl<H: HugrView> StabilizerDataflow<H> {
    fn new(hugr: &H, parent: H::Node) -> Self {
        let mut in_cols: HashMap<OutgoingPort, usize> = HashMap::default();
        let mut frontier_cols: HashMap<(H::Node, IncomingPort), usize> = HashMap::default();
        let mut n_in_qubits = 0;
        let inp = hugr.children(parent).filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_))).exactly_one().ok().unwrap();
        for (out, out_type) in hugr.out_value_types(inp) {
            if out_type == qb_t() {
                in_cols.insert(out, 2*n_in_qubits);
                let (next, next_p) = hugr.single_linked_input(inp, out).unwrap();
                frontier_cols.insert((next, next_p), 2*n_in_qubits + 1);
                n_in_qubits = n_in_qubits + 1;
            }
        }
        let tab = Tableau::new(2*n_in_qubits);
        //TODO:: Add rows to tableau
        Self{
            tab: tab,
            in_cols: in_cols,
            out_cols: HashMap::default(),
            frontier_cols: frontier_cols,
            internal_in_cols: HashMap::default(),
            internal_out_cols: HashMap::default(),
            nested_in_cols: HashMap::default(),
            nested_out_cols: HashMap::default(),
            nested_analysis: HashMap::default(),
        }
    }

    pub fn run_dfg(hugr: &H, parent: H::Node, fun_op: &FunctionOpacity) -> StabilizerDataflow<H> {
        let mut analysis = StabilizerDataflow::new(hugr, parent);
        let (region, node_map) = hugr.region_portgraph(parent);
        let mut topo = pv::Topo::new(&region);
        while let Some(pgnode) = topo.next(&region) {
            let node = node_map.from_portgraph(pgnode);
            let optype: &OpType = hugr.get_optype(node);
            match optype {
                OpType::ExtensionOp(op) => {
                    match TketOp::from_extension_op(op) {
                        Ok(tkop) => analysis.apply_quantum_gate(hugr, node, tkop),
                        Err(_) => analysis.apply_opaque(hugr, node)
                    }
                }
                OpType::Conditional(_) => {
                    let cond_analysis = StabilizerDataflow::run_conditional(hugr, node, fun_op);
                    analysis.nested_analysis.insert(node, cond_analysis);
                    analysis.apply_analysis(hugr, node);
                }
                OpType::TailLoop(_) => {
                    let loop_analysis = StabilizerDataflow::run_tail_loop(hugr, node, fun_op);
                    analysis.nested_analysis.insert(node, loop_analysis);
                    analysis.apply_analysis(hugr, node);
                }
                OpType::Call(_) => {
                    match *fun_op {
                        FunctionOpacity::Opaque => {
                            analysis.apply_opaque(hugr, node);
                        }
                        FunctionOpacity::Boundary => {
                            let call_port = optype.static_input_port().unwrap();
                            let (fun_def_node, _) = hugr.linked_outputs(node, call_port).exactly_one().ok().unwrap();
                            let fun_analysis = StabilizerDataflow::run_dfg(hugr, fun_def_node, fun_op);
                            //TODO:: Project out non-IO columns
                            analysis.nested_analysis.insert(node, fun_analysis);
                            analysis.apply_analysis(hugr, node);
                        }
                        FunctionOpacity::Inline => {
                            let call_port = optype.static_input_port().unwrap();
                            let (fun_def_node, _) = hugr.linked_outputs(node, call_port).exactly_one().ok().unwrap();
                            let fun_analysis = StabilizerDataflow::run_dfg(hugr, fun_def_node, fun_op);
                            analysis.nested_analysis.insert(node, fun_analysis);
                            analysis.apply_analysis(hugr, node);
                        }
                    }
                }
                OpType::Output(_) => {
                    // Frontier finished, move it to out_cols
                    for ((_, port), col) in analysis.frontier_cols.iter() {
                        analysis.out_cols.insert(*port, *col);
                    }
                    analysis.frontier_cols = HashMap::default();
                }
                _ => {
                    analysis.apply_opaque(hugr, node)
                }
            }
        }
        analysis
    }

    fn run_conditional(hugr: &H, node: H::Node, fun_op: &FunctionOpacity) -> StabilizerDataflow<H> {
        let mut summary: Option<StabilizerDataflow<H>> = None;
        for cond_node in hugr.children(node) {
            let analysis = StabilizerDataflow::run_dfg(hugr, cond_node, fun_op);
            let mut tab = analysis.tab.clone();
            //TODO:: Project out non-IO columns
            match summary {
                Some(_) => {
                    //TODO:: Reorder and remove columns of tab to match pfa
                    //TODO:: Compute join of tab and pfa.tab
                }
                None => {
                    // Determine consistent column indexing for inputs and outputs
                    let mut in_cols : HashMap<OutgoingPort, usize> = HashMap::default();
                    let mut out_cols : HashMap<IncomingPort, usize> = HashMap::default();
                    let mut n_qbs = 0;
                    let inp = hugr.children(cond_node).filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_))).exactly_one().ok().unwrap();
                    for (out_port, out_type) in hugr.out_value_types(inp) {
                        if out_type == qb_t() {
                            in_cols.insert(out_port, n_qbs);
                            n_qbs = n_qbs + 1;
                        }
                    }
                    let outp = hugr.children(cond_node).filter(|n| matches!(hugr.get_optype(*n), OpType::Output(_))).exactly_one().ok().unwrap();
                    for (in_port, in_type) in hugr.in_value_types(outp) {
                        if in_type == qb_t() {
                            out_cols.insert(in_port, n_qbs);
                            n_qbs = n_qbs + 1;
                        }
                    }
                    //TODO:: Reorder and remove columns of tab to match this column indexing
                    summary = Some(StabilizerDataflow {
                        tab: tab,
                        in_cols: in_cols,
                        out_cols: out_cols,
                        frontier_cols: HashMap::default(),
                        internal_in_cols: HashMap::default(),
                        internal_out_cols: HashMap::default(),
                        nested_in_cols: HashMap::default(),
                        nested_out_cols: HashMap::default(),
                        nested_analysis: HashMap::default(),
                    });
                    summary.as_mut().unwrap().nested_analysis.insert(cond_node, analysis);
                }
            }
        }
        summary.unwrap()
    }

    fn run_tail_loop(hugr: &H, node: H::Node, fun_op: &FunctionOpacity) -> StabilizerDataflow<H> {
        let child_node = hugr.children(node).exactly_one().ok().unwrap();
        let child_analysis = StabilizerDataflow::run_dfg(hugr, child_node, fun_op);
        let mut analysis = StabilizerDataflow{
            tab: Tableau::new(0),
            in_cols: HashMap::default(),
            out_cols: HashMap::default(),
            frontier_cols: HashMap::default(),
            internal_in_cols: HashMap::default(),
            internal_out_cols: HashMap::default(),
            nested_in_cols: HashMap::default(),
            nested_out_cols: HashMap::default(),
            nested_analysis: HashMap::default()
        };
        for (port, _) in child_analysis.in_cols.iter() {
            let new_col = analysis.tab.add_col();
            analysis.in_cols.insert(*port, new_col);
        }
        for (port, _) in child_analysis.out_cols.iter() {
            let new_col = analysis.tab.add_col();
            analysis.out_cols.insert(*port, new_col);
        }
        for (port, in_col) in analysis.in_cols.iter() {
            let out_port = IncomingPort::from(port.index());
            let out_col = analysis.out_cols.get(&out_port);
            //TODO:: Add rows for identity in_col--out_col
        }
        let mut tab = child_analysis.tab.clone();
        //TODO:: Project out non-IO columns
        //TODO:: Reorder and remove columns of tab to match analysis.tab
        //TODO:: Compute join of tabs
        analysis.nested_analysis.insert(child_node, child_analysis);
        analysis
    }

    fn apply_quantum_gate(&mut self, hugr : &H, node: H::Node, op: TketOp) {
        match op {
            TketOp::H => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_h(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::CX => {
                let col0: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col1: usize = self.frontier_cols.remove(&(node, IncomingPort::from(1))).unwrap();
                self.tab.append_cx(vec![col0, col1]);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
            }
            TketOp::CY => {
                let col0: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col1: usize = self.frontier_cols.remove(&(node, IncomingPort::from(1))).unwrap();
                self.tab.append_s(col1);
                self.tab.append_z(col1);
                self.tab.append_cx(vec![col0, col1]);
                self.tab.append_s(col1);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
            }
            TketOp::CZ => {
                let col0: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col1: usize = self.frontier_cols.remove(&(node, IncomingPort::from(1))).unwrap();
                self.tab.append_cz(vec![col0, col1]);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
            }
            TketOp::CRz => {
                let col_in0: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col_in1: usize = self.frontier_cols.remove(&(node, IncomingPort::from(1))).unwrap();
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
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front0);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col_front1);
            }
            TketOp::T | TketOp::Tdg | TketOp::Rz | TketOp::Measure => {
                let col_in: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_col();
                let col_front: usize = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                //TODO:: Add row for ZZ over col_in--col_out and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            TketOp::S => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_s(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::Sdg => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_s(col);
                self.tab.append_z(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::X => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_x(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::Y => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_x(col);
                self.tab.append_z(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::Z => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_z(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::Rx => {
                let col_in: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_col();
                let col_front: usize = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                //TODO:: Add row for XX over col_in--col_out and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            TketOp::Ry => {
                let col_in: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_col();
                let col_front: usize = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                //TODO:: Add row for YY over col_in--col_out and project to commuting
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            TketOp::MeasureFree => {
                let col_in: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
            }
            TketOp::QAlloc => {
                let col_front: usize = self.tab.add_col();
                //TODO:: Add row for Z over col_front
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
            }
            TketOp::QFree => {
                let col_in: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                //TODO:: Project out non-commuting rows and remove column from tableau
            }
            TketOp::Reset => {
                let col_in: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                //TODO:: Project out non-commuting rows
                // Reuse col_in for the output qubit
                //TODO:: Add row for Z over col_in
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_in);
            }
            TketOp::V => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_v(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            TketOp::Vdg => {
                let col: usize = self.frontier_cols.remove(&(node, IncomingPort::from(0))).unwrap();
                self.tab.append_v(col);
                self.tab.append_x(col);
                self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
            }
            _ => {
              self.apply_opaque(hugr, node)
            }
        }
    }

    fn apply_opaque(&mut self, hugr: &H, node: H::Node) {
        // For each Qubit input, move the column from frontier_cols to internal_in_cols
        for (p, t) in hugr.in_value_types(node) {
            if t == qb_t() {
                let col: usize = self.frontier_cols.remove(&(node, p)).unwrap();
                self.internal_in_cols.insert((node, p), col);
            }
        }
        // For each Qubit output, create a pair of columns with the identity for internal_out_cols and frontier_cols
        for (p, t) in hugr.out_value_types(node) {
            if t == qb_t() {
                let col_out = self.tab.add_col();
                let col_front = self.tab.add_col();
                //TODO:: Add rows for identity col_out--col_front
                self.internal_out_cols.insert((node, p), col_out);
                self.frontier_cols.insert(hugr.single_linked_input(node, p).unwrap(), col_front);
            }
        }
    }

    /// Suppose we have already recursively calculated a StabilizerDataflow for node and stored it in nested_analysis; performs sequential composition to append it to the appropriate qubits here
    fn apply_analysis(&mut self, hugr: &H, node: H::Node) {
        let node_analysis : &StabilizerDataflow<H> = self.nested_analysis.get(&node).unwrap();
        let old_n_qbs = self.tab.nb_qubits;
        let n_added_qbs = node_analysis.tab.nb_qubits;
        for i in 0..n_added_qbs {
            self.tab.add_col();
        }
        for (port, col) in node_analysis.in_cols.iter() {
            self.nested_in_cols.insert((node, *port), *col + old_n_qbs);
        }
        for (port, col) in node_analysis.out_cols.iter() {
            self.nested_out_cols.insert((node, *port), *col + old_n_qbs);
        }
        // It is decided at analysis construction whether or not we project to IO or keep internals, so always copy any internals that have remained (e.g. for inlining function calls)
        for (node_port, col) in node_analysis.internal_in_cols.iter() {
            self.internal_in_cols.insert(*node_port, *col + old_n_qbs);
        }
        for (node_port, col) in node_analysis.internal_out_cols.iter() {
            self.internal_out_cols.insert(*node_port, *col + old_n_qbs);
        }
        for (node_port, col) in node_analysis.nested_in_cols.iter() {
            self.nested_in_cols.insert(*node_port, *col + old_n_qbs);
        }
        for (node_port, col) in node_analysis.nested_out_cols.iter() {
            self.nested_out_cols.insert(*node_port, *col + old_n_qbs);
        }
        for i in 0..node_analysis.tab.nb_rows {
            let stab : &PauliProduct = node_analysis.tab.stabs.get(i);
            let mut new_z = BitVector::new(old_n_qbs);
            new_z.extend_vec(stab.z.get_boolean_vec(), old_n_qbs);
            let mut new_x = BitVector::new(old_n_qbs);
            new_x.extend_vec(stab.x.get_boolean_vec(), old_n_qbs);
            self.tab.add_row(PauliProduct{
                z: new_z,
                x: new_x,
                sign: stab.sign
            });
        }
        for port in hugr.node_inputs(node) {
            let out_port = OutgoingPort::from(port.index());
            let internal_col = self.frontier_cols.remove(&(node, port)).unwrap();
            self.internal_in_cols.insert((node, port), internal_col);
            let nested_col = self.nested_in_cols.get(&(node, out_port));
            //TODO:: Project ZZ and XX to compose nested_col and internal_col
        }
        for port in hugr.node_outputs(node) {
            let in_port = IncomingPort::from(port.index());
            let nested_col = self.nested_out_cols.get(&(node, in_port));
            let internal_col = self.tab.add_col();
            self.internal_out_cols.insert((node, port), internal_col);
            let front_col = self.tab.add_col();
            self.frontier_cols.insert(hugr.single_linked_input(node, port).unwrap(), front_col);
            //TODO:: Add rows for identity internal_col--front_col
            //TODO:: Project ZZ and XX to compose nested_col and internal_col
        }
    }

}