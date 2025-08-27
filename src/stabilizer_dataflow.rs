// Going to base this on Mark's python phase folding implementation instead of the hugr dataflow framework which I struggle to see how to adapt to relational values since we can't easily attribute them to individual wires

use std::collections::HashMap;
use std::hash::Hash;
use bimap::BiHashMap;
use hugr::ops::DataflowOpTrait;
use hugr::PortIndex;
use hugr_core::hugr::internal::{self, PortgraphNodeMap};
use hugr_core::{HugrView, IncomingPort, OutgoingPort};
use hugr_core::ops::OpType;
use hugr::extension::prelude::qb_t;
use itertools::{chain, Itertools};
use petgraph::visit::{self as pv, Data};
use tket::hugr::extension::simple_op::MakeExtensionOp;
use tket::TketOp;
use crate::bit_vector::BitVector;
use crate::pauli_product::PauliProduct;
use crate::symplectic_tableau::{PauliXZ, SymplecticTableau};

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

/// For each qubit in a StabilizerDataflow analysis, decides its role and attaches it to a node if relevant
#[derive(Clone, PartialEq, Eq, Hash)]
pub enum DataflowPoint<N : Copy + Eq + Hash> {
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

#[derive(Clone)]
pub struct StabilizerDataflow<H: HugrView> {
    /// Relational dataflow value captured as a set of stabilizer relations on the Choi-state of the circuit skeleton
    tab: SymplecticTableau,
    /// Maps from wires of the program to columns of the tableau. We separately need to track columns for:
    /// - Each input qubit (indexed by OutgoingPorts of the unique Input node)
    /// - Each output qubit (indexed by IncomingPorts of the unique Output node)
    /// - A frontier that moves forward through the program (eventually becoming the output qubits and being removed from here)
    /// - For any internal non-Clifford (or opaque) node, we use columns for each input and output qubit separately; for nodes with stabilizers across them (e.g. Rz has Z_i Z_o), we impose these via projections on the tableau rather than reducing the number of qubits used as this allows every node kind to be handled identically and preventing more tableau management from column elimination
    /// - For any hierarchical node, we use additional columns for each input and output port within their internal representation that we compose to "internal" columns here by projections on the tableau, again so we don't fuss with column elimination
    // in_cols: BiHashMap<OutgoingPort, usize>,
    // out_cols: BiHashMap<IncomingPort, usize>,
    // frontier_cols: BiHashMap<(H::Node, IncomingPort), usize>,
    // internal_in_cols: BiHashMap<(H::Node, IncomingPort), usize>,
    // internal_out_cols: BiHashMap<(H::Node, OutgoingPort), usize>,
    // nested_in_cols: BiHashMap<(H::Node, OutgoingPort), usize>,
    // nested_out_cols: BiHashMap<(H::Node, IncomingPort), usize>,
    q_index_map: BiHashMap<DataflowPoint<H::Node>, usize>,

    // For any control-flow region or hierarchical node, store the analysis for its internal calculations
    nested_analysis: HashMap<H::Node, StabilizerDataflow<H>>,
}

impl<H: HugrView> StabilizerDataflow<H> {
    fn new(hugr: &H, parent: H::Node) -> Self {
        // let mut in_cols: BiHashMap<OutgoingPort, usize> = BiHashMap::default();
        // let mut frontier_cols: BiHashMap<(H::Node, IncomingPort), usize> = BiHashMap::default();
        let mut q_ind_map: BiHashMap<DataflowPoint<H::Node>, usize> = BiHashMap::default();
        let mut n_in_qubits = 0;
        let inp = hugr.children(parent).filter(|n| matches!(hugr.get_optype(*n), OpType::Input(_))).exactly_one().ok().unwrap();
        for (out, out_type) in hugr.out_value_types(inp) {
            if out_type == qb_t() {
                // in_cols.insert(out, 2*n_in_qubits);
                q_ind_map.insert(DataflowPoint::Input(out), 2*n_in_qubits);
                let (next, next_p) = hugr.single_linked_input(inp, out).unwrap();
                // frontier_cols.insert((next, next_p), 2*n_in_qubits + 1);
                q_ind_map.insert(DataflowPoint::Frontier(next, next_p), 2*n_in_qubits + 1);
                n_in_qubits = n_in_qubits + 1;
            }
        }
        let mut tab = SymplecticTableau::new(2*n_in_qubits);
        // For each input wire, add an identity operation to the tableau
        for i in 0..n_in_qubits {
            let mut ii = BitVector::new(2*n_in_qubits);
            ii.xor_bit(2*i);
            ii.xor_bit(2*i + 1);
            // XX
            tab.add_stab(BitVector::new(2*n_in_qubits), ii.clone(), false);
            // ZZ
            tab.add_stab(ii, BitVector::new(2*n_in_qubits), false);
        }
        Self{
            tab: tab,
            // in_cols: in_cols,
            // out_cols: BiHashMap::default(),
            // frontier_cols: frontier_cols,
            // internal_in_cols: BiHashMap::default(),
            // internal_out_cols: BiHashMap::default(),
            // nested_in_cols: BiHashMap::default(),
            // nested_out_cols: BiHashMap::default(),
            q_index_map: q_ind_map,
            nested_analysis: HashMap::default(),
        }
    }

    fn remove_non_io_qubits(&mut self) {
        // Project out non-IO columns
        let non_ios: Vec<usize> = self.q_index_map.iter().map_while(|(dfp, q)| if matches!(dfp, DataflowPoint::Input(_)) || matches!(dfp, DataflowPoint::Output(_)) { None } else { Some(*q) }).collect_vec();
        // let non_ios = chain!(
        //     self.internal_in_cols.right_values(),
        //     self.internal_out_cols.right_values(),
        //     self.nested_in_cols.right_values(),
        //     self.nested_out_cols.right_values(),
        // );
        let project_cols = chain!(
            non_ios.iter().map(|i| (*i, PauliXZ::X)),
            non_ios.iter().map(|i| (*i, PauliXZ::Z)),
        ).collect_vec();
        self.tab.project(&project_cols);
        for i in non_ios {
            let moved_qb = self.tab.delete_qubit(i);
            match moved_qb {
                Some(mq) => {
                    // If the removal caused an IO qubit to have changed index, update it
                    let (removed_dfp, _) = self.q_index_map.remove_by_right(&mq).unwrap();
                    if matches!(removed_dfp, DataflowPoint::Input(_)) || matches!(removed_dfp, DataflowPoint::Output(_)) {
                        self.q_index_map.insert(removed_dfp, i);
                    }
                }
                None => {
                    self.q_index_map.remove_by_right(&i);
                }
            }
            // if moved_qb.is_some() {
            //     let removed_res: Option<(OutgoingPort, usize)> = self.in_cols.remove_by_right(&moved_qb.unwrap());
            //     match removed_res {
            //         Some((port, _)) => {
            //             self.in_cols.insert(port, *i);
            //         }
            //         None => {
            //             let removed_res = self.out_cols.remove_by_right(&moved_qb.unwrap());
            //             match removed_res {
            //                 Some((port, _)) => {
            //                     self.out_cols.insert(port, *i);
            //                 }
            //                 _ => {}
            //             }
            //         }
            //     }
            // }
        }
        // self.internal_in_cols.clear();
        // self.internal_out_cols.clear();
        // self.nested_in_cols.clear();
        // self.nested_out_cols.clear();
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
                            let mut fun_analysis = StabilizerDataflow::run_dfg(hugr, fun_def_node, fun_op);
                            fun_analysis.remove_non_io_qubits();
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
                    // let to_move: Vec<(DataflowPoint<H::Node>, usize)> = analysis.q_index_map.iter().filter(|(dfp, _)| matches!(dfp, DataflowPoint::Frontier(_, _))).map(|(dfp, q)| (dfp.clone(), *q)).collect_vec();
                    let to_move: Vec<(DataflowPoint<H::Node>, usize)> = analysis.q_index_map.iter().map_while(|(dfp, q)| if matches!(dfp, DataflowPoint::Frontier(_, _)) { Some((dfp.clone(), *q)) } else { None }).collect_vec();
                    for (dfp, q) in to_move {
                        if let DataflowPoint::Frontier(_node, port) = dfp {
                            analysis.q_index_map.remove_by_right(&q);
                            analysis.q_index_map.insert(DataflowPoint::Output(port), q);
                        }
                    }
                    // for ((_, port), col) in analysis.frontier_cols.iter() {
                    //     analysis.out_cols.insert(*port, *col);
                    // }
                    // analysis.frontier_cols = BiHashMap::default();
                }
                _ => {
                    analysis.apply_opaque(hugr, node)
                }
            }
        }
        analysis
    }

    fn run_conditional(hugr: &H, node: H::Node, fun_op: &FunctionOpacity) -> StabilizerDataflow<H> {
        // Assume no information is passed about Qubits within the Sum types, so our summary only incorporates the Qubits in the other args
        let cond = hugr.get_optype(node).as_conditional().unwrap();
        let sig = cond.signature();
        // Determins consistent column indexing for inputs and outputs
        // let mut unified_in_cols : BiHashMap<OutgoingPort, usize> = BiHashMap::default();
        let mut unified_q_index : BiHashMap<DataflowPoint<H::Node>, usize> = BiHashMap::default();
        let mut n_unified_qbs = 0;
        for in_port in sig.input_ports() {
            if *sig.in_port_type(in_port).unwrap() == qb_t() {
                // unified_in_cols.insert(OutgoingPort::from(in_port.index()), n_unified_qbs);
                unified_q_index.insert(DataflowPoint::Input(OutgoingPort::from(in_port.index())), n_unified_qbs);
                n_unified_qbs = n_unified_qbs + 1;
            }
        }
        // let mut unified_out_cols : BiHashMap<IncomingPort, usize> = BiHashMap::default();
        for out_port in sig.output_ports() {
            if *sig.out_port_type(out_port).unwrap() == qb_t() {
                // unified_out_cols.insert(IncomingPort::from(out_port.index()), n_unified_qbs);
                unified_q_index.insert(DataflowPoint::Output(IncomingPort::from(out_port.index())), n_unified_qbs);
                n_unified_qbs = n_unified_qbs + 1;
            }
        }
        let mut summary: Option<StabilizerDataflow<H>> = None;
        for (cond_i, cond_node) in hugr.children(node).enumerate() {
            let analysis = StabilizerDataflow::run_dfg(hugr, cond_node, fun_op);
            let mut projected_tab = analysis.tab.clone();
            // Number of ports from the condition row; given port p on input, corresponds to IncomingPort::from(p + 1 - cond_len) to the Conditional
            let cond_len = cond.sum_rows.get(cond_i).unwrap().len();
            // Project out non-IO columns (including any qubits from the condition row)
            let non_ios: Vec<usize> = analysis.q_index_map.iter().map_while(|(dfp, q)| {
                match dfp {
                    DataflowPoint::Input(port) => {
                        if port.index() < cond_len { Some(*q) }
                        else { None }
                    }
                    DataflowPoint::Output(_) => { None }
                    _ => { Some(*q) }
                }
            }).collect_vec();
            // let non_ios = chain!(
            //     analysis.in_cols.right_values().filter(|i| **i >= cond_len),
            //     analysis.internal_in_cols.right_values(),
            //     analysis.internal_out_cols.right_values(),
            //     analysis.nested_in_cols.right_values(),
            //     analysis.nested_out_cols.right_values(),
            // );
            let project_cols: Vec<(usize, PauliXZ)> = chain!(
                non_ios.iter().map(|i| (*i, PauliXZ::X)),
                non_ios.iter().map(|i| (*i, PauliXZ::Z)),
            ).collect_vec();
            projected_tab.project(&project_cols);
            // Rebuild projected_tab with the column order given by unified_X_cols
            let mut unified_order_tab = SymplecticTableau::new(n_unified_qbs);
            for i in 0..projected_tab.nb_stabs {
                let mut z = BitVector::new(n_unified_qbs);
                let mut x = BitVector::new(n_unified_qbs);
                // for (port, col) in &unified_in_cols {
                //     let old_col = analysis.in_cols.get_by_left(&OutgoingPort::from(port.index() + cond_len - 1)).unwrap();
                //     if projected_tab.z[i].get(*old_col) { z.xor_bit(*col); }
                //     if projected_tab.x[i].get(*old_col) { x.xor_bit(*col); }
                // }
                // for (port, col) in &unified_out_cols {
                //     let old_col = analysis.out_cols.get_by_left(&port).unwrap();
                //     if projected_tab.z[i].get(*old_col) { z.xor_bit(*col); }
                //     if projected_tab.x[i].get(*old_col) { x.xor_bit(*col); }
                // }
                for (dfp, col) in &unified_q_index {
                    match dfp {
                        DataflowPoint::Input(port) => {
                            let old_col = analysis.q_index_map.get_by_left(&DataflowPoint::Input(OutgoingPort::from(port.index() + cond_len - 1))).unwrap();
                            if projected_tab.z[i].get(*old_col) { z.xor_bit(*col); }
                            if projected_tab.x[i].get(*old_col) { x.xor_bit(*col); }
                        }
                        DataflowPoint::Output(port) => {
                            let old_col = analysis.q_index_map.get_by_left(&DataflowPoint::Output(*port)).unwrap();
                            if projected_tab.z[i].get(*old_col) { z.xor_bit(*col); }
                            if projected_tab.x[i].get(*old_col) { x.xor_bit(*col); }
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
                    summ.nested_analysis.insert(cond_node, analysis);
                }
                None => {
                    summary = Some(StabilizerDataflow {
                        tab: unified_order_tab,
                        // in_cols: unified_in_cols.clone(),
                        // out_cols: unified_out_cols.clone(),
                        // frontier_cols: BiHashMap::default(),
                        // internal_in_cols: BiHashMap::default(),
                        // internal_out_cols: BiHashMap::default(),
                        // nested_in_cols: BiHashMap::default(),
                        // nested_out_cols: BiHashMap::default(),
                        q_index_map: unified_q_index.clone(),
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
            tab: SymplecticTableau::new(0),
            // in_cols: BiHashMap::default(),
            // out_cols: BiHashMap::default(),
            // frontier_cols: BiHashMap::default(),
            // internal_in_cols: BiHashMap::default(),
            // internal_out_cols: BiHashMap::default(),
            // nested_in_cols: BiHashMap::default(),
            // nested_out_cols: BiHashMap::default(),
            q_index_map: BiHashMap::default(),
            nested_analysis: HashMap::default()
        };
        let tl = hugr.get_optype(node).as_tail_loop().unwrap();
        // tl.just_inputs only appear in final signature within a Sum, so qubits there will be projected away
        // tl.just_outputs do appear in the final signature, but we will not have any information about the qubits there
        for (out_port, out_type) in tl.just_outputs.iter().enumerate() {
            if *out_type == qb_t() {
                let new_col = analysis.tab.add_qubits(1);
                // analysis.out_cols.insert(IncomingPort::from(out_port), new_col);
                analysis.q_index_map.insert(DataflowPoint::Output(IncomingPort::from(out_port)), new_col);
            }
        }
        // tl.rest appear in the final input signature from port 1 onwards and in the output signature from port (tl.just_outputs.len()) onwards
        for (port_index, port_type) in tl.rest.iter().enumerate() {
            if *port_type == qb_t() {
                let in_col = analysis.tab.add_qubits(2);
                // analysis.in_cols.insert(OutgoingPort::from(port_index + 1), in_col);
                analysis.q_index_map.insert(DataflowPoint::Input(OutgoingPort::from(port_index + 1)), in_col);
                let out_col = in_col + 1;
                // analysis.out_cols.insert(IncomingPort::from(port_index + tl.just_outputs.len()), out_col);
                analysis.q_index_map.insert(DataflowPoint::Output(IncomingPort::from(port_index + tl.just_outputs.len())), out_col);
                // Add rows for identity in_col--out_col
                let mut ii = BitVector::new(analysis.tab.nb_qubits);
                ii.xor_bit(in_col);
                ii.xor_bit(out_col);
                analysis.tab.add_stab(BitVector::new(analysis.tab.nb_qubits), ii.clone(), false);
                analysis.tab.add_stab(ii, BitVector::new(analysis.tab.nb_qubits), false);
            }
        }
        let mut projected_tab = child_analysis.tab.clone();
        // Project out non-IO columns and those not shared by input and output
        let non_ios: Vec<usize> = analysis.q_index_map.iter().map_while(|(dfp, q)| {
            match dfp {
                DataflowPoint::Input(port) => {
                    if port.index() == 0 { Some(*q) }
                    else { None }
                }
                DataflowPoint::Output(port) => {
                    if port.index() < tl.just_outputs.len() { Some(*q) }
                    else { None }
                }
                _ => { Some(*q) }
            }
        }).collect_vec();
        let project_cols: Vec<(usize, PauliXZ)> = chain!(
            non_ios.iter().map(|i| (*i, PauliXZ::X)),
            non_ios.iter().map(|i| (*i, PauliXZ::Z)),
        ).collect_vec();
        projected_tab.project(&project_cols);
        // Rebuild projected_tab with the colum order given by analysis.q_index_map
        let mut unified_order_tab = SymplecticTableau::new(analysis.q_index_map.len());
        for i in 0..projected_tab.nb_stabs {
            let mut z = BitVector::new(analysis.q_index_map.len());
            let mut x = BitVector::new(analysis.q_index_map.len());
            for (dfp, new_col) in analysis.q_index_map.iter() {
                match dfp {
                    DataflowPoint::Input(port) => {
                        if port.index() != 0 {
                            let old_col = child_analysis.q_index_map.get_by_left(&DataflowPoint::Input(OutgoingPort::from(port.index() - 1))).unwrap();
                            if projected_tab.z[i].get(*old_col) { z.xor_bit(*new_col) }
                            if projected_tab.x[i].get(*old_col) { x.xor_bit(*new_col) }
                        }
                    }
                    DataflowPoint::Output(port) => {
                        if port.index() >= tl.just_outputs.len() {
                            let old_col = child_analysis.q_index_map.get_by_left(&DataflowPoint::Output(IncomingPort::from(port.index() - tl.just_outputs.len()))).unwrap();
                            if projected_tab.z[i].get(*old_col) { z.xor_bit(*new_col) }
                            if projected_tab.x[i].get(*old_col) { x.xor_bit(*new_col) }
                        }
                    }
                    _ => {}
                }
            }
            unified_order_tab.add_stab(z, x, projected_tab.signs.get(i));
        }
        // Compute join of tabs
        analysis.tab = SymplecticTableau::join(&unified_order_tab, &analysis.tab);
        analysis.nested_analysis.insert(child_node, child_analysis);
        analysis
    }

    fn apply_quantum_gate(&mut self, hugr : &H, node: H::Node, op: TketOp) {
        match op {
            TketOp::H => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_h(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::CX => {
                // let (_, col0) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                // let (_, col1) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(1))).unwrap();
                let (_, col0) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let (_, col1) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(1))).unwrap();
                self.tab.append_cx(col0, col1);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
                let (next_node0, next_port0) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node0, next_port0), col0);
                let (next_node1, next_port1) = hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node1, next_port1), col1);
            }
            TketOp::CY => {
                // let (_, col0) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                // let (_, col1) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(1))).unwrap();
                let (_, col0) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let (_, col1) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(1))).unwrap();
                self.tab.append_s(col1);
                self.tab.append_z(col1);
                self.tab.append_cx(col0, col1);
                self.tab.append_s(col1);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
                let (next_node0, next_port0) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node0, next_port0), col0);
                let (next_node1, next_port1) = hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node1, next_port1), col1);
            }
            TketOp::CZ => {
                // let (_, col0) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                // let (_, col1) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(1))).unwrap();
                let (_, col0) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let (_, col1) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(1))).unwrap();
                self.tab.append_cz(col0, col1);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col0);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col1);
                let (next_node0, next_port0) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node0, next_port0), col0);
                let (next_node1, next_port1) = hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node1, next_port1), col1);
            }
            TketOp::CRz => {
                // let (_, col_in0) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                // let (_, col_in1) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(1))).unwrap();
                let (_, col_in0) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let (_, col_in1) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(1))).unwrap();
                let col_out0: usize = self.tab.add_qubits(4);
                let col_out1: usize = col_out0 + 1;
                let col_front0: usize = col_out0 + 2;
                let col_front1: usize = col_out0 + 3;
                // Add rows for identities col_out0/1--col_front0/1
                let mut out_front_0 = BitVector::new(self.tab.nb_qubits);
                out_front_0.xor_bit(col_out0);
                out_front_0.xor_bit(col_front0);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front_0.clone(), false);
                self.tab.add_stab(out_front_0, BitVector::new(self.tab.nb_qubits), false);
                let mut out_front_1 = BitVector::new(self.tab.nb_qubits);
                out_front_1.xor_bit(col_out1);
                out_front_1.xor_bit(col_front1);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front_1.clone(), false);
                self.tab.add_stab(out_front_1, BitVector::new(self.tab.nb_qubits), false);
                // Add rows for ZZ over col_in0/1--col_out0/1 and project to commuting
                let mut in_out_0 = BitVector::new(self.tab.nb_qubits);
                in_out_0.xor_bit(col_in0);
                in_out_0.xor_bit(col_out0);
                let mut in_out_1 = BitVector::new(self.tab.nb_qubits);
                in_out_1.xor_bit(col_in1);
                in_out_1.xor_bit(col_out1);
                self.tab.project_commuting_with(&in_out_0, &BitVector::new(self.tab.nb_qubits));
                self.tab.project_commuting_with(&in_out_1, &BitVector::new(self.tab.nb_qubits));
                self.tab.add_stab(in_out_0, BitVector::new(self.tab.nb_qubits), false);
                self.tab.add_stab(in_out_1, BitVector::new(self.tab.nb_qubits), false);
                // self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in0);
                // self.internal_in_cols.insert((node, IncomingPort::from(1)), col_in1);
                // self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out0);
                // self.internal_out_cols.insert((node, OutgoingPort::from(1)), col_out1);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(0)), col_in0);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(1)), col_in1);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(0)), col_out0);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(1)), col_out1);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front0);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col_front1);
                let (next_node0, next_port0) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node0, next_port0), col_front0);
                let (next_node1, next_port1) = hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node1, next_port1), col_front1);
            }
            TketOp::T | TketOp::Tdg | TketOp::Rz | TketOp::Measure => {
                // let (_, col_in) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col_in) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_qubits(2);
                let col_front: usize = col_out + 1;
                // Add rows for identity col_out--col_front
                let mut out_front = BitVector::new(self.tab.nb_qubits);
                out_front.xor_bit(col_out);
                out_front.xor_bit(col_front);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front.clone(), false);
                self.tab.add_stab(out_front, BitVector::new(self.tab.nb_qubits), false);
                // Add row for ZZ over col_in--col_out and project to commuting
                let mut in_out = BitVector::new(self.tab.nb_qubits);
                in_out.xor_bit(col_in);
                in_out.xor_bit(col_out);
                self.tab.project_commuting_with(&in_out, &BitVector::new(self.tab.nb_qubits));
                self.tab.add_stab(in_out, BitVector::new(self.tab.nb_qubits), false);
                // self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                // self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(0)), col_in);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(0)), col_out);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col_front);
            }
            TketOp::S => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_s(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::Sdg => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_s(col);
                self.tab.append_z(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::X => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_x(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::Y => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_x(col);
                self.tab.append_z(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::Z => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_z(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::Rx => {
                // let (_, col_in) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col_in) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_qubits(2);
                let col_front: usize = col_out + 1;
                // Add rows for identity col_out--col_front
                let mut out_front = BitVector::new(self.tab.nb_qubits);
                out_front.xor_bit(col_out);
                out_front.xor_bit(col_front);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front.clone(), false);
                self.tab.add_stab(out_front, BitVector::new(self.tab.nb_qubits), false);
                // Add row for XX over col_in--col_out and project to commuting
                let mut in_out = BitVector::new(self.tab.nb_qubits);
                in_out.xor_bit(col_in);
                in_out.xor_bit(col_out);
                self.tab.project_commuting_with(&BitVector::new(self.tab.nb_qubits), &in_out);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), in_out, false);
                // self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                // self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(0)), col_in);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(0)), col_out);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col_front);
            }
            TketOp::Ry => {
                // let (_, col_in) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col_in) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let col_out: usize = self.tab.add_qubits(2);
                let col_front: usize = col_out;
                // Add rows for identity col_out--col_front
                let mut out_front = BitVector::new(self.tab.nb_qubits);
                out_front.xor_bit(col_out);
                out_front.xor_bit(col_front);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front.clone(), false);
                self.tab.add_stab(out_front, BitVector::new(self.tab.nb_qubits), false);
                // Add row for -YY (negative because of the partial transpose under the Choi isomorphism) over col_in--col_out and project to commuting
                let mut in_out = BitVector::new(self.tab.nb_qubits);
                in_out.xor_bit(col_in);
                in_out.xor_bit(col_out);
                self.tab.project_commuting_with(&in_out, &in_out);
                self.tab.add_stab(in_out.clone(), in_out, true);
                // self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                // self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(0)), col_in);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(0)), col_out);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col_front);
            }
            TketOp::Toffoli => {
                // let (_, col_in0) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                // let (_, col_in1) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(1))).unwrap();
                // let (_, col_in2) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(2))).unwrap();
                let (_, col_in0) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                let (_, col_in1) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(1))).unwrap();
                let (_, col_in2) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(2))).unwrap();
                let col_out0: usize = self.tab.add_qubits(6);
                let col_out1: usize = col_out0 + 1;
                let col_out2: usize = col_out0 + 2;
                let col_front0: usize = col_out0 + 3;
                let col_front1: usize = col_out0 + 4;
                let col_front2: usize = col_out0 + 5;
                // Add rows for identities col_out0/1/2--col_front0/1/2
                let mut out_front_0 = BitVector::new(self.tab.nb_qubits);
                out_front_0.xor_bit(col_out0);
                out_front_0.xor_bit(col_front0);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front_0.clone(), false);
                self.tab.add_stab(out_front_0, BitVector::new(self.tab.nb_qubits), false);
                let mut out_front_1 = BitVector::new(self.tab.nb_qubits);
                out_front_1.xor_bit(col_out1);
                out_front_1.xor_bit(col_front1);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front_1.clone(), false);
                self.tab.add_stab(out_front_1, BitVector::new(self.tab.nb_qubits), false);
                let mut out_front_2 = BitVector::new(self.tab.nb_qubits);
                out_front_2.xor_bit(col_out2);
                out_front_2.xor_bit(col_front2);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front_2.clone(), false);
                self.tab.add_stab(out_front_2, BitVector::new(self.tab.nb_qubits), false);
                // Add rows for ZZ/ZZ/XX over col_in0/1/2--col_out0/1/2 and project to commuting
                let mut in_out_0 = BitVector::new(self.tab.nb_qubits);
                in_out_0.xor_bit(col_in0);
                in_out_0.xor_bit(col_out0);
                let mut in_out_1 = BitVector::new(self.tab.nb_qubits);
                in_out_1.xor_bit(col_in1);
                in_out_1.xor_bit(col_out1);
                let mut in_out_2 = BitVector::new(self.tab.nb_qubits);
                in_out_2.xor_bit(col_in2);
                in_out_2.xor_bit(col_out2);
                self.tab.project_commuting_with(&in_out_0, &BitVector::new(self.tab.nb_qubits));
                self.tab.project_commuting_with(&in_out_1, &BitVector::new(self.tab.nb_qubits));
                self.tab.project_commuting_with( &BitVector::new(self.tab.nb_qubits), &in_out_2);
                self.tab.add_stab(in_out_0, BitVector::new(self.tab.nb_qubits), false);
                self.tab.add_stab(in_out_1, BitVector::new(self.tab.nb_qubits), false);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), in_out_2, false);
                // self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in0);
                // self.internal_in_cols.insert((node, IncomingPort::from(1)), col_in1);
                // self.internal_in_cols.insert((node, IncomingPort::from(2)), col_in2);
                // self.internal_out_cols.insert((node, OutgoingPort::from(0)), col_out0);
                // self.internal_out_cols.insert((node, OutgoingPort::from(1)), col_out1);
                // self.internal_out_cols.insert((node, OutgoingPort::from(2)), col_out2);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(0)), col_in0);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(1)), col_in1);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(2)), col_in2);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(0)), col_out0);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(1)), col_out1);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, OutgoingPort::from(2)), col_out2);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front0);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap(), col_front1);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(2)).unwrap(), col_front2);
                let (next_node0, next_port0) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node0, next_port0), col_front0);
                let (next_node1, next_port1) = hugr.single_linked_input(node, OutgoingPort::from(1)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node1, next_port1), col_front1);
                let (next_node2, next_port2) = hugr.single_linked_input(node, OutgoingPort::from(2)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node2, next_port2), col_front2);
            }
            TketOp::MeasureFree => {
                // let (_, col_in) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col_in) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                // self.internal_in_cols.insert((node, IncomingPort::from(0)), col_in);
                self.q_index_map.insert(DataflowPoint::InternalIn(node, IncomingPort::from(0)), col_in);
            }
            TketOp::QAlloc => {
                let col_front: usize = self.tab.add_qubits(1);
                // Add row for Z over col_front
                let mut front_bv = BitVector::new(self.tab.nb_qubits);
                front_bv.xor_bit(col_front);
                self.tab.add_stab(front_bv, BitVector::new(self.tab.nb_qubits), false);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_front);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col_front);
            }
            TketOp::QFree => {
                // let (_, col_in) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col_in) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                // Project out non-commuting rows and remove column from tableau
                self.tab.project(&vec![(col_in, PauliXZ::X), (col_in, PauliXZ::Z)]);
                let moved_qb = self.tab.delete_qubit(col_in);
                if moved_qb.is_some() {
                    let mq = moved_qb.unwrap();
                    let (removed_dfp, _) = self.q_index_map.remove_by_right(&mq).unwrap();
                    self.q_index_map.insert(removed_dfp, col_in);
                }
            }
            TketOp::Reset => {
                // let (_, col_in) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col_in) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                // Project out non-commuting rows
                self.tab.project(&vec![(col_in, PauliXZ::X), (col_in, PauliXZ::Z)]);
                // Reuse col_in for the output qubit
                // Add row for Z over col_in
                let mut bv = BitVector::new(self.tab.nb_qubits);
                bv.xor_bit(col_in);
                self.tab.add_stab(bv, BitVector::new(self.tab.nb_qubits), false);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col_in);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col_in);
            }
            TketOp::V => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_v(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
            }
            TketOp::Vdg => {
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, IncomingPort::from(0))).unwrap();
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, IncomingPort::from(0))).unwrap();
                self.tab.append_v(col);
                self.tab.append_x(col);
                // self.frontier_cols.insert(hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap(), col);
                let (next_node, next_port) = hugr.single_linked_input(node, OutgoingPort::from(0)).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col);
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
                // let (_, col) = self.frontier_cols.remove_by_left(&(node, p)).unwrap();
                // self.internal_in_cols.insert((node, p), col);
                let (_, col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, p)).unwrap();
                self.q_index_map.insert(DataflowPoint::InternalIn(node, p), col);
            }
        }
        // For each Qubit output, create a pair of columns with the identity for internal_out_cols and frontier_cols
        for (p, t) in hugr.out_value_types(node) {
            if t == qb_t() {
                let col_out = self.tab.add_qubits(2);
                let col_front = col_out + 1;
                // Add rows for identity col_out--col_front
                let mut out_front = BitVector::new(self.tab.nb_qubits);
                out_front.xor_bit(col_out);
                out_front.xor_bit(col_front);
                self.tab.add_stab(BitVector::new(self.tab.nb_qubits), out_front.clone(), false);
                self.tab.add_stab(out_front, BitVector::new(self.tab.nb_qubits), false);
                // self.internal_out_cols.insert((node, p), col_out);
                self.q_index_map.insert(DataflowPoint::InternalOut(node, p), col_out);
                // self.frontier_cols.insert(hugr.single_linked_input(node, p).unwrap(), col_front);
                let (next_node, next_port) = hugr.single_linked_input(node, p).unwrap();
                self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), col_front);
            }
        }
    }

    /// Suppose we have already recursively calculated a StabilizerDataflow for node and stored it in nested_analysis; performs sequential composition to append it to the appropriate qubits here
    fn apply_analysis(&mut self, hugr: &H, node: H::Node) {
        let node_analysis : &StabilizerDataflow<H> = self.nested_analysis.get(&node).unwrap();
        let old_n_qbs = self.tab.nb_qubits;
        let n_added_qbs = node_analysis.tab.nb_qubits;
        self.tab.add_qubits(n_added_qbs);
        for (dfp, col) in node_analysis.q_index_map.iter() {
            match dfp {
                DataflowPoint::Input(p) => {
                    self.q_index_map.insert(DataflowPoint::NestedIn(node, *p), col + old_n_qbs);
                }
                DataflowPoint::Output(p) => {
                    self.q_index_map.insert(DataflowPoint::NestedOut(node, *p), col + old_n_qbs);
                }
                _ => {
                    // It is decided at analysis construction whether or not we project to IO or keep internals, so always copy any internals that have remained (e.g. for inlining function calls)
                    self.q_index_map.insert(dfp.clone(), col + old_n_qbs);
                }
            }
        }
        // for (port, col) in node_analysis.in_cols.iter() {
        //     self.nested_in_cols.insert((node, *port), *col + old_n_qbs);
        // }
        // for (port, col) in node_analysis.out_cols.iter() {
        //     self.nested_out_cols.insert((node, *port), *col + old_n_qbs);
        // }
        // for (node_port, col) in node_analysis.internal_in_cols.iter() {
        //     self.internal_in_cols.insert(*node_port, *col + old_n_qbs);
        // }
        // for (node_port, col) in node_analysis.internal_out_cols.iter() {
        //     self.internal_out_cols.insert(*node_port, *col + old_n_qbs);
        // }
        // for (node_port, col) in node_analysis.nested_in_cols.iter() {
        //     self.nested_in_cols.insert(*node_port, *col + old_n_qbs);
        // }
        // for (node_port, col) in node_analysis.nested_out_cols.iter() {
        //     self.nested_out_cols.insert(*node_port, *col + old_n_qbs);
        // }
        for i in 0..node_analysis.tab.nb_stabs {
            let mut new_z = BitVector::new(old_n_qbs);
            new_z.extend_vec(node_analysis.tab.z[i].get_boolean_vec(), old_n_qbs);
            let mut new_x = BitVector::new(old_n_qbs);
            new_x.extend_vec(node_analysis.tab.x[i].get_boolean_vec(), old_n_qbs);
            self.tab.add_stab(new_z, new_x, node_analysis.tab.signs.get(i));
        }
        for port in hugr.node_inputs(node) {
            let out_port = OutgoingPort::from(port.index());
            // let (_, internal_col) = self.frontier_cols.remove_by_left(&(node, port)).unwrap();
            let (_, internal_col) = self.q_index_map.remove_by_left(&DataflowPoint::Frontier(node, port)).unwrap();
            // self.internal_in_cols.insert((node, port), internal_col);
            self.q_index_map.insert(DataflowPoint::InternalIn(node, port), internal_col);
            // let nested_col = self.nested_in_cols.get_by_left(&(node, out_port));
            let nested_col = self.q_index_map.get_by_left(&DataflowPoint::NestedIn(node, out_port)).unwrap();
            // Project ZZ and XX to compose nested_col and internal_col
            let mut nested_internal = BitVector::new(self.tab.nb_qubits);
            nested_internal.xor_bit(*nested_col);
            nested_internal.xor_bit(internal_col);
            self.tab.project_commuting_with(&BitVector::new(self.tab.nb_qubits), &nested_internal);
            self.tab.project_commuting_with(&nested_internal, &BitVector::new(self.tab.nb_qubits));
            self.tab.add_stab(BitVector::new(self.tab.nb_qubits), nested_internal.clone(), false);
            self.tab.add_stab(nested_internal, BitVector::new(self.tab.nb_qubits), false);
        }
        for port in hugr.node_outputs(node) {
            let in_port = IncomingPort::from(port.index());
            // let nested_col = self.nested_out_cols.get_by_left(&(node, in_port));
            let nested_col = *self.q_index_map.get_by_left(&DataflowPoint::NestedOut(node, in_port)).unwrap();
            let internal_col = self.tab.add_qubits(2);
            // self.internal_out_cols.insert((node, port), internal_col);
            self.q_index_map.insert(DataflowPoint::InternalOut(node, port), internal_col);
            let front_col = internal_col + 1;
            // self.frontier_cols.insert(hugr.single_linked_input(node, port).unwrap(), front_col);
            let (next_node, next_port) = hugr.single_linked_input(node, port).unwrap();
            self.q_index_map.insert(DataflowPoint::Frontier(next_node, next_port), front_col);
            // Add rows for identity internal_col--front_col
            let mut internal_front = BitVector::new(self.tab.nb_qubits);
            internal_front.xor_bit(internal_col);
            internal_front.xor_bit(front_col);
            self.tab.add_stab(BitVector::new(self.tab.nb_qubits), internal_front.clone(), false);
            self.tab.add_stab(internal_front, BitVector::new(self.tab.nb_qubits), false);
            // Project ZZ and XX to compose nested_col and internal_col
            let mut nested_internal = BitVector::new(self.tab.nb_qubits);
            nested_internal.xor_bit(nested_col);
            nested_internal.xor_bit(internal_col);
            self.tab.project_commuting_with(&BitVector::new(self.tab.nb_qubits), &nested_internal);
            self.tab.project_commuting_with(&nested_internal, &BitVector::new(self.tab.nb_qubits));
            self.tab.add_stab(BitVector::new(self.tab.nb_qubits), nested_internal.clone(), false);
            self.tab.add_stab(nested_internal, BitVector::new(self.tab.nb_qubits), false);
        }
    }

}

#[cfg(test)]
mod test {
    use hugr::{builder::{endo_sig, ConditionalBuilder, DFGBuilder, Dataflow, DataflowHugr, DataflowSubContainer, HugrBuilder, SubContainer}, extension::prelude::{bool_t, qb_t, usize_t}, ops::{handle::NodeHandle, OpType, OpaqueOp}, type_row, types::Signature, HugrView, IncomingPort, OutgoingPort};
    use tket::TketOp;

    use crate::{bit_vector::BitVector, pauli_product::PauliProduct, stabilizer_dataflow::{FunctionOpacity, StabilizerDataflow}};


    #[test]
    fn test_empty_analysis() {
        let builder = DFGBuilder::new(endo_sig(vec![])).unwrap();
        let hugr = builder.finish_hugr().unwrap();
        let analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
        assert_eq!(analysis.tab.nb_qubits, 0);
        assert_eq!(analysis.tab.nb_stabs, 0);
    }

    // #[test]
    // fn test_identity_analysis() {
    //     // Add an extra integer input to make sure we only track the qubits
    //     let builder = DFGBuilder::new(endo_sig(vec![usize_t(), qb_t(), qb_t()])).unwrap();
    //     let [_, qb0, qb1] = builder.input_wires_arr();
    //     let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 4);
    //     assert_eq!(analysis.tab.nb_stabs, 4);
    //     // Check the right ports are stored for tracking the qubits
    //     assert_eq!(analysis.in_cols.len(), 2);
    //     assert_eq!(analysis.out_cols.len(), 2);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(1)).unwrap(), 0);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(0)).unwrap(), 1);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(2)).unwrap(), 2);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(1)).unwrap(), 3);
    //     // Check that the rows correspond to the identity operations
    //     //TODO:: Reduce analysis.tab to row echelon form
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_boolean_vec(), vec![true, true, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_boolean_vec(), vec![false; 4]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(1).x.get_boolean_vec(), vec![false; 4]);
    //     assert_eq!(analysis.tab.stabs.get(1).z.get_boolean_vec(), vec![true, true, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(1).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(2).x.get_boolean_vec(), vec![false, false, true, true]);
    //     assert_eq!(analysis.tab.stabs.get(2).z.get_boolean_vec(), vec![false; 4]);
    //     assert_eq!(analysis.tab.stabs.get(2).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(3).x.get_boolean_vec(), vec![false; 4]);
    //     assert_eq!(analysis.tab.stabs.get(3).z.get_boolean_vec(), vec![false, false, true, true]);
    //     assert_eq!(analysis.tab.stabs.get(3).sign, false);
    // }

    // #[test]
    // fn test_bell_state() {
    //     let mut builder = DFGBuilder::new(Signature::new(vec![], vec![qb_t(), qb_t()])).unwrap();
    //     let [qb0] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::H, [qb0]).unwrap().outputs_arr();
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
    //     let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 2);
    //     assert_eq!(analysis.tab.nb_stabs, 2);
    //     // Check that the rows correspond to the Bell state stabilizers
    //     //TODO:: Reduce analysis.tab to row echelon form
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_boolean_vec(), vec![true; 2]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_boolean_vec(), vec![false; 2]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(1).x.get_boolean_vec(), vec![false; 2]);
    //     assert_eq!(analysis.tab.stabs.get(1).z.get_boolean_vec(), vec![true; 2]);
    //     assert_eq!(analysis.tab.stabs.get(1).sign, false);
    // }
    
    // #[test]
    // fn test_opaque() {
    //     let mut builder = DFGBuilder::new(Signature::new(vec![], vec![qb_t(), qb_t()])).unwrap();
    //     let [qb0] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::H, [qb0]).unwrap().outputs_arr();
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
    //     let op = OpaqueOp::new(
    //         "ext".try_into().unwrap(),
    //         "op",
    //         vec![],
    //         Signature::new_endo(vec![qb_t()])
    //     );
    //     let opaque_op = builder.add_dataflow_op(OpType::OpaqueOp(op), [qb1]).unwrap();
    //     let [qb1] = builder.add_dataflow_op(TketOp::H, [opaque_op.out_wire(0)]).unwrap().outputs_arr();
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
    //     let hugr = builder.finish_hugr_with_outputs([qb0, qb1]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 4);
    //     assert_eq!(analysis.tab.nb_stabs, 2);
    //     // Reduce analysis.tab to row echelon form with qubit ordering [out0, op_in, op_out, out1]
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(0)).unwrap(), 0);
    //     assert_eq!(*analysis.internal_in_cols.get(&(opaque_op.node(), IncomingPort::from(0))).unwrap(), 1);
    //     assert_eq!(*analysis.internal_out_cols.get(&(opaque_op.node(), OutgoingPort::from(0))).unwrap(), 2);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(1)).unwrap(), 3);
    //     //TODO:: Row echelon
    //     // Check the rows
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_boolean_vec(), vec![true, false, true, true]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_boolean_vec(), vec![false; 4]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(1).x.get_boolean_vec(), vec![false; 4]);
    //     assert_eq!(analysis.tab.stabs.get(1).z.get_boolean_vec(), vec![true, false, true, false]);
    //     assert_eq!(analysis.tab.stabs.get(1).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(2).x.get_boolean_vec(), vec![false, true, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(2).z.get_boolean_vec(), vec![false, false, true, true]);
    //     assert_eq!(analysis.tab.stabs.get(2).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(3).x.get_boolean_vec(), vec![false, false, false, true]);
    //     assert_eq!(analysis.tab.stabs.get(3).z.get_boolean_vec(), vec![false, true, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(3).sign, false);
    // }

    // #[test]
    // fn test_clifford_gates() {
    //     // Need to cover H, CX, CY, CZ, S, Sdg, X, Y, Z, V, Vdg
    //     let mut builder = DFGBuilder::new(endo_sig(vec![qb_t(), qb_t(), qb_t()])).unwrap();
    //     let [qb0, qb1, qb2] = builder.input_wires_arr();
    //     // H;S;V;S = I
    //     let [qb0] = builder.add_dataflow_op(TketOp::H, [qb0]).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::S, [qb0]).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::V, [qb0]).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::S, [qb0]).unwrap().outputs_arr();
    //     // CX;IH;CZ;IH = II
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::H, [qb1]).unwrap().outputs_arr();
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CZ, [qb0, qb1]).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::H, [qb1]).unwrap().outputs_arr();
    //     // CX;IS;CY;ISdg = II
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::S, [qb1]).unwrap().outputs_arr();
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CY, [qb0, qb1]).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::Sdg, [qb1]).unwrap().outputs_arr();
    //     // S;Sdg = II
    //     let [qb0] = builder.add_dataflow_op(TketOp::S, [qb0]).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::Sdg, [qb0]).unwrap().outputs_arr();
    //     // V;Vdh = II
    //     let [qb1] = builder.add_dataflow_op(TketOp::V, [qb1]).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::Vdg, [qb1]).unwrap().outputs_arr();
    //     // Test Paulis explicitly
    //     let [qb0] = builder.add_dataflow_op(TketOp::Z, [qb0]).unwrap().outputs_arr();
    //     let [qb1] = builder.add_dataflow_op(TketOp::X, [qb1]).unwrap().outputs_arr();
    //     let [qb2] = builder.add_dataflow_op(TketOp::Y, [qb2]).unwrap().outputs_arr();
    //     let hugr = builder.finish_hugr_with_outputs([qb0, qb1, qb2]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 6);
    //     assert_eq!(analysis.tab.nb_stabs, 6);
    //     // Reduce analysis.tab to row echelon form with qubit ordering [in0, out0, in1, out1, in2, out2]
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(0)).unwrap(), 0);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(0)).unwrap(), 1);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(1)).unwrap(), 2);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(1)).unwrap(), 3);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(2)).unwrap(), 4);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(2)).unwrap(), 5);
    //     //TODO:: Row echelon
    //     // Check the rows
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_boolean_vec(), vec![true, true, false, false, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, true);
    //     assert_eq!(analysis.tab.stabs.get(1).x.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(1).z.get_boolean_vec(), vec![true, true, false, false, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(1).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(2).x.get_boolean_vec(), vec![false, false, true, true, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(2).z.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(2).sign, false);
    //     assert_eq!(analysis.tab.stabs.get(3).x.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(3).z.get_boolean_vec(), vec![false, false, true, true, false, false]);
    //     assert_eq!(analysis.tab.stabs.get(3).sign, true);
    //     assert_eq!(analysis.tab.stabs.get(4).x.get_boolean_vec(), vec![false, false, false, false, true, true]);
    //     assert_eq!(analysis.tab.stabs.get(4).z.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(4).sign, true);
    //     assert_eq!(analysis.tab.stabs.get(5).x.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(5).z.get_boolean_vec(), vec![false, false, false, false, true, true]);
    //     assert_eq!(analysis.tab.stabs.get(5).sign, true);
    // }

    // #[test]
    // fn test_nonclifford() {
    //     // Need to cover the separate logic for CRz, T/Tdg/Rz/Measure, Rx, Ry, Toffoli
    //     let mut builder = DFGBuilder::new(endo_sig(vec![qb_t(), qb_t(), qb_t()])).unwrap();
    //     let [qb0, qb1, qb2] = builder.input_wires_arr();
    //     let t = builder.add_dataflow_op(TketOp::T, [qb0]).unwrap();
    //     let tdg = builder.add_dataflow_op(TketOp::Tdg, [t.out_wire(0)]).unwrap();
    //     let rz = builder.add_dataflow_op(TketOp::Rz, [tdg.out_wire(0)]).unwrap();
    //     let meas = builder.add_dataflow_op(TketOp::Measure, [rz.out_wire(0)]).unwrap();
    //     let ry = builder.add_dataflow_op(TketOp::Ry, [qb1]).unwrap();
    //     let rx = builder.add_dataflow_op(TketOp::Rx, [qb2]).unwrap();
    //     let crz = builder.add_dataflow_op(TketOp::CRz, [meas.out_wire(0), ry.out_wire(0)]).unwrap();
    //     let toffoli = builder.add_dataflow_op(TketOp::Toffoli, [crz.out_wire(0), crz.out_wire(1), rx.out_wire(0)]).unwrap();
    //     let hugr = builder.finish_hugr_with_outputs(toffoli.outputs_arr::<3>()).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 28);
    //     assert_eq!(analysis.tab.nb_stabs, 28);
    //     // Reduce analysis.tab to row echelon form with qubit ordering:
    //     // [in0, t.in, in1, ry.in, in2, rx.in, t.out, tdg.in, ry.out, crz.in1, rx.out, toffoli.in2,
    //     // tdg.out, rz.in, rz.out, meas.in, meas.out, crz.in0, crz.out0, crz.out1, toffoli.in0, toffoli.in1
    //     // toffoli.out0, toffoli.out1, toffoli.out2, out0, out1, out2]
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(0)).unwrap(), 0);
    //     assert_eq!(*analysis.internal_in_cols.get(&(t.node(), IncomingPort::from(0))).unwrap(), 1);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(1)).unwrap(), 2);
    //     assert_eq!(*analysis.internal_in_cols.get(&(ry.node(), IncomingPort::from(0))).unwrap(), 3);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(2)).unwrap(), 4);
    //     assert_eq!(*analysis.internal_in_cols.get(&(rx.node(), IncomingPort::from(0))).unwrap(), 5);
    //     assert_eq!(*analysis.internal_out_cols.get(&(t.node(), OutgoingPort::from(0))).unwrap(), 6);
    //     assert_eq!(*analysis.internal_in_cols.get(&(tdg.node(), IncomingPort::from(0))).unwrap(), 7);
    //     assert_eq!(*analysis.internal_out_cols.get(&(ry.node(), OutgoingPort::from(0))).unwrap(), 8);
    //     assert_eq!(*analysis.internal_in_cols.get(&(crz.node(), IncomingPort::from(1))).unwrap(), 9);
    //     assert_eq!(*analysis.internal_out_cols.get(&(rx.node(), OutgoingPort::from(0))).unwrap(), 10);
    //     assert_eq!(*analysis.internal_in_cols.get(&(toffoli.node(), IncomingPort::from(2))).unwrap(), 11);
    //     assert_eq!(*analysis.internal_out_cols.get(&(tdg.node(), OutgoingPort::from(0))).unwrap(), 12);
    //     assert_eq!(*analysis.internal_in_cols.get(&(rz.node(), IncomingPort::from(0))).unwrap(), 13);
    //     assert_eq!(*analysis.internal_out_cols.get(&(rz.node(), OutgoingPort::from(0))).unwrap(), 14);
    //     assert_eq!(*analysis.internal_in_cols.get(&(meas.node(), IncomingPort::from(0))).unwrap(), 15);
    //     assert_eq!(*analysis.internal_out_cols.get(&(meas.node(), OutgoingPort::from(0))).unwrap(), 16);
    //     assert_eq!(*analysis.internal_in_cols.get(&(crz.node(), IncomingPort::from(0))).unwrap(), 17);
    //     assert_eq!(*analysis.internal_out_cols.get(&(crz.node(), OutgoingPort::from(0))).unwrap(), 18);
    //     assert_eq!(*analysis.internal_out_cols.get(&(crz.node(), OutgoingPort::from(1))).unwrap(), 19);
    //     assert_eq!(*analysis.internal_in_cols.get(&(toffoli.node(), IncomingPort::from(0))).unwrap(), 20);
    //     assert_eq!(*analysis.internal_in_cols.get(&(toffoli.node(), IncomingPort::from(1))).unwrap(), 21);
    //     assert_eq!(*analysis.internal_out_cols.get(&(toffoli.node(), OutgoingPort::from(0))).unwrap(), 22);
    //     assert_eq!(*analysis.internal_out_cols.get(&(toffoli.node(), OutgoingPort::from(1))).unwrap(), 23);
    //     assert_eq!(*analysis.internal_out_cols.get(&(toffoli.node(), OutgoingPort::from(2))).unwrap(), 24);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(0)).unwrap(), 25);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(1)).unwrap(), 26);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(2)).unwrap(), 27);
    //     //TODO:: Row echelon
    //     // Check the rows
    //     // Xin0
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_integer_vec(), vec![0b1100001100001111111010100100i128]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, false);
    //     // Zin0
    //     assert_eq!(analysis.tab.stabs.get(1).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(1).z.get_integer_vec(), vec![0b1000000000000000000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(1).sign, false);
    //     // Zt.in
    //     assert_eq!(analysis.tab.stabs.get(2).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(2).z.get_integer_vec(), vec![0b0100000000000000000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(2).sign, false);
    //     // Xin1
    //     assert_eq!(analysis.tab.stabs.get(3).x.get_integer_vec(), vec![0b0011000011000000000101010010i128]);
    //     assert_eq!(analysis.tab.stabs.get(3).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(3).sign, false);
    //     // Zin1
    //     assert_eq!(analysis.tab.stabs.get(4).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(4).z.get_integer_vec(), vec![0b0011000010000000000000000010i128]);
    //     assert_eq!(analysis.tab.stabs.get(4).sign, false);
    //     // Yry.in
    //     assert_eq!(analysis.tab.stabs.get(5).x.get_integer_vec(), vec![0b0001000001000000000101010010i128]);
    //     assert_eq!(analysis.tab.stabs.get(5).z.get_integer_vec(), vec![0b0001000001000000000000000000i128]);
    //     assert_eq!(analysis.tab.stabs.get(5).sign, false);
    //     // Xin2
    //     assert_eq!(analysis.tab.stabs.get(6).x.get_integer_vec(), vec![0b0000100000000000000000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(6).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(6).sign, false);
    //     // Zin2
    //     assert_eq!(analysis.tab.stabs.get(7).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(7).z.get_integer_vec(), vec![0b0000110000110000000000001001i128]);
    //     assert_eq!(analysis.tab.stabs.get(7).sign, false);
    //     // Xrx.in
    //     assert_eq!(analysis.tab.stabs.get(8).x.get_integer_vec(), vec![0b0000010000000000000000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(8).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(8).sign, false);
    //     // Zt.out
    //     assert_eq!(analysis.tab.stabs.get(9).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(9).z.get_integer_vec(), vec![0b0000001000000000000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(9).sign, false);
    //     // Ztdg.in
    //     assert_eq!(analysis.tab.stabs.get(10).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(10).z.get_integer_vec(), vec![0b0000000100000000000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(10).sign, false);
    //     // Yry.out
    //     assert_eq!(analysis.tab.stabs.get(11).x.get_integer_vec(), vec![0b0000000011000000000101010010i128]);
    //     assert_eq!(analysis.tab.stabs.get(11).z.get_integer_vec(), vec![0b0000000011000000000000000000i128]);
    //     assert_eq!(analysis.tab.stabs.get(11).sign, false);
    //     // Zcrz.in1
    //     assert_eq!(analysis.tab.stabs.get(12).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(12).z.get_integer_vec(), vec![0b0000000001000000000000000010i128]);
    //     assert_eq!(analysis.tab.stabs.get(12).sign, false);
    //     // Xrx.out
    //     assert_eq!(analysis.tab.stabs.get(13).x.get_integer_vec(), vec![0b0000000000100000000000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(13).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(13).sign, false);
    //     // Xtoffoli.in2
    //     assert_eq!(analysis.tab.stabs.get(14).x.get_integer_vec(), vec![0b0000000000010000000000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(14).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(14).sign, false);
    //     // Ztdg.out
    //     assert_eq!(analysis.tab.stabs.get(15).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(15).z.get_integer_vec(), vec![0b0000000000001000000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(15).sign, false);
    //     // Zrz.in
    //     assert_eq!(analysis.tab.stabs.get(16).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(16).z.get_integer_vec(), vec![0b0000000000000100000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(16).sign, false);
    //     // Zrz.out
    //     assert_eq!(analysis.tab.stabs.get(17).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(17).z.get_integer_vec(), vec![0b0000000000000010000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(17).sign, false);
    //     // Zmeas.in
    //     assert_eq!(analysis.tab.stabs.get(18).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(18).z.get_integer_vec(), vec![0b0000000000000001000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(18).sign, false);
    //     // Zmeas.out
    //     assert_eq!(analysis.tab.stabs.get(19).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(19).z.get_integer_vec(), vec![0b0000000000000000100000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(19).sign, false);
    //     // Zcrz.in0
    //     assert_eq!(analysis.tab.stabs.get(20).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(20).z.get_integer_vec(), vec![0b0000000000000000010000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(20).sign, false);
    //     // Zcrz.out0
    //     assert_eq!(analysis.tab.stabs.get(21).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(21).z.get_integer_vec(), vec![0b0000000000000000001000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(21).sign, false);
    //     // Zcrz.out1
    //     assert_eq!(analysis.tab.stabs.get(22).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(22).z.get_integer_vec(), vec![0b0000000000000000000100000010i128]);
    //     assert_eq!(analysis.tab.stabs.get(22).sign, false);
    //     // Ztoffoli.in0
    //     assert_eq!(analysis.tab.stabs.get(23).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(23).z.get_integer_vec(), vec![0b0000000000000000000010000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(23).sign, false);
    //     // Ztoffoli.in1
    //     assert_eq!(analysis.tab.stabs.get(24).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(24).z.get_integer_vec(), vec![0b0000000000000000000001000010i128]);
    //     assert_eq!(analysis.tab.stabs.get(24).sign, false);
    //     // Ztoffoli.out0
    //     assert_eq!(analysis.tab.stabs.get(25).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(25).z.get_integer_vec(), vec![0b0000000000000000000000100100i128]);
    //     assert_eq!(analysis.tab.stabs.get(25).sign, false);
    //     // Ztoffoli.out1
    //     assert_eq!(analysis.tab.stabs.get(26).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(26).z.get_integer_vec(), vec![0b0000000000000000000000010010i128]);
    //     assert_eq!(analysis.tab.stabs.get(26).sign, false);
    //     // Xtoffoli.out2
    //     assert_eq!(analysis.tab.stabs.get(27).x.get_integer_vec(), vec![0b0000000000000000000000001001i128]);
    //     assert_eq!(analysis.tab.stabs.get(27).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(27).sign, false);
    // }

    // #[test]
    // fn test_mixed_clifford() {
    //     // Need to cover MeasureFree, QFree and Reset
    //     let mut builder = DFGBuilder::new(Signature::new(vec![qb_t(), qb_t()], vec![qb_t()])).unwrap();
    //     let [qb0, qb1] = builder.input_wires_arr();
    //     let [qb2] = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap().outputs_arr();
    //     let [qb0, qb1] = builder.add_dataflow_op(TketOp::CX, [qb0, qb1]).unwrap().outputs_arr();
    //     let [qb0] = builder.add_dataflow_op(TketOp::Reset, [qb0]).unwrap().outputs_arr();
    //     let meas = builder.add_dataflow_op(TketOp::MeasureFree, [qb1]).unwrap().node();
    //     let [qb0] = builder.add_dataflow_op(TketOp::H, [qb0]).unwrap().outputs_arr();
    //     let [qb0, qb2] = builder.add_dataflow_op(TketOp::CX, [qb0, qb2]).unwrap().outputs_arr();
    //     builder.add_dataflow_op(TketOp::QFree, [qb2]);
    //     let hugr = builder.finish_hugr_with_outputs([qb0]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 4);
    //     // Input wires, alloc, and reset-alloc give 6 qubits/stabs
    //     // Reset-free and QFree remove 2 each
    //     // MeasureFree removes 1
    //     assert_eq!(analysis.tab.nb_stabs, 1);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(0)).unwrap(), 0);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(1)).unwrap(), 1);
    //     assert_eq!(*analysis.internal_in_cols.get(&(meas, IncomingPort::from(0))).unwrap(), 2);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(0)).unwrap(), 3);
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_boolean_vec(), vec![false; 6]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_boolean_vec(), vec![true, true, true, false]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, false);
    // }

    // #[test]
    // fn test_if_simple() {
    //     let mut builder = DFGBuilder::new(endo_sig(vec![qb_t(), qb_t(), bool_t()])).unwrap();
    //     let [qb0, qb1, b] = builder.input_wires_arr();
    //     let t = builder.add_dataflow_op(TketOp::T, [qb0]).unwrap();
    //     let [qb0] = t.outputs_arr();
    //     let mut cond_builder = builder.conditional_builder(([type_row![], type_row![]], b), [(qb_t(), qb0), (qb_t(), qb1)], vec![qb_t(); 2].into()).unwrap();
    //     let mut cond0_builder = cond_builder.case_builder(0).unwrap();
    //     let [c0q0, c0q1] = cond0_builder.input_wires_arr();
    //     let [c0q0, c0q1] = cond0_builder.add_dataflow_op(TketOp::CX, [c0q0, c0q1]).unwrap().outputs_arr();
    //     cond0_builder.finish_with_outputs([c0q0, c0q1]);
    //     let cond1_builder = cond_builder.case_builder(1).unwrap();
    //     let [c1c0, c1q1] = cond1_builder.input_wires_arr();
    //     cond1_builder.finish_with_outputs([c1c0, c1q1]);
    //     let cond = cond_builder.finish_sub_container().unwrap();
    //     let [qb0, qb1] = cond.outputs_arr();
    //     let tdg = builder.add_dataflow_op(TketOp::Tdg, [qb0]).unwrap();
    //     let [qb0] = tdg.outputs_arr();
    //     let hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    //     assert_eq!(analysis.tab.nb_qubits, 16);
    //     assert_eq!(analysis.tab.nb_stabs, 14);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(0)).unwrap(), 0);
    //     assert_eq!(*analysis.internal_in_cols.get(&(t.node(), IncomingPort::from(0))).unwrap(), 1);
    //     assert_eq!(*analysis.in_cols.get(&OutgoingPort::from(1)).unwrap(), 2);
    //     assert_eq!(*analysis.internal_in_cols.get(&(cond.node(), IncomingPort::from(1))).unwrap(), 3);
    //     assert_eq!(*analysis.internal_out_cols.get(&(t.node(), OutgoingPort::from(0))).unwrap(), 4);
    //     assert_eq!(*analysis.internal_in_cols.get(&(cond.node(), IncomingPort::from(0))).unwrap(), 5);
    //     assert_eq!(*analysis.nested_in_cols.get(&(cond.node(), OutgoingPort::from(0))).unwrap(), 6);
    //     assert_eq!(*analysis.nested_in_cols.get(&(cond.node(), OutgoingPort::from(1))).unwrap(), 7);
    //     assert_eq!(*analysis.nested_out_cols.get(&(cond.node(), IncomingPort::from(0))).unwrap(), 8);
    //     assert_eq!(*analysis.nested_out_cols.get(&(cond.node(), IncomingPort::from(1))).unwrap(), 9);
    //     assert_eq!(*analysis.internal_out_cols.get(&(cond.node(), OutgoingPort::from(0))).unwrap(), 10);
    //     assert_eq!(*analysis.internal_in_cols.get(&(tdg.node(), IncomingPort::from(0))).unwrap(), 11);
    //     assert_eq!(*analysis.internal_out_cols.get(&(cond.node(), OutgoingPort::from(1))).unwrap(), 12);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(1)).unwrap(), 13);
    //     assert_eq!(*analysis.internal_out_cols.get(&(tdg.node(), OutgoingPort::from(0))).unwrap(), 14);
    //     assert_eq!(*analysis.out_cols.get(&IncomingPort::from(0)).unwrap(), 15);
    //     //TODO:: Row echelon
    //     // Zin0
    //     assert_eq!(analysis.tab.stabs.get(0).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(0).z.get_integer_vec(), vec![0b1000000000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(0).sign, false);
    //     // Zt.in
    //     assert_eq!(analysis.tab.stabs.get(1).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(1).z.get_integer_vec(), vec![0b0100000000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(1).sign, false);
    //     // Xin1
    //     assert_eq!(analysis.tab.stabs.get(2).x.get_integer_vec(), vec![0b0010000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(2).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(2).sign, false);
    //     // Xcond.in1
    //     assert_eq!(analysis.tab.stabs.get(3).x.get_integer_vec(), vec![0b0001000000000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(3).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(3).sign, false);
    //     // Zt.out
    //     assert_eq!(analysis.tab.stabs.get(4).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(4).z.get_integer_vec(), vec![0b0000100000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(4).sign, false);
    //     // Zcond.in0
    //     assert_eq!(analysis.tab.stabs.get(5).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(5).z.get_integer_vec(), vec![0b0000010000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(5).sign, false);
    //     // Zcond.nin0
    //     assert_eq!(analysis.tab.stabs.get(6).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(6).z.get_integer_vec(), vec![0b0000001000000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(6).sign, false);
    //     // Xcond.nin1
    //     assert_eq!(analysis.tab.stabs.get(7).x.get_integer_vec(), vec![0b0000000100000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(7).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(7).sign, false);
    //     // Zcond.nout0
    //     assert_eq!(analysis.tab.stabs.get(8).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(8).z.get_integer_vec(), vec![0b0000000010000001i128]);
    //     assert_eq!(analysis.tab.stabs.get(8).sign, false);
    //     // Xcond.nout1
    //     assert_eq!(analysis.tab.stabs.get(9).x.get_integer_vec(), vec![0b0000000001000100i128]);
    //     assert_eq!(analysis.tab.stabs.get(9).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(9).sign, false);
    //     // Zcond.out0
    //     assert_eq!(analysis.tab.stabs.get(10).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(10).z.get_integer_vec(), vec![0b0000000000100001i128]);
    //     assert_eq!(analysis.tab.stabs.get(10).sign, false);
    //     // Ztdg.in
    //     assert_eq!(analysis.tab.stabs.get(11).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(11).z.get_integer_vec(), vec![0b0000000000010001i128]);
    //     assert_eq!(analysis.tab.stabs.get(11).sign, false);
    //     // Xcond.out1
    //     assert_eq!(analysis.tab.stabs.get(12).x.get_integer_vec(), vec![0b0000000000001100i128]);
    //     assert_eq!(analysis.tab.stabs.get(12).z.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(12).sign, false);
    //     // Ztdg.out
    //     assert_eq!(analysis.tab.stabs.get(13).x.get_integer_vec(), vec![0i128]);
    //     assert_eq!(analysis.tab.stabs.get(13).z.get_integer_vec(), vec![0b0000000000000011i128]);
    //     assert_eq!(analysis.tab.stabs.get(13).sign, false);
    // }

    // #[test]
    // fn test_loop_null() {
    //     let mut builder = DFGBuilder::new(endo_sig(vec![qb_t(), qb_t(), bool_t()])).unwrap();
    //     let [qb0, qb1, b] = builder.input_wires_arr();
    //     let reset = builder.add_dataflow_op(TketOp::Reset, [qb0]).unwrap();
    //     let [qb0] = reset.outputs_arr();
    //     let t = builder.add_dataflow_op(TketOp::T, [qb1]).unwrap();
    //     let [qb1] = t.outputs_arr();
    //     let mut loop_builder = builder.tail_loop_builder([(bool_t(), b)], [(qb_t(), qb0), (qb_t(), qb1), (bool_t(), b)], type_row![]).unwrap();
    //     let [loop_qb0, loop_qb1, loop_b] = loop_builder.input_wires_arr();
    //     let loop_t = loop_builder.add_dataflow_op(TketOp::T, [loop_qb0]).unwrap();
    //     let [loop_qb0] = loop_t.outputs_arr();
    //     let loop_tdg = loop_builder.add_dataflow_op(TketOp::Tdg, [loop_qb1]).unwrap();
    //     let [loop_qb1] = loop_tdg.outputs_arr();
    //     // let loop_b = loop_builder.make_break(loop_builder.loop_signature().unwrap().clone(), []).unwrap();
    //     let tl = loop_builder.finish_with_outputs(loop_b, [loop_qb0, loop_qb1, loop_b]).unwrap();
    //     let [qb0, qb1, b] = tl.outputs_arr();
    //     let tdg = builder.add_dataflow_op(TketOp::Tdg, [qb1]).unwrap();
    //     let [qb1] = tdg.outputs_arr();
    //     let hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
    //     let mut analysis = StabilizerDataflow::run_dfg(&hugr, hugr.module_root(), &FunctionOpacity::Opaque);
    // }
}