use std::collections::{HashMap, HashSet};

use hugr::{
    extension::{prelude::bool_t, simple_op::MakeExtensionOp},
    hugr::hugrmut::HugrMut,
    ops::{Const, LoadConstant, OpType, Value},
    std_extensions::{arithmetic::float_ops::FloatOps, logic::LogicOp},
    HugrView, IncomingPort, OutgoingPort,
};
use hugr_core::hugr::internal::PortgraphNodeMap;
use itertools::{chain, Itertools};
use petgraph::algo::toposort;
use tket::{
    extension::{
        bool::{bool_type, BoolOp, ConstBool},
        rotation::{rotation_type, ConstRotation, RotationOp},
    },
    TketOp,
};

use crate::{
    bit_vector::BitVector,
    stabilizer_dataflow::{DataflowPoint, FunctionOpacity, SDFAnalysis},
    symplectic_tableau::PauliXZ,
};

pub struct PhaseFoldSettings {
    function_opacity: FunctionOpacity,
    preserve_measures: bool,
}

pub fn phase_fold<H: HugrMut>(hugr: &mut H, settings: &PhaseFoldSettings) {
    if settings.function_opacity == FunctionOpacity::Inline {
        unimplemented!("Inlining of functions during phase folding is not yet implemented");
    }
    let analysis = SDFAnalysis::run_hugr(hugr, &settings.function_opacity);
    for (parent, summary) in analysis.0.iter() {
        let mut echelon_tab = summary.tab.clone();
        if echelon_tab.nb_stabs == 0 {
            break;
        } // Check for empty here to save us when iterating through the stabilizers
        let topo_order_lookup: HashMap<H::Node, usize> = {
            let (region, node_map) = hugr.region_portgraph(*parent);
            HashMap::from_iter(
                toposort(&region, None)
                    .unwrap()
                    .iter()
                    .enumerate()
                    .map(|(i, ni)| (node_map.from_portgraph(*ni), i)),
            )
        };

        let mut z_rot_cols: Vec<(usize, PauliXZ)> = vec![];
        let mut other_cols: Vec<(usize, PauliXZ)> = vec![];
        let mut measure_set: HashSet<usize> = HashSet::default();
        let mut measurefree_set: HashSet<usize> = HashSet::default();
        for (dfp, col) in summary.q_index_map.iter() {
            other_cols.push((*col, PauliXZ::X));
            let mut rot_found = false;
            match dfp {
                // Only look at the InternalIn rather than including the InternalOut as well since all rotation gates have a stabilizer across them that says they act identically, so if one wire has a particular value then so will the other
                // Plus MeasureFree only has an InternalIn
                DataflowPoint::InternalIn(node, _) => {
                    if let OpType::ExtensionOp(op) = hugr.get_optype(*node) {
                        if let Ok(tkop) = TketOp::from_extension_op(op) {
                            match tkop {
                                TketOp::T | TketOp::Tdg | TketOp::Rz => {
                                    z_rot_cols.push((*col, PauliXZ::Z));
                                    rot_found = true;
                                }
                                TketOp::Measure => {
                                    measure_set.insert(*col);
                                    z_rot_cols.push((*col, PauliXZ::Z));
                                    rot_found = true;
                                }
                                TketOp::MeasureFree => {
                                    measurefree_set.insert(*col);
                                    z_rot_cols.push((*col, PauliXZ::Z));
                                    rot_found = true;
                                }
                                TketOp::Rx => {
                                    // Map X into Z basis
                                    echelon_tab.append_h(*col);
                                    z_rot_cols.push((*col, PauliXZ::Z));
                                    rot_found = true;
                                }
                                TketOp::Ry => {
                                    // Map Y into Z basis
                                    echelon_tab.append_v(*col);
                                    z_rot_cols.push((*col, PauliXZ::Z));
                                    rot_found = true;
                                }
                                _ => {}
                            }
                        }
                    }
                }
                _ => {}
            }
            if !rot_found {
                other_cols.push((*col, PauliXZ::Z));
            }
        }

        // Echelon form with rotations first will express them wrt other values of the circuit in a canonical way s.t. Za Zb is possible between two rotations iff either Za Zb is a row directly or there exists some P s.t. Za P and Zb P are rows this tableau
        // The case of Za Zb directly only occurs if everything in the group of rotations to be merged are MeasureFree ops; every other rotation type at least generates a stabilizer between its InternalIn and InternalOut, so this guarantees independent rows when we put into row echelon form
        echelon_tab.echelon(&chain!(z_rot_cols.clone(), other_cols).collect_vec());

        // For each row +- Za P corresponding to rotation a, add (a, +-) to the bucket addressed by P
        let mut fold_buckets: HashMap<(BitVector, BitVector), Vec<(usize, bool)>> =
            HashMap::default();
        let mut current_stab = 0;
        for (col, _) in z_rot_cols.iter() {
            if echelon_tab.z[current_stab].get(*col) {
                let mut pauli = (
                    echelon_tab.z[current_stab].clone(),
                    echelon_tab.x[current_stab].clone(),
                );
                pauli.0.xor_bit(*col);
                if !fold_buckets.contains_key(&pauli) {
                    fold_buckets.insert(pauli.clone(), vec![]);
                }
                fold_buckets
                    .get_mut(&pauli)
                    .unwrap()
                    .push((*col, echelon_tab.signs.get(current_stab)));
                current_stab += 1;
                if current_stab == echelon_tab.nb_stabs {
                    break;
                }
            }
        }

        for (pauli, bucket) in fold_buckets.iter_mut() {
            let mut is_const = false;
            if pauli.1.get_all_ones(echelon_tab.nb_qubits).is_empty() {
                let z_ones = pauli.0.get_all_ones(echelon_tab.nb_qubits);
                if z_ones.is_empty() {
                    // Check for the case where pauli is identity, i.e. everything in bucket is acting on a constant
                    is_const = true;
                } else if let Ok(z_col) = z_ones.iter().exactly_one() {
                    // Check for the case where pauli is still a rotation gate; necessarily it and all in the bucket must be MeasureFree ops
                    if measurefree_set.contains(z_col) {
                        // Handle this by adding the MeasureFree itself to the bucket
                        bucket.push((*z_col, false));
                    }
                }
            }

            if is_const {
                // Can just remove any non-measures
                // If settings.preserve_measures==false, we can also remove all measures
                for (col, polarity) in bucket {
                    let DataflowPoint::InternalIn(node, _) =
                        summary.q_index_map.get_by_right(&col).unwrap()
                    else {
                        unreachable!();
                    };
                    if measure_set.contains(&col) {
                        if !settings.preserve_measures {
                            let (q_pred, q_pred_port) = hugr
                                .single_linked_output(*node, IncomingPort::from(0))
                                .unwrap();
                            let (q_succ, q_succ_port) = hugr
                                .single_linked_input(*node, OutgoingPort::from(0))
                                .unwrap();
                            hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                            let measure_result = if *polarity {
                                Value::true_val()
                            } else {
                                Value::false_val()
                            };
                            let const_val_node = hugr
                                .add_node_with_parent(*parent, Into::<Const>::into(measure_result));
                            let load_const_node = hugr
                                .add_node_with_parent(*parent, LoadConstant { datatype: bool_t() });
                            hugr.connect(
                                const_val_node,
                                OutgoingPort::from(0),
                                load_const_node,
                                IncomingPort::from(0),
                            );
                            for (c_succ, c_succ_port) in hugr
                                .linked_inputs(*node, OutgoingPort::from(1))
                                .collect_vec()
                            {
                                hugr.connect(
                                    load_const_node,
                                    OutgoingPort::from(0),
                                    c_succ,
                                    c_succ_port,
                                );
                            }
                            hugr.remove_node(*node);
                        }
                    } else if measurefree_set.contains(&col) {
                        if !settings.preserve_measures {
                            let (q_pred, q_pred_port) = hugr
                                .single_linked_output(*node, IncomingPort::from(0))
                                .unwrap();
                            let qfree_node = hugr.add_node_with_parent(*parent, TketOp::QFree);
                            hugr.connect(q_pred, q_pred_port, qfree_node, IncomingPort::from(0));
                            let measure_result = ConstBool::new(*polarity);
                            let const_val_node = hugr
                                .add_node_with_parent(*parent, Const::new(measure_result.into()));
                            let load_const_node = hugr.add_node_with_parent(
                                *parent,
                                LoadConstant {
                                    datatype: bool_type(),
                                },
                            );
                            hugr.connect(
                                const_val_node,
                                OutgoingPort::from(0),
                                load_const_node,
                                IncomingPort::from(0),
                            );
                            for (c_succ, c_succ_port) in hugr
                                .linked_inputs(*node, OutgoingPort::from(0))
                                .collect_vec()
                            {
                                hugr.connect(
                                    load_const_node,
                                    OutgoingPort::from(0),
                                    c_succ,
                                    c_succ_port,
                                );
                            }
                            hugr.remove_node(*node);
                        }
                    } else {
                        // Is some kind of rotation gate
                        let (q_pred, q_pred_port) = hugr
                            .single_linked_output(*node, IncomingPort::from(0))
                            .unwrap();
                        let (q_succ, q_succ_port) = hugr
                            .single_linked_input(*node, OutgoingPort::from(0))
                            .unwrap();
                        hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                        hugr.remove_node(*node);
                    }
                }
            } else if bucket
                .iter()
                .any(|(col, _)| measure_set.contains(col) || measurefree_set.contains(col))
            {
                // Can just remove any non-measures
                // If settings.preserve_measures==false, we can also remove all but the earliest measure (keeping the earliest to make sure the outcomes are available whenever they get used)
                let (earliest_measure_node, earliest_measure_port, earliest_measure_polarity, _) =
                    bucket
                        .iter()
                        .filter(|(col, _)| {
                            measure_set.contains(col) || measurefree_set.contains(col)
                        })
                        .map(|(col, polarity)| {
                            let DataflowPoint::InternalIn(node, _) =
                                summary.q_index_map.get_by_right(col).unwrap()
                            else {
                                unreachable!()
                            };
                            let port = if measure_set.contains(col) {
                                OutgoingPort::from(1)
                            } else {
                                OutgoingPort::from(0)
                            };
                            (node, port, polarity, topo_order_lookup.get(node).unwrap())
                        })
                        .min_by(|(_, _, _, pos0), (_, _, _, pos1)| pos0.cmp(pos1))
                        .unwrap();
                let (_, earliest_measure_result_type) = hugr
                    .out_value_types(*earliest_measure_node)
                    .find(|(port, _)| *port == earliest_measure_port)
                    .unwrap();
                let mut prelude_bool_result_loc: Option<(H::Node, OutgoingPort)> =
                    if earliest_measure_result_type == bool_t() {
                        Some((*earliest_measure_node, earliest_measure_port))
                    } else {
                        None
                    };
                let mut tket_bool_result_loc: Option<(H::Node, OutgoingPort)> =
                    if earliest_measure_result_type == bool_type() {
                        Some((*earliest_measure_node, earliest_measure_port))
                    } else {
                        None
                    };
                for (col, polarity) in bucket.iter() {
                    let DataflowPoint::InternalIn(node, _) =
                        summary.q_index_map.get_by_right(col).unwrap()
                    else {
                        unreachable!();
                    };
                    if measure_set.contains(col) {
                        if !settings.preserve_measures && node != earliest_measure_node {
                            let (q_pred, q_pred_port) = hugr
                                .single_linked_output(*node, IncomingPort::from(0))
                                .unwrap();
                            let (q_succ, q_succ_port) = hugr
                                .single_linked_input(*node, OutgoingPort::from(0))
                                .unwrap();
                            hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                            if prelude_bool_result_loc.is_none() {
                                // Earliest measure was a MeasureFree
                                // Add a read op to unpack the tket bool outcome
                                let read_node = hugr.add_node_with_parent(*parent, BoolOp::read);
                                hugr.connect(
                                    *earliest_measure_node,
                                    earliest_measure_port,
                                    read_node,
                                    IncomingPort::from(0),
                                );
                                prelude_bool_result_loc = Some((read_node, OutgoingPort::from(0)));
                            }
                            let (result_node, result_port) = prelude_bool_result_loc.unwrap();
                            if polarity == earliest_measure_polarity {
                                for (c_succ, c_succ_port) in hugr
                                    .linked_inputs(*node, OutgoingPort::from(1))
                                    .collect_vec()
                                {
                                    hugr.connect(result_node, result_port, c_succ, c_succ_port);
                                }
                            } else {
                                let not_node = hugr.add_node_with_parent(*parent, LogicOp::Not);
                                hugr.connect(
                                    result_node,
                                    result_port,
                                    not_node,
                                    IncomingPort::from(0),
                                );
                                for (c_succ, c_succ_port) in hugr
                                    .linked_inputs(*node, OutgoingPort::from(1))
                                    .collect_vec()
                                {
                                    hugr.connect(
                                        not_node,
                                        OutgoingPort::from(0),
                                        c_succ,
                                        c_succ_port,
                                    );
                                }
                            }
                            hugr.remove_node(*node);
                        }
                    } else if measurefree_set.contains(col) {
                        if !settings.preserve_measures && node != earliest_measure_node {
                            let (q_pred, q_pred_port) = hugr
                                .single_linked_output(*node, IncomingPort::from(0))
                                .unwrap();
                            let qfree_node = hugr.add_node_with_parent(*parent, TketOp::QFree);
                            hugr.connect(q_pred, q_pred_port, qfree_node, IncomingPort::from(0));
                            if tket_bool_result_loc.is_none() {
                                // Earliest measure was a Measure
                                // Add a make_opaque op to pack the prelude::bool into a tket::bool
                                let make_opaque_node =
                                    hugr.add_node_with_parent(*parent, BoolOp::make_opaque);
                                hugr.connect(
                                    *earliest_measure_node,
                                    earliest_measure_port,
                                    make_opaque_node,
                                    IncomingPort::from(0),
                                );
                                tket_bool_result_loc =
                                    Some((make_opaque_node, OutgoingPort::from(0)));
                            }
                            let (result_node, result_port) = tket_bool_result_loc.unwrap();
                            if polarity == earliest_measure_polarity {
                                for (c_succ, c_succ_port) in hugr
                                    .linked_inputs(*node, OutgoingPort::from(0))
                                    .collect_vec()
                                {
                                    hugr.connect(result_node, result_port, c_succ, c_succ_port);
                                }
                            } else {
                                let not_node = hugr.add_node_with_parent(*parent, BoolOp::not);
                                hugr.connect(
                                    result_node,
                                    result_port,
                                    not_node,
                                    IncomingPort::from(0),
                                );
                                for (c_succ, c_succ_port) in hugr
                                    .linked_inputs(*node, OutgoingPort::from(0))
                                    .collect_vec()
                                {
                                    hugr.connect(
                                        not_node,
                                        OutgoingPort::from(0),
                                        c_succ,
                                        c_succ_port,
                                    );
                                }
                            }
                            hugr.remove_node(*node);
                        }
                    } else {
                        // Is some kind of rotation gate
                        let (q_pred, q_pred_port) = hugr
                            .single_linked_output(*node, IncomingPort::from(0))
                            .unwrap();
                        let (q_succ, q_succ_port) = hugr
                            .single_linked_input(*node, OutgoingPort::from(0))
                            .unwrap();
                        hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                        hugr.remove_node(*node);
                    }
                }
            } else {
                // Merge into latest rotation to make sure all angles are available by the time of the rotation
                let (latest_rotation_node, latest_rotation_polarity, _) = bucket
                    .iter()
                    .map(|(col, polarity)| {
                        let DataflowPoint::InternalIn(node, _) =
                            summary.q_index_map.get_by_right(col).unwrap()
                        else {
                            unreachable!()
                        };
                        (node, polarity, topo_order_lookup.get(node).unwrap())
                    })
                    .max_by(|(_, _, pos0), (_, _, pos1)| pos0.cmp(pos1))
                    .unwrap();
                let mut acc_static = 0.;
                let mut acc_dynamic: Option<(H::Node, OutgoingPort)> = None;
                for (col, polarity) in bucket.iter() {
                    let DataflowPoint::InternalIn(node, _) =
                        summary.q_index_map.get_by_right(col).unwrap()
                    else {
                        unreachable!();
                    };
                    let OpType::ExtensionOp(op) = hugr.get_optype(*node) else {
                        unreachable!()
                    };
                    let Ok(tkop) = TketOp::from_extension_op(op) else {
                        unreachable!()
                    };
                    match tkop {
                        TketOp::T => {
                            acc_static += if polarity == latest_rotation_polarity {
                                0.25
                            } else {
                                -0.25
                            };
                        }
                        TketOp::Tdg => {
                            acc_static += if polarity == latest_rotation_polarity {
                                -0.25
                            } else {
                                0.25
                            };
                        }
                        TketOp::Rx | TketOp::Ry | TketOp::Rz => {
                            let (source_node, source_port) = hugr
                                .single_linked_output(*node, IncomingPort::from(1))
                                .unwrap();
                            let (value_node, value_port) = if polarity == latest_rotation_polarity {
                                (source_node, source_port)
                            } else {
                                let to_float =
                                    hugr.add_node_with_parent(*parent, RotationOp::to_halfturns);
                                hugr.connect(
                                    source_node,
                                    source_port,
                                    to_float,
                                    IncomingPort::from(0),
                                );
                                let neg_node = hugr.add_node_with_parent(*parent, FloatOps::fneg);
                                hugr.connect(
                                    to_float,
                                    OutgoingPort::from(0),
                                    neg_node,
                                    IncomingPort::from(0),
                                );
                                let from_float = hugr.add_node_with_parent(
                                    *parent,
                                    RotationOp::from_halfturns_unchecked,
                                );
                                hugr.connect(
                                    neg_node,
                                    OutgoingPort::from(0),
                                    from_float,
                                    IncomingPort::from(0),
                                );
                                (from_float, OutgoingPort::from(0))
                            };
                            match acc_dynamic {
                                Some((n, p)) => {
                                    let add_node =
                                        hugr.add_node_with_parent(*parent, RotationOp::radd);
                                    hugr.connect(n, p, add_node, IncomingPort::from(0));
                                    hugr.connect(
                                        value_node,
                                        value_port,
                                        add_node,
                                        IncomingPort::from(1),
                                    );
                                    acc_dynamic = Some((add_node, OutgoingPort::from(0)));
                                }
                                None => {
                                    acc_dynamic = Some((value_node, value_port));
                                }
                            }
                        }
                        _ => {
                            // Other gate types are not considered in phase folding
                            // This includes Clifford phase gates (S, Sdg, Z) which would be abstracted away
                        }
                    }
                    if node != latest_rotation_node {
                        let (q_pred, q_pred_port) = hugr
                            .single_linked_output(*node, IncomingPort::from(0))
                            .unwrap();
                        let (q_succ, q_succ_port) = hugr
                            .single_linked_input(*node, OutgoingPort::from(0))
                            .unwrap();
                        hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                        hugr.remove_node(*node);
                    }
                }
                let OpType::ExtensionOp(op) = hugr.get_optype(*latest_rotation_node) else {
                    unreachable!()
                };
                let Ok(tkop) = TketOp::from_extension_op(op) else {
                    unreachable!()
                };
                let (q_pred, q_pred_port) = hugr
                    .single_linked_output(*latest_rotation_node, IncomingPort::from(0))
                    .unwrap();
                let (q_succ, q_succ_port) = hugr
                    .single_linked_input(*latest_rotation_node, OutgoingPort::from(0))
                    .unwrap();
                match acc_dynamic {
                    Some((mut source_node, mut source_port)) => {
                        if acc_static != 0. {
                            let const_val_node = hugr.add_node_with_parent(
                                *parent,
                                Into::<Const>::into(Value::extension(
                                    ConstRotation::new(acc_static % 2.).unwrap(),
                                )),
                            );
                            let load_const_node = hugr.add_node_with_parent(
                                *parent,
                                LoadConstant {
                                    datatype: rotation_type(),
                                },
                            );
                            hugr.connect(
                                const_val_node,
                                OutgoingPort::from(0),
                                load_const_node,
                                IncomingPort::from(0),
                            );
                            let add_node = hugr.add_node_with_parent(*parent, RotationOp::radd);
                            hugr.connect(source_node, source_port, add_node, IncomingPort::from(0));
                            hugr.connect(
                                load_const_node,
                                OutgoingPort::from(0),
                                add_node,
                                IncomingPort::from(1),
                            );
                            source_node = add_node;
                            source_port = OutgoingPort::from(0);
                        }
                        let new_optype = match tkop {
                            TketOp::T | TketOp::Tdg | TketOp::Rz => TketOp::Rz,
                            TketOp::Rx => TketOp::Rx,
                            TketOp::Ry => TketOp::Ry,
                            _ => {
                                unreachable!()
                            }
                        };
                        let new_rotation_node = hugr.add_node_with_parent(*parent, new_optype);
                        hugr.connect(
                            q_pred,
                            q_pred_port,
                            new_rotation_node,
                            IncomingPort::from(0),
                        );
                        hugr.connect(
                            source_node,
                            source_port,
                            new_rotation_node,
                            IncomingPort::from(1),
                        );
                        hugr.connect(
                            new_rotation_node,
                            OutgoingPort::from(0),
                            q_succ,
                            q_succ_port,
                        );
                    }
                    None => {
                        let cliff_t_seq: Option<Vec<TketOp>> = match acc_static % 2. {
                            0. => Some(vec![]),
                            0.25 | -1.75 => Some(vec![TketOp::T]),
                            0.5 | -1.5 => Some(vec![TketOp::S]),
                            0.75 | -1.25 => Some(vec![TketOp::Z, TketOp::Tdg]),
                            1. | -1. => Some(vec![TketOp::Z]),
                            1.25 | -0.75 => Some(vec![TketOp::Z, TketOp::T]),
                            1.5 | -0.5 => Some(vec![TketOp::Sdg]),
                            1.75 | -0.25 => Some(vec![TketOp::Tdg]),
                            _ => None,
                        };
                        match cliff_t_seq {
                            Some(mut c_t_seq) => {
                                let mut front_node = q_pred;
                                let mut front_port = q_pred_port;
                                match tkop {
                                    TketOp::Rx => {
                                        c_t_seq =
                                            chain!([TketOp::H], c_t_seq, [TketOp::H]).collect_vec();
                                    }
                                    TketOp::Ry => {
                                        c_t_seq = chain!([TketOp::V], c_t_seq, [TketOp::Vdg])
                                            .collect_vec();
                                    }
                                    _ => {}
                                }
                                for gate in c_t_seq {
                                    let gate_node = hugr.add_node_with_parent(*parent, gate);
                                    hugr.connect(
                                        front_node,
                                        front_port,
                                        gate_node,
                                        IncomingPort::from(0),
                                    );
                                    front_node = gate_node;
                                    front_port = OutgoingPort::from(0);
                                }
                                hugr.connect(front_node, front_port, q_succ, q_succ_port);
                            }
                            None => {
                                let const_val_node = hugr.add_node_with_parent(
                                    *parent,
                                    Into::<Const>::into(Value::extension(
                                        ConstRotation::new(acc_static % 2.).unwrap(),
                                    )),
                                );
                                let load_const_node = hugr.add_node_with_parent(
                                    *parent,
                                    LoadConstant {
                                        datatype: rotation_type(),
                                    },
                                );
                                hugr.connect(
                                    const_val_node,
                                    OutgoingPort::from(0),
                                    load_const_node,
                                    IncomingPort::from(0),
                                );
                                let new_optype = match tkop {
                                    TketOp::T | TketOp::Tdg | TketOp::Rz => TketOp::Rz,
                                    TketOp::Rx => TketOp::Rx,
                                    TketOp::Ry => TketOp::Ry,
                                    _ => {
                                        unreachable!()
                                    }
                                };
                                let new_rotation_node =
                                    hugr.add_node_with_parent(*parent, new_optype);
                                hugr.connect(
                                    q_pred,
                                    q_pred_port,
                                    new_rotation_node,
                                    IncomingPort::from(0),
                                );
                                hugr.connect(
                                    load_const_node,
                                    OutgoingPort::from(0),
                                    new_rotation_node,
                                    IncomingPort::from(1),
                                );
                                hugr.connect(
                                    new_rotation_node,
                                    OutgoingPort::from(0),
                                    q_succ,
                                    q_succ_port,
                                );
                            }
                        }
                    }
                }
                hugr.remove_node(*latest_rotation_node);
            }
        }
    }
}

#[cfg(test)]
mod test {
    use hugr::{
        builder::{
            endo_sig, Dataflow, DataflowHugr, DataflowSubContainer, FunctionBuilder, HugrBuilder,
            SubContainer,
        },
        extension::prelude::{bool_t, qb_t, usize_t},
        type_row,
        types::Signature,
        Hugr, HugrView,
    };
    use rstest::fixture;
    use tket::{
        extension::{bool::BoolOp, rotation::rotation_type},
        TketOp,
    };

    use crate::{
        phase_fold::{phase_fold, PhaseFoldSettings},
        stabilizer_dataflow::FunctionOpacity,
    };

    #[test]
    fn test_empty() {
        let builder = FunctionBuilder::new("empty", endo_sig(vec![])).unwrap();
        let mut hugr = builder.finish_hugr().unwrap();
        assert_eq!(hugr.num_nodes(), 4);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        assert_eq!(hugr.num_nodes(), 4);
    }

    #[test]
    fn test_identity() {
        let builder =
            FunctionBuilder::new("identity", endo_sig(vec![usize_t(), qb_t(), qb_t()])).unwrap();
        let [i, qb0, qb1] = builder.input_wires_arr();
        let mut hugr = builder.finish_hugr_with_outputs([i, qb0, qb1]).unwrap();
        assert_eq!(hugr.num_nodes(), 4);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        assert_eq!(hugr.num_nodes(), 4);
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
        let mut hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
        assert_eq!(hugr.num_nodes(), 14);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        assert_eq!(hugr.num_nodes(), 12);
        assert!(hugr.validate().is_ok());
    }

    #[test]
    fn test_loop_simple() {
        let mut builder =
            FunctionBuilder::new("loop_simple", endo_sig(vec![qb_t(), qb_t(), bool_t()])).unwrap();
        let [qb0, qb1, b] = builder.input_wires_arr();
        let t = builder.add_dataflow_op(TketOp::T, [qb0]).unwrap();
        let [qb0] = t.outputs_arr();
        let mut loop_builder = builder
            .tail_loop_builder(
                [],
                [(qb_t(), qb0), (qb_t(), qb1), (bool_t(), b)],
                type_row![],
            )
            .unwrap();
        let [loop_qb0, loop_qb1, loop_b] = loop_builder.input_wires_arr();
        let [loop_qb0, loop_qb1] = loop_builder
            .add_dataflow_op(TketOp::CX, [loop_qb0, loop_qb1])
            .unwrap()
            .outputs_arr();
        let loop_cond = loop_builder
            .make_break(loop_builder.loop_signature().unwrap().clone(), [])
            .unwrap();
        let tl = loop_builder
            .finish_with_outputs(loop_cond, [loop_qb0, loop_qb1, loop_b])
            .unwrap();
        let [qb0, qb1, b] = tl.outputs_arr();
        let tdg = builder.add_dataflow_op(TketOp::Tdg, [qb0]).unwrap();
        let [qb0] = tdg.outputs_arr();
        let mut hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
        assert_eq!(hugr.num_nodes(), 11);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        assert_eq!(hugr.num_nodes(), 9);
        assert!(hugr.validate().is_ok());
    }

    #[test]
    fn test_loop_swap() {
        let mut builder =
            FunctionBuilder::new("loop_swap", endo_sig(vec![qb_t(), qb_t(), bool_t()])).unwrap();
        let [qb0, qb1, b] = builder.input_wires_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::T, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let mut loop_builder = builder
            .tail_loop_builder(
                [],
                [(qb_t(), qb0), (qb_t(), qb1), (bool_t(), b)],
                type_row![],
            )
            .unwrap();
        let [loop_qb0, loop_qb1, loop_b] = loop_builder.input_wires_arr();
        let loop_cond = loop_builder
            .make_continue(loop_builder.loop_signature().unwrap().clone(), [])
            .unwrap();
        let tl = loop_builder
            .finish_with_outputs(loop_cond, [loop_qb1, loop_qb0, loop_b])
            .unwrap();
        let [qb0, qb1, b] = tl.outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let [qb1] = builder
            .add_dataflow_op(TketOp::Tdg, [qb1])
            .unwrap()
            .outputs_arr();
        let [qb0, qb1] = builder
            .add_dataflow_op(TketOp::CX, [qb0, qb1])
            .unwrap()
            .outputs_arr();
        let mut hugr = builder.finish_hugr_with_outputs([qb0, qb1, b]).unwrap();
        assert_eq!(hugr.num_nodes(), 14);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        assert_eq!(hugr.num_nodes(), 12);
        assert!(hugr.validate().is_ok());
    }

    #[test]
    fn test_merge_rotation_cases() {
        let mut builder = FunctionBuilder::new(
            "merge_rotation_cases",
            Signature::new(
                vec![
                    qb_t(),
                    rotation_type(),
                    rotation_type(),
                    rotation_type(),
                    rotation_type(),
                ],
                vec![qb_t()],
            ),
        )
        .unwrap();
        let [q, a0, a1, a2, a3] = builder.input_wires_arr();
        // Merge four rotation gates by addition of phases (with different original bases)
        let [q] = builder
            .add_dataflow_op(TketOp::Rz, [q, a0])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Rx, [q, a1])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Vdg, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Ry, [q, a2])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::V, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Rz, [q, a3])
            .unwrap()
            .outputs_arr();
        // Begin new chain
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        // Merge rotation gates with different polarities (final gate being Rx)
        let [q] = builder
            .add_dataflow_op(TketOp::Rz, [q, a0])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::X, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Rx, [q, a1])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        // Begin new chain
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        // Merge rotation gates with constant phases (final gate being Ry)
        let [q] = builder
            .add_dataflow_op(TketOp::T, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Vdg, [q])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::Ry, [q, a1])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::V, [q])
            .unwrap()
            .outputs_arr();
        // Begin new chain
        let [q] = builder
            .add_dataflow_op(TketOp::H, [q])
            .unwrap()
            .outputs_arr();
        // Merge rotation gates with constant phases with different parities (final gate being a constant phase)
        let [q] = builder
            .add_dataflow_op(TketOp::Rz, [q, a0])
            .unwrap()
            .outputs_arr();
        let [q] = builder
            .add_dataflow_op(TketOp::X, [q])
            .unwrap()
            .outputs_arr();
        let [mut q] = builder
            .add_dataflow_op(TketOp::T, [q])
            .unwrap()
            .outputs_arr();
        for i in 1..=8 {
            // Begin new chain
            let [new_q] = builder
                .add_dataflow_op(TketOp::H, [q])
                .unwrap()
                .outputs_arr();
            q = new_q;
            // Merge multiples of T gate phases
            for _ in 0..i {
                let [new_q] = builder
                    .add_dataflow_op(TketOp::T, [q])
                    .unwrap()
                    .outputs_arr();
                q = new_q;
            }
        }
        let mut hugr = builder.finish_hugr_with_outputs([q]).unwrap();
        assert_eq!(hugr.num_nodes(), 71);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        // Merging Rz's replaces each Rz by an Add
        // When polarities differ, add an extra node for the Neg
        // When sum includes a constant, the constant phase is replaced by a constant, load constant, and Add
        // Replace the 36 T gates in the strings by 9 gates from the LUT
        assert_eq!(hugr.num_nodes(), 54);
        assert!(hugr.validate().is_ok());
    }

    #[fixture]
    fn merge_measure_cases() -> Hugr {
        let mut builder = FunctionBuilder::new(
            "merge_measure_cases",
            Signature::new(vec![qb_t(), qb_t(), rotation_type()], vec![bool_t(); 6]),
        )
        .unwrap();
        let [q0, q1, a] = builder.input_wires_arr();
        // Merge two measures and rotations (earliest is Measure)
        let [q0] = builder
            .add_dataflow_op(TketOp::Rz, [q0, a])
            .unwrap()
            .outputs_arr();
        let [q0, b0] = builder
            .add_dataflow_op(TketOp::Measure, [q0])
            .unwrap()
            .outputs_arr();
        let [b1] = builder
            .add_dataflow_op(TketOp::MeasureFree, [q0])
            .unwrap()
            .outputs_arr();
        let [b1] = builder
            .add_dataflow_op(BoolOp::read, [b1])
            .unwrap()
            .outputs_arr(); // MeasureFree produces a tket.bool rather than hugr.prelude.bool
                            // Merge two measures with different polarity (earliest is MeasureFree)
        let [q2] = builder
            .add_dataflow_op(TketOp::QAlloc, [])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::X, [q2])
            .unwrap()
            .outputs_arr();
        let [q1, q2] = builder
            .add_dataflow_op(TketOp::CX, [q1, q2])
            .unwrap()
            .outputs_arr();
        let [b2] = builder
            .add_dataflow_op(TketOp::MeasureFree, [q1])
            .unwrap()
            .outputs_arr();
        let [b2] = builder
            .add_dataflow_op(BoolOp::read, [b2])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::H, [q2])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::Rx, [q2, a])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::H, [q2])
            .unwrap()
            .outputs_arr();
        let [q2, b3] = builder
            .add_dataflow_op(TketOp::Measure, [q2])
            .unwrap()
            .outputs_arr();
        // Merge measures and rotations for constant
        let [q2] = builder
            .add_dataflow_op(TketOp::Reset, [q2])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::T, [q2])
            .unwrap()
            .outputs_arr();
        let [q2, b4] = builder
            .add_dataflow_op(TketOp::Measure, [q2])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::X, [q2])
            .unwrap()
            .outputs_arr();
        let [q2] = builder
            .add_dataflow_op(TketOp::Tdg, [q2])
            .unwrap()
            .outputs_arr();
        let [b5] = builder
            .add_dataflow_op(TketOp::MeasureFree, [q2])
            .unwrap()
            .outputs_arr();
        let [b5] = builder
            .add_dataflow_op(BoolOp::read, [b5])
            .unwrap()
            .outputs_arr();
        builder
            .finish_hugr_with_outputs([b0, b1, b2, b3, b4, b5])
            .unwrap()
    }

    #[test]
    fn test_merge_measure_cases_preserve_measures() {
        let mut hugr = merge_measure_cases();
        assert_eq!(hugr.num_nodes(), 24);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: true,
            },
        );
        // The only gate removed should be the Rz, Rx, T, Tdg
        assert_eq!(hugr.num_nodes(), 20);
        assert!(hugr.validate().is_ok());
    }

    #[test]
    fn test_merge_measure_cases_remove_measures() {
        let mut hugr = merge_measure_cases();
        assert_eq!(hugr.num_nodes(), 24);
        phase_fold(
            &mut hugr,
            &PhaseFoldSettings {
                function_opacity: FunctionOpacity::Opaque,
                preserve_measures: false,
            },
        );
        // Should remove all 4 rotations, two Measures and two MeasureFrees
        // Each removed MeasureFree adds back a QFree
        // The two constant measures are replaced by Const and LoadConstant
        // For the other removed measures, we get a make_opaque and read for the conversions, plus a negation
        assert!(hugr.validate().is_ok());
        assert_eq!(hugr.num_nodes(), 25);
    }
}
