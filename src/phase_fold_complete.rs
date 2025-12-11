use std::collections::{BTreeMap, HashMap, HashSet};

use bimap::BiHashMap;
use hugr::{
    extension::simple_op::MakeExtensionOp,
    hugr::hugrmut::HugrMut,
    ops::{Const, LoadConstant, OpType, Tag, Value},
    std_extensions::{arithmetic::float_ops::FloatOps, logic::LogicOp},
    types::TypeRow,
    HugrView, IncomingPort, OutgoingPort,
};
use hugr_core::hugr::internal::PortgraphNodeMap;
use indexmap::{IndexMap, IndexSet};
use itertools::{chain, Itertools};
use petgraph::{algo::toposort, unionfind::UnionFind};
use tket::{
    extension::{
        bool::{bool_type, BoolOp, ConstBool},
        rotation::{rotation_type, ConstRotation, RotationOp},
    },
    TketOp,
};

use crate::{
    bit_vector::BitVector,
    stabilizer_dataflow_complete::{DataflowPoint, StabilizerDataflow},
    symplectic_tableau::{PauliXZ, SymplecticTableau},
};

/**
 * Phase folding identifies pairs of gates that can be combined with one another or those that can be replaced by a constant, primarily focussing on rotations and measurements. We will only look at the nodes for which there exist `DataflowPoint::Rotation`s in the provided analysis; therefore we require DataflowSettings::rotation_override==true for any nodes we wish to consider for phase folding.
 *
 * A pair of basic rotation gates (T, Tdg, Rx, Ry, Rz) on which we find a +-ZZ stabilizer can be merged. The sign indicates whether or not to add or subtract the angles of rotation.
 * A pair of measurement gates (Measure or MeasureFree) on which we find a +-ZZ stabilizer can be merged. The sign indicates whether the outcomes are the same or opposite.
 * A pair of one basic rotation gate and one measurement gate on which we find a +-ZZ stabilizer can be replaced by the measurement followed by a classically-conditioned global phase. The sign indicates the sign of the global phase introduced.
 * A basic rotation gate on which we find a +-Z stabilizer can be replaced with a global phase. The sign indicates the sign of the global phase introduced.
 * A measurement gate on which we find a +-Z stabilizer can be replaced with a constant classical value. The sign indicates whether the constant value is false (positive) or true (negative).
 * If there is a stabilizer whose non-identity Paulis are Z on a rotation gate and all others are on inputs to QFree or Reset gates, then the rotation gate has no effect and can be removed (we need not attempt to preserve global phase here since the discard is a CPTP-map).
 *
 * We can sort basic rotations and measurements into buckets.
 * - A special "null" bucket for discarded rotations.
 * - A special "constant" bucket for rotations that are just global phases and measurements that are just constants.
 * - A "merge" bucket for each other equivalence class of rotations and measurements under +-ZZ stabilizers.
 * Within each bucket, we fix an interpretation of the signs and record the sign of each rotation/measurement wrt the bucket's sign. E.g. in a merge bucket, fix one of the rotation/measurements A to have +, then let the sign of each other B be the sign of the stabilizer Z_A Z_B.
 *
 * When searching for the goal stabilizers above, we may also have any number of RoleControl qubits in X (flow information). Z (interface information) should only ever occur alongside a non-identity Pauli on that node's interface qubits, but should also be permitted to allow us to search for folds into QFree or Reset nodes (it does not hurt to also allow Zs in general since we will still check for identities on the interface qubits).
 * Sum controls will be present whenever this opportunity for merging is correllated with some classical data. Note that this is correllation; the nodes might be dependent on the classical data like with folding gates within the same conditional or hoisting a gate out of the conditional, or the classical data might be dependent on the node such as with measurement. We should permit non-identity values on controls that are ancestors of the rotation/measurements as well as the classical output of that particular measurement (i.e. we have to undo any inference that used classical propagations). Collectively, these may be SumInControls and SumInPhControls (for conditionals) and SumOutPhControls(for measurements) but never SumOutControls.
 *
 * We will call a qubit "permissible" if it is either a control (role or sum, not SumOutControl) or the input to a discard (QFree or Reset), since these are the qubits with variable inclusion in the stabilizer patterns we are searching for.
 * We begin by taking our entire analysis tableau and reducing it to echelon form, solving for Zs on non-permissible qubits first, then Xs on non-permissibles, then permissibles. The bottom stabs of the resulting tableau are generators for the subgroup of stabilizers over only the permissible qubits. All non-permissible Zs (i.e. the actual components of the patterns we are searching for) exist in the top stabilizers, with at least one unique Z in each stab.
 * Patterns with a single Z (i.e. for the constant and null bucket) must be formed by precisely one stabilizer from the top set, along with some number of stabilizers from the bottom set. Patterns with a pair of Zs (i.e. for merges) can be formed by either one or two stabilizers from the top set, along with some number of stabilizers from the bottom set - in the case where it is a product of two, it is only worth checking those pairs once they have been rejected from matching constant and null patterns.
 * Depending on the non-permissibles involved in the pattern, we may need to enforce identity Paulis in some of the permissibles (i.e. in sum controls for non-ancestor nodes), our current "violations". We have already removed all of those permissibles that are leading qubits in some stabilizers, so all remaining violations must not be leading qubits of any stabilizer. In attempting to resolve violations, we can only re-apply stabilizers whose leading qubits would not reintroduce violations. This hopefully gives us only a relatively small number of stabilizers to search for a combination that resolves all violations. These can be copied out to a new tableau and solved by reducing to echelon form.
 */

pub struct PhaseFoldSettings {
    // When two rotations act on the same data, allow merging them by summing the angles
    merge_rr: bool,
    // When a rotation and a measurement act on the same data, allow converting the rotation to a conditional global phase [NOTE:: global phase not currently implemented, so just removes the rotation]
    merge_rm: bool,
    // When two measurements act on the same data, allow merging them into a single measurement and a copy (possibly with negation)
    merge_mm: bool,
    // When a rotation acts on a stabilizer, allow converting the rotation to a global phase [NOTE:: global phase not currently implemented, so just removes the rotation]
    // If this is false, we still allow such rotations to be merged when multiple exist (subject to the merge_rr and merge_rm settings)
    constant_r: bool,
    // When a measurement acts on a stabilizer, allow replacing it with a constant classical value
    // If this is false, we still allow such measurements to be merged (subject to the merge_rm and merge_mm settings)
    constant_m: bool,
    // When a rotation can be discarded by a QFree or Reset, allow removing it
    // If this is false, we still allow such rotations to be merged (subject to the merge_rr setting)
    null_r: bool,
}

pub struct PhaseFold<H: HugrMut> {
    // Collect rotations into buckets using petgraph::UnionFind.
    // Each integer index corresponds to a rotation gate within the program, with two additional fake rotations to help identify the constant and null buckets.
    bucket_lookup: UnionFind<usize>,
    null_bucket: usize,
    const_bucket: usize,
    // A map between rotation nodes and indices into bucket_lookup.
    rotation_lookup: BiHashMap<H::Node, usize>,
    // For each rotation, look whether it is positive (false) or negative (true) relative to the rest of the bucket.
    // Uses the numerical indexing of bucket_lookup for the rotations.
    rotation_signs: Vec<bool>,
}

impl<H: HugrMut> PhaseFold<H> {
    pub fn new() -> Self {
        PhaseFold {
            bucket_lookup: UnionFind::new(2),
            null_bucket: 0,
            const_bucket: 1,
            rotation_lookup: BiHashMap::new(),
            rotation_signs: vec![false, false],
        }
    }

    pub fn find_folds(&mut self, hugr: &H, summary: &mut StabilizerDataflow<H>) {
        let mut rotations: Vec<usize> = vec![];
        let mut other_non_permissibles: Vec<usize> = vec![];
        let mut sums: Vec<usize> = vec![];
        let mut sum_node_lookup: HashMap<usize, H::Node> = HashMap::new();
        let mut discards: Vec<usize> = vec![];
        let mut rcs: Vec<usize> = vec![];
        for i in 0..summary.tab.nb_qubits {
            match summary.q_index_map.get_by_right(&i).unwrap() {
                DataflowPoint::NodeIn(n, _, _) => {
                    if let OpType::ExtensionOp(extop) = hugr.get_optype(*n) {
                        if let Ok(tkop) = TketOp::from_extension_op(extop) {
                            if tkop == TketOp::Reset || tkop == TketOp::QFree {
                                discards.push(i);
                            } else {
                                other_non_permissibles.push(i);
                            }
                        } else {
                            other_non_permissibles.push(i);
                        }
                    } else {
                        other_non_permissibles.push(i);
                    }
                }
                DataflowPoint::Rotation(n) => {
                    if !self.rotation_lookup.contains_left(n) {
                        self.rotation_lookup
                            .insert(*n, self.bucket_lookup.new_set());
                        self.rotation_signs.push(false);
                    }
                    rotations.push(i);
                }
                DataflowPoint::SumInControl(n, _, _, _, _)
                | DataflowPoint::SumInPhControl(n, _, _, _)
                | DataflowPoint::SumOutPhControl(n, _, _, _) => {
                    sums.push(i);
                    sum_node_lookup.insert(i, *n);
                }
                DataflowPoint::RoleControl(_, _) => {
                    rcs.push(i);
                }
                _ => {
                    other_non_permissibles.push(i);
                }
            }
        }

        let echelon_order: Vec<(usize, PauliXZ)> = chain!(
            rotations.iter().map(|q| (*q, PauliXZ::Z)),
            rotations.iter().map(|q| (*q, PauliXZ::X)),
            other_non_permissibles.iter().map(|q| (*q, PauliXZ::Z)),
            other_non_permissibles.iter().map(|q| (*q, PauliXZ::X)),
            sums.iter().map(|q| (*q, PauliXZ::Z)),
            sums.iter().map(|q| (*q, PauliXZ::X)),
            // By solving for discards, we obtain constants over nulls whenever possible.
            // This is because if a true stabilizer is passed into a discard, we get a stabilizer e.g. Zqfree.in (+ role controls) which could be added to anything and therefore convert matches for a constant pattern to matches for a null pattern.
            discards.iter().map(|q| (*q, PauliXZ::Z)),
            discards.iter().map(|q| (*q, PauliXZ::X)),
        )
        .collect_vec();
        summary.tab.echelon(&echelon_order);

        // For each qubit Q in the tableau, return the index of the stabilizer whose leading qubit is Q, if one exists
        let mut leader_to_stab: BTreeMap<(usize, PauliXZ), usize> = BTreeMap::new();
        let mut current_stab = 0;
        for (q, p) in echelon_order {
            if (p == PauliXZ::Z && summary.tab.z[current_stab].get(q))
                || (p == PauliXZ::X && summary.tab.x[current_stab].get(q))
            {
                leader_to_stab.insert((q, p), current_stab);
                current_stab += 1;
                if current_stab == summary.tab.nb_stabs {
                    break;
                }
            }
        }

        // Any matches for constant or null patterns (i.e. a single rotation with Z), we should be able to just read it off from the tableau
        let mut non_permissible_mask = BitVector::new(summary.tab.nb_qubits);
        for q in rotations.iter() {
            non_permissible_mask.xor_bit(*q);
        }
        for q in other_non_permissibles {
            non_permissible_mask.xor_bit(q);
        }

        let mut ancestors: HashMap<H::Node, HashSet<H::Node>> = HashMap::new();
        let mut unused_stabs: Vec<usize> = vec![];
        for q in rotations.iter() {
            let DataflowPoint::Rotation(r_node) = summary.q_index_map.get_by_right(q).unwrap()
            else {
                unreachable!()
            };
            let mut anc_set: HashSet<H::Node> = HashSet::new();
            let mut anc: H::Node = *r_node;
            while let Some(anc_parent) = hugr.get_parent(anc) {
                anc_set.insert(anc_parent);
                anc = anc_parent;
            }
            ancestors.insert(*r_node, anc_set.clone());
            if let Some(stab) = leader_to_stab.get(&(*q, PauliXZ::Z)) {
                // Check that no other non-permissibles are in the stabilizer
                let mut z_with_mask = summary.tab.z[*stab].clone();
                z_with_mask.xor_bit(*q);
                z_with_mask.and(&non_permissible_mask);
                let mut x_with_mask = summary.tab.x[*stab].clone();
                x_with_mask.and(&non_permissible_mask);
                if z_with_mask == x_with_mask
                    && z_with_mask == BitVector::new(summary.tab.nb_qubits)
                {
                    let mut solver = SymplecticTableau::new(summary.tab.nb_qubits);
                    let mut cols_to_solve: Vec<(usize, PauliXZ)> = vec![(*q, PauliXZ::Z)];
                    solver.add_stab(
                        summary.tab.z[*stab].clone(),
                        summary.tab.x[*stab].clone(),
                        summary.tab.signs.get(*stab),
                    );
                    // Matched the basic pattern, look up the node and check the non-identity permissibles for violations
                    for (sum_q, sum_node) in sum_node_lookup.iter() {
                        if anc_set.contains(sum_node) {
                            if let Some(sum_stab) = leader_to_stab.get(&(*sum_q, PauliXZ::Z)) {
                                // If sum_node is an ancestor (i.e. it is a conditional that contains r_node) and there exists a stabilizer lead by sum_q, we may use that stabilizer for removing other violations
                                solver.add_stab(
                                    summary.tab.z[*sum_stab].clone(),
                                    summary.tab.x[*sum_stab].clone(),
                                    summary.tab.signs.get(*sum_stab),
                                );
                            }
                            if let Some(sum_stab) = leader_to_stab.get(&(*sum_q, PauliXZ::X)) {
                                solver.add_stab(
                                    summary.tab.z[*sum_stab].clone(),
                                    summary.tab.x[*sum_stab].clone(),
                                    summary.tab.signs.get(*sum_stab),
                                );
                            }
                        } else {
                            // If sum_node is neither, then add it to the list of qubits to eliminate
                            cols_to_solve.push((*sum_q, PauliXZ::Z));
                            cols_to_solve.push((*sum_q, PauliXZ::X));
                        }
                    }
                    // Since first entry in cols_to_solve is the rotation, we can just look at the first stabilizer and check whether we were able to remove all violations
                    solver.echelon(&cols_to_solve);
                    if cols_to_solve.iter().all(|(solve_q, p)| {
                        solve_q == q
                            || match p {
                                PauliXZ::X => !solver.x[0].get(*solve_q),
                                PauliXZ::Z => !solver.z[0].get(*solve_q),
                            }
                    }) {
                        // All were identity, we have successfully matched a pattern!
                        let rotation_bucket = self.rotation_lookup.get_by_left(r_node).unwrap();
                        // This could have still been either the constant or null pattern, which we can distinguish by looking at the discard qubits
                        if discards
                            .iter()
                            .any(|dq| solver.x[0].get(*dq) || solver.z[0].get(*dq))
                        {
                            // Uses at least one discard, so this is a null pattern
                            self.bucket_lookup.union(self.null_bucket, *rotation_bucket);
                        } else {
                            // Does not use any discards, so this is a constant pattern
                            self.bucket_lookup
                                .union(self.const_bucket, *rotation_bucket);
                            // No merging has happened yet, so only this rotation's sign needs updating
                            self.rotation_signs[*rotation_bucket] = solver.signs.get(0);
                        }
                    } else {
                        unused_stabs.push(*stab);
                    }
                } else {
                    unused_stabs.push(*stab);
                }
            }
        }

        // Now to match against merge patterns, i.e. a ZZ stabilizer between two rotation qubits.
        // These can occur in the same stabilizer, where one is the leader and the other is not a leader of any stabilizer, or as the product of two stabilizers, where the rotations are the leaders of each stabilizer and all other rotations cancel.
        // If a stabilizer P is found to contain Z_leader Z_other itself, we don't need to check any combinations with it - if there is some other Q such that PQ is Z_leader Z_leader' then Q must be Z_leader Z_other so we will find the merge just by looking at Q.
        // Similarly, if PQ gives a merge Z_LP Z_LQ, then we don't need to check other combinations with Q - if QR is Z_LQ Z_LR, then PR = (PQ)QR = Z_LP Z_LR
        let mut used_in_merges: HashSet<usize> = HashSet::new();
        // Start by finding those stabilizers with ZZ on themselves
        for stab in unused_stabs.iter() {
            let mut x_with_mask = summary.tab.x[*stab].clone();
            x_with_mask.and(&non_permissible_mask);
            if x_with_mask == BitVector::new(summary.tab.nb_qubits) {
                let mut z_with_mask = summary.tab.z[*stab].clone();
                z_with_mask.and(&non_permissible_mask);
                if let [r0q, r1q] = z_with_mask.get_all_ones(summary.tab.nb_qubits).as_slice() {
                    if let DataflowPoint::Rotation(rot0) =
                        summary.q_index_map.get_by_right(r0q).unwrap()
                    {
                        if let DataflowPoint::Rotation(rot1) =
                            summary.q_index_map.get_by_right(r1q).unwrap()
                        {
                            // This stabilizer matches the pattern on its own, so check for no violations
                            // Find those nodes that are an ancestor of either rot0 or rot1
                            let mut anc_set: HashSet<H::Node> = HashSet::new();
                            let mut anc: H::Node = *rot0;
                            while let Some(anc_parent) = hugr.get_parent(anc) {
                                anc_set.insert(anc_parent);
                                anc = anc_parent;
                            }
                            anc = *rot1;
                            while let Some(anc_parent) = hugr.get_parent(anc) {
                                anc_set.insert(anc_parent);
                                anc = anc_parent;
                            }

                            // Check if we can remove any violations, same as before
                            let mut solver = SymplecticTableau::new(summary.tab.nb_qubits);
                            // Only need to solve for rot0, since no other stabilizers of solver will contain rot1
                            let mut cols_to_solve: Vec<(usize, PauliXZ)> = vec![(*r0q, PauliXZ::Z)];
                            solver.add_stab(
                                summary.tab.z[*stab].clone(),
                                summary.tab.x[*stab].clone(),
                                summary.tab.signs.get(*stab),
                            );
                            for (sum_q, sum_node) in sum_node_lookup.iter() {
                                if anc_set.contains(sum_node) {
                                    if let Some(sum_stab) =
                                        leader_to_stab.get(&(*sum_q, PauliXZ::Z))
                                    {
                                        solver.add_stab(
                                            summary.tab.z[*sum_stab].clone(),
                                            summary.tab.x[*sum_stab].clone(),
                                            summary.tab.signs.get(*sum_stab),
                                        );
                                    }
                                    if let Some(sum_stab) =
                                        leader_to_stab.get(&(*sum_q, PauliXZ::X))
                                    {
                                        solver.add_stab(
                                            summary.tab.z[*sum_stab].clone(),
                                            summary.tab.x[*sum_stab].clone(),
                                            summary.tab.signs.get(*sum_stab),
                                        );
                                    }
                                } else {
                                    cols_to_solve.push((*sum_q, PauliXZ::Z));
                                    cols_to_solve.push((*sum_q, PauliXZ::X));
                                }
                            }
                            solver.echelon(&cols_to_solve);
                            if cols_to_solve.iter().all(|(solve_q, p)| {
                                solve_q == r0q
                                    || match p {
                                        PauliXZ::X => !solver.x[0].get(*solve_q),
                                        PauliXZ::Z => !solver.z[0].get(*solve_q),
                                    }
                            }) {
                                // All were identity, we have successfully matched a pattern!
                                let r0_bucket = self.rotation_lookup.get_by_left(rot0).unwrap();
                                let r1_bucket = self.rotation_lookup.get_by_left(rot1).unwrap();
                                self.bucket_lookup.union(*r0_bucket, *r1_bucket);
                                // Whichever one of these is the leader is guaranteed to not have been merged into any other bucket yet, so we can use the sign of the stabilizer as defining its phase relative to the other one
                                if leader_to_stab.contains_key(&(*r0q, PauliXZ::Z)) {
                                    self.rotation_signs[*r0_bucket] = solver.signs.get(0);
                                } else {
                                    self.rotation_signs[*r1_bucket] = solver.signs.get(0);
                                }
                                used_in_merges.insert(*stab);
                            }
                        }
                    }
                }
            }
        }

        // Then look at pairs of the remaining stabilizers
        for (i, stab0) in unused_stabs.iter().enumerate() {
            if used_in_merges.contains(stab0) {
                continue;
            }
            for j in 0..i {
                let stab1 = unused_stabs.get(j).unwrap();
                if used_in_merges.contains(stab1) {
                    continue;
                }
                // Compute product
                let mut solver = SymplecticTableau::new(summary.tab.nb_qubits);
                solver.add_stab(
                    summary.tab.z[*stab0].clone(),
                    summary.tab.x[*stab0].clone(),
                    summary.tab.signs.get(*stab0),
                );
                solver.add_stab(
                    summary.tab.z[*stab1].clone(),
                    summary.tab.x[*stab1].clone(),
                    summary.tab.signs.get(*stab1),
                );
                solver.stab_mult(1, 0, 0);
                solver.delete_stab(1);

                let mut x_with_mask = solver.x[0].clone();
                x_with_mask.and(&non_permissible_mask);
                if x_with_mask == BitVector::new(summary.tab.nb_qubits) {
                    let mut z_with_mask = solver.z[0].clone();
                    z_with_mask.and(&non_permissible_mask);
                    if let [r0q, r1q] = z_with_mask.get_all_ones(summary.tab.nb_qubits).as_slice() {
                        // Both are guaranteed to be Rotations since they must be the leaders of stab0 and stab1
                        let DataflowPoint::Rotation(rot0) =
                            summary.q_index_map.get_by_right(r0q).unwrap()
                        else {
                            unreachable!()
                        };
                        let DataflowPoint::Rotation(rot1) =
                            summary.q_index_map.get_by_right(r1q).unwrap()
                        else {
                            unreachable!()
                        };
                        let mut anc_set: HashSet<H::Node> = HashSet::new();
                        let mut anc: H::Node = *rot0;
                        while let Some(anc_parent) = hugr.get_parent(anc) {
                            anc_set.insert(anc_parent);
                            anc = anc_parent;
                        }
                        anc = *rot1;
                        while let Some(anc_parent) = hugr.get_parent(anc) {
                            anc_set.insert(anc_parent);
                            anc = anc_parent;
                        }
                        let mut cols_to_solve: Vec<(usize, PauliXZ)> = vec![(*r0q, PauliXZ::Z)];
                        for (sum_q, sum_node) in sum_node_lookup.iter() {
                            if anc_set.contains(sum_node) {
                                if let Some(sum_stab) = leader_to_stab.get(&(*sum_q, PauliXZ::Z)) {
                                    solver.add_stab(
                                        summary.tab.z[*sum_stab].clone(),
                                        summary.tab.x[*sum_stab].clone(),
                                        summary.tab.signs.get(*sum_stab),
                                    );
                                }
                                if let Some(sum_stab) = leader_to_stab.get(&(*sum_q, PauliXZ::X)) {
                                    solver.add_stab(
                                        summary.tab.z[*sum_stab].clone(),
                                        summary.tab.x[*sum_stab].clone(),
                                        summary.tab.signs.get(*sum_stab),
                                    );
                                }
                            } else {
                                cols_to_solve.push((*sum_q, PauliXZ::Z));
                                cols_to_solve.push((*sum_q, PauliXZ::X));
                            }
                        }
                        solver.echelon(&cols_to_solve);
                        if cols_to_solve.iter().all(|(solve_q, p)| {
                            solve_q == r0q
                                || match p {
                                    PauliXZ::X => !solver.x[0].get(*solve_q),
                                    PauliXZ::Z => !solver.z[0].get(*solve_q),
                                }
                        }) {
                            // All were identity, we have successfully matched a pattern!
                            let r0_bucket = self.rotation_lookup.get_by_left(rot0).unwrap();
                            let r1_bucket = self.rotation_lookup.get_by_left(rot1).unwrap();
                            self.bucket_lookup.union(*r0_bucket, *r1_bucket);
                            // Whichever one of these is the leader of stab1 is guaranteed to not have been merged into any other bucket yet, so we can use the sign of the stabilizer as defining its phase relative to the leader of stab0
                            if leader_to_stab.get(&(*r0q, PauliXZ::Z)).unwrap() == stab1 {
                                self.rotation_signs[*r0_bucket] = solver.signs.get(0);
                            } else {
                                self.rotation_signs[*r1_bucket] = solver.signs.get(0);
                            }
                            used_in_merges.insert(*stab1);
                        }
                    }
                }
            }
        }
    }

    pub fn apply_folds(&self, hugr: &mut H, pfsettings: &PhaseFoldSettings) {
        // Split each bucket into rotations and measurements
        // IndexSet maintains the insertion order so we can iterate through and make changes to the hugr in a deterministic order (it breaks if you remove individual elements from the set, but in this routine the only removals occur when clearing the entire set)
        let mut bucket_rotations: Vec<IndexSet<H::Node>> =
            vec![IndexSet::new(); self.bucket_lookup.len()];
        let mut bucket_measures: Vec<IndexSet<H::Node>> =
            vec![IndexSet::new(); self.bucket_lookup.len()];
        let mut bucket_measure_frees: Vec<IndexSet<H::Node>> =
            vec![IndexSet::new(); self.bucket_lookup.len()];
        for (r_node, r_index) in &self.rotation_lookup {
            let r_bucket = self.bucket_lookup.find(*r_index);
            let OpType::ExtensionOp(op) = hugr.get_optype(*r_node) else {
                unreachable!()
            };
            let Ok(tkop) = TketOp::from_extension_op(op) else {
                unreachable!()
            };
            match tkop {
                TketOp::Measure => {
                    bucket_measures[r_bucket].insert(*r_node);
                }
                TketOp::MeasureFree => {
                    bucket_measure_frees[r_bucket].insert(*r_node);
                }
                _ => {
                    bucket_rotations[r_bucket].insert(*r_node);
                }
            }
        }

        if pfsettings.null_r {
            for r_node in bucket_rotations.get(self.null_bucket).unwrap() {
                // We can just remove this rotation gate
                let (q_pred, q_pred_port) = hugr
                    .single_linked_output(*r_node, IncomingPort::from(0))
                    .unwrap();
                let (q_succ, q_succ_port) = hugr
                    .single_linked_input(*r_node, OutgoingPort::from(0))
                    .unwrap();
                hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                hugr.remove_node(*r_node);
            }
            bucket_rotations[self.null_bucket].clear();
            // If pfsetting.null_r == false, the bucket will still be there for us to iterate through and merge the rotations
        }

        if pfsettings.constant_r {
            // This is currently the same as for null; this may change in future to replace the rotations with global phases
            for r_node in bucket_rotations.get(self.const_bucket).unwrap() {
                let (q_pred, q_pred_port) = hugr
                    .single_linked_output(*r_node, IncomingPort::from(0))
                    .unwrap();
                let (q_succ, q_succ_port) = hugr
                    .single_linked_input(*r_node, OutgoingPort::from(0))
                    .unwrap();
                hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                hugr.remove_node(*r_node);
            }
            bucket_rotations[self.const_bucket].clear();
        }

        if pfsettings.constant_m {
            for m_node in bucket_measures.get(self.const_bucket).unwrap() {
                // Replace the Measure with an identity wire and a constant classical value
                let (q_pred, q_pred_port) = hugr
                    .single_linked_output(*m_node, IncomingPort::from(0))
                    .unwrap();
                let (q_succ, q_succ_port) = hugr
                    .single_linked_input(*m_node, OutgoingPort::from(0))
                    .unwrap();
                hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                let polarity = match self
                    .rotation_signs
                    .get(*self.rotation_lookup.get_by_left(m_node).unwrap())
                {
                    Some(p) => *p,
                    None => false,
                };
                let parent = hugr.get_parent(*m_node).unwrap();
                let const_val_node = hugr.add_node_with_parent(
                    parent,
                    Tag::new(if polarity { 1 } else { 0 }, vec![TypeRow::new(); 2]),
                );
                for (c_succ, c_succ_port) in hugr
                    .linked_inputs(*m_node, OutgoingPort::from(1))
                    .collect_vec()
                {
                    hugr.connect(const_val_node, OutgoingPort::from(0), c_succ, c_succ_port);
                }
                hugr.remove_node(*m_node);
            }
            bucket_measures[self.const_bucket].clear();
            for mf_node in bucket_measure_frees.get(self.const_bucket).unwrap() {
                // Replace the MeasureFree with a QFree and a constant classical value
                let (q_pred, q_pred_port) = hugr
                    .single_linked_output(*mf_node, IncomingPort::from(0))
                    .unwrap();
                let parent = hugr.get_parent(*mf_node).unwrap();
                let qfree_node = hugr.add_node_with_parent(parent, TketOp::QFree);
                hugr.connect(q_pred, q_pred_port, qfree_node, IncomingPort::from(0));
                let polarity = match self
                    .rotation_signs
                    .get(*self.rotation_lookup.get_by_left(mf_node).unwrap())
                {
                    Some(p) => *p,
                    None => false,
                };
                let measure_result = ConstBool::new(polarity);
                let const_val_node =
                    hugr.add_node_with_parent(parent, Const::new(measure_result.into()));
                let load_const_node = hugr.add_node_with_parent(
                    parent,
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
                    .linked_inputs(*mf_node, OutgoingPort::from(0))
                    .collect_vec()
                {
                    hugr.connect(load_const_node, OutgoingPort::from(0), c_succ, c_succ_port);
                }
                hugr.remove_node(*mf_node);
            }
            bucket_measure_frees[self.const_bucket].clear();
            // Similarly, if pfsettings.constant_m == false, then the buckets will still be there for us to iterate through and merge the measurements with each other or with rotations
        }

        if pfsettings.merge_rm {
            for b in 0..self.bucket_lookup.len() {
                if !bucket_measures[b].is_empty() || !bucket_measure_frees[b].is_empty() {
                    for r_node in bucket_rotations.get(b).unwrap() {
                        // We can just remove this rotation gate
                        let (q_pred, q_pred_port) = hugr
                            .single_linked_output(*r_node, IncomingPort::from(0))
                            .unwrap();
                        let (q_succ, q_succ_port) = hugr
                            .single_linked_input(*r_node, OutgoingPort::from(0))
                            .unwrap();
                        hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                        hugr.remove_node(*r_node);
                    }
                    bucket_rotations[b].clear();
                }
            }
        }

        // For merging rotations, we need to merge into the latest one so all angles are available.
        // For merging measurements, we need to merge into the earliest one so the outcome is always available for all uses.
        // Both of these require knowing the relative ordering between the nodes, factoring in any classical dependencies which aren't tracked in the stabilizer analysis.
        // We use a deterministic topological sorting algorithm to provide a total ordering over all nodes with the same parent that respects all data dependencies.
        // When rotations don't have the same parent, we need to hoist one or both of them which is a more involved rewrite and is left for future implementation.
        // We will compute the topological orderings of a sibling graph lazily when we first encounter a rotation in that sibling graph.
        let mut topo_order_lookup: HashMap<H::Node, usize> = HashMap::new();

        if pfsettings.merge_rr {
            for bucket in bucket_rotations {
                // Subdivide the bucket based on parent node to merge rotations with the same parent
                let mut bucket_by_parent: IndexMap<H::Node, Vec<H::Node>> = IndexMap::new();
                for node in bucket {
                    let parent = hugr.get_parent(node).unwrap();
                    if let Some(buck) = bucket_by_parent.get_mut(&parent) {
                        buck.push(node);
                    } else {
                        bucket_by_parent.insert(parent, vec![node]);
                    }
                }
                for (parent, sibling_bucket) in bucket_by_parent {
                    // If none or only one rotation, nothing to merge
                    if sibling_bucket.len() < 2 {
                        continue;
                    }
                    let latest_rotation_node = sibling_bucket
                        .iter()
                        .max_by_key(|n| {
                            if !topo_order_lookup.contains_key(n) {
                                let (region, node_map) = hugr.region_portgraph(parent);
                                for (i, ni) in toposort(&region, None).unwrap().iter().enumerate() {
                                    topo_order_lookup.insert(node_map.from_portgraph(*ni), i);
                                }
                            }
                            topo_order_lookup.get(n).unwrap().clone()
                        })
                        .unwrap();
                    let latest_rotation_polarity = self
                        .rotation_signs
                        .get(
                            *self
                                .rotation_lookup
                                .get_by_left(latest_rotation_node)
                                .unwrap(),
                        )
                        .unwrap_or(&false);
                    let mut acc_static = 0.;
                    let mut acc_dynamic: Option<(H::Node, OutgoingPort)> = None;
                    for r_node in sibling_bucket.iter() {
                        let OpType::ExtensionOp(op) = hugr.get_optype(*r_node) else {
                            unreachable!()
                        };
                        let Ok(tkop) = TketOp::from_extension_op(op) else {
                            unreachable!()
                        };
                        let polarity = self
                            .rotation_signs
                            .get(*self.rotation_lookup.get_by_left(r_node).unwrap())
                            .unwrap_or(&false);
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
                                    .single_linked_output(*r_node, IncomingPort::from(1))
                                    .unwrap();
                                let (value_node, value_port) =
                                    if polarity == latest_rotation_polarity {
                                        (source_node, source_port)
                                    } else {
                                        let to_float = hugr
                                            .add_node_with_parent(parent, RotationOp::to_halfturns);
                                        hugr.connect(
                                            source_node,
                                            source_port,
                                            to_float,
                                            IncomingPort::from(0),
                                        );
                                        let neg_node =
                                            hugr.add_node_with_parent(parent, FloatOps::fneg);
                                        hugr.connect(
                                            to_float,
                                            OutgoingPort::from(0),
                                            neg_node,
                                            IncomingPort::from(0),
                                        );
                                        let from_float = hugr.add_node_with_parent(
                                            parent,
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
                                            hugr.add_node_with_parent(parent, RotationOp::radd);
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
                                // Other gate types are not considered in phase folding.
                                // This includes Clifford phase gates (S, Sdg, Z) which would be abstracted away
                            }
                        }
                        if r_node != latest_rotation_node {
                            let (q_pred, q_pred_port) = hugr
                                .single_linked_output(*r_node, IncomingPort::from(0))
                                .unwrap();
                            let (q_succ, q_succ_port) = hugr
                                .single_linked_input(*r_node, OutgoingPort::from(0))
                                .unwrap();
                            hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                            hugr.remove_node(*r_node);
                        }
                    }
                    // Now to merge the accumulated rotations into the latest node.
                    // We pick the rotation gate Rx/Ry/Rz to match the basis of the original gate.
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
                                    parent,
                                    Into::<Const>::into(Value::extension(
                                        ConstRotation::new(acc_static % 2.).unwrap(),
                                    )),
                                );
                                let load_const_node = hugr.add_node_with_parent(
                                    parent,
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
                                let add_node = hugr.add_node_with_parent(parent, RotationOp::radd);
                                hugr.connect(
                                    source_node,
                                    source_port,
                                    add_node,
                                    IncomingPort::from(0),
                                );
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
                            let new_rotation_node = hugr.add_node_with_parent(parent, new_optype);
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
                                            c_t_seq = chain!([TketOp::H], c_t_seq, [TketOp::H])
                                                .collect_vec();
                                        }
                                        TketOp::Ry => {
                                            c_t_seq = chain!([TketOp::V], c_t_seq, [TketOp::Vdg])
                                                .collect_vec();
                                        }
                                        _ => {}
                                    }
                                    for gate in c_t_seq {
                                        let gate_node = hugr.add_node_with_parent(parent, gate);
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
                                        parent,
                                        Into::<Const>::into(Value::extension(
                                            ConstRotation::new(acc_static % 2.).unwrap(),
                                        )),
                                    );
                                    let load_const_node = hugr.add_node_with_parent(
                                        parent,
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
                                        hugr.add_node_with_parent(parent, new_optype);
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

        if pfsettings.merge_mm {
            for (bucket_m, bucket_mf) in bucket_measures.iter().zip(bucket_measure_frees) {
                let mut measures_by_parent: IndexMap<H::Node, Vec<H::Node>> = IndexMap::new();
                for node in bucket_m {
                    let parent = hugr.get_parent(*node).unwrap();
                    if let Some(buck) = measures_by_parent.get_mut(&parent) {
                        buck.push(*node);
                    } else {
                        measures_by_parent.insert(parent, vec![*node]);
                    }
                }
                let mut measure_frees_by_parent: IndexMap<H::Node, Vec<H::Node>> = IndexMap::new();
                for node in bucket_mf {
                    let parent = hugr.get_parent(node).unwrap();
                    if let Some(buck) = measure_frees_by_parent.get_mut(&parent) {
                        buck.push(node);
                    } else {
                        measure_frees_by_parent.insert(parent, vec![node]);
                    }
                }
                let mut all_parents: HashSet<H::Node> = HashSet::new();
                all_parents.extend(measures_by_parent.keys());
                all_parents.extend(measure_frees_by_parent.keys());
                for parent in all_parents {
                    let empty: Vec<H::Node> = vec![]; // Just needed to give references below
                    let ms = measures_by_parent.get(&parent).unwrap_or(&empty);
                    let mfs = measure_frees_by_parent.get(&parent).unwrap_or(&empty);
                    if ms.len() + mfs.len() < 2 {
                        continue;
                    }
                    let earliest_measure: Option<&H::Node> = ms.iter().max_by_key(|n| {
                        if !topo_order_lookup.contains_key(n) {
                            let (region, node_map) = hugr.region_portgraph(parent);
                            for (i, ni) in toposort(&region, None).unwrap().iter().enumerate() {
                                topo_order_lookup.insert(node_map.from_portgraph(*ni), i);
                            }
                        }
                        topo_order_lookup.get(n).unwrap().clone()
                    });
                    let earliest_measure_free: Option<&H::Node> = mfs.iter().max_by_key(|n| {
                        if !topo_order_lookup.contains_key(n) {
                            let (region, node_map) = hugr.region_portgraph(parent);
                            for (i, ni) in toposort(&region, None).unwrap().iter().enumerate() {
                                topo_order_lookup.insert(node_map.from_portgraph(*ni), i);
                            }
                        }
                        topo_order_lookup.get(n).unwrap().clone()
                    });
                    let earliest_node: H::Node = match earliest_measure {
                        Some(m) => match earliest_measure_free {
                            Some(mf) => {
                                if topo_order_lookup.get(&m).unwrap()
                                    < topo_order_lookup.get(&mf).unwrap()
                                {
                                    *m
                                } else {
                                    *mf
                                }
                            }
                            None => *m,
                        },
                        None => {
                            // earliest_measure_free must contain a value
                            *earliest_measure_free.unwrap()
                        }
                    };
                    let earliest_polarity = self
                        .rotation_signs
                        .get(*self.rotation_lookup.get_by_left(&earliest_node).unwrap())
                        .unwrap_or(&false);
                    // Measure gates produce a hugr bool_t (i.e. sum of two units) whereas MeasureFree gates produce a tket bool_type (i.e. an opaque boolean)
                    // Depending on which types are needed, maintain at most one source of each type
                    let mut prelude_bool_result_loc: Option<(H::Node, OutgoingPort)> =
                        if earliest_measure == Some(&&earliest_node) {
                            Some((earliest_node, OutgoingPort::from(1)))
                        } else {
                            None
                        };
                    let mut tket_bool_result_loc: Option<(H::Node, OutgoingPort)> =
                        if earliest_measure_free == Some(&&earliest_node) {
                            Some((earliest_node, OutgoingPort::from(0)))
                        } else {
                            None
                        };
                    for m_node in ms {
                        if *m_node == earliest_node {
                            continue;
                        }
                        let (q_pred, q_pred_port) = hugr
                            .single_linked_output(*m_node, IncomingPort::from(0))
                            .unwrap();
                        let (q_succ, q_succ_port) = hugr
                            .single_linked_input(*m_node, OutgoingPort::from(0))
                            .unwrap();
                        hugr.connect(q_pred, q_pred_port, q_succ, q_succ_port);
                        if prelude_bool_result_loc.is_none() {
                            // Earliest measure was a MeasureFree
                            // Add a read op to unpack the tket bool outcome
                            let read_node = hugr.add_node_with_parent(parent, BoolOp::read);
                            hugr.connect(
                                earliest_node,
                                OutgoingPort::from(0),
                                read_node,
                                IncomingPort::from(0),
                            );
                            prelude_bool_result_loc = Some((read_node, OutgoingPort::from(0)));
                        }
                        let (result_node, result_port) = prelude_bool_result_loc.unwrap();
                        let polarity = self
                            .rotation_signs
                            .get(*self.rotation_lookup.get_by_left(m_node).unwrap())
                            .unwrap_or(&false);
                        if polarity == earliest_polarity {
                            for (c_succ, c_succ_port) in hugr
                                .linked_inputs(*m_node, OutgoingPort::from(1))
                                .collect_vec()
                            {
                                hugr.connect(result_node, result_port, c_succ, c_succ_port);
                            }
                        } else {
                            let not_node = hugr.add_node_with_parent(parent, LogicOp::Not);
                            hugr.connect(result_node, result_port, not_node, IncomingPort::from(0));
                            for (c_succ, c_succ_port) in hugr
                                .linked_inputs(*m_node, OutgoingPort::from(1))
                                .collect_vec()
                            {
                                hugr.connect(not_node, OutgoingPort::from(0), c_succ, c_succ_port);
                            }
                        }
                        hugr.remove_node(*m_node);
                    }
                    for mf_node in mfs {
                        if *mf_node == earliest_node {
                            continue;
                        }
                        let (q_pred, q_pred_port) = hugr
                            .single_linked_output(*mf_node, IncomingPort::from(0))
                            .unwrap();
                        let qfree_node = hugr.add_node_with_parent(parent, TketOp::QFree);
                        hugr.connect(q_pred, q_pred_port, qfree_node, IncomingPort::from(0));
                        if tket_bool_result_loc.is_none() {
                            // Earliest measure was a Measure
                            // Add a make_opaque op to pack the prelude::bool into a tket::bool
                            let make_opaque_node =
                                hugr.add_node_with_parent(parent, BoolOp::make_opaque);
                            hugr.connect(
                                earliest_node,
                                OutgoingPort::from(1),
                                make_opaque_node,
                                IncomingPort::from(0),
                            );
                            tket_bool_result_loc = Some((make_opaque_node, OutgoingPort::from(0)));
                        }
                        let (result_node, result_port) = tket_bool_result_loc.unwrap();
                        let polarity = self
                            .rotation_signs
                            .get(*self.rotation_lookup.get_by_left(mf_node).unwrap())
                            .unwrap_or(&false);
                        if polarity == earliest_polarity {
                            for (c_succ, c_succ_port) in hugr
                                .linked_inputs(*mf_node, OutgoingPort::from(0))
                                .collect_vec()
                            {
                                hugr.connect(result_node, result_port, c_succ, c_succ_port);
                            }
                        } else {
                            let not_node = hugr.add_node_with_parent(parent, BoolOp::not);
                            hugr.connect(result_node, result_port, not_node, IncomingPort::from(0));
                            for (c_succ, c_succ_port) in hugr
                                .linked_inputs(*mf_node, OutgoingPort::from(0))
                                .collect_vec()
                            {
                                hugr.connect(not_node, OutgoingPort::from(0), c_succ, c_succ_port);
                            }
                        }
                        hugr.remove_node(*mf_node);
                    }
                }
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
        ops::handle::NodeHandle,
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
        phase_fold_complete::{PhaseFold, PhaseFoldSettings},
        stabilizer_dataflow_complete::{DataflowFlowDetail, DataflowSettings, SDFAnalysis},
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

    const MAX_PFSETTINGS: PhaseFoldSettings = PhaseFoldSettings {
        merge_rr: true,
        merge_rm: true,
        merge_mm: true,
        constant_r: true,
        constant_m: true,
        null_r: true,
    };

    #[test]
    fn test_empty() {
        let builder = FunctionBuilder::new("empty", endo_sig(vec![])).unwrap();
        let mut hugr = builder.finish_hugr().unwrap();
        assert_eq!(hugr.num_nodes(), 4);
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        assert!(pf.rotation_lookup.is_empty());
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
        assert_eq!(hugr.num_nodes(), 4);
    }

    #[test]
    fn test_identity() {
        let builder =
            FunctionBuilder::new("identity", endo_sig(vec![usize_t(), qb_t(), qb_t()])).unwrap();
        let [i, qb0, qb1] = builder.input_wires_arr();
        let mut hugr = builder.finish_hugr_with_outputs([i, qb0, qb1]).unwrap();
        assert_eq!(hugr.num_nodes(), 4);
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        assert!(pf.rotation_lookup.is_empty());
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
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
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
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
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
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
        let loop_builder = builder
            .tail_loop_builder(
                [],
                [(qb_t(), qb0), (qb_t(), qb1), (bool_t(), b)],
                type_row![],
            )
            .unwrap();
        let [loop_qb0, loop_qb1, loop_b] = loop_builder.input_wires_arr();
        let tl = loop_builder
            .finish_with_outputs(loop_b, [loop_qb1, loop_qb0, loop_b])
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
        assert_eq!(hugr.num_nodes(), 13);
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
        assert_eq!(hugr.num_nodes(), 11);
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
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
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
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        let pfsettings = PhaseFoldSettings {
            merge_rr: true,
            merge_rm: true,
            merge_mm: false,
            constant_r: true,
            constant_m: false,
            null_r: true,
        };
        pf.apply_folds(&mut hugr, &pfsettings);
        // The only gate removed should be the Rz, Rx, T, Tdg
        assert_eq!(hugr.num_nodes(), 20);
        assert!(hugr.validate().is_ok());
    }

    #[test]
    fn test_merge_measure_cases_remove_measures() {
        let mut hugr = merge_measure_cases();
        assert_eq!(hugr.num_nodes(), 24);
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut default_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
        // Should remove all 4 rotations, two Measures and two MeasureFrees
        // Each removed MeasureFree adds back a QFree
        // The constant Measure is replaced by a Tag, and the constant MeasureFree is replaced by a Const and a LoadConstant
        // For the other removed measures, we get a make_opaque and read for the conversions, plus a negation
        assert!(hugr.validate().is_ok());
        assert_eq!(hugr.num_nodes(), 24);
    }

    #[test]
    fn test_merge_through_teleportation() {
        // Without tracking the classical data, it is impossible to commute any stabilizer flow through the qubit teleportation circuit.
        // But with our handling of boolean data via Sum types, we should be able to infer that any gate can move through.
        let mut builder = FunctionBuilder::new("teleport", endo_sig(vec![qb_t()])).unwrap();
        let [qb0] = builder.input_wires_arr();
        let t = builder.add_dataflow_op(TketOp::T, [qb0]).unwrap();
        let alloc0 = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let alloc1 = builder.add_dataflow_op(TketOp::QAlloc, []).unwrap();
        let h0 = builder
            .add_dataflow_op(TketOp::H, [alloc0.out_wire(0)])
            .unwrap();
        let cx0 = builder
            .add_dataflow_op(TketOp::CX, [h0.out_wire(0), alloc1.out_wire(0)])
            .unwrap();
        let cx1 = builder
            .add_dataflow_op(TketOp::CX, [t.out_wire(0), cx0.out_wire(0)])
            .unwrap();
        let h1 = builder
            .add_dataflow_op(TketOp::H, [cx1.out_wire(0)])
            .unwrap();
        let meas0 = builder
            .add_dataflow_op(TketOp::MeasureFree, [h1.out_wire(0)])
            .unwrap();
        let meas0_convert = builder
            .add_dataflow_op(BoolOp::read, [meas0.out_wire(0)])
            .unwrap();
        let meas1 = builder
            .add_dataflow_op(TketOp::MeasureFree, [cx1.out_wire(1)])
            .unwrap();
        let meas1_convert = builder
            .add_dataflow_op(BoolOp::read, [meas1.out_wire(0)])
            .unwrap();
        let mut cond0_builder = builder
            .conditional_builder(
                ([type_row![], type_row![]], meas1_convert.out_wire(0)),
                [(qb_t(), cx0.out_wire(1))],
                vec![qb_t()].into(),
            )
            .unwrap();
        let cond00_builder = cond0_builder.case_builder(0).unwrap();
        let [c0q] = cond00_builder.input_wires_arr();
        cond00_builder.finish_with_outputs([c0q]).ok();
        let mut cond01_builder = cond0_builder.case_builder(1).unwrap();
        let [c1q] = cond01_builder.input_wires_arr();
        let x = cond01_builder.add_dataflow_op(TketOp::X, [c1q]).unwrap();
        cond01_builder.finish_with_outputs([x.out_wire(0)]).ok();
        let cond0 = cond0_builder.finish_sub_container().unwrap();
        let mut cond1_builder = builder
            .conditional_builder(
                ([type_row![], type_row![]], meas0_convert.out_wire(0)),
                [(qb_t(), cond0.out_wire(0))],
                vec![qb_t()].into(),
            )
            .unwrap();
        let cond10_builder = cond1_builder.case_builder(0).unwrap();
        let [c0q] = cond10_builder.input_wires_arr();
        cond10_builder.finish_with_outputs([c0q]).ok();
        let mut cond11_builder = cond1_builder.case_builder(1).unwrap();
        let [c1q] = cond11_builder.input_wires_arr();
        let z = cond11_builder.add_dataflow_op(TketOp::Z, [c1q]).unwrap();
        cond11_builder.finish_with_outputs([z.out_wire(0)]).ok();
        let cond1 = cond1_builder.finish_sub_container().unwrap();
        let tdg = builder
            .add_dataflow_op(TketOp::Tdg, [cond1.out_wire(0)])
            .unwrap();
        let mut hugr = builder.finish_hugr_with_outputs([tdg.out_wire(0)]).unwrap();
        let func_node = hugr.first_child(hugr.module_root()).unwrap();
        let mut custom_settings = |n: hugr::Node| {
            if n == meas0.node() || n == meas1.node() {
                DataflowSettings {
                    flow_level: DataflowFlowDetail::All,
                    include_external_interface: true,
                    include_sum_types: true,
                    include_loop_body: false,
                    rotation_override: false,
                }
            } else {
                default_settings(n)
            }
        };
        let summary_map =
            SDFAnalysis::summarise_targets(&hugr, &vec![func_node], &mut custom_settings);
        let mut summary = summary_map.get(&func_node).unwrap().clone();
        assert_eq!(hugr.num_nodes(), 32);
        let mut pf = PhaseFold::new();
        pf.find_folds(&hugr, &mut summary);
        pf.apply_folds(&mut hugr, &MAX_PFSETTINGS);
        assert!(hugr.validate().is_ok());
        assert_eq!(hugr.num_nodes(), 30);
    }
}
