use crate::bit_vector::BitVector;
use std::fmt;
use std::{cmp::min, fmt::Display};
use std::iter::zip;
use itertools::{interleave, Itertools};

#[derive(Debug, Clone)]
pub struct SymplecticTableau {
    // Total number of qubits in the system; each may represent an input, output, or intermediary point in the original circuit, but are uniformly considered outputs of the Choi-state considered here
    pub nb_qubits: usize,

    // Number of stabilizers in the tableau; we do not impose any requirements on how this compares to nb_qubits
    pub nb_stabs: usize,

    // Binary tables; since we expect to perform a lot of row multiplications, we use a StringMajor ordering - first index for string, then index into the BitVector for qubits
    pub z: Vec<BitVector>,
    pub x: Vec<BitVector>,

    // Keep signs in a single vector
    pub signs: BitVector,

}

#[derive(Debug, Clone)]
pub enum PauliXZ {
    X,
    Z,
}

impl SymplecticTableau {
    pub fn new(nb_qubits: usize) -> Self {
        SymplecticTableau {
            nb_qubits: nb_qubits,
            nb_stabs: 0,
            z: vec![],
            x: vec![],
            signs: BitVector::new(0),
        }
    }

    pub fn add_stab(&mut self, z: BitVector, x: BitVector, sign: bool) -> usize {
        let stab_id = self.nb_stabs;
        self.z.push(z);
        self.x.push(x);
        self.signs.extend_vec(vec![sign], self.nb_stabs);
        self.nb_stabs += 1;
        stab_id
    }

    pub fn add_qubits(&mut self, nb_new_qbs: usize) -> usize {
        let qb_base = self.nb_qubits;
        self.nb_qubits += nb_new_qbs;
        for zv in &mut self.z {
            zv.extend_vec(vec![false; nb_new_qbs], nb_new_qbs);
        }
        for xv in &mut self.x {
            xv.extend_vec(vec![false; nb_new_qbs], nb_new_qbs);
        }
        qb_base
    }

    pub fn append_z(&mut self, qubit: usize) {
        for (i, xv) in self.x.iter().enumerate() {
            if xv.get(qubit) {
                self.signs.xor_bit(i);
            }
        }
    }

    pub fn append_x(&mut self, qubit: usize) {
        for (i, zv) in self.z.iter().enumerate() {
            if zv.get(qubit) {
                self.signs.xor_bit(i);
            }
        }
    }

    pub fn append_s(&mut self, qubit: usize) {
        for (i, xv) in self.x.iter().enumerate() {
            if xv.get(qubit) {
                let zv : &mut BitVector = self.z.get_mut(i).unwrap();
                if zv.get(qubit) {
                    self.signs.xor_bit(i);
                }
                zv.xor_bit(qubit);
            }
        }
    }

    pub fn append_v(&mut self, qubit: usize) {
        for (i, zv) in self.z.iter().enumerate() {
            if zv.get(qubit) {
                let xv : &mut BitVector = self.x.get_mut(i).unwrap();
                if !xv.get(qubit) {
                    self.signs.xor_bit(i);
                }
                xv.xor_bit(qubit);
            }
        }
    }

    pub fn append_h(&mut self, qubit: usize) {
        for (i, (zv, xv)) in zip(self.z.iter_mut(), self.x.iter_mut()).enumerate() {
            let z = zv.get(qubit);
            let x = xv.get(qubit);
            // Swap values
            if z ^ x {
                zv.xor_bit(qubit);
                xv.xor_bit(qubit);
            }
            // Phase flip if Y
            if z && x {
                self.signs.xor_bit(i);
            }
        }
    }

    pub fn append_cx(&mut self, ctrl: usize, trgt: usize) {
        for (i, (zv, xv)) in zip(self.z.iter_mut(), self.x.iter_mut()).enumerate() {
            let zc = zv.get(ctrl);
            let zt = zv.get(trgt);
            let xc = xv.get(ctrl);
            let xt = xv.get(trgt);
            // Phase flip if XZ or YY
            if xc && zt && !(zc ^ xt) {
                self.signs.xor_bit(i);
            }
            // If xc, flip xt
            if xc {
                xv.xor_bit(trgt);
            }
            // If zt, flip zc
            if zt {
                zv.xor_bit(ctrl);
            }
        }
    }
    
    pub fn append_cy(&mut self, ctrl: usize, trgt: usize) {
        for (i, (zv, xv)) in zip(self.z.iter_mut(), self.x.iter_mut()).enumerate() {
            let zc = zv.get(ctrl);
            let zt = zv.get(trgt);
            let xc = xv.get(ctrl);
            let xt = xv.get(trgt);
            // Phase flip if XZ or YX
            if xc && ((!zc && !xt && zt) || (zc && xt && !zt)) {
                self.signs.xor_bit(i);
            }
            // If xc, flip zt and xt
            if xc {
                xv.xor_bit(trgt);
                zv.xor_bit(trgt);
            }
            // If zt ^ xt, flip zc
            if zt ^ xt {
                zv.xor_bit(ctrl);
            }
        }
    }

    pub fn append_cz(&mut self, ctrl: usize, trgt: usize) {
        for (i, (zv, xv)) in zip(self.z.iter_mut(), self.x.iter_mut()).enumerate() {
            let zc = zv.get(ctrl);
            let zt = zv.get(trgt);
            let xc = xv.get(ctrl);
            let xt = xv.get(trgt);
            // Phase flip if XX or YY
            if xc && xt && !(zc ^ zt) {
                self.signs.xor_bit(i);
            }
            // If xc, flip zt
            if xc {
                zv.xor_bit(trgt);
            }
            // If xt, flip zc
            if xt {
                zv.xor_bit(ctrl);
            }
        }
    }

    pub fn append_swap(&mut self, q0: usize, q1: usize) {
        for (zv, xv) in zip(self.z.iter_mut(), self.x.iter_mut()) {
            // To swap, just negate both bits if they differ
            if zv.get(q0) ^ zv.get(q1) {
                zv.xor_bit(q0);
                zv.xor_bit(q1);
            }
            if xv.get(q0) ^ xv.get(q1) {
                xv.xor_bit(q0);
                xv.xor_bit(q1);
            }
        }
    }

    // Compute i^coeff stabs[sr] * stabs[sw] and store in stabs[sw]
    pub fn stab_mult(&mut self, sr: usize, sw: usize, coeff: usize) {
        let zr = self.z.get(sr).unwrap().clone();
        let xr = self.x.get(sr).unwrap().clone();
        let zw = self.z.get_mut(sw).unwrap();
        let xw = self.x.get_mut(sw).unwrap();
        let mut num_is = coeff;
        // i for each Y=iXZ in each input string
        let mut yr = xr.clone();
        yr.and(&zr);
        let mut yw = xw.clone();
        yw.and(zw);
        // num_is += yr.get_all_ones(self.nb_qubits).len() + yw.get_all_ones(self.nb_qubits).len();
        num_is += ((yr.popcount() + yw.popcount()) % 4) as usize;
        // -1 for reordering (xr, zr, xw, zw) to (xr, xw, zr, zw)
        let mut zr_and_xw = zr.clone();
        zr_and_xw.and(xw);
        // num_is += 2 * zr_and_xw.get_all_ones(self.nb_qubits).len();
        num_is += 2 * ((zr_and_xw.popcount() % 4) as usize);
        // Calculate output string
        zw.xor(&zr);
        xw.xor(&xr);
        // -i for each Y=iXZ in output string
        let mut yw = xw.clone();
        yw.and(zw);
        // num_is += 3 * yw.get_all_ones(self.nb_qubits).len();
        num_is += 3 * ((yw.popcount() % 4) as usize);
        // Adjust sign
        if self.signs.get(sr) ^ (num_is % 4 == 2) {
            self.signs.xor_bit(sw);
        }
    }

    // Removes a stabilizer from the tableau
    // If the stabilizer to be removed was the last one, returns None
    // Otherwise, in order to keep things dense, swap with the last stabilizer before removing; returns the index of the last stabilizer, i.e. the old index that is now at to_delete
    pub fn delete_stab(&mut self, to_delete: usize) -> Option<usize> {
        self.nb_stabs -= 1;
        if to_delete == self.nb_stabs {
            self.z.pop();
            self.x.pop();
            // BitVector has no record of the number of bits it contains, so we just reset the bit in case it gets reused later
            if self.signs.get(self.nb_stabs) {
                self.signs.xor_bit(self.nb_stabs);
            }
            None
        }
        else {
            // Overwrite with last stabilizer
            self.z[to_delete] = self.z.pop().unwrap();
            self.x[to_delete] = self.x.pop().unwrap();
            if self.signs.get(to_delete) != self.signs.get(self.nb_stabs) {
                self.signs.xor_bit(to_delete);
            }
            // BitVector has no record of the number of bits it contains, so we just reset the bit in case it gets reused later
            if self.signs.get(self.nb_stabs) {
                self.signs.xor_bit(self.nb_stabs);
            }
            Some(self.nb_stabs)
        }
    }
    
    // Removes a qubit from the tableau
    // If the qubit to be removed was the last one, returns None
    // Otherwise, in order to keep things dense, swap with the last qubit before removing; returns the index of the last qubit, i.e. the old index that is now at to_delete
    pub fn delete_qubit(&mut self, to_delete: usize) -> Option<usize> {
        self.nb_qubits -= 1;
        if to_delete == self.nb_qubits {
            // If any stabilizers involve the qubit, we cannot delete
            for zv in self.z.iter() {
                assert_eq!(zv.get(self.nb_qubits), false);
            }
            for xv in self.x.iter() {
                assert_eq!(xv.get(self.nb_qubits), false);
            }
            // BitVector has no record of the number of bits it contains, so we can just leave the bits set as 0
            None
        }
        else {
            // If any stabilizers involve the qubit, we cannot delete
            // Move any value from the last qubit to to_delete and reset the last element in the BitVector so it may be safely reused later
            for zv in self.z.iter_mut() {
                assert_eq!(zv.get(to_delete), false);
                if zv.get(self.nb_qubits) {
                    zv.xor_bit(to_delete);
                    zv.xor_bit(self.nb_qubits);
                }
            }
            for xv in self.x.iter_mut() {
                assert_eq!(xv.get(to_delete), false);
                if xv.get(self.nb_qubits) {
                    xv.xor_bit(to_delete);
                    xv.xor_bit(self.nb_qubits);
                }
            }
            Some(self.nb_qubits)
        }
    }

    // Reduce to row echelon form over the stabilizers
    // Given the ordering of columns (qubit, z/x), call stab_mult to achieve reduced row-echelon form
    // col_order need not include every column, in which case we terminate after solving just the columns provided
    // Feel free to suggest a better interface here
    // We may also want a version that allows us to simultaneously perform this over a pair of tableaux
    pub fn echelon(&mut self, col_order: &Vec<(usize, PauliXZ)>) {
        let mut pivot_stab = 0;
        for (qubit, p) in col_order {
            match p {
                PauliXZ::X => {
                    let mut pivot_found = false;
                    // Find next row
                    for i in pivot_stab..self.nb_stabs {
                        if self.x[i].get(*qubit) {
                            // Found new stabilizer to pivot
                            pivot_found = true;
                            if i != pivot_stab {
                                // Make sure pivot_stab contains this element
                                self.stab_mult(i, pivot_stab, 0);
                            }
                            break;
                        }
                    }
                    // Eliminate entries from all other rows
                    if pivot_found {
                        for i in 0..self.nb_stabs {
                            if (i != pivot_stab) && self.x[i].get(*qubit) {
                                self.stab_mult(pivot_stab, i, 0);
                            }
                        }
                        pivot_stab += 1;
                    }
                }
                PauliXZ::Z => {
                    let mut pivot_found = false;
                    // Find next row
                    for i in pivot_stab..self.nb_stabs {
                        if self.z[i].get(*qubit) {
                            // Found new stabilizer to pivot
                            pivot_found = true;
                            if i != pivot_stab {
                                // Make sure pivot_stab contains this element
                                self.stab_mult(i, pivot_stab, 0);
                            }
                            break;
                        }
                    }
                    // Eliminate entries from all other rows
                    if pivot_found {
                        for i in 0..self.nb_stabs {
                            if (i != pivot_stab) && self.z[i].get(*qubit) {
                                self.stab_mult(pivot_stab, i, 0);
                            }
                        }
                        pivot_stab += 1;
                    }
                }
            }
        }
    }

    pub fn all_columns(& self) -> Vec<(usize, PauliXZ)> {
        interleave(
            (0..self.nb_stabs).map(|c| (c, PauliXZ::X)),
            (0..self.nb_stabs).map(|c| (c, PauliXZ::Z))
        ).collect_vec()
    }

    // Call echelon to minimise the number of rows with non-zero components in the given columns, then remove those rows with such non-zero components
    pub fn project(&mut self, cols: &Vec<(usize, PauliXZ)>) {
        self.echelon(cols);
        let mut num_to_remove = 0;
        for (qubit, p) in cols {
            match p {
                PauliXZ::X => {
                    if self.x[num_to_remove].get(*qubit) {
                        num_to_remove += 1;
                    }
                }
                PauliXZ::Z => {
                    if self.z[num_to_remove].get(*qubit) {
                        num_to_remove += 1;
                    }
                }
            }
        }
        // Even though stabs to remove will be at the top, if num_to_remove is greater than nb_stabs/2 we can choose an order which maximises the number of stab removals from the end which are cheaper as they don't involve moving rows around
        let num_remaining = self.nb_stabs - num_to_remove;
        for i in 0..min(num_to_remove, num_remaining) {
            self.delete_stab(i);
        }
        for i in (num_remaining..num_to_remove).rev() {
            self.delete_stab(i);
        }
    }

    pub fn anticommutes_with(&mut self, stab: usize, z: &BitVector, x: &BitVector) -> bool {
        let mut sz_and_x = self.z.get(stab).unwrap().clone();
        sz_and_x.and(x);
        let mut sx_and_z = self.x.get(stab).unwrap().clone();
        sx_and_z.and(z);
        sz_and_x.xor(&sx_and_z);
        (sz_and_x.popcount() % 2) == 1
    }

    // Apply row combinations to leave at most one row anticommuting with the target Pauli string, and remove it
    pub fn project_commuting_with(&mut self, z: &BitVector, x: &BitVector) {
        let mut anticommuting_stab : Option<usize> = None;
        for i in 0..self.nb_stabs {
            if self.anticommutes_with(i, z, x) {
                match anticommuting_stab {
                    Some(ac_stab) => {
                        self.stab_mult(ac_stab, i, 0);
                    }
                    None => {
                        anticommuting_stab = Some(i);
                    }
                }
            }
        }
        if anticommuting_stab.is_some() {
            self.delete_stab(anticommuting_stab.unwrap());
        }
    }

    pub fn join(left: &SymplecticTableau, right: &SymplecticTableau) -> SymplecticTableau {
        // Initialise two tableaus forming the following stacks:
        // / L | L \
        // \ R | 0 /
        let mut lr_tab = left.clone();
        let mut l0_tab = left.clone();
        for i in 0..right.nb_stabs {
            lr_tab.add_stab(right.z[i].clone(), right.x[i].clone(), right.signs.get(i));
            l0_tab.add_stab(BitVector::new(left.nb_qubits), BitVector::new(left.nb_qubits), false);
        }
        // Perform simultaneous row combinations to convert lr_tab into row echelon form (no need for reduced row echelon)
        let mut pivot_stab = 0;
        for q in 0..lr_tab.nb_qubits {
            // Reduce X column
            let mut pivot_found = false;
            for i in pivot_stab..lr_tab.nb_stabs {
                if lr_tab.x[i].get(q) {
                    // Found new stabilizer to pivot
                    pivot_found = true;
                    if i != pivot_stab {
                        // Make sure pivot_stab contains this element
                        lr_tab.stab_mult(i, pivot_stab, 0);
                        l0_tab.stab_mult(i, pivot_stab, 0);
                    }
                    break;
                }
            }
            // Eliminate entries from lower rows
            if pivot_found {
                for i in (pivot_stab+1)..lr_tab.nb_stabs {
                    if lr_tab.x[i].get(q) {
                        lr_tab.stab_mult(pivot_stab, i, 0);
                        l0_tab.stab_mult(pivot_stab, i, 0);
                    }
                }
                pivot_stab += 1;
            }
            // Reduce Z column
            pivot_found = false;
            for i in pivot_stab..lr_tab.nb_stabs {
                if lr_tab.z[i].get(q) {
                    // Found new stabilizer to pivot
                    pivot_found = true;
                    if i != pivot_stab {
                        // Make sure pivot_stab contains this element
                        lr_tab.stab_mult(i, pivot_stab, 0);
                        l0_tab.stab_mult(i, pivot_stab, 0);
                    }
                    break;
                }
            }
            // Eliminate entries from lower rows
            if pivot_found {
                for i in (pivot_stab+1)..lr_tab.nb_stabs {
                    if lr_tab.z[i].get(q) {
                        lr_tab.stab_mult(pivot_stab, i, 0);
                        l0_tab.stab_mult(pivot_stab, i, 0);
                    }
                }
                pivot_stab += 1;
            }
        }
        // Empty rows of stacked will be at the bottom; the final value of pivot_stab is the first zero row
        // Rather than remove rows from l0_tab, since <half will remain we just copy the good ones into a new tab
        let mut res_tab = SymplecticTableau::new(l0_tab.nb_qubits);
        for i in pivot_stab..l0_tab.nb_stabs {
            res_tab.add_stab(l0_tab.z[i].clone(), l0_tab.x[i].clone(), l0_tab.signs.get(i));
        }
        for i in 0..res_tab.nb_stabs { assert_eq!(res_tab.z[i].blocks.len(), 1); }
        res_tab
    }

}

impl Display for SymplecticTableau {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..self.nb_stabs {
            let mut pauli_str = String::new();
            if self.signs.get(i) {
                pauli_str.push('-');
            }
            else {
                pauli_str.push('+');
            }
            for q in 0..self.nb_qubits {
                match (self.x[i].get(q), self.z[i].get(q)) {
                    (false, false) => {
                        pauli_str.push(' ');
                    }
                    (false, true) => {
                        pauli_str.push('Z');
                    }
                    (true, false) => {
                        pauli_str.push('X');
                    }
                    (true, true) => {
                        pauli_str.push('Y');
                    }
                }
            }
            write!(f, "{}\n", pauli_str)?;
        }
        Ok(())
    }
}