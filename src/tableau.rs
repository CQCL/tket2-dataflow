use crate::bit_vector::BitVector;
use crate::pauli_product::PauliProduct;
use tket::TketOp;

type Command = (TketOp, Vec<usize>);

#[derive(Debug, Clone)]
pub struct Tableau {
    pub nb_qubits: usize,
    pub z: Vec<BitVector>,
    pub x: Vec<BitVector>,
    pub signs: BitVector,
}

impl Tableau {
    pub fn new(nb_qubits: usize) -> Self {
        Tableau {
            nb_qubits,
            z: Tableau::init_z(nb_qubits),
            x: Tableau::init_x(nb_qubits),
            signs: BitVector::new(nb_qubits << 1),
        }
    }

     fn init_z(nb_qubits: usize) -> Vec<BitVector> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits << 1);
            bv.xor_bit(i);
            vec.push(bv);
        }
        vec
    }

     fn init_x(nb_qubits: usize) -> Vec<BitVector> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits << 1);
            bv.xor_bit(i + nb_qubits);
            vec.push(bv);
        }
        vec
    }

    pub fn append_x(&mut self, qubit: usize) {
        self.signs.xor(&self.z[qubit]);
    }

    pub fn append_z(&mut self, qubit: usize) {
        self.signs.xor(&self.x[qubit]);
    }

    pub fn append_v(&mut self, qubit: usize) {
        let mut a = self.x[qubit].clone();
        a.negate();
        a.and(&self.z[qubit]);
        self.signs.xor(&a);
        self.x[qubit].xor(&self.z[qubit]);
    }

    pub fn append_s(&mut self, qubit: usize) {
        let mut a = self.z[qubit].clone();
        a.and(&self.x[qubit]);
        self.signs.xor(&a);
        self.z[qubit].xor(&self.x[qubit]);
    }

    pub fn append_h(&mut self, qubit: usize) {
        self.append_s(qubit);
        self.append_v(qubit);
        self.append_s(qubit);
    }

    pub fn append_cx(&mut self, qubits: Vec<usize>) {
        let mut a =  self.z[qubits[0]].clone();
        a.negate();
        a.xor(&self.x[qubits[1]]);
        a.and(&self.z[qubits[1]]);
        a.and(&self.x[qubits[0]]);
        self.signs.xor(&a);
        let a = self.z[qubits[1]].clone();
        self.z[qubits[0]].xor(&a);
        let a = self.x[qubits[0]].clone();
        self.x[qubits[1]].xor(&a);
    }

    pub fn append_cz(&mut self, qubits: Vec<usize>) {
        self.append_s(qubits[0]);
        self.append_s(qubits[1]);
        self.append_cx(qubits.to_vec());
        self.append_s(qubits[1]);
        self.append_z(qubits[1]);
        self.append_cx(qubits);
    }

    pub fn to_circ(&self, inverse: bool) -> Vec<Command> {
        let mut tab = self.clone();
        let mut c = Vec::new();
        for i in 0..self.nb_qubits {
            if let Some(index) = tab.x.iter().position(|x| x.get(i)) {
                for j in (i+1)..self.nb_qubits {
                    if tab.x[j].get(i) && j != index {
                        tab.append_cx(vec![index, j]);
                        c.push((TketOp::CX, vec![index, j]));
                    }
                }
                if tab.z[index].get(i) {
                    tab.append_s(index);
                    c.push((TketOp::S, vec![index]));
                }
                tab.append_h(index);
                c.push((TketOp::H, vec![index]));
            }
            if !tab.z[i].get(i) {
                let index = tab.z.iter().position(|z| z.get(i)).unwrap();
                tab.append_cx(vec![i, index]);
                c.push((TketOp::CX, vec![i, index]));
            }
            for j in 0..self.nb_qubits {
                if tab.z[j].get(i) && j != i {
                    tab.append_cx(vec![j, i]);
                    c.push((TketOp::CX, vec![j, i]));
                }
            }
            for j in 0..self.nb_qubits {
                if tab.x[j].get(i + self.nb_qubits) && j != i {
                    tab.append_cx(vec![i, j]);
                    c.push((TketOp::CX, vec![i, j]));
                }
            }
            for j in 0..self.nb_qubits {
                if tab.z[j].get(i + self.nb_qubits) && j != i {
                    tab.append_cx(vec![i, j]);
                    c.push((TketOp::CX, vec![i, j]));
                    tab.append_s(j);
                    c.push((TketOp::S, vec![j]));
                    tab.append_cx(vec![i, j]);
                    c.push((TketOp::CX, vec![i, j]));
                }
            }
            if tab.z[i].get(i + self.nb_qubits) {
                tab.append_s(i);
                c.push((TketOp::S, vec![i]));
            }
            if tab.signs.get(i) {
                tab.append_x(i);
                c.push((TketOp::X, vec![i]));
            }
            if tab.signs.get(i + self.nb_qubits) {
                tab.append_z(i);
                c.push((TketOp::Z, vec![i]));
            }
        }
        if !inverse {
            let mut c2 = Vec::new();
            for (gate, qubits) in c.into_iter().rev() {
                c2.push((gate, qubits.to_vec()));
                if gate == TketOp::S { c2.push((TketOp::Z, qubits.to_vec())); }
            }
            return c2;
        }
        c
    }
}
#[derive(Debug, Clone)]
pub struct TableauColumnMajor {
    pub nb_qubits: usize,
    pub stabs: Vec<PauliProduct>,
    pub destabs: Vec<PauliProduct>,
}

impl TableauColumnMajor {
    pub fn new(nb_qubits: usize) -> Self {
        TableauColumnMajor {
            nb_qubits,
            stabs: TableauColumnMajor::init_stabs(nb_qubits),
            destabs: TableauColumnMajor::init_destabs(nb_qubits),
        }
    }

     fn init_stabs(nb_qubits: usize) -> Vec<PauliProduct> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits);
            bv.xor_bit(i);
            vec.push(PauliProduct::new(bv, BitVector::new(nb_qubits), false));
        }
        vec
    }

     fn init_destabs(nb_qubits: usize) -> Vec<PauliProduct> {
        let mut vec = Vec::new();
        for i in 0..nb_qubits {
            let mut bv = BitVector::new(nb_qubits);
            bv.xor_bit(i);
            vec.push(PauliProduct::new(BitVector::new(nb_qubits), bv, false));
        }
        vec
    }

    pub fn prepend_x(&mut self, qubit: usize) {
        self.stabs[qubit].sign ^= true;
    }

    pub fn prepend_z(&mut self, qubit: usize) {
        self.destabs[qubit].sign ^= true;
    }

    pub fn prepend_v(&mut self, qubit: usize) {
        self.stabs[qubit].pauli_product_mult(&self.destabs[qubit]);
    }

    pub fn prepend_s(&mut self, qubit: usize) {
        self.destabs[qubit].pauli_product_mult(&self.stabs[qubit]);
    }

    pub fn prepend_h(&mut self, qubit: usize) {
        self.prepend_s(qubit);
        self.prepend_v(qubit);
        self.prepend_s(qubit);
    }

    pub fn prepend_cx(&mut self, qubits: Vec<usize>) {
        let p = self.stabs[qubits[0]].clone();
        self.stabs[qubits[1]].pauli_product_mult(&p);
        let p = self.destabs[qubits[1]].clone();
        self.destabs[qubits[0]].pauli_product_mult(&p);
    }

    pub fn to_circ(&self, inverse: bool) -> Vec<Command> {
        let mut tab = self.clone();
        // let mut c = RestrictedSubcircuit::new(tab.nb_qubits, HashSet::new());
        let mut c = Vec::new();
        for i in 0..tab.nb_qubits {
            if let Some(index) = tab.stabs.iter().position(|p| p.x.get(i) ) {
                for j in (i+1)..tab.nb_qubits {
                    if tab.stabs[j].x.get(i) && j != index {
                        tab.prepend_cx(vec![index, j]);
                        // c.gates.push((TketOp::CX, vec![index, j]));
                        c.push((TketOp::CX, vec![index, j]));
                    }
                }
                if tab.destabs[index].x.get(i) {
                    tab.prepend_s(index);
                    // c.gates.push((TketOp::S, vec![index]));
                    c.push((TketOp::S, vec![index]));
                }
                tab.prepend_h(index);
                // c.gates.push((TketOp::H, vec![index]));
                c.push((TketOp::H, vec![index]));
            }
            if !tab.destabs[i].x.get(i) {
                let index = tab.destabs.iter().position(|p| p.x.get(i)).unwrap();
                tab.prepend_cx(vec![i, index]);
                // c.gates.push((TketOp::CX, vec![i, index]));
                c.push((TketOp::CX, vec![i, index]));
            }
            for j in 0..tab.nb_qubits {
                if tab.destabs[j].x.get(i) && j != i {
                    tab.prepend_cx(vec![j, i]);
                    // c.gates.push((TketOp::CX, vec![j, i]));
                    c.push((TketOp::CX, vec![j, i]));
                }
            }
            for j in 0..tab.nb_qubits {
                if tab.stabs[j].z.get(i) && j != i {
                    tab.prepend_cx(vec![i, j]);
                    // c.gates.push((TketOp::CX, vec![i, j]));
                    c.push((TketOp::CX, vec![i, j]));
                }
            }
            for j in 0..tab.nb_qubits {
                if tab.destabs[j].z.get(i) && j != i {
                    tab.prepend_cx(vec![i, j]);
                    // c.gates.push((TketOp::CX, vec![i, j]));
                    c.push((TketOp::CX, vec![i, j]));
                    tab.prepend_s(j);
                    // c.gates.push((TketOp::S, vec![j]));
                    c.push((TketOp::S, vec![j]));
                    tab.prepend_cx(vec![i, j]);
                    // c.gates.push((TketOp::CX, vec![i, j]));
                    c.push((TketOp::CX, vec![i, j]));
                }
            }
            if tab.destabs[i].z.get(i) {
                tab.prepend_s(i);
                // c.gates.push((TketOp::S, vec![i]));
                c.push((TketOp::S, vec![i]));
            }
            if tab.stabs[i].sign {
                tab.prepend_x(i);
                // c.gates.push((TketOp::X, vec![i]));
                c.push((TketOp::X, vec![i]));
            }
            if tab.destabs[i].sign {
                tab.prepend_z(i);
                // c.gates.push((TketOp::Z, vec![i]));
                c.push((TketOp::Z, vec![i]));
            }
        }
        // c.gates.reverse();
        c.reverse();
        if !inverse {
            // let mut c2 = RestrictedSubcircuit::new(tab.nb_qubits, HashSet::new());
            // for (gate, qubits) in c.gates.into_iter().rev() {
            //     c2.gates.push((gate, qubits.to_vec()));
            //     if gate == TketOp::S { c2.gates.push((TketOp::Z, qubits.to_vec())); }
            // }
            // return c2;
            let mut c2 = Vec::new();
            for (gate, qubits) in c.into_iter().rev() {
                c2.push((gate, qubits.to_vec()));
                if gate == TketOp::S { c2.push((TketOp::Z, qubits.to_vec())); }
            }
            return c2;       
        }
        c
    }
}
