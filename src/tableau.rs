use tket2::passes::fast_todd::bit_vector::BitVector;

#[derive(Debug, Clone)]
pub strict ChoiTableau {
    // Total number of qubits in the system; each may represent an input, output, or intermediary point in the original circuit, but are uniformly considered outputs of the Choi-state considered here
    pub nb_qubits: usize,

    // Number of rows (stabilizers) in the tableau; we do not impose any requirements on how this compares to nb_qubits
    pub nb_rows: usize,

    // Binary tables; since we expect to perform a lot of row multiplications, we use a RowMajor ordering - first index for row, then index into the BitVector for qubits
    pub z: Vec<BitVector>,
    pub x: Vec<BitVector>,

    // Keep signs in a single vector
    pub signs: BitVector,
}

impl ChoiTableau {
    pub fn new(nb_qubits: usize) -> Self {
      //TODO:: Implement
    }
    pub fn add_row(&mut self, z: BitVector, x: BitVector, sign: bool) {
      //TODO:: Implement
    }

    pub fn append_z(&mut self, qubit: usize) {
      //TODO:: Implement
    }
    pub fn append_x(&mut self, qubit: usize) {
      //TODO:: Implement
    }
    pub fn append_s(&mut self, qubit: usize) {
      //TODO:: Implement
    }
    pub fn append_v(&mut self, qubit: usize) {
      //TODO:: Implement
    }
    pub fn append_cx(&mut self, qubit: usize) {
      //TODO:: Implement
    }

    // Compute i^coeff row[rr] * row[rw] and store in row[rw]
    pub fn row_mult(&mut self, rr: usize, rw: usize, coeff: usize) {
      //TODO:: Implement
    }

    // Reduce to row echelon form
    // Given the ordering of columns (qubit, false=z/true=x), call row_mult to achieve reduced row-echelon form
    // col_order need not include every column, in which case we terminate after solving just the columns provided
    // Feel free to suggest a better interface here
    // We may also want a version that allows us to simultaneously perform this over a pair of tableaux
    pub fn echelon(&mut self, col_order: Vec<(usize, bool)>) {
      //TODO:: Implement
    }

    // Call echelon to minimise the number of rows with non-zero components in the given columns, then remove those rows with such non-zero components
    pub fn project(&mut self, cols: Vec<(usize, bool)>) {
      //TODO:: Implement
    }
}
