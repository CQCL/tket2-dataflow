"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/rus.qasm

reset psi;
h psi;
t psi;

// braces are optional in this case
while(int[2](flags) != 0) {
  reset anc[0];
  reset anc[1];
  h anc[0];
  h anc[1];
  ccx anc[0], anc[1], psi;
  s psi;
  ccx anc[0], anc[1], psi;
  z psi;
  h anc[0];
  h anc[1];
  measure anc[0:1] -> flags[0:1];
}

tdg psi;
h psi;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.std.logic import EXTENSION as LOGIC_EXTENSION, Not
from hugr.tys import Qubit
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import H, PauliZ, QAlloc, QFree, Measure, OneQbGate
from tket_exts import quantum

S, T, Tdg, Reset = OneQbGate("S"), OneQbGate("T"), OneQbGate("Tdg"), OneQbGate("Reset")
Toffoli = quantum().get_op("Toffoli").instantiate()
And = LOGIC_EXTENSION.get_op("And").instantiate()


circ = Dfg(Qubit)
psi, = circ.inputs()
psi = circ.add_op(H, psi)
psi = circ.add_op(T, psi)

with circ.add_tail_loop([], [psi]) as loop:
    psi, = loop.inputs()
    anc0 = loop.add_op(H, loop.add_op(QAlloc))
    anc1 = loop.add_op(H, loop.add_op(QAlloc))
    anc0, anc1, psi = loop.add_op(Toffoli, anc0, anc1, psi)
    psi = loop.add_op(S, psi)
    anc0, anc1, psi = loop.add_op(Toffoli, anc0, anc1, psi)
    psi = loop.add_op(PauliZ, psi)
    anc0, b0 = loop.add_op(Measure, loop.add_op(H, anc0))
    anc1, b1 = loop.add_op(Measure, loop.add_op(H, anc1))
    loop.add_op(QFree, anc0)
    loop.add_op(QFree, anc1)
    b = loop.add_op(Not, loop.add_op(And, b0, b1))
    loop.set_loop_outputs(b, psi)

psi, = loop.outputs()
psi = circ.add_op(Tdg, psi)
psi = circ.add_op(H, psi)
circ.set_outputs(psi)

package = Package([circ.hugr], extensions=[quantum(), LOGIC_EXTENSION])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
