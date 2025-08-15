"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-nonlinear.qasm

reset c;
x a;
ccx a, b, c;
t c;
ccx a, b, c;
x a;
do {
    cx a, b;
} while (x)
x a;
ccx a, b, c;
tdg c;
ccx a, b, c;
x a;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.tys import Qubit, Bool
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import H, PauliX, CX, OneQbGate
from tket_exts import quantum

T, Tdg, Reset = OneQbGate("T"), OneQbGate("Tdg"), OneQbGate("Reset")
Toffoli = quantum().get_op("Toffoli").instantiate()

circ = Dfg(Qubit, Qubit, Qubit, Bool)
a, b, c, x = circ.inputs()

c = circ.add_op(Reset, c)
a = circ.add_op(PauliX, a)
a, b, c = circ.add_op(Toffoli, a, b, c)
c = circ.add_op(T, c)
a, b, c = circ.add_op(Toffoli, a, b, c)
a = circ.add_op(PauliX, a)

with circ.add_tail_loop([], [a, b]) as loop:
    a, b = loop.inputs()
    a, b = loop.add_op(CX, a, b)
    loop.set_loop_outputs(x, a, b)

a, b = loop.outputs()
a = circ.add_op(PauliX, a)
a, b, c = circ.add_op(Toffoli, a, b, c)
c = circ.add_op(Tdg, c)
a, b, c = circ.add_op(Toffoli, a, b, c)
a = circ.add_op(PauliX, a)
circ.set_outputs(a, b, c)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
