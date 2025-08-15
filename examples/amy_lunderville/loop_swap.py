"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-swap.qasm

cx a, b;
t b;
cx a, b;
do {
    swap a, b;
} while (x)
cx a, b;
tdg b;
cx a, b;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.tys import Qubit, Bool
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import CX, OneQbGate
from tket_exts import quantum

T, Tdg = OneQbGate("T"), OneQbGate("Tdg")

circ = Dfg(Qubit, Qubit, Bool)
a, b, x = circ.inputs()
a, b = circ.add_op(CX, a, b)
b = circ.add_op(T, b)
a, b = circ.add_op(CX, a, b)

with circ.add_tail_loop([], [a, b]) as loop:
    a, b = loop.inputs()
    loop.set_loop_outputs(x, b, a)

a, b = loop.outputs()
a, b = circ.add_op(CX, a, b)
b = circ.add_op(Tdg, b)
a, b = circ.add_op(CX, a, b)
circ.set_outputs(a, b)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
