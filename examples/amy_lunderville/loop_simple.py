"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-simple.qasm

t a;
do {
    cx a, b;
} while (x)
tdg a;
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
a = circ.add_op(T, a)

with circ.add_tail_loop([], [a, b]) as loop:
    a, b = loop.inputs()
    a, b = loop.add_op(CX, a, b)
    loop.set_loop_outputs(x, a, b)

a, b = loop.outputs()
a = circ.add_op(Tdg, a)
circ.set_outputs(a, b)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
