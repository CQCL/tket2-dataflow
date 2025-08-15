"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-cycle.qasm

t b;
do {
    swap a, b;
    swap b, c;
} while (x)
tdg b;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.tys import Qubit, Bool
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import OneQbGate
from tket_exts import quantum

T, Tdg = OneQbGate("T"), OneQbGate("Tdg")

circ = Dfg(Qubit, Qubit, Qubit, Bool)
a, b, c, x = circ.inputs()
b = circ.add_op(T, b)

with circ.add_tail_loop([], [a, b, c]) as loop:
    a, b, c = loop.inputs()
    loop.set_loop_outputs(x, b, c, a)

a, b, c = loop.outputs()
b = circ.add_op(Tdg, b)
circ.set_outputs(a, b, c)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
