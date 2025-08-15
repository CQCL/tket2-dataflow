"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-null.qasm

t b;
reset a;
do {
    t b;
    t a;
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

T, Tdg, Reset = OneQbGate("T"), OneQbGate("Tdg"), OneQbGate("Reset")

circ = Dfg(Qubit, Qubit, Bool)
a, b, x = circ.inputs()
b = circ.add_op(T, b)
a = circ.add_op(Reset, a)

with circ.add_tail_loop([], [a, b]) as loop:
    a, b = loop.inputs()
    a = loop.add_op(T, a)
    b = loop.add_op(T, b)
    loop.set_loop_outputs(x, a, b)

a, b = loop.outputs()
b = circ.add_op(Tdg, b)
circ.set_outputs(a, b)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
