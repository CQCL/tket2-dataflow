"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-block.qasm

t a;
h a;
do {
    t a;
} while (x)
h a;
tdg a;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.tys import Qubit, Bool
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import H, OneQbGate
from tket_exts import quantum

T, Tdg = OneQbGate("T"), OneQbGate("Tdg")

circ = Dfg(Qubit, Bool)
a, x = circ.inputs()
a = circ.add_op(T, a)
a = circ.add_op(H, a)

with circ.add_tail_loop([], [a]) as loop:
    a, = loop.inputs()
    a = loop.add_op(T, a)
    loop.set_loop_outputs(x, a)

a, = loop.outputs()
a = circ.add_op(H, a)
a = circ.add_op(Tdg, a)
circ.set_outputs(a)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
