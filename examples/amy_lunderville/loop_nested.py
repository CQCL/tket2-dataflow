"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/loop-nested.qasm

t b;
do {
    t a;
    do {
        x b;
    } while (y)
} while(x)
tdg b;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.tys import Qubit, Bool
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import H, PauliX, OneQbGate
from tket_exts import quantum

T, Tdg = OneQbGate("T"), OneQbGate("Tdg")

circ = Dfg(Qubit, Qubit, Bool, Bool)
a, b, x, y = circ.inputs()
b = circ.add_op(T, b)

with circ.add_tail_loop([], [a, b]) as outer_loop:
    a, b = outer_loop.inputs()
    a = outer_loop.add_op(T, a)
    with outer_loop.add_tail_loop([], [b]) as inner_loop:
        b, = inner_loop.inputs()
        b = inner_loop.add_op(PauliX, b)
        inner_loop.set_loop_outputs(y, b)
    b, = inner_loop.outputs()
    outer_loop.set_loop_outputs(x, a, b)

a, b = outer_loop.outputs()
b = circ.add_op(Tdg, b)
circ.set_outputs(a, b)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
