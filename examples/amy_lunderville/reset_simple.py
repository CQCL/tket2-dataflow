"""https://github.com/meamy/feynman/blob/master/benchmarks/qasm3/reset-simple.qasm

t a;
reset a;
t a;
"""

from pathlib import Path

from hugr.build import Dfg
from hugr.package import Package
from hugr.tys import Qubit
from selene_hugr_qis_compiler import check_hugr
from tket.circuit.build import OneQbGate
from tket_exts import quantum

T, Reset = OneQbGate("T"), OneQbGate("Reset")

circ = Dfg(Qubit)
a, = circ.inputs()
a = circ.add_op(T, a)
a = circ.add_op(Reset, a)
a = circ.add_op(T, a)
circ.set_outputs(a)

package = Package([circ.hugr], extensions=[quantum()])
b = package.to_bytes()
check_hugr(b)
with open(Path(__file__).with_suffix(".hugr"), "wb") as f:
    f.write(b)
