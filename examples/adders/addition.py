"""'Gidney' Addition. This is lifted directly from guppy-algs."""

from guppylang import guppy
from guppylang.std.quantum import (
    qubit,
    x,
    h,
    cx,
    t,
    tdg,
    s,
    cz,
    project_z,
    discard_array,
    measure_array,
)
from guppylang.std.builtins import array, owned
from selene_hugr_qis_compiler import check_hugr


from hugr.package import Package
from pathlib import Path
from tket_exts import quantum


@guppy
def cxx(control: qubit, target1: qubit, target2: qubit) -> None:
    """Apply CXX."""
    cx(control, target1)
    cx(control, target2)


@guppy
def compute_logical_and(x: qubit, y: qubit, t_qubit: qubit) -> None:
    """Compute logical and operation."""
    cx(x, t_qubit)
    cx(y, t_qubit)
    cxx(t_qubit, y, x)
    tdg(x)
    tdg(y)
    t(t_qubit)
    cxx(t_qubit, y, x)
    h(t_qubit)
    s(t_qubit)


@guppy
def uncompute_logical_and(x: qubit, y: qubit, a: qubit) -> None:
    """Uncompute logical and operation."""
    h(a)
    if project_z(a):
        cz(x, y)


K = guppy.nat_var("K")
N = guppy.nat_var("N")  # Note: N = K - 1


@guppy
def add_quantum_integers(
    i_reg: array[qubit, K],
    t_reg: array[qubit, K],
    a_reg: array[qubit, N] @ owned,
) -> None:
    """Compute addition."""
    compute_logical_and(i_reg[0], t_reg[0], a_reg[0])

    for i in range(N - 1):
        cxx(a_reg[i], i_reg[i + 1], t_reg[i + 1])
        compute_logical_and(i_reg[i + 1], t_reg[i + 1], a_reg[i + 1])
        cx(a_reg[i], a_reg[i + 1])

    cx(a_reg[N - 1], t_reg[K - 1])

    for i in range(N - 1):
        cx(a_reg[N - 2 - i], a_reg[N - 1 - i])
        uncompute_logical_and(i_reg[K - 2 - i], t_reg[K - 2 - i], a_reg[N - 1 - i])
        cx(a_reg[N - 2 - i], i_reg[K - 2 - i])

    uncompute_logical_and(i_reg[0], t_reg[0], a_reg[0])

    discard_array(a_reg)

    for j in range(K):
        cx(i_reg[j], t_reg[j])




@guppy
def get_7_ts() -> array[qubit, 7]:
    """Get 7 ts."""
    arr = array(qubit() for _ in range(7))
    for i in range(7):
        h(arr[i])
        t(arr[i])

    return arr


@guppy
def encode_integer() -> array[qubit, 8]:
    """Encode integer."""
    arr = array(qubit() for _ in range(8))
    x(arr[0])
    x(arr[1])
    x(arr[2])
    x(arr[3])
    x(arr[4])
    x(arr[5])
    x(arr[6])
    x(arr[7])
    return arr


@guppy
def main() -> None:
    """Main function for addition example."""
    t_states = get_7_ts()

    i_reg = encode_integer()
    t_reg = encode_integer()

    add_quantum_integers(i_reg, t_reg, t_states)

    discard_array(i_reg)

    result("c", measure_array(t_reg))


compiled_hugr = main.compile()

b = compiled_hugr.to_bytes()

filename = "add_264.hugr"  
with open(Path(filename), "wb") as f:
    f.write(b)