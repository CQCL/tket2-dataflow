#################### COPY cnx.py FROM GUPPYALGOS #################
"""Multicontrolled x gate."""

from guppylang.decorator import guppy
from guppylang.std.builtins import array
from guppylang.std.quantum import cx, discard, h, qubit, t, tdg, toffoli, x


from typing import no_type_check

n_controls = guppy.nat_var("n_controls")


@guppy
@no_type_check
def cnx(control: array[qubit, n_controls], target: qubit) -> None:
    r"""Apply efficient cnx in terms of number of 2q gates.

    Based on https://arxiv.org/pdf/1508.03273, given $n - 1$ control qubits and the
    target, uses $\lceil{(n-3)/2}\rceil$ ancillary qubits to apply the multicontrolled-x
    gate on the target using $6n -12$ CNOT gates.

    Args:
        control (array[qubit, n_controls]): register of control qubits.
        target (qubit): target qubit.

    """
    is_odd, num_ancillas = _get_variables_cnx(n_controls)
    if n_controls == 0:
        x(target)
    elif n_controls == 1:
        cx(control[0], target)
    elif n_controls == 2:
        toffoli(control[0], control[1], target)
    elif n_controls == 3:
        ancilla = qubit()
        _tof4(control[0], control[1], control[2], target, ancilla)
        discard(ancilla)
    elif n_controls == 4:
        ancilla = qubit()
        _tof5(control[0], control[1], control[2], control[3], target, ancilla)
        discard(ancilla)
    else:
        ancilla = qubit()
        _rtof4(control[0], control[1], control[2], ancilla)
        _cnx_aux(control, target, ancilla, num_ancillas, 3, 0, is_odd)
        _irtof4(control[0], control[1], control[2], ancilla)
        discard(ancilla)


@guppy
@no_type_check
def _rtof(a: qubit, b: qubit, c: qubit) -> None:
    """Apply the relative tof block."""
    h(c)
    t(c)
    cx(b, c)
    tdg(c)
    cx(a, c)
    t(c)
    cx(b, c)
    tdg(c)
    h(c)


@guppy
@no_type_check
def _irtof(a: qubit, b: qubit, c: qubit) -> None:
    """Apply the inverse of relative tof block. rtof is self-inverse."""
    _rtof(a, b, c)


@guppy
@no_type_check
def _rtof4(a: qubit, b: qubit, c: qubit, d: qubit) -> None:
    """Apply the relative tof4 block."""
    h(d)
    t(d)
    cx(c, d)
    tdg(d)
    h(d)
    cx(a, d)
    t(d)
    cx(b, d)
    tdg(d)
    cx(a, d)
    t(d)
    cx(b, d)
    tdg(d)
    h(d)
    t(d)
    cx(c, d)
    tdg(d)
    h(d)


@guppy
@no_type_check
def _irtof4(a: qubit, b: qubit, c: qubit, d: qubit) -> None:
    """Apply the inverse of relative tof4 block."""
    h(d)
    t(d)
    cx(c, d)
    tdg(d)
    h(d)
    t(d)
    cx(b, d)
    tdg(d)
    cx(a, d)
    t(d)
    cx(b, d)
    tdg(d)
    cx(a, d)
    h(d)
    t(d)
    cx(c, d)
    tdg(d)
    h(d)


@guppy
@no_type_check
def _tof4(
    control0: qubit,
    control1: qubit,
    control2: qubit,
    target: qubit,
    ancilla: qubit,
) -> None:
    """Apply tof4 block."""
    _rtof(control0, control1, ancilla)
    toffoli(ancilla, control2, target)
    _irtof(control0, control1, ancilla)


@guppy
@no_type_check
def _tof5(
    control0: qubit,
    control1: qubit,
    control2: qubit,
    control3: qubit,
    target: qubit,
    ancilla: qubit,
) -> None:
    """Apply tof5 block."""
    _rtof4(control0, control1, control2, ancilla)
    toffoli(ancilla, control3, target)
    _irtof4(control0, control1, control2, ancilla)


@guppy.comptime
@no_type_check
def _get_variables_cnx(n_controls: int) -> tuple[int, int]:
    """Auxiliary method for cnx that computes the necessary variables."""
    is_odd = n_controls % 2
    num_ancillas = (n_controls - 2) // 2 + (n_controls - 2) % 2
    return is_odd, num_ancillas


@guppy
@no_type_check
def _cnx_aux(
    control: array[qubit, n_controls],
    target: qubit,
    ancilla: qubit,
    n_ancillas: int,
    control_counter: int,
    ancilla_counter: int,
    is_odd: int,
) -> None:
    """Auxiliary method for cnx that recursively builds circuit."""
    new_ancilla = qubit()
    if ancilla_counter == n_ancillas - 2:
        if is_odd == 1:
            _tof4(
                ancilla,
                control[n_controls - 2],
                control[n_controls - 1],
                target,
                new_ancilla,
            )
        else:
            _tof5(
                ancilla,
                control[n_controls - 3],
                control[n_controls - 2],
                control[n_controls - 1],
                target,
                new_ancilla,
            )

    else:
        _rtof4(
            ancilla, control[control_counter], control[control_counter + 1], new_ancilla
        )
        _cnx_aux(
            control,
            target,
            new_ancilla,
            n_ancillas,
            control_counter + 2,
            ancilla_counter + 1,
            is_odd,
        )
        _irtof4(
            ancilla, control[control_counter], control[control_counter + 1], new_ancilla
        )
    discard(new_ancilla)

#################### END OF cnx.py #####################

from guppylang.std.angles import angle
from guppylang.std.builtins import result
from guppylang.std.quantum import z, crz, v, vdg, cz, measure, measure_array

@guppy
def oracle(q: array[qubit, 20], target: qubit) -> None:
    r"""Oracle for Grover's search.

    For simplicity, the marked state is hard-coded, and we use a single cnx gate to implement the oracle.
    """
    x(q[0])
    x(q[3])
    x(q[4])
    x(q[7])
    x(q[9])
    x(q[10])
    x(q[13])
    x(q[17])
    x(q[19])
    cnx(q, target)
    x(q[0])
    x(q[3])
    x(q[4])
    x(q[7])
    x(q[9])
    x(q[10])
    x(q[13])
    x(q[17])
    x(q[19])

@guppy
def grover_diffusion(q: array[qubit, 20], target: qubit) -> None:
    r"""Grover diffusion operator
    """

    h(q[0])
    h(q[1])
    h(q[2])
    h(q[3])
    h(q[4])
    h(q[5])
    h(q[6])
    h(q[7])
    h(q[8])
    h(q[9])
    h(q[10])
    h(q[11])
    h(q[12])
    h(q[13])
    h(q[14])
    h(q[15])
    h(q[16])
    h(q[17])
    h(q[18])
    h(q[19])
    x(q[0])
    x(q[1])
    x(q[2])
    x(q[3])
    x(q[4])
    x(q[5])
    x(q[6])
    x(q[7])
    x(q[8])
    x(q[9])
    x(q[10])
    x(q[11])
    x(q[12])
    x(q[13])
    x(q[14])
    x(q[15])
    x(q[16])
    x(q[17])
    x(q[18])
    x(q[19])
    z(target)
    cnx(q, target)
    z(target)
    x(q[0])
    x(q[1])
    x(q[2])
    x(q[3])
    x(q[4])
    x(q[5])
    x(q[6])
    x(q[7])
    x(q[8])
    x(q[9])
    x(q[10])
    x(q[11])
    x(q[12])
    x(q[13])
    x(q[14])
    x(q[15])
    x(q[16])
    x(q[17])
    x(q[18])
    x(q[19])
    h(q[0])
    h(q[1])
    h(q[2])
    h(q[3])
    h(q[4])
    h(q[5])
    h(q[6])
    h(q[7])
    h(q[8])
    h(q[9])
    h(q[10])
    h(q[11])
    h(q[12])
    h(q[13])
    h(q[14])
    h(q[15])
    h(q[16])
    h(q[17])
    h(q[18])
    h(q[19])

@guppy
def weak_probe_rotation(q: array[qubit, 20], target: qubit, probe: qubit) -> None:
    r"""E_{k,Q} unitary to prepare a probe qubit for weak measurement
    Rotation angle .1 halfturns gives k = sin(.1)^2 ~= .01
    """

    oracle(q, target)
    cz(target, probe)
    v(probe)
    crz(target, probe, angle(.1))
    vdg(probe)
    # Oracles are always self-inverse
    oracle(q, target)

@guppy
def grover_with_weak_measurement() -> None:
    # Prepare initial state for Grover's
    q = array(qubit() for _ in range(20))
    h(q[0])
    h(q[1])
    h(q[2])
    h(q[3])
    h(q[4])
    h(q[5])
    h(q[6])
    h(q[7])
    h(q[8])
    h(q[9])
    h(q[10])
    h(q[11])
    h(q[12])
    h(q[13])
    h(q[14])
    h(q[15])
    h(q[16])
    h(q[17])
    h(q[18])
    h(q[19])
    target = qubit()
    h(target)

    # k-while loop
    probe = qubit()
    while not measure(probe):
        # Apply Grover iterate
        oracle(q, target)
        grover_diffusion(q, target)
        # Reset probe
        probe = qubit()
        # Apply weak measurement unitary
        weak_probe_rotation(q, target, probe)

    result("marked", measure_array(q))
    discard(target)

# from hugr.cli import mermaid
from hugr.hugr.render import DotRenderer

# print(mermaid(grover_with_weak_measurement.compile().to_bytes()))

grover_with_weak_measurement.compile().modules[0].render_dot().save("grover_with_weak.dot")