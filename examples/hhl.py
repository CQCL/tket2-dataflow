# from hugr import Hugr
from hugr.cli import mermaid

# fr = open("hhl.txt", "r")
# h = Hugr.from_str(fr.read())
# fw = open("hhl.mmd", "w")
# fr.write(mermaid(h.to_bytes()))

from guppylang import guppy
from guppylang.std.builtins import mem_swap, result
from guppylang.std.quantum import qubit, h, x, z, t, s, v, vdg, tdg, cx, cz, qubit, discard, measure

@guppy
def cs(ctrl: qubit, trgt: qubit) -> None:
    t(ctrl)
    t(trgt)
    cx(ctrl, trgt)
    tdg(trgt)
    cx(ctrl, trgt)

@guppy
def csdg(ctrl: qubit, trgt: qubit) -> None:
    tdg(ctrl)
    tdg(trgt)
    cx(ctrl, trgt)
    t(trgt)
    cx(ctrl, trgt)

@guppy
def qpe(phase_reg_q0: qubit, phase_reg_q1: qubit, state_reg_q0: qubit, state_reg_q1: qubit) -> None:
    # qft dagger
    mem_swap(phase_reg_q0, phase_reg_q1)
    h(phase_reg_q1)
    csdg(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q0)

    cx(phase_reg_q0, phase_reg_q1)
    h(state_reg_q0)
    h(state_reg_q1)
    cs(state_reg_q1, phase_reg_q0)
    cs(state_reg_q1, phase_reg_q1)
    cx(phase_reg_q0, phase_reg_q1)
    cz(phase_reg_q1, state_reg_q0)

    # qft
    h(phase_reg_q0)
    cs(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q1)
    mem_swap(phase_reg_q0, phase_reg_q1)

    # qft dagger
    mem_swap(state_reg_q0, state_reg_q1)
    h(state_reg_q1)
    csdg(state_reg_q1, state_reg_q0)
    h(state_reg_q0)

@guppy
def qpe_dg(phase_reg_q0: qubit, phase_reg_q1: qubit, state_reg_q0: qubit, state_reg_q1: qubit) -> None:
    # qft dagger
    mem_swap(phase_reg_q0, phase_reg_q1)
    h(phase_reg_q1)
    csdg(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q0)

    # qft
    h(state_reg_q0)
    cs(state_reg_q1, state_reg_q0)
    h(state_reg_q1)
    mem_swap(state_reg_q0, state_reg_q1)

    cz(phase_reg_q1, state_reg_q0)
    cx(phase_reg_q0, phase_reg_q1)
    csdg(state_reg_q1, phase_reg_q1)
    csdg(state_reg_q1, phase_reg_q0)
    h(state_reg_q0)
    h(state_reg_q1)
    cx(phase_reg_q0, phase_reg_q1)

    # qft
    h(phase_reg_q0)
    cs(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q1)
    mem_swap(phase_reg_q0, phase_reg_q1)

@guppy
def invert_eigenvalue(state_reg_q0: qubit, state_reg_q1: qubit, check: qubit) -> None:
    cx(state_reg_q1, check)
    x(state_reg_q0)
    cx(state_reg_q0, check)
    x(state_reg_q0)
    v(check)
    t(check)
    cx(state_reg_q0, check)
    tdg(check)
    cx(state_reg_q0, check)
    vdg(check)

@guppy
def hhl(phase_reg_q0: qubit,
        phase_reg_q1: qubit,
        state_reg_q0: qubit,
        state_reg_q1: qubit,
        check: qubit) -> None:
    h(phase_reg_q0)
    h(phase_reg_q1)
    s(phase_reg_q0)
    z(phase_reg_q1)
    # qpe
    # qft dagger
    mem_swap(phase_reg_q0, phase_reg_q1)
    h(phase_reg_q1)
    csdg(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q0)

    cx(phase_reg_q0, phase_reg_q1)
    h(state_reg_q0)
    h(state_reg_q1)
    cs(state_reg_q1, phase_reg_q0)
    cs(state_reg_q1, phase_reg_q1)
    cx(phase_reg_q0, phase_reg_q1)
    cz(phase_reg_q1, state_reg_q0)

    # qft
    h(phase_reg_q0)
    cs(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q1)
    mem_swap(phase_reg_q0, phase_reg_q1)

    # qft dagger
    mem_swap(state_reg_q0, state_reg_q1)
    h(state_reg_q1)
    csdg(state_reg_q1, state_reg_q0)
    h(state_reg_q0)
    # invert eigenvalue
    cx(state_reg_q1, check)
    x(state_reg_q0)
    cx(state_reg_q0, check)
    x(state_reg_q0)
    v(check)
    t(check)
    cx(state_reg_q0, check)
    tdg(check)
    cx(state_reg_q0, check)
    vdg(check)
    # qpe dg
    # qft dagger
    mem_swap(phase_reg_q0, phase_reg_q1)
    h(phase_reg_q1)
    csdg(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q0)

    # qft
    h(state_reg_q0)
    cs(state_reg_q1, state_reg_q0)
    h(state_reg_q1)
    mem_swap(state_reg_q0, state_reg_q1)

    cz(phase_reg_q1, state_reg_q0)
    cx(phase_reg_q0, phase_reg_q1)
    csdg(state_reg_q1, phase_reg_q1)
    csdg(state_reg_q1, phase_reg_q0)
    h(state_reg_q0)
    h(state_reg_q1)
    cx(phase_reg_q0, phase_reg_q1)

    # qft
    h(phase_reg_q0)
    cs(phase_reg_q1, phase_reg_q0)
    h(phase_reg_q1)
    mem_swap(phase_reg_q0, phase_reg_q1)

from guppylang.std.builtins import result, barrier, array
from guppylang.std.quantum import measure, measure_array, discard_array
from guppylang.std.debug import state_result

@guppy
def main_jpmc() -> None:
    q_data_0, q_data_1, q_qpe_0, q_qpe_1, check = qubit(), qubit(), qubit(), qubit(), qubit()

    h(q_data_0)
    h(q_data_1)
    s(q_data_0)
    z(q_data_1)
    # qpe
    # qft dagger
    mem_swap(q_data_0, q_data_1)
    h(q_data_1)
    csdg(q_data_1, q_data_0)
    h(q_data_0)

    cx(q_data_0, q_data_1)
    h(q_qpe_0)
    h(q_qpe_1)
    cs(q_qpe_1, q_data_0)
    cs(q_qpe_1, q_data_1)
    cx(q_data_0, q_data_1)
    cz(q_data_1, q_qpe_0)

    # qft
    h(q_data_0)
    cs(q_data_1, q_data_0)
    h(q_data_1)
    mem_swap(q_data_0, q_data_1)

    # qft dagger
    mem_swap(q_qpe_0, q_qpe_1)
    h(q_qpe_1)
    csdg(q_qpe_1, q_qpe_0)
    h(q_qpe_0)
    # invert eigenvalue
    cx(q_qpe_1, check)
    x(q_qpe_0)
    cx(q_qpe_0, check)
    x(q_qpe_0)
    v(check)
    t(check)
    cx(q_qpe_0, check)
    tdg(check)
    cx(q_qpe_0, check)
    vdg(check)
    # qpe dg
    # qft dagger
    mem_swap(q_data_0, q_data_1)
    h(q_data_1)
    csdg(q_data_1, q_data_0)
    h(q_data_0)

    # qft
    h(q_qpe_0)
    cs(q_qpe_1, q_qpe_0)
    h(q_qpe_1)
    mem_swap(q_qpe_0, q_qpe_1)

    cz(q_data_1, q_qpe_0)
    cx(q_data_0, q_data_1)
    csdg(q_qpe_1, q_data_1)
    csdg(q_qpe_1, q_data_0)
    h(q_qpe_0)
    h(q_qpe_1)
    cx(q_data_0, q_data_1)

    # qft
    h(q_data_0)
    cs(q_data_1, q_data_0)
    h(q_data_1)
    mem_swap(q_data_0, q_data_1)

    result("measure_0", measure(q_qpe_0))
    result("measure_1", measure(q_qpe_1))
    result("ancilla", measure(check))

    # Print state for debugging purposes (would fail if submitted to the device)
    state_result("data", q_data_0, q_data_1)

    # In JPMC they'd now do tomography. We skip it here, but we
    # could do it within Guppy by choosing the measurement
    # basis depending on the shot number.
    discard(q_data_0)
    discard(q_data_1)

# print(dir(main_jpmc.compile()))

print(mermaid(main_jpmc.compile().to_bytes()))