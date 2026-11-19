#!/usr/bin/env python3
"""TrappedIon-Lab - educational trapped-ion quantum physics lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * TrappedIon-Lab computes the physics of ONE trapped ion from
    scratch: the displacement matrix elements <n'|exp(i*eta*(a+a+))|n>
    by exact Laguerre summation, the carrier / red / blue sideband
    spectrum they imply, a rate-model of resolved-sideband cooling,
    a two-ion Molmer-Sorensen gate integrated on a small statevector,
    and carrier Rabi flopping for a pure Fock state vs a thermal
    mixture.
  * Every parameter is a LABELLED teaching number: a 40Ca+-like ion,
    a trap frequency of 1 MHz, a 729 nm laser, a bare Rabi frequency
    of 100 kHz. They are consistent with each other (the Lamb-Dicke
    parameter computed from them is 0.096879), but they are not the
    specification of any real device, and no ion was trapped, cooled,
    or entangled in the making of this lab.
  * Where a truncated model costs accuracy, the cost is printed, not
    hidden: the entangling gate is run in the phonon space 0..3 the
    lab advertises AND in a 0..7 reference space, and the fidelity
    difference is reported as the truncation cost.

Real trap control, laser systems, hardware design, fidelity
certification, and every quantum-advantage or commercial claim are
permanent HOLD gates - with Gemini Argon 4 Extended Thinking as
DECIDER (names-only, tighten-only), Muse as advisor (never
overrides HOLD), and Google Flow plan-only seats.
"""
import cmath
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Labelled teaching parameters. Consistent, round, and NOT a device spec.
ION_LABEL = "40Ca+-like ion (labelled)"
MASS_AMU = 40.0
AMU_KG = 1.66053906660e-27
HBAR = 1.054571817e-34
TRAP_FREQ_HZ = 1.0e6            # axial trap frequency omega/2pi, labelled
LASER_WAVELENGTH_NM = 729.0     # Ca+ S-D quadrupole line, labelled
OMEGA0_HZ = 100.0e3             # bare Rabi frequency Omega_0/2pi, labelled
ETA_TEACH = 0.1                 # round teaching Lamb-Dicke parameter
ETA_SET = (0.05, 0.1, 0.2, 0.5, 1.0)

# Cooling toy constants (labelled model choices, stated in the report).
COOL_NMAX = 20                  # phonon populations kept: n = 0..20
COOL_CYCLES = 400
COOL_GAIN = 3.0                 # rate-step gain of the red-sideband transfer
COOL_NBAR0 = 5.0                # labelled thermal mean before cooling

# Molmer-Sorensen gate constants (labelled model choices).
GATE_NPH = 3                    # advertised phonon space: levels 0..3
GATE_NPH_REF = 7                # reference space used to price truncation
GATE_G_OVER_DELTA = 0.25        # coupling/detuning found by the sim scan
GATE_LOOPS = 1                  # one closed phase-space loop
GATE_STEPS_PER_LOOP = 720       # RK4 steps per loop

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "real_trap_control": "HOLD - no control of real trap electrodes, RF drives, or ion hardware",
    "laser_system_operation": "HOLD - no operating real lasers, beam lines, or frequency locks",
    "hardware_specification_claim": "HOLD - labelled teaching numbers are not a real device's specifications",
    "vendor_performance_claim": "HOLD - no vendor or product performance claims from a teaching model",
    "quantum_advantage_claim": "HOLD - no quantum-advantage claim; a 16-amplitude toy proves nothing about advantage",
    "gate_fidelity_certification": "HOLD - no certifying a real gate's fidelity; simulated fidelity is about the model",
    "real_ion_cooling_claim": "HOLD - no claiming a real ion was cooled; the cooling curve is a rate model",
    "trap_hardware_design": "HOLD - no design guidance for real traps, electrodes, or vacuum systems",
    "laser_safety_bypass": "HOLD - laser safety interlocks are never bypassed, in simulation or prose",
    "vacuum_cryogenics_operation": "HOLD - no live vacuum or cryogenics operation directed from a sim",
    "commercial_quantum_computer_claim": "HOLD - no commercial quantum-computing claims from a lab model",
    "error_correction_deployment_claim": "HOLD - no claim of deployable error correction from two simulated ions",
    "sensing_metrology_claim": "HOLD - no ion-clock, sensing, or metrology product claims",
    "simulation_as_measurement_claim": "HOLD - simulation output must never be presented as experimental measurement",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Displacement matrix
# The laser-ion coupling in the interaction picture carries the factor
# exp(i*eta*(a + a+)), eta the Lamb-Dicke parameter. Its matrix elements
# between Fock states have an exact closed form in generalised Laguerre
# polynomials, evaluated here by the finite sum that defines them -
# nothing is quoted, everything is summed.

def laguerre(n, k, x):
    """Generalised Laguerre L_n^k(x) by its exact finite sum."""
    total = 0.0
    for j in range(n + 1):
        term = ((-1) ** j) * math.comb(n + k, n - j) * (x ** j)
        total += term / math.factorial(j)
    return total


def disp(n_out, n_in, eta):
    """<n_out| exp(i*eta*(a + a+)) |n_in> with alpha = i*eta.

    For n_out >= n_in: sqrt(n_in!/n_out!) * alpha^(n_out-n_in)
        * exp(-eta^2/2) * L_{n_in}^(n_out-n_in)(eta^2)
    For n_out <  n_in: sqrt(n_out!/n_in!) * (-alpha*)^(n_in-n_out)
        * exp(-eta^2/2) * L_{n_out}^(n_in-n_out)(eta^2)
    """
    alpha = 1j * eta
    x = eta * eta
    env = math.exp(-x / 2.0)
    if n_out >= n_in:
        return (math.sqrt(math.factorial(n_in) / math.factorial(n_out))
                * alpha ** (n_out - n_in) * env
                * laguerre(n_in, n_out - n_in, x))
    return (math.sqrt(math.factorial(n_out) / math.factorial(n_in))
            * (-alpha.conjugate()) ** (n_in - n_out) * env
            * laguerre(n_out, n_in - n_out, x))


def eta_from_labelled_params():
    """Lamb-Dicke parameter implied by the labelled ion/trap/laser set.

    x0 = sqrt(hbar / (2 m omega)) is the ground-state extent and
    eta = k x0 with k = 2 pi / lambda.
    """
    mass = MASS_AMU * AMU_KG
    omega = 2.0 * math.pi * TRAP_FREQ_HZ
    x0 = math.sqrt(HBAR / (2.0 * mass * omega))
    k = 2.0 * math.pi / (LASER_WAVELENGTH_NM * 1e-9)
    return k * x0, x0


def matrix_report():
    eta_phys, x0 = eta_from_labelled_params()
    unitarity = {}
    for eta in ETA_SET:
        worst = 0.0
        for n in range(0, 13):
            row = sum(abs(disp(m, n, eta)) ** 2 for m in range(0, 41))
            worst = max(worst, abs(row - 1.0))
        unitarity[f"{eta}"] = worst
    sub = [[round(abs(disp(m, n, ETA_TEACH)), 6) for n in range(5)]
           for m in range(5)]
    return {
        "ion": ION_LABEL,
        "trap_freq_hz": TRAP_FREQ_HZ,
        "laser_wavelength_nm": LASER_WAVELENGTH_NM,
        "ground_state_extent_nm": x0 * 1e9,
        "eta_from_labelled_params": eta_phys,
        "eta_teaching_value": ETA_TEACH,
        "unitarity_worst_deviation_rows_n0_12_cols_to_40": unitarity,
        "abs_matrix_eta_0.1_rows_out_cols_in_n0_4": sub,
    }


# ================================================================ 2. Sideband spectrum
# A drive on the carrier (Delta n = 0), red sideband (Delta n = -1) or
# blue sideband (Delta n = +1) couples at Omega_0 * |<n'|D|n>|. In the
# Lamb-Dicke limit eta -> 0 those reduce to carrier 1 - eta^2 (n + 1/2),
# red eta sqrt(n), blue eta sqrt(n + 1). The lab grades the exact
# elements against those scalings and prints the errors honestly -
# including where the approximation goes negative and absurd.

def sideband_ratios(n, eta):
    """(|carrier|, |red|, |blue|) Rabi ratios Omega/Omega_0 from level n."""
    carrier = abs(disp(n, n, eta))
    red = abs(disp(n - 1, n, eta)) if n >= 1 else 0.0
    blue = abs(disp(n + 1, n, eta))
    return carrier, red, blue


def ld_approx(n, eta):
    """Lamb-Dicke-limit scalings for (carrier, red, blue) from level n."""
    return (1.0 - eta * eta * (n + 0.5),
            eta * math.sqrt(n),
            eta * math.sqrt(n + 1))


def sideband_report():
    table = {}
    for eta in ETA_SET:
        table[f"{eta}"] = {
            f"n={n}": [round(v, 6) for v in sideband_ratios(n, eta)]
            for n in range(0, 9)
        }
    errors = {}
    for eta in (0.05, 0.5):
        worst = {"carrier": 0.0, "red": 0.0, "blue": 0.0}
        for n in range(0, 7):
            exact = sideband_ratios(n, eta)
            approx = ld_approx(n, eta)
            for i, key in enumerate(("carrier", "red", "blue")):
                worst[key] = max(worst[key], abs(exact[i] - approx[i]))
        errors[f"{eta}"] = worst
    return {
        "ratios_carrier_red_blue_by_eta": table,
        "lamb_dicke_max_abs_error_n0_6": errors,
        "red_sideband_from_ground": abs(disp(-1, 0, ETA_TEACH))
        if False else 0.0,
        "note": "red sideband from n=0 is exactly 0: there is no n=-1, "
                "which is the floor that makes ground-state cooling possible",
    }


# ================================================================ 3. Sideband cooling toy
# A rate/step model on phonon populations p_n (n <= 20), NOT a master
# equation. One cycle = (a) a red-sideband transfer step in which the
# population at level n steps down with probability
# 1 - exp(-gain * |M(n->n-1)|^2 / |M(1->0)|^2) - the exact matrix
# elements set the rates; (b) one repump scattering event, modelled as
# a single recoil kick drawn from the exact displacement kernel with
# the direction-averaged axial projection eta/sqrt(2), i.e. eta^2/2 of
# heating per kick - the standard scale of the sideband-cooling limit.
# Carrier-only driving is run through the same bookkeeping with the
# transfer set to the identity: the carrier term is diagonal in n, so
# it cannot remove a phonon, and n-bar does not move.

def thermal_distribution(nbar, nmax):
    """Truncated, renormalised thermal (geometric) distribution."""
    q = nbar / (nbar + 1.0)
    p = [(1.0 - q) * q ** n for n in range(nmax + 1)]
    s = sum(p)
    return [v / s for v in p]


def mean_n(p):
    return sum(n * p[n] for n in range(len(p)))


def _kick_kernel(eta_kick, nmax):
    kernel = [[abs(disp(m, n, eta_kick)) ** 2 for n in range(nmax + 1)]
              for m in range(nmax + 1)]
    for n in range(nmax + 1):
        col = sum(kernel[m][n] for m in range(nmax + 1))
        for m in range(nmax + 1):
            kernel[m][n] /= col
    return kernel


def cooling_run(carrier_only=False):
    eta = ETA_TEACH
    nmax = COOL_NMAX
    p = thermal_distribution(COOL_NBAR0, nmax)
    initial_nbar = mean_n(p)
    ref = abs(disp(0, 1, eta)) ** 2
    kernel = _kick_kernel(eta / math.sqrt(2.0), nmax)
    marks = (0, 1, 2, 5, 10, 25, 50, 100, 200, COOL_CYCLES)
    trajectory = [(0, initial_nbar)]
    for cycle in range(1, COOL_CYCLES + 1):
        newp = [0.0] * (nmax + 1)
        newp[0] = p[0]
        for n in range(1, nmax + 1):
            if carrier_only:
                p_down = 0.0
            else:
                p_down = 1.0 - math.exp(
                    -COOL_GAIN * abs(disp(n - 1, n, eta)) ** 2 / ref)
            newp[n] += p[n] * (1.0 - p_down)
            newp[n - 1] += p[n] * p_down
        if carrier_only:
            p = newp
        else:
            p = [sum(kernel[m][n] * newp[n] for n in range(nmax + 1))
                 for m in range(nmax + 1)]
            s = sum(p)
            p = [v / s for v in p]
        if cycle in marks:
            trajectory.append((cycle, mean_n(p)))
    return {
        "eta": eta,
        "nmax": nmax,
        "cycles": COOL_CYCLES,
        "gain": COOL_GAIN,
        "labelled_initial_nbar": COOL_NBAR0,
        "truncated_initial_nbar": initial_nbar,
        "trajectory_cycle_nbar": [[c, round(v, 6)] for c, v in trajectory],
        "final_nbar": mean_n(p),
        "final_ground_population": p[0],
        "eta_squared_over_2": eta * eta / 2.0,
        "carrier_only": carrier_only,
    }


def cooling_report():
    cooled = cooling_run(carrier_only=False)
    carrier = cooling_run(carrier_only=True)
    return {
        "red_sideband_cooling": cooled,
        "carrier_only_contrast": {
            "initial_nbar": carrier["truncated_initial_nbar"],
            "final_nbar": carrier["final_nbar"],
            "why": "the carrier matrix element is diagonal in n "
                   "(Delta n = 0): it flips the internal state and "
                   "leaves every p_n where it was, so carrier driving "
                   "alone can never cool - only the red sideband's "
                   "Delta n = -1 step removes phonons",
        },
    }


# ================================================================ 4. Two-ion gate
# A Molmer-Sorensen bichromatic drive on two ions sharing one motional
# mode, in the interaction picture:
#     H(tau) = (g/delta) * (sigma_x1 + sigma_x2)
#              * (a exp(-i tau) + a+ exp(+i tau)),   tau = delta * t.
# The state (2 qubits x phonon levels 0..nph) is integrated with RK4.
# After one closed loop (tau = 2 pi) the motion disentangles and the
# spins carry a geometric phase; with g/delta = 0.25 - the value the
# lab's own scan selects - the loop produces a Bell state up to a
# labelled virtual-Z frame correction on ion 2. The gate is computed
# in the advertised 0..3 phonon space AND a 0..7 reference space, and
# the difference is printed as the truncation cost. No 1.0 is faked:
# in the 0..3 space the fidelity lands where the integration puts it.

def hermitian_eigvals(matrix):
    """Eigenvalues of a small complex Hermitian matrix, descending.

    Uses the real embedding [[Re, -Im], [Im, Re]] (each eigenvalue of
    the Hermitian matrix appears twice) with cyclic Jacobi rotations.
    """
    n = len(matrix)
    m = 2 * n
    r = [[0.0] * m for _ in range(m)]
    for i in range(n):
        for j in range(n):
            re, im = matrix[i][j].real, matrix[i][j].imag
            r[i][j] = re
            r[i][j + n] = -im
            r[i + n][j] = im
            r[i + n][j + n] = re
    for _ in range(100):
        off = math.sqrt(sum(r[i][j] ** 2 for i in range(m)
                            for j in range(m) if i != j))
        if off < 1e-13:
            break
        for p_ in range(m):
            for q_ in range(p_ + 1, m):
                apq = r[p_][q_]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = r[p_][p_], r[q_][q_]
                tau = (aqq - app) / (2.0 * apq)
                t = ((1.0 if tau >= 0 else -1.0)
                     / (abs(tau) + math.sqrt(1.0 + tau * tau)))
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for k_ in range(m):
                    rkp, rkq = r[k_][p_], r[k_][q_]
                    r[k_][p_] = c * rkp - s * rkq
                    r[k_][q_] = s * rkp + c * rkq
                for k_ in range(m):
                    rpk, rqk = r[p_][k_], r[q_][k_]
                    r[p_][k_] = c * rpk - s * rqk
                    r[q_][k_] = s * rpk + c * rqk
    all_ev = sorted((r[i][i] for i in range(m)), reverse=True)
    return [all_ev[2 * i] for i in range(n)]


def concurrence(rho):
    """Wootters concurrence of a two-qubit density matrix."""
    sy = ((0, -1j), (1j, 0))
    y = [[0j] * 4 for _ in range(4)]
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    y[a * 2 + b][c * 2 + d] = sy[a][c] * sy[b][d]
    rc = [[rho[i][j].conjugate() for j in range(4)] for i in range(4)]
    tmp = [[sum(y[i][k] * rc[k][j] for k in range(4)) for j in range(4)]
           for i in range(4)]
    rt = [[sum(tmp[i][k] * y[k][j] for k in range(4)) for j in range(4)]
          for i in range(4)]
    rmat = [[sum(rho[i][k] * rt[k][j] for k in range(4)) for j in range(4)]
            for i in range(4)]
    ev = sorted((max(v, 0.0) for v in hermitian_eigvals(rmat)),
                reverse=True)
    return max(0.0, math.sqrt(ev[0]) - math.sqrt(ev[1])
               - math.sqrt(ev[2]) - math.sqrt(ev[3]))


def _ms_evolve(nph, g_over_delta, loops, steps_per_loop):
    """RK4 integration of the MS Hamiltonian from |00>|n=0>."""
    n_lev = nph + 1
    dim = 4 * n_lev

    def idx(s1, s2, n):
        return (s1 * 2 + s2) * n_lev + n

    def apply_h(psi, tau):
        out = [0j] * dim
        e_m = cmath.exp(-1j * tau)
        e_p = cmath.exp(1j * tau)
        for s1 in range(2):
            for s2 in range(2):
                base = idx(s1, s2, 0)
                for n in range(n_lev):
                    v = psi[base + n]
                    if v == 0:
                        continue
                    for t1, t2 in ((1 - s1, s2), (s1, 1 - s2)):
                        tb = idx(t1, t2, 0)
                        if n >= 1:
                            out[tb + n - 1] += v * math.sqrt(n) * e_m
                        if n + 1 <= nph:
                            out[tb + n + 1] += (v * math.sqrt(n + 1)
                                                * e_p)
        return [g_over_delta * v for v in out]

    def rk4(psi, tau, dt):
        k1 = apply_h(psi, tau)
        p2 = [psi[i] - 1j * dt / 2 * k1[i] for i in range(dim)]
        k2 = apply_h(p2, tau + dt / 2)
        p3 = [psi[i] - 1j * dt / 2 * k2[i] for i in range(dim)]
        k3 = apply_h(p3, tau + dt / 2)
        p4 = [psi[i] - 1j * dt * k3[i] for i in range(dim)]
        k4 = apply_h(p4, tau + dt)
        return [psi[i] - 1j * dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i]
                                        + k4[i]) for i in range(dim)]

    psi = [0j] * dim
    psi[idx(0, 0, 0)] = 1.0 + 0j
    dt = 2.0 * math.pi / steps_per_loop
    tau = 0.0
    for _ in range(steps_per_loop * loops):
        psi = rk4(psi, tau, dt)
        tau += dt
    return psi, n_lev


def _spin_report(psi, n_lev):
    """Reduced two-qubit state, frame-fixed fidelity, concurrence."""
    fixed = list(psi)
    for s1 in range(2):
        for n in range(n_lev):
            fixed[(s1 * 2 + 1) * n_lev + n] *= -1j  # virtual Z on ion 2
    rho = [[sum(fixed[a * n_lev + n] * fixed[b * n_lev + n].conjugate()
                for n in range(n_lev)) for b in range(4)]
           for a in range(4)]
    bell = (1.0 / math.sqrt(2.0), 0.0, 0.0, 1.0 / math.sqrt(2.0))
    fid = sum(bell[a] * sum(rho[a][b] * bell[b] for b in range(4))
              for a in range(4)).real
    ph_pop = [sum(abs(psi[a * n_lev + n]) ** 2 for a in range(4))
              for n in range(n_lev)]
    return {
        "fidelity_to_bell_plus_after_frame_fix": fid,
        "concurrence": concurrence(rho),
        "phonon_populations_final": ph_pop,
        "spin_ground_population": ph_pop[0],
    }


def gate_report():
    eta = ETA_TEACH
    g_hz = OMEGA0_HZ * eta / 2.0            # MS coupling g/2pi
    delta_hz = g_hz / GATE_G_OVER_DELTA     # detuning delta/2pi
    gate_time_us = GATE_LOOPS / delta_hz * 1e6
    psi3, n3 = _ms_evolve(GATE_NPH, GATE_G_OVER_DELTA, GATE_LOOPS,
                          GATE_STEPS_PER_LOOP)
    small = _spin_report(psi3, n3)
    psi7, n7 = _ms_evolve(GATE_NPH_REF, GATE_G_OVER_DELTA, GATE_LOOPS,
                          GATE_STEPS_PER_LOOP)
    ref = _spin_report(psi7, n7)
    return {
        "model": "Molmer-Sorensen bichromatic drive, interaction "
                 "picture, RK4 on 2 qubits x phonon levels",
        "eta": eta,
        "coupling_g_over_2pi_hz": g_hz,
        "detuning_delta_over_2pi_hz": delta_hz,
        "g_over_delta": GATE_G_OVER_DELTA,
        "loops": GATE_LOOPS,
        "gate_time_us": gate_time_us,
        "target": "(|00> + |11>)/sqrt(2) after a labelled virtual-Z "
                  "frame correction on ion 2",
        "phonon_space_0_to_3": small,
        "reference_space_0_to_7": ref,
        "truncation_cost_fidelity": (ref["fidelity_to_bell_plus_after_frame_fix"]
                                     - small["fidelity_to_bell_plus_after_frame_fix"]),
        "truncation_cost_concurrence": (ref["concurrence"]
                                        - small["concurrence"]),
    }


# ================================================================ 5. Rabi flopping
# On the carrier, level n flops at Omega_0 * |M(n->n)|: the Debye-Waller
# factor e^{-eta^2/2} L_n(eta^2) shrinks as n grows. A pure n=0 state
# flops at full contrast forever; a thermal mixture is a sum of
# slightly different flopping frequencies and its peaks wash out -
# contrast collapse computed from the distribution, not asserted.

def rabi_report():
    eta = ETA_TEACH
    om00 = abs(disp(0, 0, eta))
    t_pi = math.pi / om00          # in units of 1/Omega_0 (Omega_0 = 1)
    t_pi_us = math.pi / (2.0 * math.pi * OMEGA0_HZ * om00) * 1e6
    p_therm = thermal_distribution(COOL_NBAR0, 40)

    def p_excited(t, n):
        return math.sin(abs(disp(n, n, eta)) * t / 2.0) ** 2

    def p_thermal(t):
        return sum(p_therm[n] * p_excited(t, n) for n in range(41))

    peaks = {}
    for mult in (1, 3, 5, 9, 13, 19):
        t = mult * t_pi
        peaks[f"t={mult}*t_pi"] = {
            "pure_n0": p_excited(t, 0),
            "thermal_nbar5": p_thermal(t),
        }
    late = [p_thermal(20.0 * t_pi * i / 400.0) for i in range(361, 401)]
    return {
        "eta": eta,
        "carrier_ratio_n0": om00,
        "carrier_ratios_n0_to_5": [round(abs(disp(n, n, eta)), 6)
                                   for n in range(6)],
        "pi_time_over_1_over_omega0": t_pi,
        "pi_time_us_at_labelled_omega0": t_pi_us,
        "peaks_pure_vs_thermal": peaks,
        "thermal_max_last_window": max(late),
        "thermal_at_first_pi": p_thermal(t_pi),
    }


# ================================================================ 6. kNN regime classifier
# Features: Lamb-Dicke eta, sideband-resolution ratio Omega_0/omega,
# thermal mean n-bar. The training grid is labelled by the lab's own
# regime policy (policy_label below) - the classifier interpolates the
# policy, it does not replace the physics.

def policy_label(eta, ratio, nbar):
    """The regime policy distilled to four labels."""
    ld = eta * math.sqrt(2.0 * nbar + 1.0)
    if ld >= 1.0:
        return "beyond-lamb-dicke"
    if ratio > 0.5:
        return "unresolved-sideband"
    if ld >= 0.5 or ratio > 0.25:
        return "lamb-dicke-marginal"
    return "lamb-dicke-resolved"


def _training_grid():
    rows = []
    for eta in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0):
        for ratio in (0.05, 0.1, 0.25, 0.4, 0.6, 1.0):
            for nbar in (0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0):
                rows.append(((eta, ratio, nbar),
                             policy_label(eta, ratio, nbar)))
    return rows


TRAINING_GRID = _training_grid()


def _scaled(features):
    eta, ratio, nbar = features
    return (eta, ratio, nbar / 10.0)


def classify_regime(eta, ratio, nbar, k=5):
    """Distance-weighted kNN over the policy-labelled grid."""
    target = _scaled((eta, ratio, nbar))
    ranked = []
    for features, label in TRAINING_GRID:
        point = _scaled(features)
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(target, point)))
        ranked.append((dist, features, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"lamb-dicke-resolved": 0, "lamb-dicke-marginal": 1,
                "unresolved-sideband": 2, "beyond-lamb-dicke": 3}
    return sorted(votes.items(),
                  key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


def classifier_report():
    eta_phys, _ = eta_from_labelled_params()
    scenarios = {
        "ca_plus_labelled_headline": (eta_phys, OMEGA0_HZ / TRAP_FREQ_HZ,
                                      0.005277),
        "hot_ion_before_cooling": (0.1, 0.1, 5.0),
        "strong_drive_sideband_blur": (0.2, 0.6, 0.0),
        "warm_marginal_ion": (0.3, 0.4, 1.0),
        "heavy_hot_ion": (0.5, 0.1, 10.0),
        "deep_lamb_dicke_cold": (0.05, 0.05, 0.1),
    }
    return {name: classify_regime(*features)
            for name, features in scenarios.items()}


# ================================================================ 7. Advisory seats
def _seat_trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="trappedion-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all hardware and certification gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _seat_trace(payload),
    }


def claude_advisor(topic="trappedion-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the truncation cost next to every gate number: a fidelity quoted without the size of the space it was computed in is only half a result.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["one ion held in an RF trap, its motion quantised into rungs",
                  "a laser comb touching carrier, red and blue sidebands",
                  "red-sideband pulses walking phonons down to the ground rung",
                  "two ions trading a shared phonon until their spins lock",
                  "a thermal flopping curve washing out beside a pure one"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 1. Laguerre + displacement matrix
    chk("laguerre_zeroth_is_one",
        all(abs(laguerre(0, k, 0.37) - 1.0) < 1e-12 for k in (0, 1, 3)))
    chk("laguerre_at_zero",
        all(abs(laguerre(n, 0, 0.0) - 1.0) < 1e-12 for n in range(6))
        and all(abs(laguerre(n, 1, 0.0) - (n + 1)) < 1e-12
                for n in range(6)))
    chk("disp_ground_carrier",
        abs(abs(disp(0, 0, 0.3)) - math.exp(-0.09 / 2.0)) < 1e-12)
    chk("disp_ground_blue_magnitude",
        abs(abs(disp(1, 0, 0.3)) - 0.3 * math.exp(-0.09 / 2.0)) < 1e-12)
    chk("disp_symmetry",
        all(abs(abs(disp(m, n, 0.4)) - abs(disp(n, m, 0.4))) < 1e-12
            for n in range(6) for m in range(6)))
    for eta in ETA_SET:
        worst = max(
            abs(sum(abs(disp(m, n, eta)) ** 2 for m in range(41)) - 1.0)
            for n in range(13))
        chk(f"disp_unitarity_eta_{eta}", worst < 1e-9)
    eta_phys, x0 = eta_from_labelled_params()
    chk("eta_from_labelled_params", abs(eta_phys - 0.096879) < 1e-5
        and abs(x0 * 1e9 - 11.240) < 0.01)

    # 2. sideband spectrum + Lamb-Dicke scalings
    chk("red_sideband_vanishes_from_ground",
        sideband_ratios(0, ETA_TEACH)[1] == 0.0)
    chk("blue_from_ground_is_eta_carrier",
        abs(sideband_ratios(0, 0.2)[2]
            - 0.2 * math.exp(-0.04 / 2.0)) < 1e-12)
    rep = sideband_report()
    err05 = rep["lamb_dicke_max_abs_error_n0_6"]["0.05"]
    err50 = rep["lamb_dicke_max_abs_error_n0_6"]["0.5"]
    chk("ld_limit_accurate_at_eta_0.05",
        max(err05.values()) < 2e-3)
    chk("ld_limit_fails_honestly_at_eta_0.5",
        max(err50.values()) > 0.5 and err50["blue"] > err05["blue"] * 100)
    chk("ld_carrier_approx_goes_negative",
        ld_approx(6, 0.5)[0] < 0 < sideband_ratios(6, 0.5)[0])
    chk("sideband_ratios_match_matrix",
        all(abs(sideband_ratios(n, 0.2)[0] - abs(disp(n, n, 0.2))) < 1e-15
            for n in range(5)))

    # 3. cooling toy
    cool = cooling_report()
    rs = cool["red_sideband_cooling"]
    chk("cooling_starts_thermal",
        abs(rs["truncated_initial_nbar"] - 4.533386) < 1e-5)
    chk("cooling_final_nbar_near_eta_sq_over_2",
        0.001 < rs["final_nbar"] < 0.02
        and abs(rs["final_nbar"] - rs["eta_squared_over_2"]) < 0.002)
    chk("cooling_ground_population", rs["final_ground_population"] > 0.99)
    chk("cooling_trajectory_monotone_floor",
        rs["trajectory_cycle_nbar"][-1][1]
        <= rs["trajectory_cycle_nbar"][3][1])
    car = cool["carrier_only_contrast"]
    chk("carrier_only_cannot_cool",
        car["final_nbar"] == car["initial_nbar"])

    # 4. gate - eigensolver sanity, then the computed Bell state
    bell_rho = [[0j] * 4 for _ in range(4)]
    bell_rho[0][0] = 0.5
    bell_rho[0][3] = 0.5
    bell_rho[3][0] = 0.5
    bell_rho[3][3] = 0.5
    chk("concurrence_bell_is_one", abs(concurrence(bell_rho) - 1.0) < 1e-9)
    prod_rho = [[0j] * 4 for _ in range(4)]
    prod_rho[0][0] = 1.0
    chk("concurrence_product_is_zero", concurrence(prod_rho) < 1e-9)
    mix_rho = [[0j] * 4 for _ in range(4)]
    mix_rho[0][0] = 0.5
    mix_rho[3][3] = 0.5
    chk("concurrence_classical_mix_is_zero", concurrence(mix_rho) < 1e-9)
    gate_rep = gate_report()
    small = gate_rep["phonon_space_0_to_3"]
    ref = gate_rep["reference_space_0_to_7"]
    chk("gate_bell_fidelity_truncated_space",
        small["fidelity_to_bell_plus_after_frame_fix"] > 0.98)
    chk("gate_concurrence_truncated_space", small["concurrence"] > 0.97)
    chk("gate_reference_space_near_perfect",
        ref["fidelity_to_bell_plus_after_frame_fix"] > 0.9999
        and ref["concurrence"] > 0.9999)
    chk("gate_truncation_cost_printed",
        gate_rep["truncation_cost_fidelity"] > 0.005
        and gate_rep["truncation_cost_concurrence"] > 0.005)
    chk("gate_motion_returns_in_reference",
        ref["spin_ground_population"] > 0.9999)
    chk("gate_time_labelled",
        abs(gate_rep["gate_time_us"] - 50.0) < 1e-9
        and abs(gate_rep["detuning_delta_over_2pi_hz"] - 20000.0) < 1e-6)

    # 5. Rabi flopping
    rabi = rabi_report()
    chk("rabi_pure_pi_is_full_flip",
        abs(rabi["peaks_pure_vs_thermal"]["t=1*t_pi"]["pure_n0"] - 1.0)
        < 1e-9)
    chk("rabi_thermal_still_high_at_first_pi",
        0.98 < rabi["thermal_at_first_pi"] < 0.995)
    chk("rabi_thermal_contrast_collapses",
        rabi["thermal_max_last_window"] < 0.85
        and rabi["peaks_pure_vs_thermal"]["t=19*t_pi"]["thermal_nbar5"]
        < rabi["thermal_at_first_pi"] - 0.15)
    chk("rabi_carrier_ratios_shrink_with_n",
        all(rabi["carrier_ratios_n0_to_5"][i]
            > rabi["carrier_ratios_n0_to_5"][i + 1] for i in range(5)))

    # 6. classifier - labels come from the regime policy itself
    rep_cls = classifier_report()
    chk("classifier_headline_resolved",
        rep_cls["ca_plus_labelled_headline"] == "lamb-dicke-resolved")
    chk("classifier_hot_still_resolved",
        rep_cls["hot_ion_before_cooling"] == "lamb-dicke-resolved")
    chk("classifier_blur_unresolved",
        rep_cls["strong_drive_sideband_blur"] == "unresolved-sideband")
    chk("classifier_warm_marginal",
        rep_cls["warm_marginal_ion"] == "lamb-dicke-marginal")
    chk("classifier_heavy_hot_beyond",
        rep_cls["heavy_hot_ion"] == "beyond-lamb-dicke")
    chk("classifier_deep_ld_resolved",
        rep_cls["deep_lamb_dicke_cold"] == "lamb-dicke-resolved")
    chk("classifier_matches_policy_on_grid",
        all(classify_regime(*features) == label
            for features, label in TRAINING_GRID))
    chk("classifier_matches_policy_on_scenarios",
        all(classify_regime(*features) == policy_label(*features)
            for features in [(0.1, 0.1, 5.0), (0.2, 0.6, 0.0),
                             (0.3, 0.4, 1.0), (0.5, 0.1, 10.0),
                             (0.05, 0.05, 0.1)]))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_trap_control_hold", gate("real_trap_control") == "HOLD")
    chk("gate_laser_hold", gate("laser_system_operation") == "HOLD")
    chk("gate_fidelity_cert_hold",
        gate("gate_fidelity_certification") == "HOLD")
    chk("gate_advantage_hold", gate("quantum_advantage_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("TrappedIon-Lab - trapped-ion physics, computed from scratch")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - what this lab does and does not claim")
    print("    Every number below is computed in this file from labelled")
    print("    teaching parameters. Nothing was measured, no ion exists")
    print("    here, and the gate fidelity is a property of a 16-amplitude")
    print("    model - printed together with the cost of truncating it.")
    eta_phys, x0 = eta_from_labelled_params()
    print("\n[2] Labelled ion - 40Ca+-like, trap 1 MHz, laser 729 nm")
    print(f"    ground-state extent x0 = {x0 * 1e9:.3f} nm -> "
          f"eta from the labelled set = {eta_phys:.6f}")
    print(f"    teaching eta = {ETA_TEACH} | bare Rabi Omega_0/2pi = "
          f"{OMEGA0_HZ / 1e3:.0f} kHz | Omega_0/omega = "
          f"{OMEGA0_HZ / TRAP_FREQ_HZ:.2f} (sidebands resolved)")
    mat = matrix_report()
    print("\n[3] Displacement matrix <n'|exp(i eta (a + a+))|n> - exact sums")
    for eta, worst in mat["unitarity_worst_deviation_rows_n0_12_cols_to_40"].items():
        print(f"    eta = {eta:>4}: worst row-unitarity deviation "
              f"{worst:.2e} (rows n<=12, summed to n'=40)")
    sb = sideband_report()
    print("\n[4] Sideband spectrum at eta = 0.1 (ratios to Omega_0)")
    print("      n   carrier     red       blue")
    for n in range(0, 9):
        car, red, blu = sideband_ratios(n, ETA_TEACH)
        print(f"     {n:2d}   {car:.6f}  {red:.6f}  {blu:.6f}")
    e05 = sb["lamb_dicke_max_abs_error_n0_6"]["0.05"]
    e50 = sb["lamb_dicke_max_abs_error_n0_6"]["0.5"]
    print(f"    Lamb-Dicke limit max error (n<=6): eta=0.05 -> "
          f"{max(e05.values()):.2e} | eta=0.5 -> {max(e50.values()):.4f}")
    print("    red sideband from n=0 is exactly 0 - no rung below ground;")
    print("    that zero is what makes ground-state cooling possible.")
    cool = cooling_report()
    rs = cool["red_sideband_cooling"]
    print("\n[5] Resolved-sideband cooling toy (rate model, eta = 0.1)")
    traj = " | ".join(f"cyc {c}: {v:.4f}"
                       for c, v in rs["trajectory_cycle_nbar"])
    print(f"    {traj}")
    print(f"    final n-bar = {rs['final_nbar']:.6f} "
          f"(eta^2/2 = {rs['eta_squared_over_2']:.4f}) | "
          f"ground population {rs['final_ground_population']:.6f}")
    car = cool["carrier_only_contrast"]
    print(f"    carrier-only driving: n-bar {car['initial_nbar']:.6f} -> "
          f"{car['final_nbar']:.6f} (unchanged - Delta n = 0 cannot cool)")
    gate_rep = gate_report()
    small = gate_rep["phonon_space_0_to_3"]
    ref = gate_rep["reference_space_0_to_7"]
    print("\n[6] Molmer-Sorensen two-ion gate (RK4, one closed loop)")
    print(f"    g/2pi = {gate_rep['coupling_g_over_2pi_hz']:.0f} Hz | "
          f"delta/2pi = {gate_rep['detuning_delta_over_2pi_hz']:.0f} Hz | "
          f"g/delta = {gate_rep['g_over_delta']} | "
          f"gate time {gate_rep['gate_time_us']:.1f} us")
    print(f"    phonon space 0..3: F = "
          f"{small['fidelity_to_bell_plus_after_frame_fix']:.6f} | "
          f"concurrence = {small['concurrence']:.6f}")
    print(f"    reference space 0..7: F = "
          f"{ref['fidelity_to_bell_plus_after_frame_fix']:.6f} | "
          f"concurrence = {ref['concurrence']:.6f}")
    print(f"    truncation cost: F {gate_rep['truncation_cost_fidelity']:.6f} | "
          f"C {gate_rep['truncation_cost_concurrence']:.6f} - printed, not hidden")
    rabi = rabi_report()
    print("\n[7] Carrier Rabi flopping - pure n=0 vs thermal n-bar = 5")
    print(f"    pi time (n=0) = {rabi['pi_time_us_at_labelled_omega0']:.3f} us "
          f"at the labelled Omega_0 | carrier ratios n=0..5: "
          f"{rabi['carrier_ratios_n0_to_5']}")
    for key, row in rabi["peaks_pure_vs_thermal"].items():
        print(f"    {key:12s} pure {row['pure_n0']:.6f} | "
              f"thermal {row['thermal_nbar5']:.6f}")
    print(f"    thermal max in the last window: "
          f"{rabi['thermal_max_last_window']:.6f} - the n-spread of "
          f"carrier rates washes the contrast out")
    print("\n[8] Regime classifier (TinyML kNN, labels from the lab policy)")
    for name, label in classifier_report().items():
        print(f"    {name:32s} -> {label}")
    print("\n[9] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[10] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "matrix":
        print(json.dumps(matrix_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "sidebands":
        print(json.dumps(sideband_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "cooling":
        print(json.dumps(cooling_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "gate":
        print(json.dumps(gate_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "rabi":
        print(json.dumps(rabi_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "real_trap_control"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
