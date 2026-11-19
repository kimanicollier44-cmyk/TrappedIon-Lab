"""Tests for TrappedIon-Lab: displacement matrix, sidebands, cooling,
Molmer-Sorensen gate, Rabi flopping, classifier, HOLD gates."""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import main as lab


# ---------------- Laguerre + displacement matrix ----------------

def test_laguerre_at_zero():
    for n in range(6):
        assert abs(lab.laguerre(n, 0, 0.0) - 1.0) < 1e-12
        assert abs(lab.laguerre(n, 1, 0.0) - (n + 1)) < 1e-12


def test_laguerre_known_values():
    # L_1(x) = 1 - x ; L_2(x) = 1 - 2x + x^2/2
    assert abs(lab.laguerre(1, 0, 0.4) - 0.6) < 1e-12
    assert abs(lab.laguerre(2, 0, 0.4) - (1 - 0.8 + 0.08)) < 1e-12


def test_disp_ground_carrier_debye_waller():
    for eta in (0.1, 0.3, 0.7):
        assert abs(abs(lab.disp(0, 0, eta))
                   - math.exp(-eta * eta / 2.0)) < 1e-12


def test_disp_ground_blue():
    assert abs(abs(lab.disp(1, 0, 0.25))
               - 0.25 * math.exp(-0.25 ** 2 / 2.0)) < 1e-12


def test_disp_unitarity_all_etas():
    for eta in lab.ETA_SET:
        for n in range(0, 13):
            row = sum(abs(lab.disp(m, n, eta)) ** 2 for m in range(41))
            assert abs(row - 1.0) < 1e-9


def test_disp_magnitude_symmetry():
    for n in range(6):
        for m in range(6):
            assert abs(abs(lab.disp(m, n, 0.4))
                       - abs(lab.disp(n, m, 0.4))) < 1e-12


def test_eta_from_labelled_params():
    eta, x0 = lab.eta_from_labelled_params()
    assert abs(eta - 0.096879) < 1e-5
    assert abs(x0 * 1e9 - 11.240) < 0.01


# ---------------- sideband spectrum ----------------

def test_red_sideband_vanishes_from_ground():
    for eta in lab.ETA_SET:
        assert lab.sideband_ratios(0, eta)[1] == 0.0


def test_sideband_ld_limit_small_eta():
    for n in range(0, 7):
        exact = lab.sideband_ratios(n, 0.05)
        approx = lab.ld_approx(n, 0.05)
        for e, a in zip(exact, approx):
            assert abs(e - a) < 2e-3


def test_sideband_ld_limit_breaks_at_large_eta():
    exact = lab.sideband_ratios(6, 0.5)
    approx = lab.ld_approx(6, 0.5)
    assert abs(exact[2] - approx[2]) > 0.5
    assert approx[0] < 0 < exact[0]  # LD carrier goes negative: absurd


def test_sideband_report_errors_grow():
    rep = lab.sideband_report()
    e05 = max(rep["lamb_dicke_max_abs_error_n0_6"]["0.05"].values())
    e50 = max(rep["lamb_dicke_max_abs_error_n0_6"]["0.5"].values())
    assert e05 < 2e-3 and e50 > 0.5


# ---------------- cooling ----------------

def test_thermal_distribution_normalised():
    p = lab.thermal_distribution(5.0, 20)
    assert abs(sum(p) - 1.0) < 1e-12
    assert abs(lab.mean_n(p) - 4.533386) < 1e-5


def test_cooling_reaches_eta_sq_floor():
    rep = lab.cooling_report()["red_sideband_cooling"]
    assert rep["final_nbar"] < 0.02
    assert abs(rep["final_nbar"] - 0.005) < 0.002
    assert rep["final_ground_population"] > 0.99


def test_cooling_beats_initial_by_orders():
    rep = lab.cooling_report()["red_sideband_cooling"]
    assert rep["final_nbar"] < rep["truncated_initial_nbar"] / 100.0


def test_carrier_only_cannot_cool():
    rep = lab.cooling_report()["carrier_only_contrast"]
    assert rep["final_nbar"] == rep["initial_nbar"]


# ---------------- gate ----------------

def test_concurrence_reference_states():
    bell = [[0j] * 4 for _ in range(4)]
    bell[0][0] = 0.5
    bell[0][3] = 0.5
    bell[3][0] = 0.5
    bell[3][3] = 0.5
    assert abs(lab.concurrence(bell) - 1.0) < 1e-9
    prod = [[0j] * 4 for _ in range(4)]
    prod[0][0] = 1.0
    assert lab.concurrence(prod) < 1e-9
    mix = [[0j] * 4 for _ in range(4)]
    mix[0][0] = 0.5
    mix[3][3] = 0.5
    assert lab.concurrence(mix) < 1e-9


def test_gate_truncated_space_honest_numbers():
    small = lab.gate_report()["phonon_space_0_to_3"]
    assert small["fidelity_to_bell_plus_after_frame_fix"] > 0.98
    assert small["concurrence"] > 0.97


def test_gate_reference_space_near_perfect():
    ref = lab.gate_report()["reference_space_0_to_7"]
    assert ref["fidelity_to_bell_plus_after_frame_fix"] > 0.9999
    assert ref["concurrence"] > 0.9999
    assert ref["spin_ground_population"] > 0.9999


def test_gate_truncation_cost_positive_and_printed():
    rep = lab.gate_report()
    assert rep["truncation_cost_fidelity"] > 0.005
    assert rep["truncation_cost_concurrence"] > 0.005


def test_gate_labelled_timing():
    rep = lab.gate_report()
    assert abs(rep["gate_time_us"] - 50.0) < 1e-9
    assert abs(rep["detuning_delta_over_2pi_hz"] - 20000.0) < 1e-6
    assert abs(rep["coupling_g_over_2pi_hz"] - 5000.0) < 1e-6


# ---------------- Rabi flopping ----------------

def test_rabi_pure_full_contrast():
    rep = lab.rabi_report()
    for row in rep["peaks_pure_vs_thermal"].values():
        assert abs(row["pure_n0"] - 1.0) < 1e-9


def test_rabi_thermal_collapse():
    rep = lab.rabi_report()
    assert 0.98 < rep["thermal_at_first_pi"] < 0.995
    assert rep["thermal_max_last_window"] < 0.85
    peaks = rep["peaks_pure_vs_thermal"]
    assert peaks["t=19*t_pi"]["thermal_nbar5"] < 0.65


def test_rabi_carrier_ratios_shrink():
    ratios = lab.rabi_report()["carrier_ratios_n0_to_5"]
    assert all(ratios[i] > ratios[i + 1] for i in range(5))
    assert abs(ratios[0] - 0.995012) < 1e-6


# ---------------- classifier ----------------

def test_classifier_scenarios():
    rep = lab.classifier_report()
    assert rep["ca_plus_labelled_headline"] == "lamb-dicke-resolved"
    assert rep["hot_ion_before_cooling"] == "lamb-dicke-resolved"
    assert rep["strong_drive_sideband_blur"] == "unresolved-sideband"
    assert rep["warm_marginal_ion"] == "lamb-dicke-marginal"
    assert rep["heavy_hot_ion"] == "beyond-lamb-dicke"
    assert rep["deep_lamb_dicke_cold"] == "lamb-dicke-resolved"


def test_classifier_matches_policy_on_grid():
    for features, label in lab.TRAINING_GRID:
        assert lab.classify_regime(*features) == label


def test_policy_label_boundaries():
    assert lab.policy_label(1.0, 0.1, 0.0) == "beyond-lamb-dicke"
    assert lab.policy_label(0.1, 0.6, 0.0) == "unresolved-sideband"
    assert lab.policy_label(0.1, 0.1, 0.0) == "lamb-dicke-resolved"


# ---------------- gates + seats + CLI ----------------

def test_hold_gates_14_all_hold():
    assert len(lab.HOLD_GATES) == 14
    assert all(v.startswith("HOLD") for v in lab.HOLD_GATES.values())
    assert lab.gate("real_trap_control") == "HOLD"
    assert lab.gate("laser_system_operation") == "HOLD"
    assert lab.gate("gate_fidelity_certification") == "HOLD"
    assert lab.gate("run_simulation") == "PROCEED"


def test_hold_demo_exits_2():
    assert lab.main(["hold-demo"]) == 2
    assert lab.main(["hold-demo", "quantum_advantage_claim"]) == 2
    assert lab.main(["hold-demo", "run_simulation"]) == 0


def test_seats_never_override():
    brains = lab.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["seat"] == "Gemini Argon 4 Extended Thinking"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = lab.run_checks()
    assert all(ok for _, ok in checks)
    assert len(checks) == 52
