import math
import numpy as np
import pytest
from scipy.integrate import quad
from quirksim.radiation_annihilation import (
    ALPHA_ZERO, lienard_power_gev_per_ps, dirac_singlet_partial_sigma_v,
    wkb_crossing_exponent, OrbitAveragedModel, evolve_orbit_average,
    linear_cooling_time_ps, radiation_recoil_scale, coulomb_ground_reference,
    evolve_quantum_levels)
from quirksim.physics import HBAR_GEV_PS


def test_lienard_matches_angular_larmor_integral_and_longitudinal_force_limit():
    # Integrate dP/dOmega for a nonrelativistic instantaneous rest frame.
    a = .004
    angular = quad(lambda u: ALPHA_ZERO*HBAR_GEV_PS*a*a/(4*math.pi)
                   *(1-u*u)*2*math.pi, -1, 1)[0]
    assert lienard_power_gev_per_ps([0, 0, 0], [a, 0, 0]) == pytest.approx(angular, rel=1e-12)
    for beta in (0., .4, .8):
        gamma = 1/math.sqrt(1-beta*beta)
        assert lienard_power_gev_per_ps([beta, 0, 0], [a/gamma**3, 0, 0]) == pytest.approx(angular, rel=1e-12)


def test_threshold_photon_kernel_reproduces_parapositronium_width():
    m, alpha = .00051099895, ALPHA_ZERO
    # Remove IC multiplicity from the NIC=2 photon coefficient; project the
    # unpolarized hard cross section onto its spin-singlet channel (factor 4).
    coefficients = dirac_singlet_partial_sigma_v(m, 2, .01)
    sigma_para = 4*coefficients['sigma_v_gev_minus2']['gamma_gamma']/2
    psi2 = coulomb_ground_reference(m, m, alpha)['psi_squared_gev3']
    assert psi2*sigma_para == pytest.approx(m*alpha**5/2, rel=1e-13)


def test_ic_color_trace_factor_and_partial_coverage():
    for n in (2, 3, 4):
        generators = []
        for i in range(n):
            for j in range(i+1, n):
                sym, antisym = np.zeros((n,n),complex), np.zeros((n,n),complex)
                sym[i,j] = sym[j,i] = .5
                antisym[i,j], antisym[j,i] = -.5j, .5j
                generators += [sym, antisym]
        for kdiag in range(1, n):
            diagonal = np.zeros(n)
            diagonal[:kdiag] = 1.
            diagonal[kdiag] = -kdiag
            generators.append(np.diag(diagonal)/math.sqrt(2*kdiag*(kdiag+1)))
        gluon_amplitude_squared = sum(abs(np.trace(a@b+b@a)/math.sqrt(n))**2
                                     for a in generators for b in generators)
        photon_amplitude_squared = (2*math.sqrt(n))**2
        k = dirac_singlet_partial_sigma_v(225., n, .04, alpha_em=.04)
        coeff = k['sigma_v_gev_minus2']
        assert coeff['ic_gluons']/coeff['gamma_gamma'] == pytest.approx(
            gluon_amplitude_squared/photon_amplitude_squared)
        assert k['coverage'] == 'PARTIAL_QED_AND_IC_ONLY'


def test_zero_radiation_competing_risks_match_exact_exponential():
    model = OrbitAveragedModel(225., 774., 0., 774e-9, 1e-7, .4, 258.5,
                               'test coefficients', include_photons=False)
    r = model.rates(40.)
    result = evolve_orbit_average(model, 40., 1000., .001)
    a, w = r['annihilation_rate_per_ps'], r['first_weak_rate_per_ps']
    survival = math.exp(-(a+w)*1000.)
    assert result['survival_before_first_weak_or_annihilation'] == pytest.approx(survival, abs=1e-11)
    assert result['annihilation_before_first_weak_probability'] == pytest.approx(a/(a+w)*(1-survival), abs=1e-11)
    assert result['conditional_alive_trajectory_energy_loss_gev'] == 0
    assert abs(result['probability_residual']) < 1e-12


def test_zero_overlap_never_infers_annihilation_from_radial_orbit():
    model = OrbitAveragedModel(225., 774., .1, 774e-9, 1., 0., 258.5, 'test')
    result = evolve_orbit_average(model, 40., 1000., .001)
    assert result['annihilation_before_first_weak_probability'] == 0
    assert result['physical_prediction_status'] == 'UNKNOWN'
    assert result['conditional_alive_trajectory_energy_loss_gev'] > 0


def test_cooling_feedback_matches_independent_energy_quadrature():
    # Deliberately accelerated diagnostic parameters make the matching event
    # reachable quickly. No benchmark physics is claimed by this test.
    model = OrbitAveragedModel(225., 1e6, 1., .1, 0., 0., 1e30,
                               'diagnostic', include_photons=False)
    cooling = linear_cooling_time_ps(model, .2, .01)
    result = evolve_orbit_average(model, .2, 2*cooling, .01, max_step_ps=cooling/30)
    assert result['termination'] == 'LINEAR_TO_QUANTUM_MATCHING'
    assert result['time_ps'] == pytest.approx(cooling, rel=1e-8)
    assert result['conditional_alive_trajectory_energy_loss_gev'] == pytest.approx(.19, abs=1e-10)


def test_two_stage_quantum_radiation_then_annihilation_exact_chain_and_energy():
    t = 4.
    r = evolve_quantum_levels([451., 450.], [[0., 2.], [0., 0.]],
                             [0., 1.], [0., 0.], [1., 0.], t, 'analytic test chain')
    p0 = math.exp(-2*t)
    p1 = 2*(math.exp(-t)-math.exp(-2*t))
    np.testing.assert_allclose(r['level_probabilities'], [p0, p1], atol=1e-13)
    assert r['annihilation_probability'] == pytest.approx(1-p0-p1, abs=1e-13)
    assert r['expected_radiated_energy_gev'] == pytest.approx(1-p0, abs=1e-13)
    assert abs(r['probability_residual']) < 1e-13
    assert abs(r['energy_residual_gev']) < 1e-10


def test_quantum_weak_competition_and_no_upward_spontaneous_transition():
    r = evolve_quantum_levels([450.], [[0.]], [1.], [2.], [1.], 10., 'test')
    assert r['annihilation_probability'] == pytest.approx((1-math.exp(-30))/3)
    assert r['first_weak_probability'] == pytest.approx(2*(1-math.exp(-30))/3)
    with pytest.raises(ValueError):
        evolve_quantum_levels([450., 451.], [[0., 1.], [0., 0.]], [1., 1.], [0., 0.], [1., 0.], 1., 'test')


def test_radiation_small_energy_can_still_have_a_large_short_distance_recoil_scale():
    r = radiation_recoil_scale(4*774e-9, 225., 1000., 774.)
    assert r['scale_over_confinement_reference'] > 10
    assert r['worldsheet_geometry_error_bound'] is None
    assert wkb_crossing_exponent(225., .4, 1e-7, 0.) == 0
