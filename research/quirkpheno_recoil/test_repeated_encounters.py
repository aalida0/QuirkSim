from dataclasses import replace
import math
import unittest
from quirksim.repeated_encounters import EncounterPrescription,evolve,pericenter_gev_inverse
from quirksim.radiation_annihilation import OrbitAveragedModel,linear_cooling_time_ps


class RepeatedEncounterTest(unittest.TestCase):
    def test_turning_point_satisfies_hamiltonian(self):
        m=225.;e=40.;ell=1e7;lam=774.
        r=pericenter_gev_inverse(m,e,ell,lam)
        h=2*(math.hypot(m,ell/r)-m)+(lam*1e-9)**2*r
        self.assertAlmostEqual(h/e,1.,places=11)
        self.assertEqual(pericenter_gev_inverse(m,e,0.,lam),0.)
        with self.assertRaises(ValueError):pericenter_gev_inverse(m,1e-6,1e12,lam)

    def test_zero_L_cooling_reproduces_existing_integral_and_ledgers(self):
        model=EncounterPrescription(captured_singlet_fraction=.25)
        result=evolve(40.,0.,1e13,model)
        clock=OrbitAveragedModel(225.,774.,.1,774e-9,0.,0.,1.,'control',charge_magnitude=0.,include_photons=False)
        expected=linear_cooling_time_ps(clock,40.,8*774e-9)
        self.assertEqual(result['status'],'CAPTURED')
        self.assertAlmostEqual(result['capture_time_ps']/expected,1.,places=10)
        self.assertAlmostEqual(result['annihilation_probability_in_declared_model'],.25)
        self.assertLess(abs(result['energy_ledger_residual_gev']),1e-10)
        self.assertEqual(result['angular_ledger_residual_hbar'],0.)

    def test_gate_closed_preserves_actions_over_repeated_passes(self):
        result=evolve(40.,1e12,1e5,EncounterPrescription(angular_loss_per_effective_emission=1.,annihilation_optical_depth_per_eligible_encounter=.01))
        self.assertEqual(result['status'],'EMISSION_GATE_CLOSED')
        self.assertGreater(result['radial_encounters_in_secular_clock'],100.)
        self.assertEqual(result['eligible_encounters_in_secular_clock'],0.)
        self.assertEqual(result['final_L'],1e12)
        self.assertEqual(result['radiated_cooling_energy_gev'],0.)
        self.assertEqual(result['annihilation_probability_in_declared_model'],0.)

    def test_angular_loss_changes_capture_outcome(self):
        conserved=evolve(40.,1e7,1e13,EncounterPrescription())
        drained=evolve(40.,1e7,1e13,EncounterPrescription(angular_loss_per_effective_emission=1.))
        self.assertEqual(conserved['status'],'EMISSION_GATE_CLOSED')
        self.assertEqual(drained['status'],'CAPTURED')
        self.assertAlmostEqual(drained['orbital_L_removed_in_cooling'],1e7,places=5)

    def test_no_radiation_exact_competing_exponential(self):
        model=EncounterPrescription(epsilon_per_full_period=0.,annihilation_optical_depth_per_eligible_encounter=.01)
        r=evolve(40.,0.,1000.,model)
        self.assertEqual(r['status'],'NO_RADIATION')
        self.assertAlmostEqual(r['annihilation_probability_in_declared_model'], -math.expm1(-.01*r['eligible_encounters_in_secular_clock']),places=14)

    def test_finite_horizon_and_capture_spin_prescription(self):
        model=EncounterPrescription(captured_singlet_fraction=0.)
        r=evolve(40.,0.,1000.,model)
        self.assertEqual(r['status'],'HORIZON_REACHED');self.assertEqual(r['capture_probability_in_declared_model'],0.)
        self.assertLess(r['final_energy_gev'],40.)
        r=evolve(40.,0.,1e13,model)
        self.assertEqual(r['capture_probability_in_declared_model'],1.)
        self.assertEqual(r['annihilation_probability_in_declared_model'],0.)
        self.assertEqual(r['captured_other_spin_sector_probability'],1.)

    def test_step_refinement(self):
        a=evolve(40.,1e7,1e13,EncounterPrescription())
        b=evolve(40.,1e7,1e13,EncounterPrescription(energy_step_fraction=.01))
        self.assertAlmostEqual(a['final_energy_gev']/b['final_energy_gev'],1.,places=10)
        self.assertAlmostEqual(a['time_ps']/b['time_ps'],1.,places=10)


if __name__=='__main__':unittest.main()
