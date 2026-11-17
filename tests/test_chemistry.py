import math
import unittest

from thermodynamic_emergence.chemistry import (
    Reaction,
    chemical_affinity_j_mol,
    chemical_power_density_w_m3,
    evaluate_proton_potential_profile,
    evaluate_reaction_profile,
    proton_transfer_delta_mu,
    reaction_gibbs_energy,
    redox_reaction_delta_g,
)


class LocalChemicalDisequilibriumTests(unittest.TestCase):
    def setUp(self):
        self.reaction = Reaction("A to B", {"A": -1, "B": 1})

    def test_reaction_gibbs_energy_uses_activity_quotient(self):
        temperature = 298.15
        delta_g0 = 5_000.0
        result = reaction_gibbs_energy(
            self.reaction,
            delta_g0_j_mol=delta_g0,
            temperature_K=temperature,
            activities={"A": 0.1, "B": 0.2},
        )
        expected = delta_g0 + 8.31446261815324 * temperature * math.log(2.0)
        self.assertAlmostEqual(result, expected, places=10)
        self.assertAlmostEqual(chemical_affinity_j_mol(result), -result, places=10)

    def test_equilibrium_quotient_gives_zero_reaction_gibbs_energy(self):
        temperature = 310.0
        delta_g0 = -12_000.0
        q_equilibrium = math.exp(-delta_g0 / (8.31446261815324 * temperature))
        result = reaction_gibbs_energy(
            self.reaction,
            delta_g0_j_mol=delta_g0,
            temperature_K=temperature,
            activities={"A": 1.0, "B": q_equilibrium},
        )
        self.assertAlmostEqual(result, 0.0, places=9)

    def test_three_unit_pH_gradient_has_expected_chemical_component(self):
        result = proton_transfer_delta_mu(
            source_pH=5.0,
            destination_pH=8.0,
            temperature_K=298.15,
        )
        self.assertAlmostEqual(result, -17_124.0, delta=5.0)

    def test_electrical_potential_is_included_in_proton_transfer(self):
        no_voltage = proton_transfer_delta_mu(
            source_pH=5.0, destination_pH=8.0, temperature_K=298.15
        )
        with_voltage = proton_transfer_delta_mu(
            source_pH=5.0,
            destination_pH=8.0,
            temperature_K=298.15,
            destination_minus_source_potential_V=0.1,
        )
        self.assertAlmostEqual(with_voltage - no_voltage, 9648.533212, places=5)

    def test_redox_delta_g_uses_acceptor_minus_donor_potential(self):
        result = redox_reaction_delta_g(
            electrons_transferred=2,
            acceptor_reduction_potential_V=0.2,
            donor_reduction_potential_V=-0.1,
        )
        self.assertAlmostEqual(result, -57_891.199272, places=5)

    def test_chemical_power_density_is_rate_times_affinity(self):
        result = chemical_power_density_w_m3(
            net_reaction_rate_mol_m3_s=1e-4,
            delta_g_j_mol=-10_000.0,
        )
        self.assertAlmostEqual(result, 1.0, places=12)

    def test_spatial_reaction_profile(self):
        result = evaluate_reaction_profile(
            self.reaction,
            positions_m=[0.0, 0.5, 1.0],
            temperatures_K=[280.0, 290.0, 300.0],
            delta_g0_j_mol=[-10_000.0, -8_000.0, -6_000.0],
            activities_by_species={"A": [1.0, 1.0, 1.0], "B": [0.1, 1.0, 10.0]},
        )
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0]["direction"], "forward_exergonic")
        self.assertEqual(result[2]["direction"], "reverse_exergonic")

    def test_proton_potential_profile_reference_is_zero(self):
        result = evaluate_proton_potential_profile(
            positions_m=[0, 1, 2],
            pH_profile=[5, 6, 8],
            potential_profile_V=[0.0, 0.0, 0.0],
            temperature_K=298.15,
        )
        self.assertAlmostEqual(result[0]["delta_mu_proton_from_reference_j_mol"], 0.0)
        self.assertLess(result[2]["delta_mu_proton_from_reference_j_mol"], 0.0)

    def test_missing_species_activity_is_rejected(self):
        with self.assertRaises(KeyError):
            reaction_gibbs_energy(
                self.reaction,
                delta_g0_j_mol=0.0,
                temperature_K=298.15,
                activities={"A": 1.0},
            )

    def test_zero_activity_is_rejected(self):
        with self.assertRaises(ValueError):
            reaction_gibbs_energy(
                self.reaction,
                delta_g0_j_mol=0.0,
                temperature_K=298.15,
                activities={"A": 1.0, "B": 0.0},
            )

    def test_profile_requires_monotonic_positions(self):
        with self.assertRaises(ValueError):
            evaluate_proton_potential_profile(
                positions_m=[0, 1, 0.5],
                pH_profile=[5, 6, 7],
                potential_profile_V=[0, 0, 0],
                temperature_K=298.15,
            )


if __name__ == "__main__":
    unittest.main()
