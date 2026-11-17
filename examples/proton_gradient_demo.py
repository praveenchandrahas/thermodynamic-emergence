"""Illustrative proton-potential profile; values are synthetic, not field data."""

from thermodynamic_emergence.chemistry import evaluate_proton_potential_profile

profile = evaluate_proton_potential_profile(
    positions_m=[0.0, 0.01, 0.02, 0.03],
    pH_profile=[5.0, 6.0, 7.0, 8.0],
    potential_profile_V=[0.0, 0.0, 0.0, 0.0],
    temperature_K=298.15,
    reference_index=0,
)

print("Synthetic pH-gradient demonstration")
print("Position (m) | delta-pH | delta-mu_H+ (kJ/mol)")
for point in profile:
    print(
        f"{point['position_m']:>12.3f} | "
        f"{point['delta_pH_from_reference']:>8.2f} | "
        f"{point['delta_mu_proton_from_reference_j_mol'] / 1000:>20.3f}"
    )
