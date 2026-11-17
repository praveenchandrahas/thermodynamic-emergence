"""Local chemical-disequilibrium calculations.

This module evaluates thermodynamic driving forces from supplied, dimensionless
species activities and environmental profiles. It does not calculate aqueous
speciation, reaction kinetics, transport, or abiogenesis.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence

GAS_CONSTANT_J_MOL_K = 8.31446261815324
FARADAY_CONSTANT_C_MOL = 96485.33212


@dataclass(frozen=True)
class Reaction:
    """A written chemical reaction represented by stoichiometric numbers.

    Negative coefficients are reactants; positive coefficients are products.
    The reaction must be atom- and charge-balanced by the caller. This module
    cannot verify elemental balance without a chemical formula/speciation
    database.
    """
    name: str
    stoichiometry: Mapping[str, float]

    def validate(self) -> None:
        if not self.name.strip():
            raise ValueError("reaction name cannot be empty")
        if len(self.stoichiometry) < 2:
            raise ValueError("a reaction needs at least two species")
        values = []
        for species, coefficient in self.stoichiometry.items():
            if not species.strip():
                raise ValueError("species names cannot be empty")
            if not math.isfinite(coefficient) or coefficient == 0:
                raise ValueError("stoichiometric coefficients must be finite and non-zero")
            values.append(coefficient)
        if not any(v < 0 for v in values) or not any(v > 0 for v in values):
            raise ValueError("reaction must contain reactants and products")


def _positive_finite(value: float, name: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def reaction_log_quotient(
    reaction: Reaction,
    activities: Mapping[str, float],
) -> float:
    """Return ln(Q), with Q = product(a_i ** nu_i).

    Activities are dimensionless relative to their specified standard states.
    Concentrations may only substitute for activities under a stated ideal-
    solution approximation and consistent standard-state convention.
    """
    reaction.validate()
    log_q = 0.0
    for species, coefficient in reaction.stoichiometry.items():
        if species not in activities:
            raise KeyError(f"missing activity for species {species!r}")
        activity = activities[species]
        _positive_finite(activity, f"activity[{species}]")
        log_q += coefficient * math.log(activity)
    return log_q


def reaction_gibbs_energy(
    reaction: Reaction,
    *,
    delta_g0_j_mol: float,
    temperature_K: float,
    activities: Mapping[str, float],
) -> float:
    """Calculate reaction Gibbs energy in J/mol: ΔrG = ΔrG°(T,P) + RT ln Q.

    The supplied standard reaction Gibbs energy must correspond to the local
    temperature, pressure, standard states and chosen chemical species.
    """
    _positive_finite(temperature_K, "temperature_K")
    if not math.isfinite(delta_g0_j_mol):
        raise ValueError("delta_g0_j_mol must be finite")
    log_q = reaction_log_quotient(reaction, activities)
    return delta_g0_j_mol + GAS_CONSTANT_J_MOL_K * temperature_K * log_q


def chemical_affinity_j_mol(delta_g_j_mol: float) -> float:
    """Return chemical affinity A = -ΔrG in J/mol."""
    if not math.isfinite(delta_g_j_mol):
        raise ValueError("delta_g_j_mol must be finite")
    return -delta_g_j_mol


def proton_transfer_delta_mu(
    *,
    source_pH: float,
    destination_pH: float,
    temperature_K: float,
    destination_minus_source_potential_V: float = 0.0,
) -> float:
    """Electrochemical potential change for one mole of H+ moving source→destination.

    Returns mu_destination - mu_source in J/mol:
      delta_mu_H+ = -2.303 R T (pH_destination - pH_source) + F delta_phi

    A negative result means this transfer is downhill in electrochemical
    potential. The pH term alone is not proof that a natural membrane exists
    or that a prebiotic molecular mechanism can harvest the gradient.
    """
    _positive_finite(temperature_K, "temperature_K")
    if not math.isfinite(source_pH) or not math.isfinite(destination_pH):
        raise ValueError("pH values must be finite")
    if not math.isfinite(destination_minus_source_potential_V):
        raise ValueError("potential difference must be finite")
    chemical = (
        -math.log(10.0)
        * GAS_CONSTANT_J_MOL_K
        * temperature_K
        * (destination_pH - source_pH)
    )
    electrical = FARADAY_CONSTANT_C_MOL * destination_minus_source_potential_V
    return chemical + electrical


def redox_reaction_delta_g(
    *,
    electrons_transferred: int,
    acceptor_reduction_potential_V: float,
    donor_reduction_potential_V: float,
) -> float:
    """Return delta_G = -n F delta_E for donor→acceptor electron transfer.

    Potentials must use the same reference electrode and must already be
    corrected (e.g. by Nernst/speciation calculations) for the stated local
    temperature, pressure, pH and activities. delta_E is E_acceptor - E_donor.
    """
    if isinstance(electrons_transferred, bool) or not isinstance(electrons_transferred, int):
        raise ValueError("electrons_transferred must be an integer")
    if electrons_transferred < 1:
        raise ValueError("electrons_transferred must be a positive integer")
    if not math.isfinite(acceptor_reduction_potential_V):
        raise ValueError("acceptor_reduction_potential_V must be finite")
    if not math.isfinite(donor_reduction_potential_V):
        raise ValueError("donor_reduction_potential_V must be finite")
    delta_e = acceptor_reduction_potential_V - donor_reduction_potential_V
    return -electrons_transferred * FARADAY_CONSTANT_C_MOL * delta_e


def chemical_power_density_w_m3(
    *,
    net_reaction_rate_mol_m3_s: float,
    delta_g_j_mol: float,
) -> float:
    """Return -r delta_G in W/m3 for a signed net rate in the written direction.

    A positive value indicates downhill Gibbs-energy release for the specified
    net reaction flux. A negative value signals that the imposed net flux is
    uphill; it may require coupling to other processes and must not be counted
    as positive dissipative power on its own.
    """
    if not math.isfinite(net_reaction_rate_mol_m3_s):
        raise ValueError("net_reaction_rate_mol_m3_s must be finite")
    if not math.isfinite(delta_g_j_mol):
        raise ValueError("delta_g_j_mol must be finite")
    return -net_reaction_rate_mol_m3_s * delta_g_j_mol


def evaluate_reaction_profile(
    reaction: Reaction,
    *,
    positions_m: Sequence[float],
    temperatures_K: Sequence[float],
    delta_g0_j_mol: Sequence[float],
    activities_by_species: Mapping[str, Sequence[float]],
) -> list[dict[str, float | str]]:
    """Evaluate reaction thermodynamic driving force along a supplied 1-D profile.

    Every activity and standard reaction Gibbs-energy value is supplied by
    the caller. In a real geochemical system these should come from
    measurements or a speciation/thermodynamic database, not from linear
    interpolation of pH or redox potential.
    """
    reaction.validate()
    n = len(positions_m)
    if n == 0:
        raise ValueError("profile cannot be empty")
    if len(temperatures_K) != n or len(delta_g0_j_mol) != n:
        raise ValueError("positions, temperatures and delta_g0 profiles must have equal length")
    if any(not math.isfinite(x) for x in positions_m):
        raise ValueError("positions must be finite")
    if any(positions_m[i] >= positions_m[i + 1] for i in range(n - 1)):
        raise ValueError("positions must be strictly increasing")
    for species in reaction.stoichiometry:
        if species not in activities_by_species:
            raise KeyError(f"missing activity profile for species {species!r}")
        if len(activities_by_species[species]) != n:
            raise ValueError(f"activity profile for {species!r} has incorrect length")

    output: list[dict[str, float | str]] = []
    for i, position in enumerate(positions_m):
        activities = {
            species: activities_by_species[species][i]
            for species in reaction.stoichiometry
        }
        delta_g = reaction_gibbs_energy(
            reaction,
            delta_g0_j_mol=delta_g0_j_mol[i],
            temperature_K=temperatures_K[i],
            activities=activities,
        )
        output.append({
            "position_m": float(position),
            "temperature_K": float(temperatures_K[i]),
            "delta_g_j_mol": delta_g,
            "affinity_j_mol": -delta_g,
            "ln_reaction_quotient": reaction_log_quotient(reaction, activities),
            "direction": (
                "forward_exergonic" if delta_g < 0
                else "reverse_exergonic" if delta_g > 0
                else "equilibrium"
            ),
        })
    return output


def evaluate_proton_potential_profile(
    *,
    positions_m: Sequence[float],
    pH_profile: Sequence[float],
    potential_profile_V: Sequence[float],
    temperature_K: float,
    reference_index: int = 0,
) -> list[dict[str, float]]:
    """Map H+ electrochemical potential differences from one reference location.

    The electric potential and pH profiles must be independently supplied.
    This is a thermodynamic diagnostic; it does not simulate proton transport
    or demonstrate that the environment contains a selective barrier.
    """
    n = len(positions_m)
    if n == 0 or len(pH_profile) != n or len(potential_profile_V) != n:
        raise ValueError("positions, pH and potential profiles must have equal non-zero length")
    if any(not math.isfinite(x) for x in positions_m):
        raise ValueError("positions must be finite")
    if any(positions_m[i] >= positions_m[i + 1] for i in range(n - 1)):
        raise ValueError("positions must be strictly increasing")
    if isinstance(reference_index, bool) or not isinstance(reference_index, int):
        raise ValueError("reference_index must be an integer")
    if not 0 <= reference_index < n:
        raise ValueError("reference_index is outside the profile")
    source_pH = pH_profile[reference_index]
    source_phi = potential_profile_V[reference_index]
    output = []
    for x, pH, phi in zip(positions_m, pH_profile, potential_profile_V):
        delta_mu = proton_transfer_delta_mu(
            source_pH=source_pH,
            destination_pH=pH,
            temperature_K=temperature_K,
            destination_minus_source_potential_V=phi - source_phi,
        )
        output.append({
            "position_m": float(x),
            "delta_pH_from_reference": float(pH - source_pH),
            "delta_phi_from_reference_V": float(phi - source_phi),
            "delta_mu_proton_from_reference_j_mol": delta_mu,
        })
    return output
