# 02 — Local chemical disequilibria
## From planetary energy flows to geochemical opportunities for organization

**Status:** Research framework, version 0.1  
**Relationship to the first model:** This note defines the next scale down from the global stellar–planetary energy ledger.

---

## 1. Why the global energy budget is not enough

The first model establishes the planetary-scale radiative setting: a star supplies radiation, the planet absorbs and transforms part of it, and the planet emits energy back to space. That global balance is necessary context, but it does not tell us which local environments can support chemical organization.

A global-mean value such as absorbed solar power per square metre does not identify:
- Which chemical species coexist locally.
- Which reactants and products are far from chemical equilibrium.
- Which reactions have a negative Gibbs energy under the local temperature, pressure and composition.
- Whether a mineral surface, pore, membrane, phase boundary or other compartment can bring reactants together.
- Whether reaction rates are fast enough to compete with diffusion, mixing, dilution and decomposition.
- Whether the gradient is replenished faster than it is consumed or dissipated.

The relevant next step is therefore not to distribute the planetary energy budget over a grid and call each cell a chemical environment. We need local state variables and explicit geochemical boundary conditions.

**The central bridge is from energy flux at planetary scale to chemical-potential differences at local scale, and then from thermodynamic availability to actual reaction rates and maintained throughput.**

## 2. The fundamental thermodynamic quantity is chemical potential

The chemical potential of species i is its partial molar Gibbs energy:

[
mu_i = left(rac{partial G}{partial n_i}ight)_{T,P,n_{j
e i}}.
]

For a species in a specified standard state, its chemical potential can be written

[
mu_i(T,P,a_i)=mu_i^circ(T,P)+RTln a_i,
]

where a_i is its dimensionless activity, R is the gas constant and T is absolute temperature. Activity is not generally identical to concentration. In dilute ideal solutions, concentration-based approximations may be acceptable; in saline, mineral-rich or otherwise non-ideal fluids, activity coefficients and aqueous speciation matter [1].

For a reaction written with stoichiometric coefficients ν_i (negative for reactants and positive for products),

[
Delta_r G = sum_i 
u_imu_i
= Delta_r G^circ(T,P)+RTln Q,
]

[
Q=prod_i a_i^{
u_i}.
]

The reaction affinity is

[
mathcal{A}=-Delta_r G.
]

With this sign convention:
- ΔrG<0: the written forward direction is thermodynamically downhill.
- ΔrG>0: the reverse direction is downhill under the specified conditions.
- ΔrG=0: the reaction is at equilibrium.

The value of ΔrG° is not a universal constant independent of environment. It depends on temperature, pressure, standard-state conventions and the species included in the reaction. A serious geochemical model must calculate aqueous speciation and activities consistently before evaluating Q.

This equation turns the broad idea of chemical disequilibrium into something calculable: the environment supports a reaction opportunity only relative to a specified reaction, species set and local state.

## 3. Different gradients are related, but they are not interchangeable

### 3.1 Redox disequilibrium

Redox reactions transfer electrons between donor and acceptor couples. When reduction potentials are expressed on the same reference scale and corrected for the local temperature, pH and activities, the Gibbs energy for electron transfer can be written

[
Delta G=-nFDelta E_{mathrm{cell}},
]

where n is the number of electrons transferred, F is Faraday's constant and ΔE_cell=E_acceptor−E_donor.

A more positive acceptor reduction potential than donor reduction potential gives a positive cell potential and a negative Gibbs energy for the coupled electron-transfer reaction.

Redox potential is not itself a concentration or a conserved field. It is a conditional property of a redox couple and its composition. For real geochemical systems it must be calculated from the couples and species present, with the appropriate Nernst corrections.

### 3.2 pH and proton electrochemical disequilibrium

The chemical contribution to the electrochemical potential of a proton depends on proton activity. If a proton moves from a source compartment to a destination compartment, then

[
Deltawidetilde{mu}_{H^+}
=
widetilde{mu}_{mathrm{destination}}
-
widetilde{mu}_{mathrm{source}}
]

[
=
-2.303RT(pH_{mathrm{destination}}-pH_{mathrm{source}})
+F(phi_{mathrm{destination}}-phi_{mathrm{source}}).
]

The first term is the chemical contribution and the second is the electrical contribution. A calculation using only ΔpH is incomplete when an electrical potential difference is present.

At 298.15 K, a three-unit pH difference corresponds to about 17.1 kJ per mole of protons in the chemical component alone. A 100 mV electrical potential difference contributes about 9.65 kJ per mole of monovalent charge. These are illustrative thermodynamic magnitudes, not evidence that a natural environment has a membrane capable of preserving or harvesting the gradient.

That distinction matters. A pH contrast between two fluids is not automatically a proton-motive force available to a molecular machine. A usable trans-compartment gradient needs a physical separation with suitable permeability and a mechanism that couples proton transfer to work.

### 3.3 Thermal disequilibrium

A temperature gradient can drive heat flow. For conductive transport, Fourier's law is

[
mathbf{q}=-k
abla T,
]

where q is heat flux and k is thermal conductivity. Thermal gradients can influence reaction rates, equilibria, convection and thermophoresis. They are not automatically equivalent to a chemical reaction's Gibbs energy. The thermal driving force and its coupling to chemical transformations must be represented explicitly.

### 3.4 Concentration, water-activity and wet–dry cycling

Evaporation, freezing, adsorption and phase changes can concentrate reactants, change ionic strength and water activity, and alter reaction equilibria. Dehydration can favour some condensation reactions that are unfavourable in bulk water. Rewetting can redistribute products and reactants.

These processes may create non-equilibrium sequences rather than one steady gradient. The energy source, water flux, material transport and cycle timing must all be included. A wet–dry cycle is not free energy without an account of the evaporation, heating, cooling and transport that maintain it.

### 3.5 Photochemical disequilibrium

Ultraviolet and visible photons can drive reactions that are not spontaneous under local dark equilibrium conditions. Photochemical models need wavelength-resolved photon flux, absorption cross-sections, quantum yields and product pathways. The stellar radiation ledger can eventually constrain the available photon supply, but it cannot by itself establish the chemistry or yield of a local photochemical network.

## 4. Candidate geochemical settings — a comparison, not a selection

Several environments can create combinations of chemical, thermal and physical non-equilibrium. None should be declared the origin-of-life setting by assumption.

| Environment | Potentially relevant non-equilibria | Important constraints |
|---|---|---|
| Deep-sea hydrothermal mixing interfaces | Reduced hydrothermal fluids mixing with ocean water; redox, temperature, pH and mineral gradients | High temperature can damage products; dilution, flow and mineral precipitation alter gradients |
| Serpentinizing systems | Water–rock reactions can generate H₂, alkaline fluids, reduced minerals and carbon-containing products | The early Earth's geology, ocean chemistry and in-situ pH are uncertain; proposed proton-gradient harvesting mechanisms are contested |
| Shallow alkaline hydrothermal fields | Temperature variation, wet–dry cycling, seawater/freshwater mixing and minerals | Conditions vary in space and time; present-day analogues are not exact early-Earth replicas |
| Subaerial hot springs and volcanic pools | Evaporation, concentration, UV exposure, thermal cycling and mineral surfaces | Reactants must be supplied; salinity, hydrolysis and product degradation matter |
| Mineral pores and sediment interfaces | Adsorption, confinement, redox and pH microgradients, concentration and retention | Pore geometry, diffusion, surface chemistry and mineral identity determine actual behaviour |
| Ice and freeze–thaw compartments | Solute concentration in brines, changing water activity and temperature | Low temperatures slow many reaction rates; freeze concentration and phase behaviour need explicit treatment |
| Atmosphere–ocean interfaces | Photochemistry, gas exchange and multi-phase chemical disequilibrium | Gas-only estimates can miss aqueous chemistry; abiotic and biological sources can produce disequilibria |

The literature supports investigating a range of settings. Ianeselli and colleagues review heated rock pores, hydrothermal vents, evaporating ponds, freeze–thaw compartments, UV irradiation and high-CO₂ conditions as possible non-equilibrium contexts for prebiotic nucleic-acid chemistry [2]. Barge and Price argue that shallow alkaline hydrothermal systems should be considered because they combine wet–dry cycles, temperature variation and changing freshwater/saltwater inputs [3]. A geological review by Rodriguez and colleagues discusses multiple early-Earth environments and the uncertain geochemical transitions between them [4].

Serpentinization is a particularly interesting candidate because water–rock reactions can generate hydrogen and strongly reducing fluids, alongside mineral and pH contrasts [5,6]. But that does not settle the origin question. Jackson's review challenges whether natural alkaline-vent pH gradients could have been held across the inorganic membranes required by a proton-harvesting scenario [7]. More recently, Tutolo's 2026 PNAS article argues that the geological and chemical conditions of ancient serpentinizing systems may differ significantly from modern Lost City analogues, including whether strong pH gradients existed at relevant hydrothermal temperatures [8]. This is a critique of a specific vent-based mechanism, not a refutation of all hydrothermal or geochemical origin scenarios.

The correct modelling response is to represent competing environmental scenarios and uncertain parameter ranges, not to encode one popular hypothesis as the default reality.

## 5. Thermodynamic availability is not kinetic accessibility

A negative reaction Gibbs energy means a specified forward reaction can lower Gibbs energy under the stated conditions. It does not tell us:
- How quickly the reaction proceeds.
- Whether an activation barrier makes the rate negligible.
- Whether the required catalyst or mineral surface exists.
- Whether products are stable long enough to participate in further reactions.
- Whether reactants remain in contact at useful concentrations.
- Whether the environment is replenished before the gradient is exhausted.

If a net reaction rate r is specified in mol m⁻³ s⁻¹, a local rate of Gibbs-energy release can be written

[
p_{mathrm{chem}}=-rDelta_r G,
]

with units W m⁻³. For a network of reactions,

[
p_{mathrm{chem}}=-sum_r r_rDelta_rG_r.
]

This is a chemical power density only when the rates are defined as net reaction-extent rates in the written directions and the Gibbs energies use the same local conditions. An uphill step may occur when coupled to another process; one must account for the coupled network rather than count the uphill step as an independent source of power.

A favourable ΔrG is therefore an energy-availability screen, not a prediction of chemical production or a demonstration of a prebiotic pathway. Power requires both a thermodynamic driving force and an actual rate.

## 6. The spatial model we should build next

The first local model should be a one-dimensional diagnostic of a geochemical interface, not yet a full fluid-dynamics or origin-of-life simulation.

### 6.1 Geometry and environmental profiles

Represent a transect from one fluid or environment into another. For a hydrothermal test case, this might be a mixing interface between a hydrothermal endmember and ambient water. The same model interface should later accept other scenarios, such as a hot-spring wet–dry boundary or a mineral pore.

At each position x, the input state should include:
- Temperature and pressure.
- Dimensionless activities or a complete analytical composition sufficient for speciation.
- pH and electrical potential where they are independently constrained.
- Redox couples and their measured or calculated potentials.
- Water activity and ionic strength.
- Relevant mineral phases and surface areas.
- Fluid velocity, diffusion coefficients and residence time.

Do not linearly interpolate pH or redox potential between endmembers. pH is logarithmic in proton activity; redox potentials depend on the redox couples and their activities. Instead, mix analytical concentrations under a declared transport/mixing model, solve chemical speciation, then calculate pH, activities and redox potentials from the resulting state.

### 6.2 Reaction thermodynamics

Select a small set of atom- and charge-balanced candidate reactions. For every position, calculate:
- The reaction quotient Q from species activities.
- The standard reaction Gibbs energy at local temperature and pressure.
- The actual reaction Gibbs energy.
- The reaction affinity.
- The maximum reversible work per reaction extent for downhill processes.

Do not choose a candidate reaction solely because its products are biologically familiar. Selection should be based on geochemical plausibility, reactant availability, catalyst evidence and sufficiently reliable thermodynamic data.

### 6.3 Transport and kinetics

A general species balance is

[
rac{partial c_i}{partial t}
=
-
ablacdotmathbf{N}_i
+sum_r 
u_{ir}r_r+s_i,
]

where N_i is the transport flux, ν_ir is the stoichiometric coefficient and s_i represents declared sources or sinks.

For neutral dilute species, an initial transport approximation could be

[
mathbf{N}_i=mathbf{u}c_i-D_i
abla c_i.
]

Charged species require a more complete treatment, such as Nernst–Planck transport with electric fields, activity corrections and electroneutrality constraints. A simple Fickian term is not enough for a proton/electrochemical-gradient model.

Reaction kinetics must be supplied independently from thermodynamics. Rate laws should account for the relevant catalysts, minerals, surface area, temperature dependence, inhibition and competing reactions. We should add these only after choosing a specific chemistry for which defensible data exist.

### 6.4 Gradient maintenance

A gradient is not a permanent reservoir. It is maintained by external boundary fluxes and is eroded by reaction, diffusion, mixing and leakage.

The model must compare:
- Supply rates of reactants or reducing/oxidizing equivalents.
- Reaction consumption rates.
- Diffusive and advective transport.
- Residence and replenishment times.
- Degradation or precipitation of potentially useful products.

A large initial ΔG that disappears quickly may be less consequential than a smaller gradient maintained for a long period. The relevant quantity for sustained organization is not merely instantaneous availability but the time-integrated, kinetically accessible chemical power under replenished conditions.

## 7. What the initial software module does — and does not do

The companion module at src/thermodynamic_emergence/chemistry.py implements a limited thermodynamic layer:

- Reaction Gibbs energy from supplied standard reaction Gibbs energy and dimensionless activities.
- Chemical affinity.
- The pH and electrical-potential contributions to proton electrochemical potential.
- Redox reaction Gibbs energy from a donor/acceptor potential difference.
- Chemical power density when a net reaction rate is supplied.
- Spatial profile evaluation when environmental activities and standard reaction Gibbs energies are supplied at each position.

It does not:
- Perform aqueous speciation or compute activity coefficients.
- Calculate redox potentials from full chemical composition.
- Interpolate fluid mixing or solve transport.
- Predict reaction rates or mineral catalysis.
- Model membrane formation, protocells, replication or life.

This distinction is intentional. A calculation based on supplied activities is a transparent thermodynamic tool; inventing the local chemistry would create false precision.

## 8. What counts as a useful local disequilibrium?

For this project, call a local disequilibrium a candidate usable chemical opportunity only after four separate questions have been addressed:

1. **Thermodynamic drive:** Is there a specified, balanced reaction or coupled process with available Gibbs energy at local conditions?
2. **Kinetic pathway:** Is there a plausible pathway and catalyst with a rate that is meaningful on environmental timescales?
3. **Gradient maintenance:** Are reactants replenished and products/heat removed so the disequilibrium can persist?
4. **Organizational compatibility:** Can products survive, concentrate, interact or be retained long enough to contribute to more elaborate chemical networks?

Only the first item is answered by equilibrium thermodynamics alone. A chemical opportunity is not yet a dissipative structure; a dissipative structure is not automatically alive.

## 9. Research guardrails

- Do not add pH, redox potential, temperature difference and photon flux into a single arbitrary energy score. They are different thermodynamic forces and can be coupled; double-counting is possible.
- Do not add proton-gradient energy to a reaction Gibbs energy if the proton chemical contribution is already represented through the reaction quotient and stoichiometry.
- Do not treat a negative ΔrG as proof of a reaction occurring at a useful rate.
- Do not treat present-day hydrothermal sites as direct replicas of Hadean environments.
- Do not assume hydrothermal vents are the unique or most likely origin site.
- Do not use planetary entropy production as an evolutionary objective.
- Do not claim abiogenesis has been modelled while the system starts with an imposed reaction network, catalyst, compartment or replicator that already assumes the key transition.

## 10. Immediate next research decisions

Before adding full transport or kinetics, select a benchmark geochemical scenario and assemble defensible endmember data. The benchmark should be chosen for data quality and testability, not because it best supports the research premise.

The next dataset should record, with provenance and uncertainty:
- Endmember temperature, pressure and fluid composition.
- Major and trace dissolved species, gas composition and mineral phases.
- Activity/speciation assumptions and redox reference conventions.
- Spatial mixing or flow measurements, where available.
- Candidate reactions and standard thermodynamic data over the relevant temperature-pressure range.

The first result should be a spatial map of calculated reaction Gibbs energies with uncertainty bounds, not a claim about life or the origin of life.

---

## References

[1] IUPAC. “Chemical potential.” Compendium of Chemical Terminology (Gold Book), 5th ed., 2025. https://doi.org/10.1351/goldbook.C01032

[2] Ianeselli, A. et al. (2023). “Physical non-equilibria for prebiotic nucleic acid chemistry.” Nature Reviews Physics, 5, 185–195. https://doi.org/10.1038/s42254-022-00550-3

[3] Barge, L. M., & Price, R. E. (2022). “Diverse geochemical conditions for prebiotic chemistry in shallow-sea alkaline hydrothermal vents.” Nature Geoscience, 15, 976–981. https://doi.org/10.1038/s41561-022-01067-1

[4] Rodriguez, L. E. et al. (2024). “A Geological and Chemical Context for the Origins of Life on Early Earth.” Astrobiology, 24(1 Suppl.). https://doi.org/10.1089/ast.2021.0139

[5] Russell, M. J. et al. (2010). “Serpentinization as a source of energy at the origin of life.” International Journal of Astrobiology. https://doi.org/10.1111/j.1472-4669.2010.00249.x

[6] “Serpentinization as the source of energy, electrons, organics, catalysts, nutrients and pH gradients for the origin of LUCA and life.” (2023). Open-access review. https://pmc.ncbi.nlm.nih.gov/articles/PMC10581274/

[7] Jackson, J. B. (2016). “Natural pH Gradients in Hydrothermal Alkali Vents Were Unlikely to Have Played a Role in the Origin of Life.” Journal of Molecular Evolution, 83, 1–11. https://doi.org/10.1007/s00239-016-9756-6

[8] Tutolo, B. M. (2026). “Rethinking the origin of life at seafloor hydrothermal vents.” Proceedings of the National Academy of Sciences, 123(36), e2613108123. https://doi.org/10.1073/pnas.2613108123

[9] Moldogazieva, N. T. et al. (2026). “Redox chemistry of early Earth and the origin of life.” Communications Chemistry, 9, 143. https://doi.org/10.1038/s42004-026-01969-w

[10] Krissansen-Totton, J., Bergsman, D. S., & Catling, D. C. (2016). “On Detecting Biospheres from Chemical Thermodynamic Disequilibrium in Planetary Atmospheres.” Astrobiology, 16(1), 39–67. https://doi.org/10.1089/ast.2015.1327
