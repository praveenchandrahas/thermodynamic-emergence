# Local chemical-disequilibrium toolkit (v0.1)

This package adds a thermodynamic calculator for local chemistry. It is a diagnostic layer over supplied chemical activities and environmental profiles, not a fluid simulation or an origin-of-life model.

## What it calculates

- Reaction Gibbs energy from a balanced reaction, standard reaction Gibbs energy at local T/P, and dimensionless activities.
- Reaction affinity.
- Proton electrochemical potential difference from pH and electric potential.
- Redox reaction Gibbs energy from Nernst-corrected donor and acceptor reduction potentials.
- Chemical power density when a net reaction-extent rate is supplied.
- Spatial profiles of reaction Gibbs energy when activities and standard reaction Gibbs energy are supplied at every point.

## Important input requirement

The module does not infer chemical activities from concentrations. A geochemical workflow must first perform aqueous speciation and activity corrections, or clearly document an ideal-solution approximation. Standard reaction Gibbs energy must correspond to local temperature, pressure and standard-state conventions.

Do not linearly interpolate pH or redox potential between endmembers. Mix analytical composition under a specified mixing model, calculate speciation and activities, then derive pH and redox potentials.

## Run

Install the project from its root:

    python -m pip install -e .

Run the illustrative pH-gradient calculation:

    python examples/proton_gradient_demo.py

Run the unit tests:

    python -m unittest discover -s tests -v

## Demonstration

The example uses a synthetic pH profile only to demonstrate the electrochemical-potential calculation. It is not a measurement or a proposed reconstruction of a particular hydrothermal vent.

At 298.15 K, a three-unit pH contrast contributes approximately 17.1 kJ/mol per mole of protons to the chemical component of the electrochemical potential difference. An electrical potential difference changes this value by F times the potential difference. Neither number establishes that a naturally occurring membrane or prebiotic molecular machine can access the gradient.

## What is deliberately not included

- Aqueous speciation or activity coefficients.
- Mineral equilibria, adsorption and precipitation.
- Redox-potential prediction from composition.
- Advection, diffusion, electromigration or fluid dynamics.
- Reaction kinetics or catalytic mechanisms.
- Autocatalysis, compartment growth, replication or evolution.

The next model requires a benchmark geochemical scenario with traceable endmember data and thermodynamic constants. The software should not invent such data in order to produce an attractive map.
