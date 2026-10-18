# Adiabatic Flame Temperature & Emissions Analysis (C3H8 - H2 Blends)

[![MATLAB](https://img.shields.io/badge/Language-MATLAB-orange.svg)](https://www.mathworks.com/)
[![Status](https://img.shields.io/badge/Status-Completed-success.svg)]()

Numerical investigation of adiabatic flame temperatures (T_ad) and specific carbon/pollutant emissions for varying propane (C3H8) and hydrogen (H2) fuel blends across a wide range of equivalence ratios (phi = 0.1 to 2.0).

---

## Overview
This project models the combustion thermochemistry of hydrocarbon-hydrogen fuel blends reacting with air at constant pressure (p_0 = 1 atm) and standard inlet temperature (T_0 = 298 K). 

The primary objective is to evaluate how blending molecular hydrogen into a propane stream affects:
* Peak adiabatic flame temperatures (T_ad).
* High-temperature product dissociation mass fractions (w_CO, w_H2).
* Specific carbon monoxide (CO) and carbon dioxide (CO2) emissions under fuel-rich conditions (phi > 1.0).

---

## Key Features & Methodology
* **NASA 7-Coefficient Polynomials:** Employs NASA thermodynamic property tables to calculate temperature-dependent molar heat capacities (c_p), sensible enthalpies, and entropies.
* **Equilibrium Modeling:** 
  * Complete oxidation assumed for lean and stoichiometric conditions (phi <= 1.0).
  * Water-gas shift equilibrium (CO + H2O = CO2 + H2) solved numerically using K_p(T) for rich conditions (phi > 1.0).
* **Root-Finding Solver:** Iteratively solves the enthalpy balance equation (H_reactants(T_0) = H_products(T_ad)) via MATLAB's `fzero` routine using Brent's method.
* **Diversified Fuel Blends:** Evaluates 5 distinct volumetric compositions:
  1. 100% C3H8 (Pure Propane)
  2. 80% C3H8 + 20% H2
  3. 50% C3H8 + 50% H2
  4. 20% C3H8 + 80% H2
  5. 100% H2 (Pure Hydrogen)

---

## Project Structure
```text
├── lab1CombustionFINAL.m     # Main MATLAB script containing model, solver, and plotting routines
├── Lab_1_Combustion.pdf      # Detailed lab report and visualization outputs
└── README.md                 # Project documentation
