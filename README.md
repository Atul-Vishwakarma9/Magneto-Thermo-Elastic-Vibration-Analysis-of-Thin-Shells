# Nonlinear Vibration Analysis of Magneto–Thermo–Elastic Coupled Cylindrical Shells

This repository contains the analytical formulation, MATLAB implementation, and numerical results of a **nonlinear free vibration analysis of isotropic and functionally graded (FGM) cylindrical shells** subjected to **combined magnetic and thermal fields**.

The work is carried out as a **B.Tech Mechanical Engineering final year project** at **NIT Raipur**.

---

## Problem Overview

Cylindrical shell structures are widely used in aerospace, nuclear, marine, and high-temperature engineering applications. When subjected to **thermal gradients and magnetic fields**, their dynamic behavior becomes strongly nonlinear.

This project investigates:
- **Nonlinear vibration behavior**
- **Frequency shifts due to magneto-thermo-elastic coupling**
- **Influence of material gradation laws**

---

##  Objectives

- Study **nonlinear free vibration** of cylindrical shells
- Incorporate **magnetic (Lorentz) forces** and **thermal stresses**
- Compare **Isotropic vs Functionally Graded Materials**
- Analyze different FGM laws:
  - Power-law
  - Sigmoidal
  - Hyperbolic
- Solve governing equations using **Galerkin and Multi-Scale methods**

---

## Theory & Methodology

### ✔ Shell Theory
- Classical Shell Theory (Kirchhoff–Love assumptions)
- Thin shell formulation

### ✔ Energy Formulation
- Hamilton’s Principle
- Kinetic Energy
- Strain Energy
- Work done by Lorentz forces

### ✔ Numerical Techniques
- Galerkin Method (PDE → ODE)
- Multi-Scale Perturbation Method (Nonlinear solution)

---

##  Material Modeling

### Functionally Graded Materials (FGM)
Material properties vary continuously through thickness using:

- **Power Law**
- **Sigmoidal Function**
- **Hyperbolic Function**

Properties considered:
- Young’s Modulus
- Density
- Thermal expansion coefficient
- Thermal conductivity (temperature-dependent)

### Isotropic Materials
- Aluminum
- Nickel
- Silicon Carbide (SiC)

---

##  MATLAB Code Structure

| File | Description |
|-----|------------|
| `main.m` | Main execution file |
| `material_properties.m` | Temperature-dependent material properties |
| `stiffness_matrices.m` | A, B, D stiffness matrices |
| `magnetic_force.m` | Lorentz force computation |
| `thermal_load.m` | Thermal stress calculation |
| `galerkin_reduction.m` | Governing ODE derivation |
| `multiscale_solution.m` | Nonlinear frequency solution |
| `plot_results.m` | Frequency response plots |

---

## Key Results

- **Natural frequency decreases** with:
  - Increasing temperature
  - Increasing magnetic field
  - Higher power-law index
- **Frequency increases** with:
  - Higher thickness ratio (h/R)
  - Higher circumferential mode number (n)
- FGM shells show **better stiffness retention** than isotropic shells
- Sigmoidal FGMs provide **more stable nonlinear response**

---

## Model Validation

The formulation is validated by:
- Reducing FGM to isotropic case
- Neglecting magnetic & thermal fields
- Linearizing nonlinear equations
- Comparing results with published literature

Deviation observed: **< 5%**

---

## Repository Contents

- `matlab_code/` – Complete simulation codes
- `results/` – Parametric study plots
- `report/` – Full B.Tech thesis (PDF)
- `figures/` – Frequency response graphs
- `validation/` – Benchmark comparison

---

##  Requirements

- MATLAB R2021a or later
- Symbolic Math Toolbox (optional)

---

##  Authors

- **Atul Kumar Vishwakarma**  
 

Department of Mechanical Engineering  
National Institute of Technology Raipur

---

## Supervisor

**Dr. Ankur Gupta**  
Assistant Professor, NIT Raipur

---

##  License

This project is licensed for **academic and research use only**.
