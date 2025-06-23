# DLRA_Feedback_Control

This repository provides MATLAB scripts to solve a PDE-constrained optimal control problem for a one-dimensional Burgers'-type equation. It implements both a **Full-Order Model (FOM)** solver and a **Dynamic Low-Rank Approximation (DLRA)** solver, compares their performance, and reports error metrics.

## Repository Structure

- **main1D.m** 
  - Entry-point script that orchestrates the full-order and low-rank simulations.
  - Sets up spatial and temporal discretization, physical parameters, control/observation operators, and initial condition ensembles.
  - Runs the **Full-Order Model** (via `control_full.m`) across parameter realizations and time steps, collects state trajectories and controls.
  - Initializes low-rank basis via SVD, optionally enriched by information on the Riccati solutions, and runs the **DLRA** simulation using `control_proj_T.m`.
  - Computes and prints error metrics between FOM and DLRA results, and reports runtimes.

- **control_full.m**
  - Implements the full-dimensional optimal control solver for a single state snapshot `y`.
  - Solves the continuous-time algebraic Riccati equation (CARE) using a Newton–Kleinman iteration (`newton_kleinman`) with optional line search.
  - Computes the optimal feedback control \(u = -R^{-1} B^T P y\) and returns the control vector, Riccati solution \(P\), iteration count, and CPU time.

- **control_proj_T.m**
  - Constructs a reduced-order projected control solver for the DLRA framework.
  - For each parameter realization, solves a reduced CARE via Newton–Kleinman, reusing previous solutions in blocks of size `nk`.
  - Computes the reduces controls.


## Notes

- **Adaptive Rank**: When `isAdaptive = true`, the low-rank basis is enriched with information from the Riccati solution to capture control dynamics.
- **Riccati Solver**: The Newton–Kleinman routine supports line search; you can toggle `linesearch` and adjust tolerances for convergence behavior.


