# DLRA_Feedback_Control

This repository provides MATLAB scripts to solve a PDE-constrained optimal control problem for a one-dimensional Burgers'-type equation. It implements both a **Full-Order Model (FOM)** solver and a **Dynamic Low-Rank Approximation (DLRA)** solver, compares their performance, and reports error metrics.

## Repository Structure

- **main1D.m**
  - Entry-point script that orchestrates the full-order and low-rank simulations.
  - Sets up spatial and temporal discretization, physical parameters, control/observation operators, and initial condition ensembles.
  - Runs the **Full-Order Model** (via `control_full.m`) across parameter realizations and time steps, collects state trajectories, controls, and costs.
  - Initializes low-rank basis via SVD, optionally augmented by Riccati directions, and runs the **DLRA** simulation using `control_proj_T.m`.
  - Computes and prints error metrics between FOM and DLRA results, and reports runtimes.

- **control_full.m**
  - Implements the full-dimensional optimal control solver for a single state snapshot `y`.
  - Linearizes the system around the current state, forming matrix \(A_y = A + \text{diag}(D y)\).
  - Solves the continuous-time algebraic Riccati equation (CARE) using a Newton–Kleinman iteration (`newton_kleinman`) with optional line search.
  - Computes the optimal feedback control \(u = -R^{-1} B^T P y\) and returns the control vector, Riccati solution \(P\), iteration count, and CPU time.

- **control_proj_T.m**
  - Constructs a reduced-order projected control solver for the DLRA framework.
  - Projects full-order operators \(A, B, Q, B_y\) onto a low-rank basis \(U\), forming reduced matrices \(A_p, Q_p, B_p\).
  - Builds a tensor contraction for the nonlinear Burgers' term in the reduced space.
  - For each parameter realization, solves a reduced CARE via Newton–Kleinman, reusing previous solutions in blocks of size `nk`.
  - Computes both the reduced dual control contributions (`umatUT`) and the full-order control signals (`Ufull`), as well as a feedforward correction term.

## Dependencies

- MATLAB R2020a or later (for sparse matrices and `svd('econ')`)
- No additional toolboxes are strictly required; the code relies on base MATLAB functions and a custom `newton_kleinman` implementation.

## Usage

1. Clone the repository and add its directory to the MATLAB path.
2. Adjust parameters (grid size, time step, regularization weight, adaptive rank settings) at the top of **main1D.m**.
3. Run `main1D.m` in MATLAB:
   ```matlab
   >> main1D
   ```
4. Review printed outputs for simulation timings, mean errors, and optionally inspect stored trajectories in the workspace.

## Notes

- **Adaptive Rank**: When `isAdaptive = true`, the low-rank basis is enriched with directions from the Riccati solution to capture control dynamics.
- **Riccati Solver**: The Newton–Kleinman routine supports line search; you can toggle `linesearch` and adjust tolerances for convergence behavior.
- **Extensibility**: The modular structure allows swapping in different time integrators (e.g., higher-order RK) or control solvers.

---
*For questions or contributions, please open an issue or submit a pull request.*

