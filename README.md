# DLRA-Feedback-Control

This repository implements dynamical low-rank approximation (DLRA) techniques for optimal control of the nonlinear Burgers' equation in 1D and 2D using feedback laws obtained from the solution of matrix Riccati equations.

For more detail we refer to
L. Saluzzi and M. Strazzullo, "Dynamical Low-Rank Approximation Strategies for Nonlinear Feedback Control Problems" arXiv preprint arXiv:2501.07439 (2025).

## Features

- Approximation of the solution via dynamical low-rank method.
- Feedback control using Riccati-based Newton-Kleinman iterations and the cascade approach.
- Rank enrichment strategy.
  
## Files Overview

### Main Scripts

- `burgers_control_main.m`: Main script to run the control simulation for the Burgers' equation.
- `setting_Burgers1D.m`: Defines spatial/temporal discretization and operators for the 1D case.
- `setting_Burgers2D.m`: Defines spatial/temporal discretization and operators for the 2D case.

### Auxiliary Functions

- `control_full.m`: Solves the full-order control problem using Newton-Kleinman.
- `newton_kleinman.m`: Newton-Kleinman iterative solver for the continuous-time algebraic Riccati equation.
- `second_uflow.m`: Second-order Runge-Kutta-Munthe-Kaas method to evolve the low-rank factor.
- `retract.m`: Retraction operator ensuring orthonormality of updated basis on the Stiefel manifold.  
- `invDretract.m`: Inverse differential of the retraction map for accurate tangent updates.

## Notes

- **Adaptive Rank**: When `isAdaptive = true`, the low-rank basis is enriched with information from the Riccati solution to capture control dynamics.
- **Riccati Solver**: The Newton–Kleinman routine supports line search; you can toggle `linesearch` and adjust tolerances for convergence behavior.


