# DLRA_Feedback_Control

This repository provides MATLAB scripts to compute feedback control with PDE constrants via **Dynamic Low-Rank Approximation (DLRA)** solver.

## Repository Structure

- **main_Burgers1D.m**: Compute feedback control for the 1D Burgers' equation, comparing FOM and DLRA strategy.  
- **control_full.m**: Full-order Riccati solver and control computation.  
- **control_proj.m**: Projected control solver for reduced-order model.  
- **newton_kleinman.m**: Newton–Kleinman iterative solver for the continuous-time algebraic Riccati equation.  
- **second_uflow.m**: 2nd-order Runge–Kutta integrator on the Stiefel manifold for low-rank basis evolution.  
- **retract.m**: Retraction operator ensuring orthonormality of updated basis on the Stiefel manifold.  
- **invDretract.m**: Inverse differential of the retraction map for accurate tangent updates.  

## Notes

- **Adaptive Rank**: When `isAdaptive = true`, the low-rank basis is enriched with information from the Riccati solution to capture control dynamics.
- **Riccati Solver**: The Newton–Kleinman routine supports line search; you can toggle `linesearch` and adjust tolerances for convergence behavior.


