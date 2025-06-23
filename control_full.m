function [control, P, nkIters, cpuTime] = control_full(A, Q, By, D, y, Pguess, tol_trunc, doLineSearch, maxLineSearchIters, maxOuterIters, coeffBurger, invRB, resA)

% ------------------------------------------------------------------------
%   [control, P, nkIters, cpuTime] = CONTROL_FULL(...)
%   computes the optimal control for a single state realization y by
%   solving the continuous‐time algebraic Riccati equation (CARE)
%   using a Newton–Kleinman iteration, then forming the feedback law.
%
%   Inputs:
%     A             - system matrix (Nh×Nh)
%     B             - control input matrix (Nh×mB)
%     Q             - state weighting matrix (Nh×Nh)
%     R             - control weighting matrix (mB×mB)
%     By            - B * inv(R) * B' (Nh×Nh)
%     D             - discretized differential operator for nonlinear term (Nh×Nh)
%     y             - current state vector (Nh×1)
%     Pguess        - initial Riccati guess (Nh×Nh)
%     tol_trunc     - tolerance for truncation in Riccati solver
%     doLineSearch  - flag to enable line search in Newton–Kleinman
%     maxLineSearchIters - maximum inner line search iterations
%     maxOuterIters - maximum Newton–Kleinman outer iterations
%     coeffBurger   - Burgers’ convection coefficient (scalar)
%     invRB         - map inv(R)*B' (mB×Nh)
%     resA          - function handle for Riccati residual: resA(Ay, X)
%
%   Outputs:
%     control       - optimal control vector for current y (mB×1)
%     P             - Riccati solution matrix (Nh×Nh)
%     nkIters       - number of Newton–Kleinman outer iterations
%     cpuTime       - CPU time for Newton–Kleinman solve (s)
% ------------------------------------------------------------------------

% 1) Form linearized system matrix with nonlinear term diag(D*y)
Ay = A + coeffBurger * diag(D * y);

% 2) Residual function for current Ay
residual = @(X) resA(Ay, X);

% 3) Solve CARE via Newton–Kleinman: P solves
%    Ay' P + P Ay - P By P + Q = 0
[P, nkIters, ~, cpuTime] = newton_kleinman(Ay, By, Q, Pguess, tol_trunc, maxOuterIters, ...
                                       residual, doLineSearch, maxLineSearchIters);

% 4) Compute feedback control: u = - R^{-1} B' P y = - invRB * P * y
control = -invRB * (P * y);
end