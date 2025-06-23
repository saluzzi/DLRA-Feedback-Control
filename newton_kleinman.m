function [X, numIter, resHistory, cpuTime] = newton_kleinman(A, B, C, X0, tol, maxIter, resFun, doLineSearch, maxLineSearch)
% NEWTON_KLEINMAN  Solve continuous-time algebraic Riccati equation (CARE)
%   A'*X + X*A - X*B*X + C = 0
% using the Newton–Kleinman iterative method with optional line search.
%
% INPUTS:
%   A               - System matrix (n x n)
%   B               - Control weight matrix in CARE (n x n)
%   C               - State cost matrix (n x n)
%   X0              - Initial guess for X (n x n)
%   tol             - Convergence tolerance for residual norm
%   maxIter         - Maximum Newton–Kleinman iterations
%   resFun          - Function handle: resFun(X) returns normalized
%                     CARE residual norm
%   doLineSearch    - Boolean flag to enable Armijo-like line search
%   maxLineSearch   - Maximum number of line-search iterations
%
% OUTPUTS:
%   X               - Approximate solution to CARE
%   numIter         - Number of Newton iterations performed
%   resHistory      - Vector of residual norms per iteration
%   cpuTime         - Total CPU time for the iterative solver (seconds)

    % Initialize
    numIter    = 0;
    X          = X0;
    res_old    = resFun(X);
    resHistory = [];

    % Predefine residual operator for line search (if used)
    if doLineSearch
        Rfunc = @(Acur, Xcur) (Acur' * Xcur + Xcur * Acur - Xcur * B * Xcur + C);
        % quartic model
        modelPhi = @(t, a, b, c) a*(1 - t)^2 - 2*b*(1 - t)*t^2 + c*t^4;
    end

    tic;  % Start timing
    % Main Newton–Kleinman loop
    while res_old > max(tol, 1e-10) && numIter < maxIter
        numIter = numIter + 1;

        % Form linear Lyapunov operator: A1 = A - B*X
        A1 = A - B * X;

        % Compute right-hand side
        Qrhs = X * B * X + C;

        % Solve Sylvester equation: A1' * Xnew + Xnew * A1 + Qrhs = 0
        % Equivalent to sylvester(A1', A1, -Qrhs)
        Xnew = lyap(A1', Qrhs);
        Xnew = sparse(Xnew);

        % Optional line search to improve convergence
        if doLineSearch && numIter <= maxLineSearch
            Z = Xnew - X;
            V = Z * B * Z;
            Rx = Rfunc(A, X);
            a  = trace(Rx * Rx);
            b  = trace(Rx * V);
            c  = trace(V * V);
            % Minimize phi(t) over t in [0,2]
            phi = @(t) modelPhi(t, a, b, c);
            t_opt = fminbnd(phi, 0, 2);
            X = X + t_opt * Z;
        else
            % Standard Newton step
            X = Xnew;
        end

        % Record residual and CPU time
        res_new      = resFun(X);
        resHistory(end+1) = res_new;
        res_old      = res_new;
    end
    cpuTime = toc;

end