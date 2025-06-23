function output = control_proj(A,B,Q,By,U,Ztilde,tol_trunc,linesearch,it_linesearch,outer_it,np,nk,m_B,invRB,T,l,Guess)

%   output = CONTROL_PROJ_T(..)
%   computes the control actions in reduced space using the Newton–Kleinman
%   iteration on the projected Riccati equation.
%
%   Inputs:
%     A          - full-order system matrix (Nh×Nh)
%     B          - full-order control input matrix (Nh×mB)
%     Q          - full-order state weighting matrix (Nh×Nh)
%     By         - B*inv(R)*B' for full-order (Nh×Nh)
%     U          - projection basis (Nh×r)
%     Ztilde     - reduced coordinates (r×np)
%     tol_trunc  - truncation tolerance for Riccati updates
%     linesearch - flag for line search in Newton–Kleinman
%     it_linesearch - max iterations for line search
%     outer_it - max Newton–Kleinman outer iterations
%     np         - number of parameter realizations
%     nk         - block size for reusing Riccati guesses
%     m_B         - number of actuators
%     invRB      - inv(R)*B' full-order map (m_B×Nh)
%     T          - tensor operator for nonlinear term (Nh×Nh^2)
%     l          - reduced dimension (rank)
%     Guess    - initial guess for reduced Riccati (l×l)
%
%   Output (cell array):
%     {1} umatUT  - projected control contributions for U (np×l)
%     {2} umatZT  - projected control contributions for Z (Nh×np)
%     {3} Ufull   - full-order control signals (m_B×np)
%     {4} Pguess0 - initial Riccati guess for first block

% Precompute projected operators
Z = Ztilde';
Tp = U'*T*(kron(U,U));
UAU = U'*A*U;
umatUT = zeros(np,l);
Byp = U'*By*U;
Qp = U'*Q*U;
invRBp = invRB*U;
inv_B_U = invRBp'*B'*U;

% Preallocate outputs
control = zeros(m_B,np);
Pmat = zeros(l,np);

% Riccati residual function handle
res_A = @(Ap,X) norm(Ap'*X+X*Ap-X*Byp*X+Qp,'fro')/norm(Qp,'fro');


% --- Process first parameter ---
Xguess = Guess;
Ayp = zeros(l);
for j = 1:l
    Ayp(:,j) = Tp(:,(j-1)*l+(1:l))*Ztilde(:,1);
end
Ayp = Ayp+UAU;

res = @(X) res_A(Ayp,X);
% Solve Riccati
P = newton_kleinman(Ayp,Byp,Qp,Xguess,tol_trunc,outer_it,res,linesearch,it_linesearch);

Pguess = P;
P_guess_mu0 = P;
Pmat(:,1) = P*Ztilde(:,1);
control(:,1) = -invRBp*(P*Ztilde(:,1));
umatUT(1,:) = -Z(1,:)*(P*inv_B_U);

% --- Loop over remaining parameterss ---
for j = 2:np

    Ayp = zeros(l);

    for k = 1:l

        Ayp(:,k) = Tp(:,(k-1)*l+(1:l))*Ztilde(:,j);
    end

    Ayp = Ayp+UAU;

    res = @(X) res_A(Ayp,X);

    if mod(j-1,nk)==0
        Xguess = Pguess;
    else
        Xguess = P;
    end

    P = newton_kleinman(Ayp,Byp,Qp,Xguess,tol_trunc,outer_it,res,linesearch,it_linesearch);

    if mod(j-1,nk) == 0

        Pguess = P;
    end

    Pmat(:,j) = P*Ztilde(:,j);
    control(:,j) = -invRBp*(P*Ztilde(:,j));
    umatUT(j,:) = -Z(j,:)*(P*inv_B_U);

end

umatZT = -B*(invRBp*(Pmat*Z));

% Outputs
output{1} = umatUT;
output{2} = umatZT;
output{3} = control;
output{4} = P_guess_mu0;

end