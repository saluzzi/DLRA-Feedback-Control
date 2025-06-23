%% SETTING BURGERS1 EQUATION IN 1D

% Spatial discretization
Nh         = 100;            % number of grid points
a          = -5;             % domain start
b          =  30;            % domain end
x          = linspace(a, b, Nh)';
dx         = x(2) - x(1);

% Temporal discretization
T          = 1;           % final time
dt         = 1e-3;           % time step size
nt         = floor(T / dt);  % number of time steps

% Physical parameters
coeff_burger = 1;            % nonlinear convection coefficient
c            = 20;           % linear transport speed
alpha        = 0.1;          % regularization weight

% Basis selection options

tol_ratio = 0.9999;   % Energy cutoff for rank truncation
isRankEnrich = true;           % If true: use the enriching rank strategy
if isRankEnrich
    isP0       = true;           % If true: enrich basis from the Riccati solution;
end                              % else from the feedback matrix

% Newton-Kleinman settings
tol_trunc      = 1e-9;
linesearch     = true;
it_linesearch  = 1;
outer_it       = 30;

% Control domain indicators
% B-support intervals
a_B1 = -1; b_B1 = 1;
a_B2 =  5; b_B2 = 7;
vec_B = (x>=a_B1 & x<=b_B1) | (x>=a_B2 & x<=b_B2);

% Q-support intervals
a_Q1 = 2;  b_Q1 = 4;
a_Q2 = 8;  b_Q2 = 10;
vec_Q = (x>=a_Q1 & x<=b_Q1) | (x>=a_Q2 & x<=b_Q2);

% Build Control and Observation Operators
m_B = sum(vec_B);
m_Q = sum(vec_Q);
B    = sparse(find(vec_B), 1:m_B, 1, Nh, m_B);
Q    = sparse(find(vec_Q), find(vec_Q), dx, Nh, Nh);
R    = alpha * speye(m_B) * dx;
invRB = R\ B';
By   = B * invRB;

% Spatial operators
% Upwind difference for advection
A = sparse(-eye(Nh) + diag(ones(Nh-1,1), -1));
A = c * A / dx;
% Forward difference for nonlinear term
D = sparse(-eye(Nh) + diag(ones(Nh-1,1), +1));
D = D / dx;

% Tensor for nonlinear term y.*(Dy)
TT = sparse(Nh, Nh^2);
for k = 1:Nh-1
    TT(k, (k-1)*Nh + k)     = -1;
    TT(k, (k-1)*Nh + k + 1) =  1;
end
TT(Nh, end) = -1;
TT = (coeff_burger / dx) * TT;

%% Initial Condition Ensemble
np_k  = 20;  % number of k-parameters
np_kk = 20;  % number of kk-parameters

y0fun = @(x,k, kk) k*exp(-kk*(x').^2);

k_vec = linspace(0.1,0.5,np_k);
kk_vec = linspace(0.2,1.5,np_kk);
y0 = [];
for i = 1:np_k
    for j = 1:np_kk
        y0 = [y0 y0fun(x,k_vec(i),kk_vec(j))'];
    end
end
np = np_k*np_kk;

