%% SETTING BURGERS’ EQUATION IN 2D

% Spatial discretization parameters
Nh1           = 20;              % Number of grid points per dimension
Nh            = Nh1^2;           % Total number of spatial degrees of freedom
a             = -10;             % Domain start (in both x and y)
b             =  10;             % Domain end (in both x and y)
x             = linspace(a, b, Nh1);  % 1D grid vector
dx            = x(2) - x(1);     % Grid spacing

% Construct 2D derivative operator using Kronecker products
D1            = sparse((-eye(Nh1) + diag(ones(Nh1-1,1), -1)) / dx);  
                                   % 1D backward‐difference matrix
D             = kron(D1, eye(Nh1)) + kron(eye(Nh1), D1);  
                                   % 2D Laplacian‐like operator for advective terms

% Create 2D mesh coordinates (flattened into vectors)
[X1, X2]      = meshgrid(x, x);  % 2D grid for plotting or ICs
xx            = X1(:);           % Flattened x‐coordinates
yy            = X2(:);           % Flattened y‐coordinates

%% Temporal parameters
T             = 1;            % Final time
dt            = 1e-3;            % Time step
nt            = floor(T / dt);   % Number of time steps

%% Physical parameters

c             = 20;              % Combined transport speed
alpha         = 0.1;             % Control regularization weight
coeff_burger  = 1;               % Nonlinear convection coefficient

%% Parametric sampling for initial conditions
np_k          = 20;               % Number of amplitude parameters
np_kk         = 20;               % Number of width parameters

%% Rank‐enrichment strategy flags
tol_ratio     = 0.9999;          % Energy cutoff for basis truncation
isRankEnrich  = true;            % Enable adaptive rank enrichment
if isRankEnrich
    isP0       = true;           % Use Riccati solution for enrichment if true
end

%% Newton–Kleinman solver settings
tol_trunc     = 1e-9;            % Riccati truncation tolerance
linesearch    = true;            % Enable line search in Newton–Kleinman
it_linesearch = 1;               % Max inner line‐search iterations
outer_it      = 30;              % Max Newton–Kleinman iterations

%% Define control (B) and observation (Q) support regions
a_B = [-5, -5];                  % Lower‐left corner of control rectangle
b_B = [ 0,  0];                  % Upper‐right corner
a_Q = [ 0,  0];                  % Lower‐left corner of observation
b_Q = [ 5,  5];                  % Upper‐right corner

vec_B         = (xx>=a_B(1) & xx<=b_B(1)) & (yy>=a_B(2) & yy<=b_B(2));
                                   % Logical mask for control region
vec_Q         = (xx>=a_Q(1) & xx<=b_Q(1)) & (yy>=a_Q(2) & yy<=b_Q(2));
                                   % Logical mask for observation region

% Build sparse control and observation matrices
m_B           = sum(vec_B);      % Number of control inputs
m_Q           = sum(vec_Q);      % Number of observed states
B             = sparse(Nh, m_B);
B(vec_B, :)   = speye(m_B);      % Select rows corresponding to control region
Q             = sparse(Nh, Nh);
Q(vec_Q, vec_Q) = speye(m_Q) * dx; % Weighted observation on Q‐region

% Control cost matrix and its inverse times B'
R             = alpha * speye(m_B) * dx;
invRB         = R \ B';
By            = B * invRB;       % Precomputed B R⁻¹ B'

%% Spatial advection operators
A             = c * sparse(D);   % Combined advection operator
                                   % (here D applies backward‐differences in both dims)

%% Nonlinear convection tensor for y.*(D y)
TT            = sparse(Nh, Nh^2);
% Interior points: second‐order finite difference stencil
for k = 1:Nh1
    idx = k;  % first row
    TT(idx, (idx-1)*Nh + idx)     = -2;
    if k > 1
        TT(idx, (idx-1)*Nh + idx - 1) = 1;
    end
end
% Remaining rows
for k = Nh1+1:Nh
    TT(k, (k-1)*Nh + k)     = -2;
    if mod(k-1, Nh1) ~= 0
        TT(k, (k-1)*Nh + k - 1) = 1;
    end
    TT(k, (k-1)*Nh + k - Nh1) = 1;
end
TT = (coeff_burger / dx) * TT;   % Scale by coefficient and grid spacing

%% Initial condition ensemble (Gaussian bumps)
y0fun = @(x, y, k, kk) k * exp(-kk * ((x+2).^2 + (y+2).^2));
k_vec = linspace(0.01, 0.05, np_k);
kk_vec= linspace(0.1,  0.3,  np_kk);
y0    = zeros(Nh, np_k*np_kk);
col   = 1;
for i = 1:np_k
    for j = 1:np_kk
        y0(:, col) = y0fun(xx, yy, k_vec(i), kk_vec(j));
        col = col + 1;
    end
end
np = np_k * np_kk;                % Total number of parameter samples
