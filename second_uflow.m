function U = second_uflow(U,Z,dt,F)
% Solve the equation \dot{U} = F(U,Z) with Z known
% see E. Celledoni and B. Owren. A class of intrinsic schemes for orthogonal integration.
% SIAM Journal on Numerical Analysis, 40(6):2069–2084, 2002. doi:10.1137/S0036142901385143.

% Tangent method RK-MK 2nd order
V = dt/2*F(U,Z);
R = retract(U,V);
FRV = F(R,Z);
Vtang = dt*invDretract(U,V,FRV,R);
U = retract(U,Vtang);

end
