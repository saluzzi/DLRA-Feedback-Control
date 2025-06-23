function U = second_uflow(U,Z,dt,F)
% Solve the equation
% \dot{U} = F(U,Z) with Z known

% Tangent method RK-MK 2nd order
V = dt/2*F(U,Z);
R = retract(U,V);
FRV = F(R,Z);
Vtang = dt*invDretract(U,V,FRV,R);
U = retract(U,Vtang);

end
