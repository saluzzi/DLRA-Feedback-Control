function Vtilde = invDretract(Q,V,W,R)
% Compute the inverse of the tangent map (dR_Q at V)^(-1) W
% Input: Q   Nxn unitary matrix
%        V   Nxn such that Q'V is skew-Hermitian
%        W   Nxn such that R'W is skew-Hermitian
%        R   is the retraction retract(Q,V)
% Output: Vtilde = (dR_Q at V)^(-1) W   Nxn such that Q'Vtilde is skew-Hermitian

gV = V - Q*(1/2*Q'*V);
A = [gV,-Q];
B = [Q,gV];

In = eye(size(Q,2));
T2 = (2*W-A*(B'*W))/(Q'*R+In);
T1 = - 1/2*(Q'*T2+T2'*Q)-1/2*(R'*Q+In)\((R+Q)'*T2);
gVtilde = Q*T1 + T2;
Vtilde = gVtilde + Q*(Q'*gVtilde);

end
