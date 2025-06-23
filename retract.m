function R = retract(Q,V)
% Compute the retraction R_Q(V)
% Input: Q   Nxn unitary matrix
%        V   Nxn such that Q'V is skew-hermitian
% Output: R = R_Q(V)   Nxn unitary matrix

[~, n] = size(Q);
% Compute g_Q^S(V)
gV = V - Q*(1/2*Q'*V);
% Compute the low-rank factors of V
A = [gV,-Q]; B = [Q,gV];
D = B'*A;

R = Q + A*((eye(2*n)-D./2)\(B'*Q));

end
