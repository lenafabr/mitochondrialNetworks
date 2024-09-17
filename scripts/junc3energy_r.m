function [E3,E3all] = junc3energy_r(rj,r1,r2,r3)
% 3-way junction energy, rj = junciton position, 
% r1,r2,r3 are the other nodes
% r3 is the incoming node
% E3 only includes angles containing the incoming segment
% E3all includes all the angles

v1 = r1-rj; nv1 = norm(v1);
v2 = r2-rj; nv2 = norm(v2);
v3 = r3-rj; nv3 = norm(v3);

cb(1) = -v1*v3'/nv1/nv3;
cb(2) = -v2*v3'/nv2/nv3;
cb(3) = -v2*v1'/nv2/nv1;

cb= max(cb,-1);
cb = min(cb,1);

sb = sqrt(1-cb.^2);

% cos(pi/3 - beta)
cshift = 0.5*cb + sqrt(3)/2*sb;

E3 = 2 - sum(cshift(1:2));
E3all = 3-sum(cshift);
end