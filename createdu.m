function [DU] = createdu(J,up,um,vm,vp,phi)
%Creates a dual unitary via the standard parametrisation for the
%given parameters.
%
%   DU = exp(i*phi) * kron(up,um) * V[J] * kron(vm,vp)
%   V[J] = exp(-i * (pi/4 * kron(sigx,sigx) + pi/4 * kron(sigy,sigy) + J
%   kron(sigz,sigz) ) )
%
% J, phi: real numbers
% up,um,vm,vp: SU(2) matrices. If any is given as [], it will be treated as
% the identity.
%
% For J = phi = pi/4, UD is the SWAP-gate
%
%For more details see DOI: 10.1103/PhysRevLett.123.210601
%Written by R. Milbradt

if isempty(up)
    up = eye(2);
end
if isempty(um)
    um = eye(2);
end
if isempty(vm)
    vm = eye(2);
end
if isempty(vp)
    vp = eye(2);
end

sigx = [[0,1];[1,0]];
sigy = [[0,-1i];[1i,0]];
sigz = [[1,0];[0,-1]];

V = -1i * (pi/4 * kron(sigx,sigx) + pi/4 * kron(sigy,sigy) + J * kron(sigz,sigz));
V = expm(V);
DU = exp(1i * phi) * kron(up,um) * V * kron(vm,vp);

end

