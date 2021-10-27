function [N] = solvMPS(K,u,v)
% Creates solvable MPS two-site tensors for a dual unitary circuit
% according to the given parameters. This is limited to the case of qubits,
% ie. physical dimension d=2 and virtual dimension chi=d=2.
% INPUTS:
% K: A vector of three real numbers [K1,K2,K3]. Which are the angles used
% in our construction.
% u,v: Two SU(2) matrices. Putting [] for either of them will lead to it
% being the identity.
%
% OUTPUT:
% N: A four-legged tensor. Indexing order leftbond - phys1 - phys2 - 
%  rightbond. This is the two site tensor of our solvable MPS.
%
% Written by R. Milbradt
% For further details see DOI: 10.1103/PhysRevB.101.094304

if size(K) ~= 3
    error("Size of vector K is not 3!");
end

if isempty(u)
    u = eye(2);
end
if isempty(v)
    v = eye(2);
end

N = zeros(2,2,2,2);

%To get the right normalization in the end.
v = v/sqrt(2);

%All parts to be used are calculated here
expp = exp(1i * K(3));
expn = exp(-1i * K(3));
cosp = cos(K(1));
cosn = cos(K(2));
sinp = sin(K(1));
sinn = sin(K(2));

%And part by part put together here.
%N^(1,1)
nt = [[expn*cosn,0];[0,expp*cosp]];
nt = u*nt*v;
N(:,1,1,:) = nt;

%N^(1,2)
nt = [[0,-1i*expn*sinn];[-1i*expp*sinp,0]];
nt = u*nt*v;
N(:,1,2,:) = nt;

%N^(2,1)
nt = [[0,-1i*expp*sinp];[-1i*expn*sinn,0]];
nt = u*nt*v;
N(:,2,1,:) = nt;

%N^(2,2)
nt = [[expp*cosp,0];[0,expn*cosn]];
nt = u*nt*v;
N(:,2,2,:) = nt;
end

