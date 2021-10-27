function [ng] = unitarisation(N)
%UNITARISATION Getting a unitary gate ng corresponding to N in a certain
%way.
%   N: The 4-tensor making up a solvable MPS. It us asszned
%   bond-dim.>=phys.-dim.
%   ng: The gate that is output. It will be unitary and is defined based on
%   N.
%
% Written by R. Milbradt
% For details refer to the Thesis section 8.2

% Matrix dimensions
d = size(N,2); %physical
chi = size(N,1); %bond/virtual

nchi = ceil(log(chi) / log(d)); % Number of qudits required to implement bond-dim.
zerod = d^nchi-chi; % Number of matrix elements to just be zero later.

ncol = d^nchi*d^2; % Total number of columns in our final gate

% Getting column vector of ng defined by N
v = zeros(ncol,chi); % Each column is one of the desired vectors.

tN1 = permute(N,[1,4,3,2]);
tN = reshape(tN1,[chi,ncol-zerod]);

for it = (1:chi)
    v(1:ncol-zerod,it) = tN(it,:);
end

% Putting the column vectors in the correct place
ng = zeros(ncol); % ng has to be a square matrix

%The left bond index alpha will be the first one, ie. of highest value,
%when considering the numerical system of base d
colvec = zeros(chi+1,1);

for it = (1:chi)
    bs = todbase(it-1,d);
    s = size(bs);
    col = 0;    
    
    for it2 = (1:s(2))
        col = col + d^(it2+1) * bs(it2);
    end
    colvec(it) = col+1; %Will be usefull later
    
    ng(:,col+1) = v(:,it);
end

% Completing the gate to unitarity
W = null(v'); %Find a basis orthonormal to the complex column-vectors in v.

count = 1;
for it = (1:ncol)
    %If the column was not assigned something before, assign it a vector
    %that is orthonormal to the rest.
    if ~(it == colvec(count))
        ng(:,it) = W(:,it-count+1);
    else
        count = count + 1;
    end
end

end

function [bs] = todbase(j,d)
% Converts the decimal number j to a number in base d, the different terms
% are given as a vector where j = sum_i (d^(i-1)*bs_i).

i = 0;

while j > d
    i = i+1;
    bs(i) = mod(j,d);
    j = floor(j/d);
end

bs(i+1) = j;

end