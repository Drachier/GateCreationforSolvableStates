function [U] = randSU2()
%RANDSU2 Creates a random matrix in the SU(2) group.
% Written by R. Milbradt

%Create random parameter values
r1 = rand;
r2 = sqrt(1- r1^2);
theta1 = 2*pi*rand;
theta2 = 2*pi*rand;

%Use parametrisation to create U
U = zeros(2);

U(1,1) = r1 * exp(1i * theta1);
U(1,2) = r2 * -exp(-1i * theta2);
U(2,1) = r2 * exp(1i * theta2);
U(2,2) = r1 * exp(-1i * theta1);

end

