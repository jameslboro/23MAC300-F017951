% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%% P2 functions phi_i(x,y) N=6

% Basis functions expressed as:
% phi(xi,yi) = c1*xi^2 + c2*yi^2 + c3*xi*yi + c4*xi + c5*yi + c6
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)
% 4 (1/2,0)
% 5 (1/2,1/2)
% 6 (0,1/2)

xi = [0 1 0 1/2 1/2 0];
yi = [0 0 1 0 1/2 1/2];

A1 = sym(zeros(6,6));

for i = 1:6
    A1(i,:) = [xi(i)^2, yi(i)^2, xi(i)*yi(i), xi(i), yi(i), 1];
end

b1 = sym(eye(6,6));

% Coefficients of the basis functions
c1 = linsolve(A1,b1);


syms x y

phi = sym(zeros(6,1)); % expression of the basis function phi

dphidx = sym(zeros(6,1));
dphidy = sym(zeros(6,1));
for j=1:6
    phi(j,1) = c1(1,j)*x*x + ...
               c1(2,j)*y*y + ...
               c1(3,j)*x*y + ...
               c1(4,j)*x + ...
               c1(5,j)*y + ...
               c1(6,j);
    dphidx(j,1) = diff(phi(j,1),x);
    dphidy(j,1) = diff(phi(j,1),y);
end

%% P2 functions psi_j(x,y) M=6

% NOTE:
% This 2nd section is not necessary for cases P1P1, P2P2, P3P3.
% The final matrices can be computed using only the first section as
% phi = psi, hence dphidx = dpsidx and dphidy = dpsidy.
% As a result we also find: mat0X = matX0' and mat0Y = matY0'.

% This 2nd section has been included for consistency with all other
% combinations of basisfunctionsP#P#.m codes.

% Basis functions expressed as:
% psi(xi,yi) = c1*xi^2 + c2*yi^2 + c3*xi*yi + c4*xi + c5*yi + c6
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)
% 4 (1/2,0)
% 5 (1/2,1/2)
% 6 (0,1/2)

xi = [0 1 0 1/2 1/2 0];
yi = [0 0 1 0 1/2 1/2];

A2 = sym(zeros(6,6));

for i = 1:6
    A2(i,:) = [xi(i)^2, yi(i)^2, xi(i)*yi(i), xi(i), yi(i), 1];
end

b2 = sym(eye(6,6));

% Coefficients of the basis functions
c2 = linsolve(A2,b2);


syms x y

psi = sym(zeros(6,1)); % expression of the basis function psi

dpsidx = sym(zeros(6,1));
dpsidy = sym(zeros(6,1));
for j=1:6
    psi(j,1) = c2(1,j)*x*x + ...
               c2(2,j)*y*y + ...
               c2(3,j)*x*y + ...
               c2(4,j)*x + ...
               c2(5,j)*y + ...
               c2(6,j);
    dpsidx(j,1) = diff(psi(j,1),x);
    dpsidy(j,1) = diff(psi(j,1),y);
end

%% 

mat0X = sym(zeros(6,6));
mat0Y = sym(zeros(6,6));
matX0 = sym(zeros(6,6));
matY0 = sym(zeros(6,6));

for i=1:6
    for j=1:6
        mat0X(i,j) = int(int(phi(i,1)*dpsidx(j,1),y,0,-x+1),x,0,1);
        mat0Y(i,j) = int(int(phi(i,1)*dpsidy(j,1),y,0,-x+1),x,0,1);
        matX0(i,j) = int(int(dphidx(i,1)*psi(j,1),y,0,-x+1),x,0,1);
        matY0(i,j) = int(int(dphidy(i,1)*psi(j,1),y,0,-x+1),x,0,1);
    end
end
