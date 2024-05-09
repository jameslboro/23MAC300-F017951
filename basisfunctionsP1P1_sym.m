% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%% P1 functions phi_i(x,y) N=3

% Basis functions expressed as:
% phi(xi,yi) = c1*xi + c2*yi + c3
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)

xi = [0 1 0];
yi = [0 0 1];

A1 = sym(zeros(3,3));

for i = 1:3
    A1(i,:) = [xi(i), yi(i), 1];
end

b1 = sym(eye(3,3));

% Coefficients of the basis functions
c1 = linsolve(A1,b1);


syms x y

phi = sym(zeros(3,1)); % expression of the basis function phi

dphidx = sym(zeros(3,1));
dphidy = sym(zeros(3,1));
for j = 1:3
    phi(j,1) = c1(1,j)*x + ...
               c1(2,j)*y + ...
               c1(3,j);
    dphidx(j,1) = diff(phi(j,1),x);
    dphidy(j,1) = diff(phi(j,1),y);
end

%% P1 functions psi_j(x,y) M=3 

% NOTE:
% This 2nd section is not necessary for cases P1P1, P2P2, P3P3.
% The final matrices can be computed using only the first section as
% phi = psi, hence dphidx = dpsidx and dphidy = dpsidy.
% As a result we also find: mat0X = matX0' and mat0Y = matY0'.

% This 2nd section has been included for consistency with all other
% combinations of basisfunctionsP#P#.m codes.

% Basis functions expressed as:
% psi(xi,yi) = c1*xi + c2*yi + c3
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)

xi = [0 1 0];
yi = [0 0 1];

A2 = sym(zeros(3,3));

for i = 1:3
    A2(i,:) = [xi(i), yi(i), 1];
end

b2 = sym(eye(3,3));

% Coefficients of the basis functions
c2 = linsolve(A2,b2);


syms x y

psi = sym(zeros(3,1)); % expression of the basis function psi

dpsidx = sym(zeros(3,1));
dpsidy = sym(zeros(3,1));
for j=1:3
    psi(j,1) = c2(1,j)*x + ...
               c2(2,j)*y + ...
               c2(3,j);
    dpsidx(j,1) = diff(psi(j,1),x);
    dpsidy(j,1) = diff(psi(j,1),y);
end

%% 

mat0X = sym(zeros(3,3));
mat0Y = sym(zeros(3,3));
matX0 = sym(zeros(3,3));
matY0 = sym(zeros(3,3));

for i = 1:3
    for j = 1:3
        mat0X(i,j) = int(int(phi(i,1)*dpsidx(j,1),y,0,-x+1),x,0,1);
        mat0Y(i,j) = int(int(phi(i,1)*dpsidy(j,1),y,0,-x+1),x,0,1);
        matX0(i,j) = int(int(dphidx(i,1)*psi(j,1),y,0,-x+1),x,0,1);
        matY0(i,j) = int(int(dphidy(i,1)*psi(j,1),y,0,-x+1),x,0,1);
    end
end
