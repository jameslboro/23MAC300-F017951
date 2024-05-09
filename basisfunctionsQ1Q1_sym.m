% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%% Q1 functions phi_i(x,y) N=4

% Basis functions expressed as:
% phi(xi,yi) = c1*xi*yi + c2*xi + c3*yi + c4
% Numbering of the nodes:
% 1 (-1,-1)
% 2 ( 1,-1)
% 3 ( 1, 1)
% 4 (-1, 1)

xi = [-1  1 1 -1];
yi = [-1 -1 1  1];

A1 = sym(zeros(4,4));

for i = 1:4
    A1(i,:) = [xi(i)*yi(i), ...
                     xi(i), ...
                     yi(i), ...
                         1];
end

b1 = sym(eye(4,4));
  
% Coefficients of the basis functions
c1 = linsolve(A1,b1);


syms x y

phi = sym(zeros(4,1)); % expression of the basis function phi

dphidx = sym(zeros(4,1));
dphidy = sym(zeros(4,1));
for j = 1:4
    phi(j,1) = c1(1,j)*x*y + ...
               c1(2,j)*x + ...
               c1(3,j)*y + ...
               c1(4,j);
    dphidx(j,1) = diff(phi(j,1),x);
    dphidy(j,1) = diff(phi(j,1),y);
end

%% Q1 functions psi_j(x,y) N=4

% NOTE:
% This 2nd section is not necessary for cases P1P1, P2P2, P3P3.
% The final matrices can be computed using only the first section as
% phi = psi, hence dphidx = dpsidx and dphidy = dpsidy.
% As a result we also find: mat0X = matX0' and mat0Y = matY0'.

% This 2nd section has been included for consistency with all other
% combinations of basisfunctionsQ#Q#.m codes.

% Basis functions expressed as:
% psi(xi,yi) = c1*xi*yi + c2*xi + c3*yi + c4
% Numbering of the nodes:
% 1 (-1,-1)
% 2 ( 1,-1)
% 3 ( 1, 1)
% 4 (-1, 1)

xi = [-1  1 1 -1];
yi = [-1 -1 1  1];

A2 = sym(zeros(4,4));

for i = 1:4
    A2(i,:) = [xi(i)*yi(i), ...
                     xi(i), ...
                     yi(i), ...
                         1];
end

b2 = sym(eye(4,4));
  
% Coefficients of the basis functions
c2 = linsolve(A2,b2);


syms x y

psi = sym(zeros(4,1)); % expression of the basis function psi

dpsidx = sym(zeros(4,1));
dpsidy = sym(zeros(4,1));
for j = 1:4
    psi(j,1) = c2(1,j)*x*y + ...
               c2(2,j)*x + ...
               c2(3,j)*y + ...
               c2(4,j);
    dpsidx(j,1) = diff(psi(j,1),x);
    dpsidy(j,1) = diff(psi(j,1),y);
end

%%

mat0X = sym(zeros(4,4));
mat0Y = sym(zeros(4,4));
matX0 = sym(zeros(4,4));
matY0 = sym(zeros(4,4));

for i = 1:4
    for j = 1:4
        mat0X(i,j) = int(int(phi(i,1)*dpsidx(j,1),x,-1,1),y,-1,1);
        mat0Y(i,j) = int(int(phi(i,1)*dpsidy(j,1),x,-1,1),y,-1,1);
        matX0(i,j) = int(int(dphidx(i,1)*psi(j,1),x,-1,1),y,-1,1);
        matY0(i,j) = int(int(dphidy(i,1)*psi(j,1),x,-1,1),y,-1,1);
    end
end
