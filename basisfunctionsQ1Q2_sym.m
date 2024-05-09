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

%% Q2 functions psi_j(x,y) N=9

% Basis functions expressed as:
% psi(xi,yi) = c1*xi^2*yi^2 + c2*xi^2*yi + c3*xi*yi^2 + c4*xi^2 + ...
%              c5*xi*yi + c6*yi^2 + c7*xi + c8*yi + c9
% Numbering of the nodes:
% 1 (-1,-1)
% 2 ( 1,-1)
% 3 ( 1, 1)
% 4 (-1, 1)
% 5 ( 0,-1)
% 6 ( 1, 0)
% 7 ( 0, 1)
% 8 (-1, 0)
% 9 ( 0, 0)

xi = [-1  1 1 -1  0 1 0 -1 0];
yi = [-1 -1 1  1 -1 0 1  0 0];

A2 = sym(zeros(9,9));

for i = 1:9
    A2(i,:) = [xi(i)^2*yi(i)^2, ...
                xi(i)^2*yi(i), ...
                xi(i)*yi(i)^2, ...
                      xi(i)^2, ...
                  xi(i)*yi(i), ...
                      yi(i)^2, ...
                        xi(i), ...
                        yi(i), ...
                            1];
end

b2 = sym(eye(9,9));
  
% Coefficients of the basis functions
c2 = linsolve(A2,b2);


syms x y

psi = sym(zeros(9,1)); % expression of the basis function psi

dpsidx = sym(zeros(9,1));
dpsidy = sym(zeros(9,1));
for j = 1:9
    psi(j,1) = c2(1,j)*x^2*y^2 + ...
               c2(2,j)*x^2*y + ...
               c2(3,j)*x*y^2 + ...
               c2(4,j)*x^2 + ...
               c2(5,j)*x*y + ...
               c2(6,j)*y^2 + ...
               c2(7,j)*x + ...
               c2(8,j)*y + ...
               c2(9,j);
    dpsidx(j,1) = diff(psi(j,1),x);
    dpsidy(j,1) = diff(psi(j,1),y);
end

%% 

mat0X = sym(zeros(4,9));
mat0Y = sym(zeros(4,9));
matX0 = sym(zeros(4,9));
matY0 = sym(zeros(4,9));

for i = 1:4
    for j = 1:9
        mat0X(i,j) = int(int(phi(i,1)*dpsidx(j,1),x,-1,1),y,-1,1);
        mat0Y(i,j) = int(int(phi(i,1)*dpsidy(j,1),x,-1,1),y,-1,1);
        matX0(i,j) = int(int(dphidx(i,1)*psi(j,1),x,-1,1),y,-1,1);
        matY0(i,j) = int(int(dphidy(i,1)*psi(j,1),x,-1,1),y,-1,1);
    end
end
