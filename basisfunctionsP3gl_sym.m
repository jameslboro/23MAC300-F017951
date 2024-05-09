% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

% Basis functions expressed as:
% chi(xi,yi) = c1*xi^3 + c2*yi^3 + c3*xi^2*yi + c4*xi*yi^2 + ...
%              c5*xi^2 + c6*yi^2 + c7*xi*yi + c8*xi + c9*yi + c10
% Numbering of the nodes:
%  1 (0,0)
%  2 (1,0)
%  3 (0,1)
%  4 ((5-sqrt(5))/10,0)
%  5 ((5+sqrt(5))/10,0)
%  6 (1-((5-sqrt(5))/10),5-sqrt(5)/10)
%  7 (1-((5+sqrt(5))/10),5+sqrt(5)/10)
%  8 (0,(5+sqrt(5))/10)
%  9 (0,(5-sqrt(5))/10)
% 10 (1/3,1/3)

sqrt5 = sym(sqrt(5)); % create symbolic variable for sqrt(5)

xi = [0, 1, 0, (sqrt5-1)/(2*sqrt5), (sqrt5+1)/(2*sqrt5), ...
    (sqrt5+1)/(2*sqrt5), (sqrt5-1)/(2*sqrt5), 0, 0, 1/3];
yi = [0, 0, 1, 0, 0, (sqrt5-1)/(2*sqrt5), (sqrt5+1)/(2*sqrt5), ...
    (sqrt5+1)/(2*sqrt5), (sqrt5-1)/(2*sqrt5), 1/3];

A = sym(zeros(10,10));

for i = 1:10
    A(i,:) = [xi(i)^3, yi(i)^3, xi(i)^2*yi(i), xi(i)*yi(i)^2, xi(i)^2, yi(i)^2, xi(i)*yi(i), xi(i), yi(i), 1];
end

b = sym(eye(10,10));

% Coefficients of the basis functions
c = linsolve(A,b);


syms x y

chi = sym(zeros(10,1)); % expression of the basis function chi

dchidx = sym(zeros(10,1));
dchidy = sym(zeros(10,1));
for j = 1:10
    chi(j,1) = c(1,j)*x^3 + ...
               c(2,j)*y^3 + ...
               c(3,j)*x^2*y + ...
               c(4,j)*x*y^2 + ...
               c(5,j)*x^2 + ...
               c(6,j)*y^2 + ...
               c(7,j)*x*y + ...
               c(8,j)*x + ...
               c(9,j)*y + ...
               c(10,j);
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end


symDERXX = sym(zeros(10,10));
symDERYY = sym(zeros(10,10));
symDERXY = sym(zeros(10,10));
symMASS = sym(zeros(10,10));

for i = 1:10
    for j = 1:10
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),y,0,-x+1),x,0,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),y,0,-x+1),x,0,1);
    end
end
