% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

% Basis functions expressed as:
% chi(xi,yi) = c1*xi^2 + c2*yi^2 + c3*xi*yi + c4*xi + c5*yi + c6
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)
% 4 (1/2,0)
% 5 (1/2,1/2)
% 6 (0,1/2)

xi = [0 1 0 1/2 1/2 0];
yi = [0 0 1 0 1/2 1/2];

A = sym(zeros(6,6));

for i = 1:6
    A(i,:) = [xi(i)^2, yi(i)^2, xi(i)*yi(i), xi(i), yi(i), 1];
end

b = sym(eye(6,6));

% Coefficients of the basis functions
c = linsolve(A,b);


syms x y

chi = sym(zeros(6,1)); % expression of the basis function chi

dchidx = sym(zeros(6,1));
dchidy = sym(zeros(6,1));
for j=1:6
    chi(j,1) = c(1,j)*x*x + ...
               c(2,j)*y*y + ...
               c(3,j)*x*y + ...
               c(4,j)*x + ...
               c(5,j)*y + ...
               c(6,j);
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end


symDERXX = sym(zeros(6,6));
symDERYY = sym(zeros(6,6));
symDERXY = sym(zeros(6,6));
symMASS = sym(zeros(6,6));

for i = 1:6
    for j = 1:6
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),y,0,-x+1),x,0,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),y,0,-x+1),x,0,1);
    end
end
