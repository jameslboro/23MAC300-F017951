% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

% Basis functions expressed as:
% chi(xi,yi) = c1*xi + c2*yi + c3
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)

xi = [0 1 0];
yi = [0 0 1];

A = sym(zeros(3,3));

for i = 1:3
    A(i,:) = [xi(i), yi(i), 1];
end

b = sym(eye(3,3));

% Coefficients of the basis functions
c = linsolve(A,b);


syms x y

chi = sym(zeros(3,1)); % expression of the basis function chi

dchidx = sym(zeros(3,1));
dchidy = sym(zeros(3,1));
for j=1:3
    chi(j,1) = c(1,j)*x + ...
               c(2,j)*y + ...
               c(3,j);
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end


symDERXX = sym(zeros(3,3));
symDERYY = sym(zeros(3,3));
symDERXY = sym(zeros(3,3));
symMASS = sym(zeros(3,3));

for i = 1:3
    for j = 1:3
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),y,0,-x+1),x,0,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),y,0,-x+1),x,0,1);
    end
end
