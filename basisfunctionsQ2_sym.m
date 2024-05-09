% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

% Basis functions expressed as:
% chi(xi,yi) = c1*xi^2*yi^2 + c2*xi^2*yi + c3*xi*yi^2 + c4*xi^2 + ...
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

A = sym(zeros(9,9));

for i = 1:9
    j = i;
    A(i,:) = [xi(i)^2*yi(j)^2, ...
                xi(i)^2*yi(j), ...
                xi(i)*yi(j)^2, ...
                      xi(i)^2, ...
                  xi(i)*yi(j), ...
                      yi(j)^2, ...
                        xi(i), ...
                        yi(j), ...
                            1];
end

b = sym(eye(9,9));
  
% Coefficients of the basis functions
c = linsolve(A,b);


syms x y

chi = sym(zeros(9,1)); % expression of the basis function chi

dchidx = sym(zeros(9,1));
dchidy = sym(zeros(9,1));
for j = 1:9
    chi(j,1) = c(1,j)*x^2*y^2 + ...
               c(2,j)*x^2*y + ...
               c(3,j)*x*y^2 + ...
               c(4,j)*x^2 + ...
               c(5,j)*x*y + ...
               c(6,j)*y^2 + ...
               c(7,j)*x + ...
               c(8,j)*y + ...
               c(9,j);
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end


symDERXX = sym(zeros(9,9));
symDERYY = sym(zeros(9,9));
symDERXY = sym(zeros(9,9));
symMASS = sym(zeros(9,9));

for i = 1:9
    for j = 1:9
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),x,-1,1),y,-1,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),x,-1,1),y,-1,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),x,-1,1),y,-1,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),x,-1,1),y,-1,1);
    end
end
