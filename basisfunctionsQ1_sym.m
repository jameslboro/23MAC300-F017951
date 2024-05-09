% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

% Basis functions expressed as:
% chi(xi,yi) = c1*xi*yi + c2*xi + c3*yi + c4
% Numbering of the nodes:
% 1 (-1,-1)
% 2 ( 1,-1)
% 3 ( 1, 1)
% 4 (-1, 1)

xi = [-1  1 1 -1];
yi = [-1 -1 1  1];

A = sym(zeros(4,4));

for i = 1:4
    j = i;
    A(i,:) = [xi(i)*yi(j), ...
                    xi(i), ...
                    yi(j), ...
                        1];
end

b = sym(eye(4,4));
  
% Coefficients of the basis functions
c = linsolve(A,b);


syms x y

chi = sym(zeros(4,1)); % expression of the basis function chi

dchidx = sym(zeros(4,1));
dchidy = sym(zeros(4,1));
for j=1:4
    chi(j,1) = c(1,j)*x*y + ...
        c(2,j)*x + ...
        c(3,j)*y + ...
        c(4,j);
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end


symDERXX = sym(zeros(4,4));
symDERYY = sym(zeros(4,4));
symDERXY = sym(zeros(4,4));
symMASS = sym(zeros(4,4));

for i = 1:4
    for j = 1:4
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),x,-1,1),y,-1,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),x,-1,1),y,-1,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),x,-1,1),y,-1,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),x,-1,1),y,-1,1);
    end
end
