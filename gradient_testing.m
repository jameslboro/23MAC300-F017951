% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

% x = 0:0.2:1;
% y = 0:0.2:1;
% y = y';
% u = 1/(2*pi^2)*sin(pi*x).*sin(pi*y);
% 
% % surf(x,y,u)
% % xlabel('x')
% % ylabel('y')
% % zlabel('z')
% 
% 
% [ux,uy] = gradient(u,0.2);
% 
% sigma1 = 1/(2*pi)*cos(pi*x).*sin(pi*y);
% 
% e = norm(ux-sigma1,Inf);

%% all points of mesh
clear;

numberOfMeshes = 4;

ex = zeros(numberOfMeshes,1);
ey = zeros(numberOfMeshes,1);
e = zeros(numberOfMeshes,1);

for i = 0:1:numberOfMeshes-1
    hx = 0.5;
    hx = hx*(1/2).^i;
    hy = 0.5;
    hy = hy*(1/2).^i;
    x = 0:hx:1;
    y = 0:hy:1;
    y = y';
    % Exact solution and gradient
    u = 1/(2*pi^2)*sin(pi*x).*sin(pi*y);
    sigma1 = 1/(2*pi)*cos(pi*x).*sin(pi*y);
    sigma2 = 1/(2*pi)*sin(pi*x).*cos(pi*y);

    % Numerical gradient
    [ux,uy] = gradient(u,hx,hy);

    % Error estimates
    xDifference = ux - sigma1;
    yDifference = uy - sigma2;

    xError = zeros(size(xDifference,1),1);
    for k = 1:size(xDifference,1)
        xError(k,1) = norm(xDifference(k,:),Inf);
    end

    yError = zeros(size(yDifference,2),1);
    for k = 1:size(yDifference,2)
        yError(k,1) = norm(yDifference(:,k),Inf);
    end
    
    Difference = [xDifference,yDifference];
    Error = zeros (size(Difference,2),1);
    for k = 1:size(Difference,2)
        Error(k,1) = norm(Difference(:,k),Inf);
    end

    ex(i+1,1) = norm(xError,Inf);
    ey(i+1,1) = norm(yError,Inf);
    e(i+1,1) = norm(Error,Inf);
end

exconv = zeros(length(ex)-1,1);
eyconv = zeros(length(ey)-1,1);

for i = 1:1:size(exconv)
    exconv(i,1) = log(ex(i+1)/ex(i))/(log(1/2));
    eyconv(i,1) = log(ey(i+1)/ey(i))/(log(1/2));
end

%% interior points of mesh (removing 2 rows and columns)

% ex_int = zeros(4,1);
% ey_int = zeros(4,1);
% 
% for i = 0:1:3
%     hx = 0.1;
%     hx = hx*(1/2).^i;
%     hy = 0.1;
%     hy = hy*(1/2).^i;
%     x = 0:hx:1;
%     y = 0:hy:1;
%     y = y';
% 
%     u = 1/(2*pi^2)*sin(pi*x).*sin(pi*y);
% 
%     % figure
%     % surf(x,y,u)
%     % xlabel('x')
%     % ylabel('y')
%     % zlabel('z')
% 
%     [ux,uy] = gradient(u,hx,hy);
% 
%     ux_int = ux(3:length(ux)-2,3:length(ux)-2);
%     uy_int = uy(3:length(uy)-2,3:length(uy)-2);
% 
%     sigma1 = 1/(2*pi)*cos(pi*x).*sin(pi*y);
%     sigma2 = 1/(2*pi)*sin(pi*x).*cos(pi*y);
% 
%     sigma1int = sigma1(3:length(ux)-2,3:length(ux)-2);
%     sigma2int = sigma2(3:length(uy)-2,3:length(uy)-2);
% 
%     ex_int(i+1,1) = norm(ux_int-sigma1int,Inf);
%     ey_int(i+1,1) = norm(uy_int-sigma2int,Inf);
% end
% 
% ex_intconv = zeros(length(ex_int)-1,1);
% ey_intconv = zeros(length(ey_int)-1,1);
% 
% for i = 1:1:size(ex_intconv)
%     ex_intconv(i,1) = ex_int(i)/ex_int(i+1);
%     ey_intconv(i,1) = ey_int(i)/ey_int(i+1);
% end

% %% exterior points of mesh
% 
% % ex_int = zeros(4,1);
% ex_ext = zeros(4,1);
% 
% % ey_int = zeros(4,1);
% ey_ext = zeros(4,1);
% 
% for i = 0:1:3
%     hx = 0.1;
%     hx = hx*(1/2).^i;
%     hy = 0.1;
%     hy = hy*(1/2).^i;
%     x = 0:hx:1;
%     y = 0:hy:1;
%     y = y';
% 
%     u = 1/(2*pi^2)*sin(pi*x).*sin(pi*y);
% 
%     % figure
%     % surf(x,y,u)
%     % xlabel('x')
%     % ylabel('y')
%     % zlabel('z')
% 
%     [ux,uy] = gradient(u,hx,hy);
% 
%     ux_ext = zeros(size(ux));
%     ux_ext(1,:) = ux(1,:); ux_ext(:,length(ux)) = ux(:,length(ux));
%     ux_ext(:,1) = ux(:,1); ux_ext(:,length(ux)) = ux(:,length(ux));
%     % ux_int = uy(3:length(ux)-2,3:length(ux)-2);
%     uy_ext = zeros(size(uy));
%     uy_ext(1,:) = uy(1,:); uy_ext(:,length(uy)) = uy(:,length(uy));
%     uy_ext(:,1) = uy(:,1); uy_ext(:,length(uy)) = uy(:,length(uy));
%     % uy_int = uy(3:length(uy)-2,3:length(uy)-2);
% 
%     sigma1 = 1/(2*pi)*cos(pi*x).*sin(pi*y);
%     sigma2 = 1/(2*pi)*sin(pi*x).*cos(pi*y);
% 
%     sigma1_ext = zeros(size(ux));
%     sigma1_ext(1,:) = sigma1(1,:); sigma1_ext(:,length(sigma1)) = sigma1(:,length(sigma1));
%     sigma1_ext(:,1) = sigma1(:,1); sigma1_ext(:,length(sigma1)) = sigma1(:,length(sigma1));
% 
%     sigma2_ext = zeros(size(uy));
%     sigma2_ext(1,:) = sigma2(1,:); sigma2_ext(:,length(sigma2)) = sigma2(:,length(sigma2));
%     sigma2_ext(:,1) = sigma2(:,1); sigma2_ext(:,length(sigma2)) = sigma2(:,length(sigma2));
% 
% 
%     % sigma1int = sigma1(3:length(ux)-2,3:length(ux)-2);
%     % sigma2int = sigma2(3:length(uy)-2,3:length(uy)-2);
% 
%     % ex_int(i+1,1) = norm(ux_int-sigma1int,Inf);
%     ex_ext(i+1,1) = norm(ux_ext-sigma1_ext,Inf);
% 
%     % ey_int(i+1,1) = norm(uy_int-sigma2int,Inf);
%     ey_ext(i+1,1) = norm(uy_ext-sigma2_ext,Inf);
% end
% 
% % ex_intconv = zeros(length(ex_int)-1,1);
% ex_extconv = zeros(length(ex_ext)-1,1);
% % ey_intconv = zeros(length(ey_int)-1,1);
% ey_extconv = zeros(length(ey_ext)-1,1);
% 
% for i = 1:1:size(ex_extconv)
%     % ex_intconv(i,1) = ex_int(i)/ex_int(i+1);
%     ex_extconv(i,1) = ex_ext(i)/ex_ext(i+1);
%     % ey_intconv(i,1) = ey_int(i)/ey_int(i+1);
%     ey_extconv(i,1) = ey_ext(i)/ey_ext(i+1);
% end
