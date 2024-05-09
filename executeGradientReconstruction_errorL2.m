% James Jarman (j.jarman-20@lboro.ac.uk)
% Loughborough University
% 2024

clear; close all;

%%

numberOfMeshes = 5; % number of meshes
numberOfSettings = 3; % number of combinations of parameters
errorL2 = zeros(numberOfMeshes,numberOfSettings);
hSize = zeros(numberOfMeshes,numberOfSettings);

deltaSettings = [0,0.1,1]; % set the combinations of parameter delta 
alphaSettings = [1,1,2]; % set the combinations of alpha parameter

for i = 1:1:numberOfMeshes
    % Load mesh
    meshElements = 2^i;
    meshIdentifier = strcat(polyCombi,'_',num2str(meshElements),'_',num2str(meshElements));
    filename = strcat('data/gradient/loula/mesh_',meshIdentifier,'.mat');
    load(filename);
    clear filename meshIdentifier;

    for j = 1:1:3
        % Settings for the problem
        settings.delta     = deltaSettings(j); % value of parameter delta
        settings.alpha     = alphaSettings(j); % value of parameter alpha
        settings.plot      = 0;
        settings.plotExact = 0;
        settings.error     = 1;
        
        % Execute main code
        problemDataFile = @gradientDataLoula;
        % problemDataFile = @gradientDataExp;
        [u,error] = gradientReconstruction (mesh,dof,settings,problemDataFile,{});
        
        hSize(i,j) = 1/meshElements;
        errorL2(i,j) = mean(error.L2);
    end
end


% figure
% hold on
% % for i = 1:1:numberOfSettings
% %     plot(-log10(hSize(:,i)),log10(errorL2(:,i)));
% % end

plot(-log10(hSize(:,1)),log10(errorL2(:,1)),'-^', ...
     -log10(hSize(:,2)),log10(errorL2(:,2)),'-diamond', ...
     -log10(hSize(:,3)),log10(errorL2(:,3)),'-pentagram')

title('Q3Q3')
legend('\delta = 0','\delta = 0.1, \alpha = 1','\delta = 1, \alpha = 2')
xlabel('-log_1_0(h)')
ylabel('log_1_0||e||')

errorL2conv = zeros(numberOfMeshes-1,numberOfSettings);

for i = 1:1:numberOfMeshes-1
    for j = 1:1:numberOfSettings
        errorL2conv(i,j) = log(errorL2(i+1,j)/errorL2(i,j))/(log(1/2));
    end
end

% disp(polyCombi)
% errorL2conv
