% One Layer Gaussian fitting
close all
dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror/pi:3/output';
cd(dir)

angle = 60;
trainnum = 1e6;
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
range = 1:length(alphavec);
figure
for k = 6
    alpha = alphavec(k);
    filename = [num2str(angle), 'outputx_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C1 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = [num2str(angle),'outputy_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C2 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = [num2str(angle),'outputz_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C3 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = [num2str(angle),'outputweight_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C4 = textscan(fileID,'%f');
    fclose(fileID);
    
    x = C1{1};
    y = C2{1};
    z = C3{1};
    weight = C4{1};
    observe = 6000;
    x = x/observe;
    y = y/observe;
    z = z/observe;
    epsilon = 1e-6;
    
    if sum(z<0)/length(z)>0.01
        msg = 'More than 1% rays go down, not mirror data.';
        error(msg)
    end
    
    testafter = 1e6;
    generatenum = length(x) - testafter;
    
    cd('/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting');
    
%     % fitting using phi mu
%     gaussiannumvec = 2:5;
%     phinum = 400;
%     munum = 100;
%     fitting_phimu(dir,alpha,angle,x,y,z,weight,epsilon,testafter,trainnum,...
%       generatenum, gaussiannumvec, phinum, munum);
    
%     % fitting using tan(halfvector)
%     gaussiannumvec = 1:5;
%     incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
%     thetanum = 1e3;
%     thetarange = 3;
%     fitting_tanhalfvector(dir,alpha,angle,x,y,z,weight,testafter,trainnum,...
%         generatenum, gaussiannumvec, incident, thetarange, thetanum);
    
%     % fitting using phi(halfvector), theta(halfvector)
%     gaussiannumvec = 1:5;
%     incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
%     phinum = 400;
%     thetanum = 200;
%     fitting_halfvector(dir,alpha,angle,x,y,z,weight,epsilon, testafter, trainnum, ...
%     generatenum, gaussiannumvec, incident, phinum, thetanum);

%     % fitting using x,y of halfvector
%     gaussiannumvec = 1:5;
%     incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
%     xnum = 100;
%     ynum = 100;
%     fitting_halfvector_projected(dir,alpha,angle,x,y,z,weight,testafter, trainnum, ...
%     generatenum, gaussiannumvec, incident, xnum, ynum);

% fitting using x,y of halfvector
gaussiannumvec = 1:5;
incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
xnum = 100;
ynum = 100;
fitting_halfvector_projected_compensate(dir,alpha,angle,x,y,z,weight,...
    testafter, trainnum, generatenum, gaussiannumvec, incident, xnum, ynum)
end

