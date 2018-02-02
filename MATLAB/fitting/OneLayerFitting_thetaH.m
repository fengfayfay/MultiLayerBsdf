%
% One Layer Gaussian fitting for varying incident angle heightfield data
%
% convert 3d heightfield data to tanTheta 
% fit using (accelerated) EM algorithm
% plot the input and fitted result and calculate the relative l2 error
%

close all
clear
clc

mirror = true;

if mirror
    datadir = '/Users/fengxie/work/Github/GaussianData/HeightFieldData/singleLayerUnistack07';
else
    datadir = '/Users/fengxie/Github/GaussianData/GaussianHeightField/SingleLayer/pi:3/output/';
end
fundir = '/Users/fengxie/work/Github/GaussianClean/MATLAB/fitting/';
addpath(datadir,fundir)

trainnum = 1e5;
%trainnum = 1e4;
generatenum = 1e5;
%generatenum = 1e4;
gaussiannumvec = 10; % number of gaussians vector
accelerated = true; % if true uses accelerated em, otherwise uses customized gmcluster
reflectdata = false;
maxiter = 600;
tol = 1e-5;
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
alpharange = 1:length(alphavec);
for k = 5
    alpha = alphavec(k);
    filename = ['3d_outputx_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C1 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputy_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C2 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputz_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C3 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputdepth_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C5 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputangle_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C4 = textscan(fileID,'%f');
    fclose(fileID);
    
    x = C1{1};
    y = C2{1};
    z = C3{1};
    angle = C4{1};
    observe = 6000;
    x = x/observe;
    y = y/observe;
    z = z/observe;

    depth = C5{1};
    disp(length(x));
    
    % plotting parameters
    xnum = 100;
    ynum = 100;
    znum = 90;
    
    
    if mirror
        % for mirror
        if sum(z<0)/length(z)>0.01
            msg = 'More than 1% rays go down, not mirror data.';
            error(msg)
        end
        % only take z>= data
        %angle = angle(depth>1);
        %x = x(depth>1);
        %y = y(depth>1);
        %z = z(depth>1);

        angle = angle(z>=0);
        x = x(z>=0);
        y = y(z>=0);
        z = z(z>=0);
        
        input = [x,y,z,angle];
        disp(length(input));
        
        % randomly permute input data
        input = input(randperm(length(input)),:);
      
        xnum = 200;
        znum = 90;
        zmin = 0;
        zmax = pi * .5;
        trainnum = 1000;
        generatenum = 4000;
        boundary_ratio = .999;
         
        [paramMid, angleValues, rayleighPoly,  range ] = fitting_thetaH_3D(input,...
         xnum, znum, trainnum,generatenum); 
        
        %tan2brdf(rayleighPoly);

    end% check mirror or glass end
    
end


