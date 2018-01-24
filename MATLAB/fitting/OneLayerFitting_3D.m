%
% One Layer Gaussian fitting for varying incident angle heightfield data
%
% convert 3d heightfield data to slope domain
% fit using (accelerated) EM algorithm
% plot the input and fitted result and calculate the relative l2 error
%

close all
clear
clc

mirror = true;

if mirror
    %     datadir = '/Users/mandy/Github/MultiLayerBsdf/build_clang';
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
else
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output/';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting/';
addpath(datadir,fundir)

trainnum = 1e6;
generatenum = 1e7;
gaussiannumvec = 50; % number of gaussians vector
accelerated = true; % if true uses accelerated em, otherwise uses customized gmcluster
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
alpharange = 1:length(alphavec);
for k = 4
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
        angle = angle(z>=0);
        x = x(z>=0);
        y = y(z>=0);
        z = z(z>=0);
        
        input = [x,y,z,angle];
        
        fitting_halfvector_z13D(datadir,fundir,alpha,input,...
            trainnum, generatenum, gaussiannumvec,xnum, ynum, znum,accelerated);
        
    else% glass case
        
    end% check mirror or glass end
    
end

