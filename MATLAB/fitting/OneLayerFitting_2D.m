%
% One Layer Gaussian fitting for fixed incident angle heightfield data
%
% convert 2d heightfield data to slope domain 
% fit using (accelerated) EM algorithm 
% plot the input and fitted result and calculate the relative l2 error
%

close all
clear
clc

mirror = false;

if mirror
    % datadir = '/Users/mandy/Github/MultiLayerBsdf/build_clang/';
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_2d/angle60/output/';
else
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerGlass_2d/angle60/output/';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting/';

addpath(datadir,fundir)

trainnum = 1e6; % number of data for training
generatenum = 1e6;  % number of data for testing
gaussiannumvec = 5; % number of gaussians vector
accelerated = true; % if true uses accelerated em, otherwise uses customized gmcluster
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
anglevec = [0, 10, 20, 30, 40, 50, 60, 70, 80 ,89];
alpharange = 1:length(alphavec);
anglerange = 1:length(anglevec);

for j = 7
    angle = anglevec(j);
    for k = 4
        cd(datadir)
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
        
        if mirror
            if sum(z<0)/length(z)>0.01
                msg = 'More than 1% rays go down, not mirror data.';
                error(msg)
            end
            % only take z>=0 data
            x = x(z>=0);
            y = y(z>=0);
            z = z(z>=0);
        end
        
        input = [x,y,z];
        
        xnum = 100;
        ynum = 100;

        if mirror
            
            fitting_halfvector_z1(datadir,fundir,alpha,angle,input,...
                trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated);
            
        else
            
            ior = 1.5;
            Glass_fitting_halfvector_z1(datadir,fundir,alpha,angle,input,...
                trainnum, generatenum, gaussiannumvec, xnum, ynum,ior,accelerated);
            
        end
    end
end

