function [obj,W,M,R] = fit2D(angle,softinit,W,M,R)
mirror = true;

if mirror
    datadir = '/Users/mandy/Github/MultiLayerBsdf/build/';
    %     datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_2d/angle60/output/';
else
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerGlass_2d/angle60/output/';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting/';

addpath(datadir,fundir)

trainnum = 1e6; % number of data for training
generatenum = 1e6;  % number of data for testing
gaussiannumvec = 5; % number of gaussians vector
accelerated = false; % if true uses accelerated em, othe
maxiter = 500;
tol = 1e-5;
alpha = 0.5;

% close all
cd(datadir)
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
    
    [obj,W,M,R] = fitting_halfvector_z1(datadir,fundir,alpha,angle,input,...
        trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated,maxiter,tol,softinit, W, M, R);
    
    %             thetanum = 100;
    %             angle = 60;
    %             fitting_tanthetah(datadir,fundir,alpha,angle,input,trainnum, ...
    %                 generatenum, gaussiannumvec, thetanum,accelerated,maxiter,tol);
    
else
    
    ior = 1.5;
    Glass_fitting_halfvector_z1(datadir,fundir,alpha,angle,input,...
        trainnum, generatenum, gaussiannumvec, xnum, ynum,ior,accelerated,maxiter,tol);
    
end

end