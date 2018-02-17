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

owner = 'Feng';
gaussiannumvec = 20;

if mirror
    if (strcmp(owner, 'Mandy'))
       datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
    else
       datadir = '/Users/fengxie/work/GitHub/GaussianData/HeightfieldData/singleLayer08/';
    end
else
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output/';
end


if (strcmp(owner, 'Mandy'))
    fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting/';
else
    fundir = '/Users/fengxie/work/GitHub/GaussianClean/MATLAB/fitting/';
end

addpath(datadir,fundir)
disp(datadir);

if (strcmp(owner, 'Mandy'))
    pbrtbuild = '/Users/mandy/Github/MultiLayerBsdf/build/';
else
    pbrtbuild = strcat(datadir, num2str(gaussiannumvec), '/');
    status = mkdir(pbrtbuild);
    disp(pbrtbuild);
end

trainnum = 1e6 * 2;
generatenum = 1e7;
accelerated = true; % if true uses accelerated em, otherwise uses customized gmcluster
extenddata = 1;
extendratio = .5;   % extend ratio on both ends (0 means no reflection)
maxiter = 2000;
tol = 1e-5;
alphavec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
alpharange = 1:length(alphavec);

W = [];
M = [];
R = [];

for k = 3
    alpha = .8;
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
%     znum = (reflectdata + 1) * 90;
    znum = round((extendratio * 2 + 1) * 90);

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
        
        % randomly permute input data
        input = input(randperm(length(input)),:);
        
%         obj = fitting_halfvector_z13D(datadir,fundir,alpha,input,...
%             trainnum, generatenum, gaussiannumvec,xnum, ynum, znum,accelerated,reflectdata,maxiter,tol,false,W,M,R);
        obj = fitting_halfvector_z13D_extend(datadir,fundir,alpha,input,...
            trainnum, generatenum, gaussiannumvec,xnum, ynum, znum,accelerated,extendratio,maxiter,tol,false,W,M,R);
        
        gm2pbrtinput(pbrtbuild,obj);
        
        % check brdf*cos plot and energy conservation
        gm2brdf(obj,3,alpha,extendratio, 89);

        plotGMM(obj, 10);
        
    else% glass case
        
    end% check mirror or glass end
    
end
