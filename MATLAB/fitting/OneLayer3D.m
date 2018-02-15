function [obj] = OneLayer3D(alpha, gaussiannum, extend)
%
% One Layer Gaussian fitting for varying incident angle heightfield data
%
% convert 3d heightfield data to slope domain
% fit using (accelerated) EM algorithm
% plot the input and fitted result and calculate the relative l2 error
%

close all

mirror = true;

owner = 'Feng';
gaussiannumvec = gaussiannum;

if (strcmp(owner, 'Mandy'))
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
    fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting/';
else
    datadir = '/Users/fengxie/work/GitHub/GaussianData/HeightfieldData/singleLayer/';
    fundir = '/Users/fengxie/work/GitHub/GaussianClean/MATLAB/fitting/';
end

datadir = strcat(datadir, 'alpha', num2str(alpha), '/');
addpath(datadir,fundir)
disp(datadir);

pbrtbuild = strcat(datadir, num2str(gaussiannumvec), '/extend', num2str(extend), '/');
status = mkdir(pbrtbuild);
disp(pbrtbuild);

trainnum = 1e6 * 2;
generatenum = 1e7;
accelerated = true; % if true uses accelerated em, otherwise uses customized gmcluster
extenddata = extend > 0;
extendratio = extend;   % extend ratio on both ends (0 means no reflection)

maxiter = 2000;
tol = 1e-5;

W = [];
M = [];
R = [];

for k = 3
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
        
        gm2pbrtinput(pbrtbuild,obj,extenddata);
        
        % check brdf*cos plot and energy conservation
        for j = 1:8:88
            gm2brdf(obj,3,alpha,extendratio,j);
        end;

        plotGMM(obj, gaussiannum);
        
    else% glass case
        
    end% check mirror or glass end
    
end
