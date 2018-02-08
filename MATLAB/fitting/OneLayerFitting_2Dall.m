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

mirror = true;

owner = 'Feng';

if mirror
    if (strcmp(owner,'Mandy'))
        datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_2d/angle60/output/';
    else
        datadir = '/Users/fengxie/work/Github/GaussianData/HeightfieldData/singleLayerUnistack07/';
    end
else
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerGlass_2d/angle60/output/';
end

if (strcmp(owner, 'Mandy'))
    fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting/';
else
    fundir = '/Users/fengxie/work/Github/GaussianClean/MATLAB/fitting/';
end

addpath(datadir,fundir)

trainnum = 1e4; % number of data for training
generatenum = 1e4;  % number of data for testing
gaussiannumvec = 5; % number of gaussians vector
accelerated = true; % if true uses accelerated em, othe
maxiter = 1000;
tol = 1e-5;
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
anglevec = [0:1:90];
alpharange = 1:length(alphavec);
anglerange = 1:length(anglevec);

runcount = 0;
W = [];
M = [];
R = [];
isigma = 0;

for k = 5
        close all
        cd(datadir)
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
        
        filename = ['3d_outputweight_', num2str(alpha),'.txt'];
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
        weight = C5{1};
        angle = C4{1};
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
            angle = angle(z>=0);
            x = x(z>=0);
            y = y(z>=0);
            z = z(z>=0);
        end

        anglevalues = unique(angle);
        anglecount = length(anglevalues);
        disp(anglecount)
        step = floor(anglecount/10)
        disp(step)

    for j = anglecount:-step:1
        
        iangle = anglevalues(j);
        ix = x(abs(angle - iangle) < .0001);
        iy = y(abs(angle - iangle) < .0001);
        iz = z(abs(angle - iangle) < .0001);

        disp(length(ix));
            
        input = [ix,iy,iz];
        
        xnum = 100;
        ynum = 100;
        runcount = runcount + 1; 
        
        if mirror
             [W, M, R, isigma] = fitting_halfvector_z1(datadir,fundir,alpha,iangle,input,...
                 trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated,maxiter,tol, runcount > 1, W, M, R, isigma);
        else
            
            ior = 1.5;
            Glass_fitting_halfvector_z1(datadir,fundir,alpha,iangle,input,...
                trainnum, generatenum, gaussiannumvec, xnum, ynum,ior,accelerated,maxiter,tol);
            
        end
    end
 end

