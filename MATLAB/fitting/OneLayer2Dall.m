%
% One Layer Gaussian fitting for fixed incident angle heightfield data
%
% convert 2d heightfield data to slope domain
% fit using (accelerated) EM algorithm
% plot the input and fitted result and calculate the relative l2 error
%

function [fitting_data, anglesforfit] =  OneLayer2Dall(alpha, gaussiannum)
%close all
%clear
%clc

mirror = true;

owner = 'Feng';

if mirror
    datadir = '/Users/fengxie/work/Github/GaussianData/HeightfieldData/singleLayerSlice/';
    fundir = '/Users/fengxie/work/Github/GaussianClean/MATLAB/fitting/';
end

datadir = strcat(datadir, 'alpha', num2str(alpha), '/');
addpath(datadir,fundir)
disp(datadir);
resultdir = strcat(datadir, num2str(gaussiannum), '/');
status = mkdir(resultdir);
disp(resultdir);


trainnum = 100000; % number of data for training
generatenum = 100000;  % number of data for testing
gaussiannumvec = gaussiannum; % number of gaussians vector
accelerated = true; % if true uses accelerated em, othe
maxiter = 1000;
tol = 1e-5;

%angle = 80;

runcount = 0;
W = [];
M = [];
R = [];
fitting_data = [];

for k = 4
        [x, y, z, angle] = read_data(datadir, alpha,mirror);
        anglevalues = unique(angle);
        anglesforfit = anglevalues;
        anglecount = length(anglevalues);
        %step = floor(anglecount/10)
        step = 1

    %for j = anglecount:-step:1
    for j = 1:step:anglecount
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
            [obj, W, M, R] = fitting_halfvector_z1(resultdir,fundir,alpha,iangle,input,...
                 trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated,maxiter,tol, runcount > 1, W, M, R);
            fitting_data = prep_for_fitting(obj, runcount == 1, j, anglecount, fitting_data);
        else
            
            ior = 1.5;
            Glass_fitting_halfvector_z1(datadir,fundir,alpha,iangle,input,...
                trainnum, generatenum, gaussiannumvec, xnum, ynum,ior,accelerated,maxiter,tol);
            
        end
    end
end
end %function end


function [x, y, z, angle] = read_data(datadir, alpha, mirror)
        close all
        cd(datadir)
        prefix = '3d_output';
        filename = [prefix, 'x_',  num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C1 = textscan(fileID,'%f');
        fclose(fileID);
        
        filename = [prefix, 'y_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C2 = textscan(fileID,'%f');
        fclose(fileID);
        
        filename = [prefix, 'z_',  num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C3 = textscan(fileID,'%f');
        fclose(fileID);
        
        filename = [prefix, 'weight_',  num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C4 = textscan(fileID,'%f');
        fclose(fileID);

        filename = [prefix, 'angle_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C5 = textscan(fileID,'%f');
        fclose(fileID);
        
        x = C1{1};
        y = C2{1};
        z = C3{1};
        weight = C4{1};
        angle = C5{1};
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
end

function [fitting_data] = prep_for_fitting(obj, start, angleindex, anglecount, fitting_data)
    weights = obj.ComponentProportion;
    means = obj.mu;
    cov = obj.Sigma;
    wsize = prod(size(weights));
    msize = prod(size(means));
    if (start)
        s = wsize+msize;
        %fitting_data = zeros(s, anglecount);
        fitting_data = zeros(anglecount, s);
    end
    
    i = angleindex; 
    for j = 1:wsize
        %fitting_data(j, i) = weights(j);
        fitting_data(i, j) = weights(j);
    end
    for j = 1:msize
        %fitting_data(j+wsize, i) = means(j)
        fitting_data(i, j+wsize) = means(j);
    end
end

        
