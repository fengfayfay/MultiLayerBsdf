%
% One Layer Gaussian fitting for fixed incident angle heightfield data
%
% convert 2d heightfield data to slope domain
% fit using (accelerated) EM algorithm
% plot the input and fitted result and calculate the relative l2 error
% fit the series of 2D gaussians over incident angle using polynomials
%

function [mean_0, mean_1, cv_0, cv_1, energyRatios, anglevalues] =  OneLayer2Dall(alpha, gaussiannum, bounce)
    %close all
    %clear
    %clc

    datadir = '/Users/fengxie/work/Github/GaussianData/HeightfieldData/singleLayerSlice/';
    fundir = '/Users/fengxie/work/Github/GaussianClean/MATLAB/fitting/';

    datadir = strcat(datadir, 'alpha', num2str(alpha), '/');
    addpath(datadir,fundir)
    disp(datadir);
    resultdir = strcat(datadir, num2str(gaussiannum), '/');
    status = mkdir(resultdir);
    resultdir = strcat(resultdir, 'depth', num2str(bounce), '/');
    status = mkdir(resultdir);
    disp(resultdir);


    trainnum = 20000; % number of data for training
    generatenum = 20000;  % number of data for testing
    gaussiannumvec = gaussiannum; % number of gaussians vector
    accelerated = true; % if true uses accelerated em, othe
    maxiter = 1000;
    tol = 1e-5;


    runcount = 0;
    W = [];
    M = [];
    R = [];
    fitting_data = [];

    [ox, oangle, x, y, z, angle] = read_data(datadir, alpha, bounce);
    [anglevalues, energyRatios] = selectValidAnglesForFitting(ox, oangle, x, angle);
    anglecount = length(anglevalues);
    step = 1;
    for j = 1:step:anglecount
        iangle = anglevalues(j);
        ix = x(abs(angle - iangle) < .0001);
        iy = y(abs(angle - iangle) < .0001);
        iz = z(abs(angle - iangle) < .0001);

        maxSampleNum = floor(length(ix)/2)
        trainnum = min(trainnum, maxSampleNum)
        generatenum = trainnum
        input = [ix,iy,iz];
        
        xnum = 100;
        ynum = 100;
        runcount = runcount + 1; 
        
        [obj, W, M, R] = fitting_halfvector_z1(resultdir,fundir,alpha,iangle,input,...
                 trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated,maxiter,tol, runcount > 1, W, M, R, energyRatios(j));
        fitting_data = prep_for_fitting(obj, runcount == 1, j, anglecount, fitting_data);
    end

    mean_0 = fitting_data(:,1)
    assert(length(mean_0) == anglecount);
    mean_1 = fitting_data(:,2)
    assert(length(mean_1) == anglecount);
    cv_0 = fitting_data(:,3)
    assert(length(cv_0) == anglecount);
    cv_1 = fitting_data(:,4)
    assert(length(cv_1) == anglecount);
   
    filename = [resultdir, 'scatter_', num2str(alpha), '.txt'];
    file = fopen(filename, 'w');
    coefCount = 5;
    fprintf(file, "%d\n", coefCount);
    p = polyfit(anglevalues, energyRatios, coefCount - 1);
    exportPolyFit(file, p);    
    disp(p);
    for i = 1:4
        f = fitting_data(:,i);
        assert(length(f) == anglecount);
        p = polyfit(anglevalues, f, coefCount - 1);
        exportPolyFit(file, p);    
        disp(p);
    end
    fclose(file);
end %function end


function exportPolyFit(file, p)
    l = length(p);
    for i = 1:l
        fprintf(file, '%4.6f ', p(i));
    end
    fprintf(file, '\n');
end

function [anglevalues, energyRatios] = selectValidAnglesForFitting(ox, oangle, x, angle)
    anglecandidates = unique(angle);
    anglecount = length(anglecandidates);

    anglevalues = zeros(anglecount,1);
    energyRatios = zeros(anglecount,1);
    validCount = 0;
    for j = 1:anglecount
        iangle = anglecandidates(j);
        iox = ox(abs(oangle - iangle) < .0001);
        ix = x(abs(angle - iangle) < .0001);
        lenIX = length(ix)
        lenIOX = length(iox)
        energy = lenIX/lenIOX
        %if energy > .02
        if energy > 0.001 
            validCount = validCount + 1;
            energyRatios(validCount) = energy;
            anglevalues(validCount) = iangle;
        end 
    end
    energyRatios = energyRatios(1:validCount)
    anglevalues = anglevalues(1:validCount)
end
    
       
        
        
    


function [ox, oangle, x, y, z, angle] = read_data(datadir, alpha, bounce)
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
        
        filename = [prefix, 'depth_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C6 = textscan(fileID,'%d');
        fclose(fileID);
        
        x = C1{1};
        y = C2{1};
        z = C3{1};
        weight = C4{1};
        angle = C5{1};
        depth = C6{1};
        observe = 6000;
        
        ox = x(z>=0);
        oangle = angle(z>=0);

        x = x/observe;
        y = y/observe;
        z = z/observe;

        if bounce == 1
            x = x(depth == 1);
            y = y(depth == 1);
            z = z(depth == 1);
            angle = angle(depth == 1);
        else  
            x = x(depth > 1);
            y = y(depth > 1);
            z = z(depth > 1);
            angle = angle(depth > 1);
        end
        angle = angle(z>=0);
        x = x(z>=0);
        y = y(z>=0);
        z = z(z>=0);
end

function [fitting_data] = prep_for_fitting(obj, start, angleindex, anglecount, fitting_data)
    weights = obj.ComponentProportion;
    means = obj.mu;
    cov = obj.Sigma;
    %disp(cov);


    csize = 2;
    msize = prod(size(means));
    if (start)
        s = 4;
        fitting_data = zeros(anglecount, s);
    end
    
    i = angleindex; 
    fitting_data(i, 1) = means(1);
    fitting_data(i, 2) = means(2);
    fitting_data(i, 3) = cov(1, 1);
    fitting_data(i, 4) = cov(2, 2);
end

        
