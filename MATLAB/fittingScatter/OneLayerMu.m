%
% One Layer Gaussian fitting for fixed incident angle heightfield data
%
% convert 2d heightfield data to slope domain
% fit using (accelerated) EM algorithm
% plot the input and fitted result and calculate the relative l2 error
% fit the series of 2D gaussians over incident angle using polynomials
% 
%   alpha: roughness
%   gaussiannum:  number of gaussians, 1 is recommended
%   bounce:  depth of 1 or greater than 1
%

function [mean_0, mean_1, cv_0, cv_1, energyRatios, anglevalues, errs, mus] =  OneLayerMu(ni, no, alpha, polyDegree, bounce, clampvalue)
    %close all
    %clear
    %clc

    %datadir = '/Users/fengxie/work/Github/GaussianData/HeightfieldData/singleLayerGlassSlice/';
    datadir = '/Users/fengxie/work/Github/GaussianData/singleLayerMirrorSlice/';
    %datadir = '/Users/fengxie/work/Github/GaussianData/singleLayerGlassWaterSlice/';
    %
    %datadir = '/Users/fengxie/work/Github/GaussianData/singleLayerGlassFresnel/';
    %datadir = '/Users/fengxie/work/Github/GaussianData/singleLayerMirrorSlice/';
    
    fundir = '/Users/fengxie/work/Github/MultiLayerBsdf/MATLAB/fittingScatter/';

    alphastr = sprintf('%.2f', alpha)
    disp(alphastr)

    resultdir = strcat(datadir, 'results/alpha', alphastr, '/', 'poly', num2str(polyDegree), '/');
    resultdir = strcat(resultdir, 'depth', num2str(bounce), '/');
    status = mkdir(resultdir);
    disp(resultdir);

    datadir = strcat(datadir, 'alpha', alphastr, '/');
    addpath(datadir,fundir)
    disp(datadir);

    trainnum = 20000; % number of data for training
    generatenum = 20000;  % number of data for testing
    gaussiannumvec = 1; % number of gaussians vector
    accelerated = true; % if true uses accelerated em, othe
    maxiter = 1000;
    tol = 1e-5;

    targetNames = ["energy ratio", "mean_x", "mean_y", "\Sigma_x", "\Sigma_y"]; 
    targetFileNames = ["energyRatio", "mean_x", "mean_y", "cov_x", "cov_y"]; 
    %targetName = targetNames{1}
    %filename = [resultdir, 'polyfit_', targetName]


    runcount = 0;
    W = [];
    M = [];
    R = [];
    fitting_data = [];

    [ox, oangle, x, y, z, angle] = read_data(datadir, alpha, bounce);
    [anglevalues, energyRatios, mus] = selectValidAnglesForFitting(ox, oangle, x, angle);
    disp(anglevalues);
    anglecount = length(anglevalues);
    disp(anglecount);

    step = 1;

    errs = zeros(anglecount, 1);
    for j = 1:step:anglecount
        iangle = anglevalues(j);
        ix = x(abs(angle - iangle) < .0001);
        %disp(size(ix))
        iy = y(abs(angle - iangle) < .0001);
        iz = z(abs(angle - iangle) < .0001);

        maxSampleNum = floor(length(ix)/2)
        trainnum = min(trainnum, maxSampleNum)
        generatenum = trainnum
        input = [ix,iy,iz];
        disp(size(input))
        
        xnum = 100;
        ynum = 100;
        runcount = runcount + 1; 
         
        [obj, err, W, M, R] = fitting_slopedomain(resultdir,fundir, ni, no, clampvalue, alpha,iangle,input,...
                 trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated,maxiter,tol, runcount > 1, W, M, R, energyRatios(j));
        fitting_data = prep_for_fitting(obj, runcount == 1, j, anglecount, fitting_data);

        errs(j) = err;
    end

    plot_error(rad2deg(anglevalues), errs, resultdir, false);
    plot_error(mus, errs, resultdir, true);

    mean_0 = fitting_data(:,1)
    assert(length(mean_0) == anglecount);
    mean_1 = fitting_data(:,2)
    assert(length(mean_1) == anglecount);
    cv_0 = fitting_data(:,3)
    assert(length(cv_0) == anglecount);
    cv_1 = fitting_data(:,4)
    assert(length(cv_1) == anglecount);
    
    createFit(alpha, polyDegree, rad2deg(anglevalues), energyRatios, fitting_data, targetNames, targetFileNames, resultdir, false);
    createFit(alpha, polyDegree, mus, energyRatios, fitting_data, targetNames, targetFileNames, resultdir, true);
        
end %function end


function exportPolyFit(file, p)
    l = length(p);
    for i = 1:l
        fprintf(file, '%4.6f ', p(i));
    end
    fprintf(file, '\n');
end

function [anglevalues, energyRatios, mus] = selectValidAnglesForFitting(ox, oangle, x, angle)
    anglecandidates = unique(angle);
    anglecount = length(anglecandidates);

    anglevalues = zeros(anglecount,1);
    mus  = zeros(anglecount,1);
    energyRatios = zeros(anglecount,1);
    validCount = 0;
    for j = 1:anglecount
        iangle = anglecandidates(j);
        iox = ox(abs(oangle - iangle) < .0001);
        
        ix = x(abs(angle - iangle) < .0001);
        lenIX = length(ix);
        lenIOX = length(iox);
        energy = lenIX/lenIOX;
        if energy > .01
        %if energy > 0.00 
            validCount = validCount + 1;
            energyRatios(validCount) = energy;
            anglevalues(validCount) = iangle;
            mus(validCount) = cos(iangle);
        end 
    end
    energyRatios = energyRatios(1:validCount);
    anglevalues = anglevalues(1:validCount);
    mus = mus(1:validCount);
    assert(length(energyRatios) == length(anglevalues));
    assert(length(mus) == length(energyRatios));

end
    
       
function createFit(alpha, polyDegree, anglevalues, energyRatios, fitting_data, targetNames, targetFileNames, resultdir, isMu)
    filename = [resultdir, 'scatter_', num2str(alpha), '.txt'];
    file = fopen(filename, 'w');
    coefCount = polyDegree + 1;
    fprintf(file, "%d\n", coefCount);
       
    p = polyfit(anglevalues, energyRatios, polyDegree)

    exportPolyFit(file, p);    
    plotGaussianFit(polyDegree, anglevalues, energyRatios, targetNames{1}, targetFileNames{1}, resultdir, isMu); 
    anglecount = length(anglevalues);
    %assert(length(energyRatios) == anglecount);

    for i = 1:4
        f = fitting_data(:,i);
        %assert(length(f) == anglecount);
        p = polyfit(anglevalues, f, polyDegree)

        exportPolyFit(file, p);    
        plotGaussianFit(polyDegree, anglevalues, f, targetNames{i+1}, targetFileNames{i+1}, resultdir, isMu); 
    end
    fclose(file);
end
        
    


function [ox, oangle, x, y, z, angle] = read_data(datadir, alpha, bounce)
        close all
        cd(datadir)
        prefix = 'reflection_';
        %prefix = 'refraction_';
        prefix = strcat('refraction_0_', '3d_output');
        prefix = strcat('reflection_0_', '3d_output');
        filename = [prefix, 'x_',  num2str(alpha),'.p'];
        fileID = fopen(filename);
        C1 = fread(fileID, 'float');
        fclose(fileID);
        disp('x array size');
        disp(size(C1));
        
        filename = [prefix, 'y_', num2str(alpha),'.p'];
        fileID = fopen(filename);
        C2 = fread(fileID, 'float');
        fclose(fileID);
        
        filename = [prefix, 'z_',  num2str(alpha),'.p'];
        fileID = fopen(filename);
        C3 = fread(fileID, 'float');
        fclose(fileID);
        
        filename = [prefix, 'weight_',  num2str(alpha),'.p'];
        fileID = fopen(filename);
        C4 = fread(fileID, 'float');
        fclose(fileID);

        filename = [prefix, 'angle_', num2str(alpha),'.p'];
        fileID = fopen(filename);
        C5 = fread(fileID, 'float');
        fclose(fileID);
        
        filename = [prefix, 'depth_', num2str(alpha),'.p'];
        fileID = fopen(filename);
        %C6 = fread(fileID, 'float');
        C6 = fread(fileID, 'int');
        fclose(fileID);
        
        %{
        x = C1{1};
        y = C2{1};
        z = C3{1};
        weight = C4{1};
        angle = C5{1};
        depth = C6{1};
        observe = 6000;
        x = x/observe;
        y = y/observe;
        z = z/observe;
        %}
        
        x = C1;
        y = C2;
        z = C3;
        weight = C4;
        angle = C5;
        depth = C6;

        ox = x(z>=0);
        oangle = angle(z>=0);


        if bounce == 1
            x = x(depth >= 1);
            y = y(depth >= 1);
            z = z(depth >= 1);
            angle = angle(depth >= 1);
        else  
            x = x(depth >= bounce);
            y = y(depth >= bounce);
            z = z(depth >= bounce);
            angle = angle(depth >=bounce);
        end
        angle = angle(z>=0);
        x = x(z>=0);
        y = y(z>=0);
        z = z(z>=0);
        
        disp('whx');
        disp(x(1:10, 1));
        disp('why');
        disp(y(1:10, 1));
        disp('whz');
        disp(z(1:10, 1));
        
        disp('angle'); 
        disp(angle(1:10, 1));
        disp('depth'); 
        disp(depth(1:10, 1));
        
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
