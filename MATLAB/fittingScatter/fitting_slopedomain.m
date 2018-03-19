function [obj, W, M, R] = fitting_slopedomain(dir,fundir,alpha,angle,input,...
    trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated,maxiter,tol,softinit, W, M, R, energyRatio)

% fitting 2d mirror heightfield data in slope domain using a mixture of gaussians
%   
%
% Input:    dir             - inputdata directory
%           fundir          - funcion location
%           alpha           - roughness
%           angle           - incident angle
%           input           - heightfield data
%           trainnum        - number of training data
%           generatenum     - number of testing data
%           gaussiannumvec  - number of gaussians vector
%           xnum            - plotting resolution in x direction
%           ynum            - plotting resolution in y direction
%           accelerated     - if true use accelerated em, otherwise use
%                               cutomized gmcluster
%

% initialize error and count vector
errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% preprocess data
incident = [sin(angle), 0, cos(angle)];
boundary_ratio = 99/100;
[train, test,range] = preprocess(input,incident,...
    boundary_ratio,xnum,ynum,trainnum,generatenum);

% plot the input data
titlestring = ['Gaussian heightfield slope distribution, alpha=', num2str(alpha),' angle=',num2str(rad2deg(angle))];
filename = [dir,'slopedomain_',num2str(rad2deg(angle)),'_alpha_',num2str(alpha), '_heightfield'];

[~,result] = plotbygrid(xnum,ynum,test,range,titlestring,filename, energyRatio);
result = result/sum(result(:));
filename = [dir,'angle_',num2str(rad2deg(angle)),'_alpha_',num2str(alpha), '_brdfsim'];
[~,brdfsimulated] = plotgrid(input, xnum, ynum, titlestring, filename, energyRatio);

hold();

cd(fundir)

for j = 1:length(gaussiannumvec)
    
    numGaussian = gaussiannumvec(j);
    
    if accelerated
        aemdir = [fundir,'accelerated_greedy_EM'];
        addpath(aemdir);
        [obj, W, M, R] = accelerated_em(train,trainnum,numGaussian,maxiter,tol, softinit, W, M, R);
    else
        [obj, W, M, R] = customized_em(train,numGaussian,maxiter,tol, softinit, W, M, R);
    end

    %plotGMM(obj, 0);
    
    % save gm result
    filename = [dir,'slopedomain_',num2str(rad2deg(angle)),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    %save(filename,'obj')
    filename = [dir,'angle_',num2str(rad2deg(angle)),'_alpha_',num2str(alpha), '_brdfgmm'];
    gm2brdf(obj, 2, rad2deg(angle),alpha, 0, filename, energyRatio);
     
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    % plot gmm generated data
    
    titlestring = ['Slope distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(rad2deg(angle)),' #G=',num2str(numGaussian)];
    filename = [dir,'slopedomain_',num2str(rad2deg(angle)),'_alpha_',num2str(alpha), '_gmm',num2str(numGaussian)];
    [count,predict] = plotbygrid(xnum,ynum,Y,range, titlestring,filename, energyRatio);
    predict = predict/sum(predict(:));
    
    % calculate relative L2 error
    [err,~,~,~] = relativel2err(result,predict);
    fprintf('err: %4.6f\n', err)
    
    errvec(j) = err;
    countvec(j) = count;
    
end

% save error and count file
errvec_filename = [dir,'slopedomain_',num2str(rad2deg(angle)),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = [dir,'slopedomain_',num2str(rad2deg(angle)),'_angle_',num2str(alpha),'_badcount.mat'];
%save(errvec_filename,'errvec')
%save(countvec_filename,'count')

end

function [train, test,range] = preprocess(input,incident,...
    boundary_ratio,xnum,ynum,trainnum,generatenum)

    h = (input + incident)/2;
    hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
    h = h./hnorm;
    xdividez = h(:,1)./h(:,3);
    ydividez = h(:,2)./h(:,3);
    sortedxz = sort(xdividez);
    xrange = sortedxz(ceil(boundary_ratio*length(xdividez)));
    disp([sortedxz(1), xrange]);
    xrange = xrange - sortedxz(1);
    sortedyz = sort(ydividez);
    yrange = sortedyz(ceil(boundary_ratio*length(xdividez)));
    disp([sortedyz(1), yrange]);
    yrange = yrange - sortedyz(1);
    range = max(xrange, yrange) * 1.25;
    fprintf('plotting range is %4.2f\n',range)

    % training data
    xtrain = xdividez(1:trainnum);
    ytrain = ydividez(1:trainnum);
    xtrain_new = xtrain(abs(xtrain)<=range/2 & abs(ytrain)<=range/2);
    ytrain_new = ytrain(abs(xtrain)<=range/2 & abs(ytrain)<=range/2);
    train = [xtrain_new, ytrain_new];

    % testing data
    indexrange = trainnum+1:trainnum+generatenum;
    test = [xdividez(indexrange),ydividez(indexrange)];
end
