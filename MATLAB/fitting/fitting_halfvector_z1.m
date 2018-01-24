function fitting_halfvector_z1(dir,fundir,alpha,angle,input,...
    trainnum, generatenum, gaussiannumvec, xnum, ynum,accelerated)

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
incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
boundary_ratio = 99/100;
[train, test,range,x_unit,y_unit] = preprocess(input,incident,...
    boundary_ratio,xnum,ynum,trainnum,generatenum);

% plot the input data
titlestring = ['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),' angle=',num2str(angle)];
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];

cd(fundir)
[~,result] = plotbygrid(xnum,ynum,test,range,x_unit,y_unit,titlestring,filename);
result = result/sum(result(:));

for j = 1:length(gaussiannumvec)
    
    numGaussian = gaussiannumvec(j);
    
    if accelerated
        aemdir = [fundir,'accelerated_greedy_EM'];
        addpath(aemdir);
        obj = accelerated_em(train,trainnum,numGaussian);
    else
        maxiter = 500;
        tol = 1e-5;
        obj = customized_em(maxiter,tol,train,numGaussian);
    end
    
    % save gm result
    filename = [dir,'half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
    % plot gmm generated data
    titlestring = ['Mirror distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)];
    filename = [dir,'half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    [count,predict] = plotbygrid(xnum,ynum,Y,range,x_unit,y_unit,titlestring,filename);
    predict = predict/sum(predict(:));
    
    % calculate relative L2 error
    [err,~,~,~] = relativel2err(result,predict);
    fprintf('err: %4.6f\n', err)
    
    errvec(j) = err;
    countvec(j) = count;
    
    
end

% save error and count file
errvec_filename = [dir,'half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = [dir,'half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')

end

function [train, test,range,x_unit,y_unit] = preprocess(input,incident,...
    boundary_ratio,xnum,ynum,trainnum,generatenum)

h = (input + incident)/2;
hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
h = h./hnorm;
xdividez = h(:,1)./h(:,3);
ydividez = h(:,2)./h(:,3);
sortedxz = sort(abs(xdividez));
xrange = sortedxz(ceil(boundary_ratio*length(xdividez)));
sortedyz = sort(abs(ydividez));
yrange = sortedyz(ceil(boundary_ratio*length(xdividez)));
range = 2*max(xrange, yrange);
fprintf('plotting range is %4.2f\n',range)
x_unit = 2*range/xnum;
y_unit = 2*range/ynum;

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
