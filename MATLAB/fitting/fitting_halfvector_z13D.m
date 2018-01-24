function fitting_halfvector_z13D(dir,fundir,alpha,input,...
    trainnum, generatenum, gaussiannumvec,xnum, ynum, znum,accelerated)

% fitting 3d mirror heightfield data in slope domain using a mixture of gaussians
%
%
% Input:    dir             - inputdata directory
%           fundir          - funcion location
%           alpha           - roughness
%           input           - heightfield data
%           trainnum        - number of training data
%           generatenum     - number of testing data
%           gaussiannumvec  - number of gaussians vector
%           xnum            - plotting resolution in x direction
%           ynum            - plotting resolution in y direction
%           znum            - plotting resolution in z direction
%           accelerated     - if true use accelerated em, otherwise use
%                               cutomized gmcluster
%

% initialize error and count vector
errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

boundary_ratio = 99/100;
zmax = pi/2;
zmin = 0;
[train, test,range,x_unit,y_unit,z_unit] = preprocess(input,...
    boundary_ratio,xnum,ynum,znum,zmax,zmin,trainnum,generatenum);

[~,result] = bygrid3d(test,xnum,ynum,znum,range,zmax,zmin,x_unit,y_unit,z_unit);
result = result/sum(result(:));

% 3d point distribution visualization
figure
scatter3(test(:,1), test(:,2), test(:,3),'filled')
xlabel('hx/hz')
ylabel('hy/hz')
zlabel('incident angle')

for i = 1:ceil(znum/9):znum
    
    titlestring = ['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~',num2str(i)];
    filename = [dir,'3d_mirror_reflect_alpha_',num2str(alpha), '~',num2str(i)];
    plot_near_angle(result,i,titlestring,filename);
    
end

for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
    fprintf('using %4d Gaussians\n',numGaussian)
    
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
    filename = [dir,'3dhalf_projected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,length(test));
    
    % calculate predict matrix
    [count,predict] = bygrid3d(Y,xnum,ynum,znum,range,zmax,zmin,x_unit,y_unit,z_unit);
    predict = predict/sum(predict(:));
    
    % calculate relative L2 error
    err = relativel2err(result,predict);
    fprintf('error is %4.6f\n', err)
    errvec(j) = err;
    countvec(j) = count;
    
    for i = 1:znum
        e(i) = relativel2err(result(:,:,i),predict(:,:,i));
    end
    figure
    plot(1:znum*90/znum,e,'linewidth',2)
    title('relative l2 error for each incident angle bin')
    xlabel('angle bin')
    ylabel('relative l2 error')
    grid on
    filename = [dir,'3drelativel2error_alpha',num2str(alpha)];
    saveas(gcf,[filename,'.jpeg'])
    
    
    for i = 1:ceil(znum/9):znum
        
        titlestring = ['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~',num2str(i)];
        filename = [dir,'3d_mirror_predict_alpha_',num2str(alpha), '~',num2str(i)];
        plot_near_angle(predict,i,titlestring,filename);
        
    end
    
end

% save error and count file
errvec_filename = [dir,'3dhalf_projected_z1',num2str(alpha),'_err.mat'];
countvec_filename = [dir,'half_projected_z1',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')

end

function [train, test,range,x_unit,y_unit,z_unit] = preprocess(input,...
    boundary_ratio,xnum,ynum,znum,zmax,zmin,trainnum,generatenum)
% preprocess raw data to eliminate extreme hx/hz hy/hz values for
% each incident angle
angle = input(:,4);
incident = [sin(angle), zeros(length(angle),1), cos(angle)];
h = (input(:,1:3) + incident)/2;
hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
h = h./hnorm;
xdividez = h(:,1)./h(:,3);
ydividez = h(:,2)./h(:,3);
sortedxz = sort(abs(xdividez));
xrange = sortedxz(ceil(boundary_ratio*length(xdividez)));
sortedyz = sort(abs(ydividez));
yrange = sortedyz(ceil(boundary_ratio*length(xdividez)));
cutoff = max(xrange, yrange);
range = 2*cutoff;
x_unit = 2 * range/xnum;
y_unit = 2 * range/ynum;
z_unit = (zmax-zmin)/znum;

xdividez_train = xdividez(1:trainnum);
ydividez_train = ydividez(1:trainnum);

angle_train = angle(1:trainnum);

trainindex = abs(xdividez_train)<=cutoff&abs(ydividez_train)<=cutoff;
xdividez_train_new = xdividez_train(trainindex);
ydividez_train_new = ydividez_train(trainindex);
angle_train_new = angle_train(trainindex);

indexrange = trainnum+1:trainnum+generatenum;
xdividez_test = xdividez(indexrange);
ydividez_test = ydividez(indexrange);
angle_test = angle(indexrange);

% % relfection around normal incidence
% xdividez_train_new = [xdividez_train_new;xdividez_train_new];
% ydividez_train_new = [ydividez_train_new;ydividez_train_new];
% angle_train_new = [angle_train_new;-angle_train_new];
% xdividez_test = [xdividez_test;xdividez_test];
% ydividez_test = [ydividez_test;ydividez_test];
% angle_test = [angle_test;-angle_test];

train = [xdividez_train_new,ydividez_train_new,angle_train_new];
test = [xdividez_test,ydividez_test,angle_test];

end

function [count,result] = bygrid3d(A,xnum,ynum,znum,range,zmax,zmin,x_unit,y_unit,z_unit)

count = 0;
result = zeros(xnum,ynum,znum);
for i = 1:length(A)
    if abs(A(i,1))<range && abs(A(i,2))<range && A(i,3)<zmax && A(i,3)>zmin
        result(ceil((A(i,1)+range)/x_unit),ceil((A(i,2)+range)/y_unit),ceil((A(i,3)-zmin)/z_unit)) = ...
            result(ceil((A(i,1)+range)/x_unit),ceil((A(i,2)+range)/y_unit),ceil((A(i,3)-zmin)/z_unit)) + 1;
    else
        count = count + 1;
    end
end

end

function plot_near_angle(A,angle,titlestring,filename)

figure
imagesc(A(:,:,angle))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(titlestring)
saveas(gcf,[filename,'.jpeg'])

end