function fitting_halfvector_z13D(dir,fundir,alpha,input,...
    trainnum, generatenum, gaussiannumvec,xnum, znum,accelerated,reflectdata,maxiter,tol)

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

[train, test,test2,range,x_unit,z_unit] = preprocess(input,...
    boundary_ratio,xnum,znum,zmax,zmin,trainnum,generatenum,reflectdata);

[invalidCount,trainresult] = bygrid2d(train,xnum,znum,range,zmax,zmin,x_unit,z_unit);
disp(invalidCount);

[invalidTestCount,result] = bygrid2d(test,xnum,znum,range,zmax,zmin,x_unit,z_unit);
result = result/sum(result(:));
disp(invalidTestCount);
plot_A(result, "GaussianHeightField", "HeightField");


[~,result2] = bygrid2d(test2,xnum,znum,range,zmax,zmin,x_unit,z_unit);
result2 = result2/sum(result2(:));

% for i = 1:ceil(znum/9):znum
%     
%     titlestring = ['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~',num2str(i)];
%     filename = [dir,'3d_mirror_reflect_alpha_',num2str(alpha), '~',num2str(i)];
%     plot_near_angle(result,i,titlestring,filename);
%     
% end

for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
    fprintf('using %4d Gaussians\n',numGaussian)
    
    if accelerated
        aemdir = [fundir,'accelerated_greedy_EM'];
        addpath(aemdir);
        obj = accelerated_em(train,trainnum,numGaussian,maxiter,tol);
    else
        obj = customized_em(train,numGaussian,maxiter,tol);
    end
    
    % save gm result
    filename = [dir,'3dhalf_projected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_reflect_',num2str(reflectdata),'.mat'];
    save(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,length(test));
    
    % calculate predict matrix
    [count,predict] = bygrid2d(Y,xnum,znum,range,zmax,zmin,x_unit,z_unit);
    predict = predict/sum(predict(:));
    plot_A(predict, "GMM", "GMM");
    
    % calculate relative L2 error
    err = relativel2err(result,predict);
    fprintf('error is %4.6f\n', err)
    errvec(j) = err;
    countvec(j) = count;
    
    % plot error by incident angle bin
    %filename = [dir,'3drelativel2error_alpha',num2str(alpha),'_reflect_',num2str(reflectdata)];
    %plot_error_by_angle(result,result2,predict,znum,reflectdata,filename)
    
    
%     for i = 1:ceil(znum/9):znum
%         
%         titlestring = ['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~',num2str(i)];
%         filename = [dir,'3d_mirror_predict_alpha_',num2str(alpha), '~',num2str(i)];
%         plot_near_angle(predict,i,titlestring,filename);
%         
%     end
    
end

% save error and count file
errvec_filename = [dir,'3dhalf_projected_z1',num2str(alpha),'_err.mat'];
countvec_filename = [dir,'half_projected_z1',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec');
save(countvec_filename,'count')

end

function [train, test,test2, range,x_unit,z_unit] = preprocess(input,...
    boundary_ratio,xnum, znum, zmax,zmin,trainnum,generatenum,reflectdata)
% preprocess raw data to eliminate extreme hx/hz hy/hz values for
% each incident angle
angle = input(:,4);
incident = [sin(angle), zeros(length(angle),1), cos(angle)];
h = input(:,1:3) + incident;
%h = (input(:,1:3) + incident)/2;
%hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
%h = h./hnorm;
xdividez = h(:,1)./h(:,3);
ydividez = h(:,2)./h(:,3);
xsq = times(xdividez, xdividez);
ysq = times(ydividez, ydividez);
tanThetaHSQ = xsq + ysq;
tanThetaH = sqrt(tanThetaHSQ);
disp(length(tanThetaH));

 
sortedtan = sort(tanThetaH);
tanrange = sortedtan(ceil(boundary_ratio*length(tanThetaH)));
disp(tanrange);
cutoff = tanrange;
range = cutoff;
disp(range);
x_unit = range/xnum;
disp(x_unit);
z_unit = (zmax-zmin)/znum;


tanThetaH_train = tanThetaH(1:trainnum);
angle_train = angle(1:trainnum);

trainindex = tanThetaH_train<=cutoff;
tanThetaH_train_new = tanThetaH_train(trainindex);
angle_train_new = angle_train(trainindex);

indexrange = trainnum+1:trainnum+generatenum;
tanThetaH_test = tanThetaH(indexrange);
angle_test = angle(indexrange);

%indexrange2 = trainnum+generatenum+1:trainnum+2*generatenum;
%tanThetaH_test2 = tanThetaH(indexrange2);
%angle_test2 = angle(indexrange2);

train = [tanThetaH_train_new,angle_train_new];
test = [tanThetaH_test,angle_test];
test2 = [tanThetaH_test,angle_test];
%test2 = [tanThetaH_test2,angle_test2];

end

function [count,result] = bygrid2d(A,xnum,znum,range,zmax,zmin,x_unit,z_unit)

count = 0;
result = zeros(xnum,znum);
for i = 1:length(A)
    if A(i,1)<range && A(i, 1) > 0 && A(i,2)<zmax && A(i,2)>zmin
        result(ceil(A(i,1)/x_unit),ceil((A(i,2)-zmin)/z_unit)) = result(ceil(A(i,1)/x_unit), ceil((A(i,2)-zmin)/z_unit)) + 1;
    else
        count = count + 1;
    end
end

end

function plot_A(A,titlestring,filename)

figure
imagesc(A)
ylabel('tanThetaH')
xlabel('incidentAngle')
colorbar()
title(titlestring)
%saveas(gcf,[filename,'.jpeg'])

end

function plot_error_by_angle(result,result2,predict,znum,reflectdata,filename)

if reflectdata
    e = zeros(znum/2,1);
    eself = zeros(znum/2,1);
    for i = znum/2+1:znum
        e(i-znum/2) = relativel2err(result(:,i),predict(:,i));
        eself(i-znum/2) = relativel2err(result(:,i),result2(:,i));
    end
else
    e = zeros(znum,1);
    eself = zeros(znum,1);
    for i = 1:znum
        e(i) = relativel2err(result(:,i),predict(:,i));
        eself(i) = relativel2err(result(:,i),result2(:,i));
    end
end

figure
plot(e,'linewidth',2)
hold on
plot(eself,'linewidth',2)
title('relative l2 error for each incident angle bin')
xlabel('angle bin')
ylabel('relative l2 error')
grid on
legend('fitting error','self error')
saveas(gcf,[filename,'.jpeg'])
end
