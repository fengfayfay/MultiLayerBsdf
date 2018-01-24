function Glass_fitting_halfvector_z1(dir, fundir, alpha,angle,input,...
    trainnum, generatenum, gaussiannumvec, xnum, ynum,ior,accelerated,maxiter,tol)

% fitting 2d galss heightfield data in slope domain using a mixture of gaussians
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
%           ior             - refractive index
%           accelerated     - if true use accelerated em, otherwise use
%                               cutomized gmcluster
%

% initialize error and count vector
errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% preprocess data
incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
boundary_ratio = 99/100;
[reflect_train, reflect_test,reflect_range,reflect_x_unit,reflect_y_unit,...
    transmit_train, transmit_test, transmit_range,transmit_x_unit,...
    transmit_y_unit] = preprocess(input,incident,boundary_ratio,xnum,ynum,...
                        trainnum,generatenum,ior);
                    
fprintf('reflect plotting range is %4.2f\n',reflect_range)
fprintf('transmit plotting range is %4.2f\n',transmit_range)

% plot the input data
titlestring = ['Gaussian Heightfiled reflect ray distribution, alpha=', num2str(alpha),' angle=',num2str(angle)];
filename = [dir,'reflect_halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
[~,reflect] = plotbygrid(xnum,ynum,reflect_test,...
    reflect_range,reflect_x_unit,reflect_y_unit,titlestring,filename);
reflect = reflect/sum(reflect(:));

titlestring = ['Gaussian Heightfiled transmit ray distribution, alpha=', num2str(alpha),' angle=',num2str(angle)];
filename = [dir,'transmit_halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
[~,transmit] = plotbygrid(xnum,ynum,transmit_test,...
    transmit_range,transmit_x_unit,transmit_y_unit,titlestring,filename);
transmit = transmit/sum(transmit(:));

for j = 1:length(gaussiannumvec)

    numGaussian = gaussiannumvec(j);

    if accelerated
        aemdir = [fundir,'accelerated_greedy_EM'];
        addpath(aemdir);
        reflect_obj = accelerated_em(reflect_train,trainnum,numGaussian,maxiter,tol);
        transmit_obj = accelerated_em(transmit_train,trainnum,numGaussian,maxiter,tol);
    else
        reflect_obj = customized_em(reflect_train,numGaussian,maxiter,tol);
        transmit_obj = customized_em(transmit_train,numGaussian,maxiter,tol);
    end
    
    % generate points from fitted model
    test_reflect_num = length(reflect_test);
    test_transmit_num = length(transmit_test);
    
    reflect_Y = random(reflect_obj,test_reflect_num);
    transmit_Y = random(transmit_obj,test_transmit_num);
    
    % plot gmm generated data
    titlestring = ['reflect distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)];
    filename = [dir,'reflect_half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    [reflect_count,reflect_predict] = plotbygrid(xnum,ynum,...
        reflect_Y,reflect_range,reflect_x_unit,reflect_y_unit,titlestring,filename);
    
    titlestring = ['transmit distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)];
    filename = [dir,'transmit_half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    [transmit_count,transmit_predict] = plotbygrid(xnum,ynum,...
        transmit_Y,transmit_range,transmit_x_unit,transmit_y_unit,titlestring,filename);
    
    % calculate error
    reflect_predict = reflect_predict/sum(reflect_predict(:));
    transmit_predict = transmit_predict/sum(transmit_predict(:));
    reflect_err = relativel2err(reflect,reflect_predict);
    transmit_err = relativel2err(transmit,transmit_predict);
    fprintf('reflect err: %4.6f\n', reflect_err)
    fprintf('transmit err: %4.6f\n', transmit_err)
    reflect_errvec(j) = reflect_err;
    transmit_errvec(j) = transmit_err;
    countvec(j) = reflect_count + transmit_count;
    
    
end

reflect_errvec_filename = [dir,'reflect_half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
transmit_errvec_filename = [dir,'transmit_half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = [dir,'glass_half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(reflect_errvec_filename,'reflect_errvec')
save(transmit_errvec_filename,'transmit_errvec')
save(countvec_filename,'countvec')

end

function [reflect_train,reflect_test,reflect_range,reflect_x_unit,reflect_y_unit,...
    transmit_train, transmit_test, transmit_range,transmit_x_unit,...
    transmit_y_unit] = preprocess(input,incident,boundary_ratio,xnum,ynum,...
    trainnum,generatenum,ior)

% reflect
reflect_index = find(input(:,3)>=0);

reflect_h = (input(reflect_index,:)+incident)/2;
reflect_hnorm = repmat(sqrt(sum(reflect_h.^2,2)),1,3);
reflect_h = reflect_h./reflect_hnorm;

reflect_xdividez = reflect_h(:,1)./reflect_h(:,3);
reflect_ydividez = reflect_h(:,2)./reflect_h(:,3);
reflect_sortedxz = sort(abs(reflect_xdividez));
reflect_xrange = reflect_sortedxz(ceil(boundary_ratio*length(reflect_xdividez)));
reflect_sortedyz = sort(abs(reflect_ydividez));
reflect_yrange = reflect_sortedyz(ceil(boundary_ratio*length(reflect_xdividez)));
cutoff = max(reflect_xrange, reflect_yrange);
reflect_range = 2*cutoff;
reflect_x_unit = 2*reflect_range/xnum;
reflect_y_unit = 2*reflect_range/ynum;

reflect_xdividez_train = reflect_xdividez(reflect_index<=trainnum);
reflect_ydividez_train = reflect_ydividez(reflect_index<=trainnum);
reflect_xdividez_new = reflect_xdividez_train(abs(reflect_xdividez_train)<=cutoff&abs(reflect_ydividez_train)<=cutoff);
reflect_ydividez_new = reflect_ydividez_train(abs(reflect_xdividez_train)<=cutoff&abs(reflect_ydividez_train)<=cutoff);
reflect_xdividez_test = reflect_xdividez(reflect_index>trainnum & reflect_index<=generatenum+trainnum);
reflect_ydividez_test = reflect_ydividez(reflect_index>trainnum & reflect_index<=generatenum+trainnum);

reflect_train = [reflect_xdividez_new, reflect_ydividez_new];
reflect_test = [reflect_xdividez_test,reflect_ydividez_test];

% transmit
transmit_index = find(input(:,3)<0);
transmit_h = ior*input(transmit_index,:)+ incident;
transmit_hnorm = repmat(sqrt(sum(transmit_h.^2,2)),1,3);
transmit_h = transmit_h./transmit_hnorm;

transmit_xdividez = transmit_h(:,1)./transmit_h(:,3);
transmit_ydividez = transmit_h(:,2)./transmit_h(:,3);
transmit_sortedxz = sort(abs(transmit_xdividez));
transmit_xrange = transmit_sortedxz(ceil(boundary_ratio*length(transmit_xdividez)));
transmit_sortedyz = sort(abs(transmit_ydividez));
transmit_yrange = transmit_sortedyz(ceil(boundary_ratio*length(transmit_xdividez)));
cutoff = max(transmit_xrange, transmit_yrange);
transmit_range = 2*cutoff;
transmit_x_unit = 2*transmit_range/xnum;
transmit_y_unit = 2*transmit_range/ynum;

transmit_xdividez_train = transmit_xdividez(transmit_index<=trainnum);
transmit_ydividez_train = transmit_ydividez(transmit_index<=trainnum);
transmit_xdividez_new = transmit_xdividez_train(abs(transmit_xdividez_train)<=cutoff&abs(transmit_ydividez_train)<=cutoff);
transmit_ydividez_new = transmit_ydividez_train(abs(transmit_xdividez_train)<=cutoff&abs(transmit_ydividez_train)<=cutoff);
transmit_xdividez_test = transmit_xdividez(transmit_index>trainnum & transmit_index<=generatenum+trainnum);
transmit_ydividez_test = transmit_ydividez(transmit_index>trainnum & transmit_index<=generatenum+trainnum);

transmit_train = [transmit_xdividez_new, transmit_ydividez_new];
transmit_test = [transmit_xdividez_test, transmit_ydividez_test];


end