function Glass_fitting_halfvector_z1(dir,alpha,angle,x,y,z,weight,...
    testafter, trainnum, generatenum, gaussiannumvec, incident, xnum, ynum,ior)

cd(dir)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% plot the input data
% reflect
reflect_index = find(z>=0);
reflect_h = ([x(z>=0), y(z>=0), z(z>=0)] + incident)/2;
reflect_hnorm = repmat(sqrt(sum(reflect_h.^2,2)),1,3);
reflect_h = reflect_h./reflect_hnorm;
reflect_xdividez = reflect_h(:,1)./reflect_h(:,3);
reflect_ydividez = reflect_h(:,2)./reflect_h(:,3);
reflect_xrange = max(max(reflect_xdividez),-min(reflect_xdividez));
reflect_yrange = max(max(reflect_ydividez),-min(reflect_ydividez));
reflect_range = 2*max(reflect_xrange, reflect_yrange);
reflect_x_unit = 2*reflect_range/xnum;
reflect_y_unit = 2*reflect_range/ynum;
reflect = zeros(xnum,ynum);
reflect_test_index = reflect_index(reflect_index>trainnum);
reflect_xdividez_train = reflect_xdividez(reflect_index<=trainnum);
reflect_ydividez_train = reflect_ydividez(reflect_index<=trainnum);
reflect_xdividez_test = reflect_xdividez(reflect_index>trainnum);
reflect_ydividez_test = reflect_ydividez(reflect_index>trainnum);

for i = 1:length(reflect_xdividez_test)
    reflect(ceil((reflect_xdividez_test(i) + reflect_range)/reflect_x_unit), ceil((reflect_ydividez_test(i) + reflect_range)/reflect_y_unit))...
            = reflect(ceil((reflect_xdividez_test(i) + reflect_range)/reflect_x_unit), ceil((reflect_ydividez_test(i) + reflect_range)/reflect_y_unit)) + ...
            weight(reflect_test_index(i));
end

% close all
% figure
% imagesc(reflect/generatenum)
% ylabel('x/z')
% xlabel('y/z')
% colorbar()
% title(['Gaussian Heightfiled reflect ray distribution, alpha=', num2str(alpha),' angle=',num2str(angle)])
% filename = ['reflect_halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
% saveas(gcf,[filename,'.jpeg'])

% transmit
transmit_index = find(z<0);
transmit_h = ior*[x(z<0), y(z<0), z(z<0)] + incident;
transmit_hnorm = repmat(sqrt(sum(transmit_h.^2,2)),1,3);
transmit_h = transmit_h./transmit_hnorm;
transmit_xdividez = transmit_h(:,1)./transmit_h(:,3);
transmit_ydividez = transmit_h(:,2)./transmit_h(:,3);

% set the plotting range
transmit_range = alpha*4;
transmit_x_unit = 2*transmit_range/xnum;
transmit_y_unit = 2*transmit_range/ynum;
transmit = zeros(xnum,ynum);
transmit_test_index = transmit_index(transmit_index>trainnum);
transmit_xdividez_train = transmit_xdividez(transmit_index<=trainnum);
transmit_ydividez_train = transmit_ydividez(transmit_index<=trainnum);

% remove extreme data from training
transmit_xdividez_new = transmit_xdividez_train(abs(transmit_xdividez_train)<=10&abs(transmit_xdividez_train)<=10);
transmit_ydividez_new = transmit_ydividez_train(abs(transmit_xdividez_train)<=10&abs(transmit_xdividez_train)<=10);

transmit_xdividez_test = transmit_xdividez(transmit_index>trainnum);
transmit_ydividez_test = transmit_ydividez(transmit_index>trainnum);

for i = 1:length(transmit_xdividez_test)
    if abs(transmit_xdividez_test(i))<transmit_range && abs(transmit_ydividez_test(i))<transmit_range
        transmit(ceil((transmit_xdividez_test(i) + transmit_range)/transmit_x_unit), ceil((transmit_ydividez_test(i) + transmit_range)/transmit_y_unit))...
            = transmit(ceil((transmit_xdividez_test(i) + transmit_range)/transmit_x_unit), ceil((transmit_ydividez_test(i) + transmit_range)/transmit_y_unit)) + ...
            weight(transmit_test_index(i));
    end
end

close all
figure
imagesc(transmit/generatenum)
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled transmit ray distribution, alpha=', num2str(alpha),' angle=',num2str(angle)])
filename = ['transmit_halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
saveas(gcf,[filename,'.jpeg'])

%training data
reflect_train = [reflect_xdividez_train, reflect_ydividez_train];
transmit_train = [transmit_xdividez_new, transmit_ydividez_new];
% 
% close all
% figure
% plot(result(xnum/2,:)/generatenum, 'linewidth', 2)
% hold on
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);

    options = statset('MaxIter',500, 'Display','final','TolFun',1e-5);
    try
        reflect_obj = fitgmdist(reflect_train,numGaussian,'Options',options);
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message
    end
    
    try
        transmit_obj = fitgmdist(transmit_train,numGaussian,'Options',options,'Start','customize');
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message
    end
    
%     disp(obj.mu)
%     
%     filename = ['reflect_half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
%     save(filename,'reflect_obj')
%     
%     filename = ['transmit_half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
%     save(filename,'transmit_obj')
    
    
%     %% visualize result
%     figure
%     scatter(train(:,1),train(:,2),10,'.')
%     hold on
%     h = ezcontour(@(x,y)pdf(obj,[x y]),[0 2*pi],[0 1]);
    
    %% generate points from fitted model
    test_reflect_num = length(reflect_xdividez_test);
    test_transmit_num = length(transmit_xdividez_test);
     
    if test_reflect_num + test_transmit_num ~= generatenum
        msg = 'number of test points is incorrect';
        error(msg)
    end
    
    reflect_Y = random(reflect_obj,test_reflect_num);
    transmit_Y = random(transmit_obj,test_transmit_num);
    
    %% calculate energy
    reflect_predict= zeros(xnum,ynum);
    reflect_count = 0;
    for i = 1:test_reflect_num
        if abs(reflect_Y(i,1))<=reflect_range && abs(reflect_Y(i,2))<=reflect_range
              reflect_predict(ceil((reflect_Y(i,1)+reflect_range)/reflect_x_unit),ceil((reflect_Y(i,2)+reflect_range)/reflect_y_unit)) = ...
                reflect_predict(ceil((reflect_Y(i,1)+reflect_range)/reflect_x_unit),ceil((reflect_Y(i,2)+reflect_range)/reflect_y_unit)) + 1;
        else
            reflect_count = reflect_count +1;
        end
    end
    
    transmit_predict= zeros(xnum,ynum);
    transmit_count = 0;
    for i = 1:test_transmit_num
        if abs(transmit_Y(i,1))<=transmit_range && abs(transmit_Y(i,2))<=transmit_range
              transmit_predict(ceil((transmit_Y(i,1)+transmit_range)/transmit_x_unit),ceil((transmit_Y(i,2)+transmit_range)/transmit_y_unit)) = ...
                transmit_predict(ceil((transmit_Y(i,1)+transmit_range)/transmit_x_unit),ceil((transmit_Y(i,2)+transmit_range)/transmit_y_unit)) + 1;
        else
            transmit_count = transmit_count +1;
        end
    end
    
%     plot(predict(xnum/2,:)/(generatenum-count), 'linewidth', 2)
    count = transmit_count + reflect_count;
%     figure
%     imagesc(reflect_predict/(generatenum-count))
%     colorbar()
%     title(['Reflect distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
%     filename = ['reflect_half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
%     saveas(gcf,[filename,'.jpeg'])
    
    figure
    imagesc(transmit_predict/(generatenum-count))
    colorbar()
    title(['Transmit distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename = ['transmit_half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    saveas(gcf,[filename,'.jpeg'])
    
%     figure
%     plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
%     title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
%     filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
%     saveas(gcf,[filename,'.jpeg'])
    
    %% calculate error
    reflect_diff = reflect_predict/(generatenum-count)-reflect/generatenum;
    transmit_diff = transmit_predict/(generatenum-count)-transmit/generatenum;
%     diff_new = sort(abs(reshape(diff,10000,1)),'descend');
%     figure
%     plot(diff_new, 'linewidth',2)
%     title(['error from large to small, 5 gaussian, alpha=',num2str(alpha)])
    reflect_err = sqrt(sum(sum(reflect_diff.*reflect_diff)));
    transmit_err = sqrt(sum(sum(transmit_diff.*transmit_diff)));
    
    reflect_errvec(j) = reflect_err;
    transmit_errvec(j) = transmit_err;
    countvec(j) = reflect_count + transmit_count;
    
    
end
% title(['Energy generated using GMM on projected h, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])

reflect_errvec_filename = ['reflect_half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
transmit_errvec_filename = ['transmit_half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = ['glass_half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(reflect_errvec_filename,'reflect_errvec')
save(transmit_errvec_filename,'transmit_errvec')
save(countvec_filename,'count')
end