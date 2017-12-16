function fitting_halfvector_disk(dir,alpha,angle,x,y,z,weight,...
    testafter, trainnum, generatenum, gaussiannumvec, incident, xynum, znum)

cd(dir)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% Energy value from data
xy_unit = 2*sqrt(2)/xynum;
z_unit = 2/znum;
result = zeros(xnum,ynum);
for i = testafter+1:length(x)
    if z(i)>=0
        h = ([x(i),y(i),z(i)] + incident)/2;
        h = h/norm(h);
        
        dim1 = sqrt(h(1)^2 + h(2)^2);
        dim2 = h(3);
        result(ceil((-dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) = result(ceil((-dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) + weight(i);
        result(ceil((dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) = result(ceil((dim1+sqrt(2))/xy_unit),ceil((dim2+1)/z_unit)) + weight(i);
    end
end

close all
figure
imagesc(result/generatenum)
colorbar()
title(['Energy of Gaussian Heightfiled, alpha=', num2str(alpha),' angle=',num2str(angle)])
filename = ['halfprojected',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
saveas(gcf,[filename,'.jpeg'])

xvec = zeros(trainnum,1);
yvec = zeros(trainnum,1);
for i = 1:trainnum
    if z(i)>=0
        h = ([x(i),y(i),z(i)] + incident)/2;
        h = h/norm(h);
        xvec(i) = h(1);
        yvec(i) = h(2);
    end
end
train = [xvec, yvec];


% close all
% figure
% plot(result(xnum/2,:)/generatenum, 'linewidth', 2)
% hold on
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
    
    options = statset('MaxIter',10000, 'Display','final');
    obj = gmdistribution.fit(train,numGaussian,'Options',options);
    
    filename = ['half_projected',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
%     %% visualize result
%     figure
%     scatter(train(:,1),train(:,2),10,'.')
%     hold on
%     h = ezcontour(@(x,y)pdf(obj,[x y]),[0 2*pi],[0 1]);
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
    %% calculate energy
    predict= zeros(xnum,ynum);
    count = 0;
    for i = 1:generatenum
        if Y(i,1)>=-1 && Y(i,1)<=1 && Y(i,2)>=-1 && Y(i,2)<=1
              predict(ceil((Y(i,1)+1)/x_unit),ceil((Y(i,2)+1)/y_unit)) = ...
                predict(ceil((Y(i,1)+1)/x_unit),ceil((Y(i,2)+1)/y_unit)) + 1;
        else
            count = count +1;
        end
    end
    
%     plot(predict(xnum/2,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict/(generatenum-count))
    colorbar()
    title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename = ['half_projected',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    saveas(gcf,[filename,'.jpeg'])
    
%     figure
%     plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
%     title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
%     filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
%     saveas(gcf,[filename,'.jpeg'])
    
    %% calculate error
    diff = predict/(generatenum-count)-result/generatenum;
    err = sqrt(sum(sum(diff.*diff)));
    
    errvec(j) = err;
    countvec(j) = count;
    
    
end
% title(['Energy generated using GMM on projected h, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=['half_projected',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])

errvec_filename = ['half',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = ['half',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')
end