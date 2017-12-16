function fitting_halfvector_z1(dir,alpha,angle,x,y,z,weight,...
    testafter, trainnum, generatenum, gaussiannumvec, incident, xnum, ynum)

cd(dir)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% plot the input data
h = ([x, y, z] + incident)/2;
hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
h = h./hnorm;
xdividez = h(:,1)./h(:,3);
ydividez = h(:,2)./h(:,3);
xrange = max(max(xdividez),-min(xdividez));
yrange = max(max(ydividez),-min(ydividez));
range = 2*max(xrange, yrange);
x_unit = 2*range/xnum;
y_unit = 2*range/ynum;
result = zeros(xnum,ynum);

for i = testafter+1:length(x)
    result(ceil((xdividez(i)+range)/x_unit),ceil((ydividez(i)+range)/y_unit)) = ...
        result(ceil((xdividez(i)+range)/x_unit),ceil((ydividez(i)+range)/y_unit)) + weight(i);
end

close all
figure
imagesc(result/generatenum)
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Energy of Gaussian Heightfiled, alpha=', num2str(alpha),' angle=',num2str(angle)])
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
saveas(gcf,[filename,'.jpeg'])

% training data
train = [xdividez(1:trainnum), ydividez(1:trainnum)];

% close all
% figure
% plot(result(xnum/2,:)/generatenum, 'linewidth', 2)
% hold on
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);

    options = statset('MaxIter',500, 'Display','final');
    try
        obj = fitgmdist(train,numGaussian,'Options',options);
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message
    end
    
    filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
%     filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
%     load(filename,'obj')
    
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
        if Y(i,1)>=-range && Y(i,1)<=range && Y(i,2)>=-range && Y(i,2)<=range
              predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) = ...
                predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) + 1;
        else
            count = count +1;
        end
    end
    
%     plot(predict(xnum/2,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict/(generatenum-count))
    colorbar()
    title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    saveas(gcf,[filename,'.jpeg'])
    
%     figure
%     plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
%     title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
%     filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
%     saveas(gcf,[filename,'.jpeg'])
    
    %% calculate error
%     diff = predict/(generatenum-count)-result/generatenum;
%     err = sqrt(sum(sum(diff.*diff)));
%     
%     errvec(j) = err;
%     countvec(j) = count;
    
    
end
% title(['Energy generated using GMM on projected h, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])

errvec_filename = ['half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = ['half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')
end