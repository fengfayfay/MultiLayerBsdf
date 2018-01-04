function fitting_halfvector_z1(dir,alpha,angle,x,y,z,...
    trainnum, generatenum, gaussiannumvec, incident, xnum, ynum)

cd(dir)
errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% plot the input data
h = ([x, y, z] + incident)/2;
hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
h = h./hnorm;
xdividez = h(:,1)./h(:,3);
ydividez = h(:,2)./h(:,3);
sortedxz = sort(abs(xdividez));
% xrange = sortedxz(ceil(9999/10000*length(xdividez)));
xrange = max(xdividez);
sortedyz = sort(abs(ydividez));
% yrange = sortedyz(ceil(9999/10000*length(xdividez)));
yrange = max(ydividez);
cutoff = max(xrange, yrange);
range = 2*max(xrange, yrange);
fprintf('range is %4.2f\n',range)
x_unit = 2*range/xnum;
y_unit = 2*range/ynum;
result = zeros(xnum,ynum);

for i = trainnum+1:trainnum+generatenum
    if abs(xdividez(i))<range && abs(ydividez(i))<range
        result(ceil((xdividez(i)+range)/x_unit),ceil((ydividez(i)+range)/y_unit)) = ...
            result(ceil((xdividez(i)+range)/x_unit),ceil((ydividez(i)+range)/y_unit)) + 1;
    end
end

close all
figure
imagesc(result/sum(sum(result)))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),' angle=',num2str(angle)])
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
% saveas(gcf,[filename,'.jpeg'])

%training data
xtrain = xdividez(1:trainnum);
ytrain = ydividez(1:trainnum);
xtrain_new = xtrain(abs(xtrain)<=cutoff & abs(ytrain)<=cutoff);
ytrain_new = ytrain(abs(xtrain)<=cutoff & abs(ytrain)<=cutoff);
train = [xtrain_new, ytrain_new];

% % % use accelearted em for fitting
% t=cputime;
% tree = buildtree(train, 0, 0, 3, 1000);
% [W,M,R,ff,Ws,Ms,Rs] = em(train,[],5,0,1,tree);
% fprintf('\nRuntime: %.2f seconds\n', cputime-t);
% Rnew = reshape(R', 2,2,5);
% obj = gmdistribution(M,Rnew,W');
% Y = random(obj,generatenum);

% close all
% figure
% plot(result(xnum/2,:)/generatenum, 'linewidth', 2)
% hold on
for j = 1:length(gaussiannumvec)
    t=cputime;
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
    
    options = statset('MaxIter',500, 'Display','final','TolFun',1e-5);
    try
        obj = fitgmdist(train,numGaussian,'Options',options,'Start','customize');
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message
    end
    fprintf('\nRuntime: %.2f seconds\n', cputime-t);
    
    %     disp(obj.mu)
    %
    %     filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    %     save(filename,'obj')
    
    %     filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    %     load(filename,'obj')
    
    %     %% visualize result
    %     figure
    %     scatter(train(:,1),train(:,2),10,'.')
    %     hold on
    %     h = ezcontour(@(x,y)pdf(obj,[x y]),[-range range],[-range range]);
    %     fprintf('\nRuntime: %.2f seconds\n', cputime-t);
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
%    calculate predict matrix
    predict= zeros(xnum,ynum);
    count = 0;
    for i = 1:generatenum
        if abs(Y(i,1))<range && abs(Y(i,2))<range
            predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) = ...
                predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit)) + 1;
        else
            count = count +1;
        end
    end
    
    %     plot(predict(xnum/2,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict/sum(sum(predict)))
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    %     saveas(gcf,[filename,'.jpeg'])
    %
    %     figure
    %     plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
    %     title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    %     filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
    %     saveas(gcf,[filename,'.jpeg'])
    
    % calculate error
    predict = predict/sum(sum(predict));
    result = result/sum(sum(result));
    diff = (predict-result);
    % L2 error
    err = sqrt(sum(sum(diff.*diff)))/sqrt(sum(sum(result.*result)));
%     % L1 error
%     err = sum(sum(abs(diff)));
    disp(err)
    errvec(j) = err;
    countvec(j) = count;
    % elementwise error
    error_el = abs(diff)./(result+0.0001);
    figure
    imagesc(error_el)
    title(['relative error at each grid, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    colorbar
    
end
% title(['Energy generated using GMM on projected h, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])

% errvec_filename = ['half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
% countvec_filename = ['half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
% save(errvec_filename,'errvec')
% save(countvec_filename,'count')
end