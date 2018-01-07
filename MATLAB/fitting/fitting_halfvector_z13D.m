function fitting_halfvector_z13D(dir,fundir,alpha,trainx,trainy,trainangle,...
    testx,testy,testangle,range,xnum, ynum,znum,gaussiannumvec)

cd(dir)

x_unit = 2 * range/xnum;
y_unit = 2 * range/ynum;
z_unit = pi/2/znum;
result = zeros(xnum,ynum,znum);
for i = 1:length(testangle)
    if testangle(i)<pi/2 && testangle(i)>0
        result(ceil((testx(i)+range)/x_unit),ceil((testy(i)+range)/y_unit),ceil(testangle(i)/z_unit)) = ...
            result(ceil((testx(i)+range)/x_unit),ceil((testy(i)+range)/y_unit),ceil(testangle(i)/z_unit)) + 1;
    end
end

figure
scatter3(testx(1:100000), testy(1:100000), testangle(1:100000),'filled')
xlabel('hx/hz')
ylabel('hy/hz')
zlabel('incident angle')

figure
imagesc(result(:,:,80)/sum(result(:)))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'xnum/2'])


% close all
figure
A = result(:,ynum/2,:)/sum(result(:));
C = permute(A,[1 3 2]);
imagesc(C)
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'xnum/2'])
filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), 'heightfield_yz'];
saveas(gcf,[filename,'.jpeg'])

figure
A = result(xnum/2,:,:)/sum(result(:));
C = permute(A,[2 3 1]);
imagesc(C)
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'ynum/2'])
filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), 'heightfield_xz'];
saveas(gcf,[filename,'.jpeg'])


%training data
train = [trainx, trainy,trainangle];

% use accelearted em for fitting
t=cputime;
tree = buildtree(train, 0, 0, 3, 10000);
[W,M,R,ff,Ws,Ms,Rs] = em(train,[],20,1,1,tree);
fprintf('\nRuntime: %.2f seconds\n', cputime-t);
Rnew = reshape(R', 3,3,20);
obj = gmdistribution(M,Rnew,W');

for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
%     fprintf('using %4d Gaussians\n',numGaussian)
%     
%     options = statset('MaxIter',500, 'Display','final','TolFun',1e-4);
%     try
%         obj = fitgmdist(train,numGaussian,'Options',options,'start','customize');
%     catch exception
%         disp('There was an error fitting the Gaussian mixture model')
%         error = exception.message
%     end
%     
%     %     disp(obj.mu)
%     %
%     filename = ['3dhalf_projected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
%     save(filename,'obj')
    
    %     filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    %     load(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,length(testangle));
    
    % calculate predict matrix
    predict = zeros(xnum,ynum,znum);
    for i = 1:length(Y)
        if abs(Y(i,1))<range && abs(Y(i,2))<range && Y(i,3)<pi/2 && Y(i,3)>0
              predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit),ceil(Y(i,3)/z_unit)) = ...
                predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit),ceil(Y(i,3)/z_unit)) + 1;
        end
    end
    
    figure
    A = predict(xnum/2,:,:)/sum(predict(:));
    C = permute(A,[2 3 1]);
    imagesc(C)
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' #G=',num2str(numGaussian),'xnum/2'])
    filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'yz'];
    saveas(gcf,[filename,'.jpeg'])
    
    figure
    A = predict(:,ynum/2,:)/sum(predict(:));
    C = permute(A,[1 3 2]);
    imagesc(C)
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' #G=',num2str(numGaussian),'ynum/2'])
    filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'xz'];
    saveas(gcf,[filename,'.jpeg'])
    
    % calculate error
    predict = predict/sum(predict(:));
    result = result/sum(result(:));
    % L2 error
    cd(fundir)
    err = relativel2err(result,predict);
    fprintf('error is %4.6f\n', err)
    errvec(j) = err;
    
    
end
% title(['Energy generated using GMM on projected h, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])

% errvec_filename = ['3dhalf_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
% countvec_filename = ['half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
% save(errvec_filename,'errvec')
% save(countvec_filename,'count')
end