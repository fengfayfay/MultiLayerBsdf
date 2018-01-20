function fitting_halfvector_z13D(dir,fundir,alpha,trainx,trainy,trainangle,...
    testx,testy,testangle,range,xnum, ynum,znum,zmax,zmin,gaussiannumvec)

cd(dir)

x_unit = 2 * range/xnum;
y_unit = 2 * range/ynum;
z_unit = (zmax-zmin)/znum;
result = zeros(xnum,ynum,znum);
for i = 1:length(testangle)
    if testangle(i)<zmax && testangle(i)>zmin
        result(ceil((testx(i)+range)/x_unit),ceil((testy(i)+range)/y_unit),ceil((testangle(i)-zmin)/z_unit)) = ...
            result(ceil((testx(i)+range)/x_unit),ceil((testy(i)+range)/y_unit),ceil((testangle(i)-zmin)/z_unit)) + 1;
    end
end

result = result/sum(result(:));

% figure
% scatter3(testx, testy, testangle,'filled')
% xlabel('hx/hz')
% ylabel('hy/hz')
% zlabel('incident angle')

figure
imagesc(result(:,:,1))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~10'])
% filename = ['3d_mirror_reflect_alpha_',num2str(alpha), '~10'];
% saveas(gcf,[filename,'.jpeg'])

figure
imagesc(result(:,:,30))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~30'])
filename = ['3d_mirror_reflect_alpha_',num2str(alpha), '~30'];
% saveas(gcf,[filename,'.jpeg'])

figure
imagesc(result(:,:,50))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~50'])
filename = ['3d_mirror_reflect_alpha_',num2str(alpha), '~50'];
% saveas(gcf,[filename,'.jpeg'])

figure
imagesc(result(:,:,70))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~70'])
filename = ['3d_mirror_reflect_alpha_',num2str(alpha), '~70'];
% saveas(gcf,[filename,'.jpeg'])

figure
imagesc(result(:,:,90))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'angle~90'])
filename = ['3d_mirror_reflect_alpha_',num2str(alpha), '~90'];
% saveas(gcf,[filename,'.jpeg'])

% figure
% A = result(:,ynum/2,:)/sum(result(:));
% C = permute(A,[1 3 2]);
% imagesc(C)
% ylabel('x/z')
% xlabel('y/z')
% colorbar()
% title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'xnum/2'])
% filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), 'heightfield_yz'];
% saveas(gcf,[filename,'.jpeg'])
% 
% figure
% A = result(xnum/2,:,:)/sum(result(:));
% C = permute(A,[2 3 1]);
% imagesc(C)
% ylabel('x/z')
% xlabel('y/z')
% colorbar()
% title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),'ynum/2'])
% filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), 'heightfield_xz'];
% saveas(gcf,[filename,'.jpeg'])


%training data
train = [trainx, trainy,trainangle];

for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
    fprintf('using %4d Gaussians\n',numGaussian)
    
    % use accelearted em for fitting
    cd('accelerated_greedy_EM')
    tic
    tree = buildtree(train, 0, 0, 3, 1000);
    [W,M,R,ff,Ws,Ms,Rs] = em(train,[],numGaussian,0,0,tree);
    fprintf('\nRuntime: %.2f seconds\n', toc);
    Rnew = reshape(R', 3,3,numGaussian);
    obj = gmdistribution(M,Rnew,W');
   
%     options = statset('MaxIter',500, 'Display','final','TolFun',1e-4);
%     try
%         obj = fitgmdist(train,numGaussian,'Options',options,'start','customize');
%     catch exception
%         disp('There was an error fitting the Gaussian mixture model')
%         error = exception.message
%     end
%     
%     filename = ['3dhalf_projected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
%     save(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,length(testangle));
    
    % calculate predict matrix
    predict = zeros(xnum,ynum,znum);
    count = 0;
    for i = 1:length(Y)
        if abs(Y(i,1))<range && abs(Y(i,2))<range && Y(i,3)<zmax && Y(i,3)>zmin
              predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit),ceil((Y(i,3)-zmin)/z_unit)) = ...
                predict(ceil((Y(i,1)+range)/x_unit),ceil((Y(i,2)+range)/y_unit),ceil((Y(i,3)-zmin)/z_unit)) + 1;
        else
            count = count+1;
        end
    end
    
    predict = predict/sum(predict(:));
    
    cd(fundir)
    for i = 1:180
        e(i) = relativel2err(result(:,:,i),result1(:,:,i));
    end
    figure
    plot(e,'linewidth',2)
    title('relative l2 error for each incident angle bin')
    xlabel('bin')
    ylabel('relative l2 error')
    grid on
    
    % L2 error
    err = relativel2err(result,predict);
    fprintf('error is %4.6f\n', err)
    errvec(j) = err;
    
%     figure
%     A = predict(xnum/2,:,:);
%     C = permute(A,[2 3 1]);
%     imagesc(C)
%     colorbar()
%     title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' #G=',num2str(numGaussian),'xnum/2'])
%     filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'yz'];
%     saveas(gcf,[filename,'.jpeg'])
%     
%     figure
%     A = predict(:,ynum/2,:);
%     C = permute(A,[1 3 2]);
%     imagesc(C)
%     colorbar()
%     title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' #G=',num2str(numGaussian),'ynum/2'])
%     filename = ['3d_halfprojected_z1_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'xz'];
%     saveas(gcf,[filename,'.jpeg'])
%     
    figure
    imagesc(predict(:,:,10))
    ylabel('x/z')
    xlabel('y/z')
    colorbar()
    title(['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~10'])
    filename = ['3d_mirror_predict_predict_alpha_',num2str(alpha), 'angle~10'];
%     saveas(gcf,[filename,'.jpeg'])
    
    figure
    imagesc(predict(:,:,30))
    ylabel('x/z')
    xlabel('y/z')
    colorbar()
    title(['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~30'])
    filename = ['3d_mirror_predict_alpha_',num2str(alpha), 'angle~30'];
%     saveas(gcf,[filename,'.jpeg'])
    
    figure
    imagesc(predict(:,:,50))
    ylabel('x/z')
    xlabel('y/z')
    colorbar()
    title(['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~50'])
    filename = ['3d_mirror_predict_predict_alpha_',num2str(alpha), 'angle~50'];
%     saveas(gcf,[filename,'.jpeg'])
    
    figure
    imagesc(predict(:,:,70))
    ylabel('x/z')
    xlabel('y/z')
    colorbar()
    title(['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~70'])
    filename = ['3d_mirror_predict_predict_alpha_',num2str(alpha), 'angle~70'];
%     saveas(gcf,[filename,'.jpeg'])
    
    figure
    imagesc(predict(:,:,90))
    ylabel('x/z')
    xlabel('y/z')
    colorbar()
    title(['GMM mirror ray distribution, alpha=', num2str(alpha),'angle~90'])
    filename = ['3d_mirror_predict_predict_alpha_',num2str(alpha), 'angle~90'];
%     saveas(gcf,[filename,'.jpeg'])
    
    
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