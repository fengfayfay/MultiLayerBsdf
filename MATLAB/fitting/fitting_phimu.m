function fitting_phimu(dir,alpha,angle,x,y,z,weight,epsilon,...
    testafter, trainnum, generatenum, gaussiannumvec, phinum, munum)

cd(dir)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% Energy value from data
phi_unit = 2*pi/phinum;
mu_unit = 1/munum;
result = zeros(phinum,munum);
for i = testafter+1:length(x)
    if z(i)>=0
        if abs(x(i))<epsilon && y(i)>=0
            phi = pi/2;
        elseif (abs(x(i))<epsilon && y(i)<0)
            phi = 3*pi/2;
        else
            phi = atan(y(i)/x(i));
            if x(i)< 0
                phi = phi + pi;
            else
                if y(i) < 0
                    phi = phi + 2*pi;
                end
            end
        end
        result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) = result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) + weight(i);
    end
end

% training data
phivec = zeros(trainnum,1);
muvec = zeros(trainnum,1);
for i = 1:trainnum
    if z(i)>=0
        if abs(x(i))<epsilon && y(i)>=0
            phi = pi/2;
        elseif (abs(x(i))<epsilon && y(i)<0)
            phi = 3*pi/2;
        else
            phi = atan(y(i)/x(i));
            if x(i)< 0
                phi = phi + pi;
            else
                if y(i) < 0
                    phi = phi + 2*pi;
                end
            end
        end
        phivec(i) = phi;
        muvec(i) = z(i);
    end
end
train = [phivec, muvec];


% close all
% figure
% plot(result(200,:)/generatenum, 'linewidth', 2)
% hold on
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians
    numGaussian = gaussiannumvec(j);
    
    options = statset('MaxIter',500, 'Display','final');
    obj = gmdistribution.fit(train,numGaussian,'Options',options);
    
    filename = [num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
    %% visualize result
    %     figure
    %     scatter(train(:,1),train(:,2),10,'.')
    %     hold on
    %     h = ezcontour(@(x,y)pdf(obj,[x y]),[0 2*pi],[0 1]);
    
    
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
    %% calculate energy
    predict= zeros(phinum,munum);
    count = 0;
    for i = 1:generatenum
        if Y(i,1)>=0 && Y(i,1)<=2*pi && Y(i,2)>=0 && Y(i,2)<=1
            predict(ceil(Y(i,1)/phi_unit),ceil(Y(i,2)/mu_unit)) = ...
                predict(ceil(Y(i,1)/phi_unit),ceil(Y(i,2)/mu_unit)) + 1;
        else
            count = count +1;
        end
    end
    
%     plot(predict(200,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict/generatenum)
    colorbar()
    title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename = [num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    saveas(gcf,[filename,'.jpeg'])
    figure
    imagesc(result/generatenum)
    colorbar()
    title(['Energy of Gaussian Heightfiled, alpha=', num2str(alpha),' angle=',num2str(angle)])
    filename = [num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
    saveas(gcf,[filename,'.jpeg'])
    
    figure
    plot(predict(200,:)/generatenum, 'linewidth', 2)
    title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
    saveas(gcf,[filename,'.jpeg'])
    figure
    plot(result(200,:)/generatenum, 'linewidth', 2)
    title(['Energy of Gaussian Heightfiled, alpha=', num2str(alpha),' angle=',num2str(angle)])
    filename = [num2str(angle),'_alpha_',num2str(alpha), 'heightfield_2D'];
    saveas(gcf,[filename,'.jpeg'])
    
    %% calculate error
    diff = predict/(generatenum-count)-result/generatenum;
    err = sqrt(sum(sum(diff.*diff)));
    
    errvec(j) = err;
    countvec(j) = count;
    
%     % visualize pdf
%     figure
%     fsurf(@(x,y)pdf(obj,[x y]),[0 2*pi 0 1])
%     colorbar()
    
end
% title(['Energy generated using GMM on phimu, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(guassiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=[num2str(angle),'_angle_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])

errvec_filename = [num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
countvec_filename = [num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
save(countvec_filename,'count')
end