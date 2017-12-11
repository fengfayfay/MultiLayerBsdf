function fitting_halfvector(dir,alpha,angle,x,y,z,weight,epsilon,...
    testafter, trainnum, generatenum, gaussiannumvec, incident, phinum, thetanum)

cd(dir)

errvec = zeros(1,length(gaussiannumvec));
countvec = zeros(1,length(gaussiannumvec));

% Energy value from data
phi_unit = 2*pi/phinum;
theta_unit = pi/2/thetanum;
result = zeros(phinum,thetanum);
for i = testafter+1:length(x)
    if z(i)>=0
        h = ([x(i),y(i),z(i)] + incident)/2;
        h = h/norm(h);
        
        if abs(h(1))<epsilon && h(2)>=0
            phi = pi/2;
        elseif (abs(h(1))<epsilon && h(2)<0)
            phi = 3*pi/2;
        else
            phi = atan(h(2)/h(1));
            if h(1)< 0
                phi = phi + pi;
            else
                if h(2) < 0
                    phi = phi + 2*pi;
                end
            end
        end
        result(ceil(phi/phi_unit),ceil(acos(h(3))/theta_unit)) = result(ceil(phi/phi_unit),ceil(acos(h(3))/theta_unit)) + weight(i);
    end
end

figure
imagesc(result/generatenum)
colorbar()
title(['Energy of Gaussian Heightfiled, alpha=', num2str(alpha),' angle=',num2str(angle)])
filename = [num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
saveas(gcf,[filename,'.jpeg'])

figure
plot(result(200,:)/generatenum, 'linewidth', 2)
title(['Energy of Gaussian Heightfiled, alpha=', num2str(alpha),' angle=',num2str(angle)])
filename = [num2str(angle),'_alpha_',num2str(alpha), 'heightfield_2D'];
saveas(gcf,[filename,'.jpeg'])


% close all
% figure
% plot(result(phinum/2,:)/generatenum, 'linewidth', 2)
% hold on
for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);

    phivec = zeros(trainnum,1);
    thetavec = zeros(trainnum,1);
    for i = 1:trainnum
        if z(i)>=0
            h = ([x(i),y(i),z(i)] + incident)/2;
            h = h/norm(h);
            
            if abs(h(1))<epsilon && h(2)>=0
                phi = pi/2;
            elseif (abs(h(1))<epsilon && h(2)<0)
                phi = 3*pi/2;
            else
                phi = atan(h(2)/h(1));
                if h(1)< 0
                    phi = phi + pi;
                else
                    if h(2) < 0
                        phi = phi + 2*pi;
                    end
                end
            end
            
            theta = acos(h(3));
            phivec(i) = phi;
            thetavec(i) = theta;
        end
    end
    train = [phivec, thetavec];
    options = statset('MaxIter',500, 'Display','final');
    obj = gmdistribution.fit(train,numGaussian,'Options',options);
    
    filename = ['half',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
%     %% visualize result
%     figure
%     scatter(train(:,1),train(:,2),10,'.')
%     hold on
%     h = ezcontour(@(x,y)pdf(obj,[x y]),[0 2*pi],[0 1]);
    
    %% generate points from fitted model
    Y = random(obj,generatenum);
    
    %% calculate energy
    predict= zeros(phinum,thetanum);
    count = 0;
    for i = 1:generatenum
        if Y(i,1)>=0 && Y(i,1)<=2*pi && Y(i,2)>=0 && Y(i,2)<=pi/2
            predict(ceil(Y(i,1)/phi_unit),ceil(Y(i,2)/theta_unit)) = ...
                predict(ceil(Y(i,1)/phi_unit),ceil(Y(i,2)/theta_unit)) + 1;
        else
            count = count +1;
        end
    end
    
%     plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict/(generatenum-count))
    colorbar()
    title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
%     filename = [num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
%     saveas(gcf,[filename,'.jpeg'])
    
    figure
    plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
    title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
    saveas(gcf,[filename,'.jpeg'])
    
    %% calculate error
    diff = predict/(generatenum-count)-result/generatenum;
    err = sqrt(sum(sum(diff.*diff)));
    
    errvec(j) = err;
    countvec(j) = count;
    
    
end
% title(['Energy generated using GMM on htheta, alpha=', num2str(alpha),' angle=',num2str(angle),'compare'])
% legendInfo{1} = ['heightfield'];
% for i = 1:length(gaussiannumvec)
%     legendInfo{i+1} = [num2str(gaussiannumvec(i)),'G'];
% end
% legend(legendInfo)
% filename=['half',num2str(angle),'_alpha_',num2str(alpha),'_2Dcompare'];
% saveas(gcf,[filename,'.jpeg'])
% 
% errvec_filename = ['half',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
% countvec_filename = ['half',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
% save(errvec_filename,'errvec')
% save(countvec_filename,'count')
end