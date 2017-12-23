function fitting_halfvector_z13D(dir,alpha,angle,trainangle,trainx,trainy,...
    testangle,testx,testy,range,xnum, ynum,gaussiannumvec)

cd(dir)

x_unit = 2 * range/xnum;
y_unit = 2 * range/ynum;
result1 = zeros(xnum,ynum);
result2 = zeros(xnum,ynum);
result3 = zeros(xnum,ynum);
result4 = zeros(xnum,ynum);
for i = 1:length(testangle)
    if testangle(i) == 20
        result1 = addpoint(result1,range,testx(i),testy(i),x_unit,y_unit);
    elseif testangle(i) == 40
        result2 = addpoint(result2,range,testx(i),testy(i),x_unit,y_unit);
    elseif testangle(i) == 60
        result3 = addpoint(result3,range,testx(i),testy(i),x_unit,y_unit);
    elseif testangle(i) == 80
        result4 = addpoint(result4,range,testx(i),testy(i),x_unit,y_unit);
    else
    end
end

% close all
figure
imagesc(result1/sum(sum(result1)))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),' angle=',num2str(20)])
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];
% saveas(gcf,[filename,'.jpeg'])

figure
imagesc(result2/sum(sum(result2)))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),' angle=',num2str(40)])
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];

figure
imagesc(result3/sum(sum(result3)))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),' angle=',num2str(60)])
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];

figure
imagesc(result4/sum(sum(result4)))
ylabel('x/z')
xlabel('y/z')
colorbar()
title(['Gaussian Heightfiled mirror ray distribution, alpha=', num2str(alpha),' angle=',num2str(80)])
filename = ['halfprojected_z1',num2str(angle),'_alpha_',num2str(alpha), 'heightfield'];

%training data
train = [trainangle, trainx, trainy];

for j = 1:length(gaussiannumvec)
    %% fit mixture of Gaussians using half vector
    numGaussian = gaussiannumvec(j);
    fprintf('using %4d Gaussians\n',numGaussian)
    
    options = statset('MaxIter',500, 'Display','final','TolFun',1e-5);
    try
        obj = fitgmdist(train,numGaussian,'Options',options);
    catch exception
        disp('There was an error fitting the Gaussian mixture model')
        error = exception.message
    end
    
    %     disp(obj.mu)
    %
    filename = ['3dhalf_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    save(filename,'obj')
    
    %     filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'.mat'];
    %     load(filename,'obj')
    
    %% generate points from fitted model
    Y = random(obj,length(testangle));
    
    %% calculate energy
    predict1= zeros(xnum,ynum);
    predict2= zeros(xnum,ynum);
    predict3= zeros(xnum,ynum);
    predict4= zeros(xnum,ynum);
    for i = 1:length(Y)
        if abs(Y(i,1)-20)<=5
            predict1 = addpoint(predict1, range, Y(i,2), Y(i,3), x_unit, y_unit);
        elseif abs(Y(i,1)-40)<=5
            predict2 = addpoint(predict2, range, Y(i,2), Y(i,3), x_unit, y_unit);
        elseif abs(Y(i,1)-60)<=5
            predict3 = addpoint(predict3, range, Y(i,2), Y(i,3), x_unit, y_unit);
        elseif abs(Y(i,1)-80)<=5
            predict4 = addpoint(predict4, range, Y(i,2), Y(i,3), x_unit, y_unit);
        else
            
        end
    end
    
    %     plot(predict(xnum/2,:)/(generatenum-count), 'linewidth', 2)
    
    figure
    imagesc(predict1/sum(sum(predict1)))
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' angle~=',num2str(20),' #G=',num2str(numGaussian)])
    %     filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian)];
    %     saveas(gcf,[filename,'.jpeg'])
    figure
    imagesc(predict2/sum(sum(predict2)))
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' angle~=',num2str(40),' #G=',num2str(numGaussian)])
    
    figure
    imagesc(predict3/sum(sum(predict3)))
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' angle~=',num2str(60),' #G=',num2str(numGaussian)])
    
    figure
    imagesc(predict4/sum(sum(predict4)))
    colorbar()
    title(['Mirror distribution generated using GMM, alpha=', num2str(alpha),' angle~=',num2str(80),' #G=',num2str(numGaussian)])
    %     figure
    %     plot(predict(phinum/2,:)/(generatenum-count), 'linewidth', 2)
    %     title(['Energy generated using GMM, alpha=', num2str(alpha),' angle=',num2str(angle),' #G=',num2str(numGaussian)])
    %     filename=[num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numGaussian),'_2D'];
    %     saveas(gcf,[filename,'.jpeg'])
    
    %% calculate error
    diff1 = predict1/sum(sum(predict1))-result1/sum(sum(result1));
    diff2 = predict2/sum(sum(predict2))-result2/sum(sum(result2));
    diff3 = predict3/sum(sum(predict3))-result3/sum(sum(result3));
    diff4 = predict4/sum(sum(predict4))-result4/sum(sum(result4));
    err = sqrt(sum(sum(diff1.*diff1 + diff2.*diff2 + diff3.*diff3 + diff4.*diff4)));
    value = sqrt(sum(sum(result1/sum(sum(result1)).*result1/sum(sum(result1)) + ...
        result2/sum(sum(result2)).*result2/sum(sum(result2)) + ...
        result3/sum(sum(result3)).*result3/sum(sum(result3)) + ...
        result4/sum(sum(result4)).*result4/sum(sum(result4)))));
    err = err/value;
    fprintf('error is %4.6f\n', err)
    
%     diff_new = sort(abs(reshape(diff,10000,1)),'descend');
    %     figure
    %     plot(diff_new, 'linewidth',2)
    %     title(['error from large to small, 5 gaussian, alpha=',num2str(alpha)])
    
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

errvec_filename = ['3dhalf_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_err.mat'];
% countvec_filename = ['half_projected_z1',num2str(angle),'_angle_',num2str(alpha),'_badcount.mat'];
save(errvec_filename,'errvec')
% save(countvec_filename,'count')
end

function predict = addpoint(predict,range,x,y,x_unit,y_unit)
if abs(x)<range && abs(y)<range
    predict(ceil((x+range)/x_unit),ceil((y+range)/y_unit)) = ...
        predict(ceil((x+range)/x_unit),ceil((y+range)/y_unit)) + 1;
end
end