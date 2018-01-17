% One Layer Gaussian fitting
close all

mirror = true;

if mirror
    % mirror dir
%     datadir = '/Users/mandy/Github/MultiLayerBsdf/build_clang';
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
else
    % glass dir
    % datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting';

cd(datadir)

trainnum = 1e6;
generatenum = 1e7;
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
alpharange = 1:length(alphavec);
trainangle = [];
trainx = [];
trainy = [];
trainz = [];
trainweight = [];
testangle = [];
testx = [];
testy = [];
testz = [];
testweight = [];
plotrange = 0;
for k = 4
    alpha = alphavec(k);
    filename = ['3d_outputx_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C1 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputy_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C2 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputz_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C3 = textscan(fileID,'%f');
    fclose(fileID);
    
    filename = ['3d_outputangle_', num2str(alpha),'.txt'];
    fileID = fopen(filename);
    C4 = textscan(fileID,'%f');
    fclose(fileID);
    
    x = C1{1};
    y = C2{1};
    z = C3{1};
    angle = C4{1};
    observe = 6000;
    x = x/observe;
    y = y/observe;
    z = z/observe;
    epsilon = 1e-6;
    
    if mirror
        % for mirror
        if sum(z<0)/length(z)>0.01
            msg = 'More than 1% rays go down, not mirror data.';
            error(msg)
        end
        % only take z>= data
        x = x(z>=0);
        y = y(z>=0);
        z = z(z>=0);
        angle = angle(z>=0);
        
%         % use tan angle
%         angle = tan(angle);
        
        % use inverr function
%         angle = sqrt(2)*erfinv(4*angle/pi - 1);
        
        % preprocess raw data to eliminate extreme hx/hz hy/hz values for
        % each incident angle
        incident = [sin(angle), zeros(length(x),1), cos(angle)];
        h = ([x, y, z] + incident)/2;
        
        hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
        h = h./hnorm;
        xdividez = h(:,1)./h(:,3);
        ydividez = h(:,2)./h(:,3);
        sortedxz = sort(abs(xdividez));
        xrange = sortedxz(ceil(99/100*length(xdividez)));
        sortedyz = sort(abs(ydividez));
        yrange = sortedyz(ceil(99/100*length(xdividez)));
        cutoff = max(xrange, yrange);
        range = 2*cutoff;
        disp(range)
        plotrange = max(plotrange, range);
        
        xdividez_train = xdividez(1:trainnum);
        ydividez_train = ydividez(1:trainnum);
        angle_train = angle(1:trainnum);
        
        
        zmax = pi/2;
        zmin = 0;
        % use tangent theta
%         zmax = ceil(max(angle));
%         zmin = 0;
%         fprintf('max tan angle is is %4.6f\n', max(tan(angle)))
        % use inverse of error function
%         zmax = ceil(max(angle));
%         zmin = floor(min(angle));
        trainindex = abs(xdividez_train)<=cutoff&abs(ydividez_train)<=cutoff;
        xdividez_train_new = xdividez_train(trainindex);
        ydividez_train_new = ydividez_train(trainindex);
        angle_train_new = angle_train(trainindex);
        
%         % relfection around normal incidence 
%         xdividez_train_new = [xdividez_train_new;xdividez_train_new];
%         ydividez_train_new = [ydividez_train_new;ydividez_train_new];
%         angle_train_new = [angle_train_new;-angle_train_new];
        
        xdividez_test = xdividez(trainnum+1:trainnum+generatenum);
        ydividez_test = ydividez(trainnum+1:trainnum+generatenum);
        testindex = abs(xdividez_test)<=range & abs(ydividez_test)<=range;
        angle_test = angle(trainnum+1:trainnum+generatenum);
        xdividez_test_new = xdividez_test(testindex);
        ydividez_test_new = ydividez_test(testindex);
        angle_test_new = angle_test(testindex);
        
%         xdividez_test_new = [xdividez_test_new;xdividez_test_new];
%         ydividez_test_new = [ydividez_test_new;ydividez_test_new];
%         angle_test_new = [angle_test_new;-angle_test_new];
               
        
    else% glass case
        
    end% check mirror or glass end
    
end

cd(fundir);
% fitting using x/z,y/z of halfvector
gaussiannumvec = 50;
xnum = 100;
ynum = 100;
znum = 90;
fitting_halfvector_z13D(datadir,fundir,alpha,xdividez_train_new,ydividez_train_new,angle_train_new,...
    xdividez_test_new,ydividez_test_new,angle_test_new,plotrange, xnum, ynum, znum, zmax,zmin,gaussiannumvec);

%     %% glass fitting
%     gaussiannumvec = 5;=
%     incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
%     xnum = 100;
%     ynum = 100;
%     ior = 1.5;
%     Glass_fitting_halfvector_z1(datadir,alpha,angle,x,y,z,...
%         testafter, trainnum, generatenum, gaussiannumvec, incident, xnum, ynum,ior)