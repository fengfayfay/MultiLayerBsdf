% One Layer Gaussian fitting
close all

mirror = true;

if mirror
    % mirror dir
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror/pi:3/output';
else
    % glass dir
    % datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting';

cd(datadir)

trainnum = 1e5;
generatenum = 1e6;
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
anglevec = [0, 10, 20, 30, 40, 50, 60, 70, 80 ,89];
alpharange = 1:length(alphavec);
anglerange = 1:length(anglevec);
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
    for j = 1:2:9
        angle = anglevec(j);
        filename = [num2str(angle), 'outputx_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C1 = textscan(fileID,'%f');
        fclose(fileID);
        
        filename = [num2str(angle),'outputy_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C2 = textscan(fileID,'%f');
        fclose(fileID);
        
        filename = [num2str(angle),'outputz_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C3 = textscan(fileID,'%f');
        fclose(fileID);
        
        filename = [num2str(angle),'outputweight_', num2str(alpha),'.txt'];
        fileID = fopen(filename);
        C4 = textscan(fileID,'%f');
        fclose(fileID);
        
        x = C1{1};
        y = C2{1};
        z = C3{1};
        weight = C4{1};
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
            
            % preprocess raw data to eliminate extreme hx/hz hy/hz values for
            % each incident angle
            incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
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
            plotrange = max(plotrange, range);
            
            xdividez_train = xdividez(1:trainnum);
            ydividez_train = ydividez(1:trainnum);
            xdividez_train_new = xdividez_train(abs(xdividez_train)<=cutoff&abs(ydividez_train)<=cutoff);
            ydividez_train_new = ydividez_train(abs(xdividez_train)<=cutoff&abs(ydividez_train)<=cutoff);
            
            xdividez_test = xdividez(trainnum+1:trainnum+generatenum);
            ydividez_test = ydividez(trainnum+1:trainnum+generatenum);
            xdividez_test_new = xdividez_test(abs(xdividez_test)<=range & abs(ydividez_test)<=range);
            ydividez_test_new = ydividez_test(abs(xdividez_test)<=range & abs(ydividez_test)<=range);
            
            trainangle = [trainangle; zeros(length(xdividez_train_new),1) + angle];
            testangle = [testangle; zeros(length(xdividez_test_new),1) + angle];
            
            trainx = [trainx;xdividez_train_new];
            testx = [testx; xdividez_test_new];
            
            trainy = [trainy;ydividez_train_new];
            testy = [testy; ydividez_test_new];
            
        else% glass case
            
        end% check mirror or glass end
        
    end
    
    cd(fundir);
    % fitting using x/z,y/z of halfvector
    gaussiannumvec = 1:20;
    xnum = 100;
    ynum = 100;
    fitting_halfvector_z13D(datadir,alpha,angle,trainangle,trainx,trainy,...
        testangle,testx,testy,plotrange, xnum, ynum,gaussiannumvec);
    
    %     %% glass fitting
    %     gaussiannumvec = 5;
    %     incident = [sin(angle*pi/180), 0, cos(angle*pi/180)];
    %     xnum = 100;
    %     ynum = 100;
    %     ior = 1.5;
    %     Glass_fitting_halfvector_z1(datadir,alpha,angle,x,y,z,weight,...
    %         testafter, trainnum, generatenum, gaussiannumvec, incident, xnum, ynum,ior)
end

