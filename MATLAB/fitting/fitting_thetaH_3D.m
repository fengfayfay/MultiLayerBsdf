function [paramMid, angleValues, rayleighPoly, range]  =  fitting_thetaH_3D(input,...
    xnum, znum, trainnum, generatenum)

% fitting 3d mirror heightfield data in slope domain using a mixture of gaussians
%
%
% Input:    dir             - inputdata directory
%           fundir          - funcion location
%           alpha           - roughness
%           input           - heightfield data
%           trainnum        - number of training data
%           generatenum     - number of testing data
%           gaussiannumvec  - number of gaussians vector
%           xnum            - plotting resolution in x direction
%           ynum            - plotting resolution in y direction
%           znum            - plotting resolution in z direction
%           accelerated     - if true use accelerated em, otherwise use
%                               cutomized gmcluster
%

% initialize error and count vector

boundary_ratio = 999/1000;
zmax = pi/2;
zmin = 0;

[xdividez, ydividez, tanThetaH, angle, range] = preprocess(input,boundary_ratio);
%[xparam, xdistPoly] = rayleighPolyFit(xdividez, angle);
%[yparam, ydistPoly] = rayleighPolyFit(ydividez, angle);
[tanparam, tanDist] = rayleighPolyFit(tanThetaH, angle);
paramMid = tanparam;
rayleighPoly = tanDist;
angleValues = unique(angle);
    

%{
generatenum = 10000;
[invalidTestCount,result] = bygrid2d(finalTest,xnum,znum,range,zmax,zmin);
result = result/sum(result(:));
disp(invalidTestCount);
plot_A(result, "GaussianHeightField", "HeightField");

generatenum = 10000;
angleSample = zeros(generatenum,1);
tanSample = zeros(generatenum,1);
%tsample = 1;
%angleUnit = pi/2/length(angleValues);
disp(rayleighPoly);
for i = 1:generatenum
    ia = rand() * (zmax-zmin) + zmin; 
    b = rayleighPoly(ia);
    tandist = makedist('Rayleigh', 'b', b);
    y = random(tandist, 1);
    angleSample(i) = ia;
    tanSample(i) = y;
end
generateTest = [tanSample,angleSample];


[invalidCount,gresult] = bygrid2d(generateTest,xnum/2,znum/2,range,zmax,zmin);
gresult = gresult/sum(gresult(:));
disp(invalidCount);
plot_A(gresult, "RayleighInterp", "HeightField");

%}

%tan2brdf(rayleighPoly);

end


function [xdividez, ydividez, tanThetaH, angle, range ] = preprocess(input, boundary_ratio)

% preprocess raw data to eliminate extreme hx/hz hy/hz values for
% each incident angle
angle = input(:,4);
disp(length(angle));

incident = [sin(angle), zeros(length(angle),1), cos(angle)];
h = input(:,1:3) + incident;
%h = (input(:,1:3) + incident)/2;
%hnorm = repmat(sqrt(sum(h.^2,2)),1,3);
%h = h./hnorm;
xdividez = h(:,1)./h(:,3);
ydividez = h(:,2)./h(:,3);
xsq = times(xdividez, xdividez);
ysq = times(ydividez, ydividez);
disp(length(xsq));
tanThetaHSQ = xsq + ysq;
tanThetaH = sqrt(tanThetaHSQ);
%tanThetaH = xdividez;
disp(length(tanThetaH));
 
sortedtan = sort(tanThetaH);
tanrange = sortedtan(ceil(boundary_ratio*length(tanThetaH)) - 1);
disp(tanrange);
range = tanrange;
disp(range);
end

function [paramMid, rayleighPoly] = rayleighPolyFit(ydata, angle)

angleSamples = 300;
angleValues = zeros(angleSamples, 1);
angleUnit = (pi/2) /length(angleValues);
for i = 1:angleSamples
    angleValues(i) = 0 + (i-1) * angleUnit;
end
angleValues = unique(angle);
angleUnit = (pi/2) /length(angleValues) * .5;
    
paramLower = zeros(length(angleValues), 1);
paramHigher = zeros(length(angleValues), 1);

trainnum = floor(length(ydata) / 3);
trainnum = trainnum * 2;
disp(trainnum);
dataTrain = ydata(1:trainnum);
angleTrain = angle(1:trainnum);

for i = 1:length(angleValues)
    train = dataTrain(abs(angleTrain - angleValues(i)) < angleUnit);
    disp(length(train));
    pd = fitdist(train, 'Rayleigh'); 
    disp(pd);
    tan2brdf(pd, angleValues(i));
    ci = paramci(pd, 'Alpha', .05);
    paramLower(i) = ci(1);
    paramHigher(i) = ci(2);
end

paramMid = paramLower + paramHigher;
paramMid = paramMid * .5;
rayleighPoly = fit(angleValues, paramMid, 'poly4');

end

function [count,result] = bygrid2d(A,xnum,znum,range,zmax,zmin)
x_unit = range/xnum;
z_unit = (zmax - zmin)/znum;
count = 0;
result = zeros(xnum, znum);
for i = 1:length(A)
    if A(i,1)<range && A(i, 1) > 0 && A(i,2)<zmax && A(i,2)>zmin
        xc = ceil((A(i,2)-zmin)/z_unit);
        yc = ceil(A(i,1)/x_unit);
        result(yc, xc) = result(yc, xc) + 1;
    else
        count = count + 1;
    end
end

end

function plot_A(A,titlestring,filename)

figure
imagesc(A)
ylabel('tanThetaH')
xlabel('incidentAngle')
colorbar()
title(titlestring)
%saveas(gcf,[filename,'.jpeg'])

end

function plot_error_by_angle(result,result2,predict,znum,filename)

    e = zeros(znum,1);
    eself = zeros(znum,1);
    for i = 1:znum
        e(i) = relativel2err(result(:,i),predict(:,i));
        eself(i) = relativel2err(result(:,i),result2(:,i));
    end

figure
plot(e,'linewidth',2)
hold on
plot(eself,'linewidth',2)
title('relative l2 error for each incident angle bin')
xlabel('angle bin')
ylabel('relative l2 error')
grid on
legend('fitting error','self error')
saveas(gcf,[filename,'.jpeg'])
end
