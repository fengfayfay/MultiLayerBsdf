% test reciprocity on heightfield simulation data

datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
alpha = 0.5;
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

% put data point into grids
phinum = 100;
munum = 100;
anglenum = 90;

phi_unit = 2*pi/phinum;
mu_unit = 1/munum;
angle_unit = pi/2/anglenum;

result = zeros(phinum,munum,anglenum);
for i = 1:length(x)
    if z(i)>0
        phi = atan2(y(i),x(i));
        if phi<0
            phi = phi + 2*pi;
        end
        result(ceil(phi/phi_unit),ceil(z(i)/mu_unit),ceil(angle(i)/angle_unit)) = ...
            result(ceil(phi/phi_unit),ceil(z(i)/mu_unit),ceil(angle(i)/angle_unit)) + 1;
    end
end


%% checking by random phi,mu,theta values

numsample = 100;
brdf = zeros(numsample,1);
brdf1 = zeros(numsample,1);

for i = 1:numsample
    % wi as incident approximate brdf value
    
    theta = pi/2*rand;
    phi = 2*pi*rand;
    mu = rand;
    
    angleindex = ceil(theta/angle_unit);
    phiindex = ceil(phi/phi_unit);
    muindex = ceil(mu/mu_unit);
    
    brdfcos = result(phiindex,muindex,angleindex)/sum(sum(result(:,:,angleindex)));
    brdf(i) = brdfcos / mu;
    
    % switch wi and wo
    newangleindex = ceil(acos(mu)/angle_unit);
    newphi = 2*pi-phi;
    newmu = cos(theta);
    
    brdfcos1 = result(ceil(newphi/phi_unit),ceil(newmu/mu_unit),newangleindex)/...
        sum(sum(result(:,:,newangleindex)));
    brdf1(i) = brdfcos1 / newmu;
    
end

figure
plot(abs(brdf-brdf1),'linewidth',2);
title('abs(brdf-brdf1), check reciprocity')
grid on

figure
plot(abs(brdf-brdf1)./brdf,'linewidth',2);
title('abs(brdf-brdf1)/brdf, check reciprocity')
grid on


%% checking by fix phi,mu,theta possible values
numsample = 100;

mu = (linspace(0,munum,munum+1)+0.5)/munum;
mu = mu(1:munum);
phi = 2*pi*(linspace(0,phinum,phinum+1)+0.5)/phinum;
phi = phi(1:phinum);
theta = pi/2*(linspace(0,anglenum,anglenum+1)+0.5)/anglenum;
theta = theta(1:anglenum);

theta_index_sample = unidrnd(anglenum,numsample,1);
phi_index_sample = unidrnd(phinum,numsample,1);
mu_index_sample = unidrnd(munum,numsample,1);

theta_sample = theta(theta_index_sample);
phi_sample = phi(phi_index_sample);
mu_sample = mu(mu_index_sample);

brdf = zeros(numsample,1);
brdf1 = zeros(numsample,1);
for i = 1:numsample
    % wi as incident, approximate brdf value
    brdfcos = result(phi_index_sample(i),mu_index_sample(i),theta_index_sample(i))/...
        sum(sum(result(:,:,theta_index_sample(i))));
    brdf(i) = brdfcos / cos(theta_sample(i));
    
    % wo as incident, approximate brdf value
    newangleindex = ceil(acos(mu_sample(i))/angle_unit);
    newphi = 2*pi-phi_sample(i);
    newmu = cos(theta_sample(i));
    
    brdfcos1 = result(ceil(newphi/phi_unit),ceil(newmu/mu_unit),newangleindex)/...
        sum(sum(result(:,:,newangleindex)));
    brdf1(i) = brdfcos1 / newmu;
    
    if isnan(brdf1(i))
        fprintf('brdfcos1 value %4.4f, newmu value %4.4f\n',brdfcos1, newmu)
        fprintf('%4d %4d %4d\n',ceil(newphi/phi_unit), ceil(newmu/mu_unit),newangleindex)
    end
    
end

figure
plot(abs(brdf-brdf1)./brdf, 'linewidth',2);