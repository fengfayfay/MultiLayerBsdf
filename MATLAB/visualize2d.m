%
%   Visualize fixed incident angle heightfield data from pbrt heightfield experiment
%

close all
clear
clc

mirror = true;

if mirror
    dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror/pi:3/output';
else
    dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output';
end
cd(dir)

alpha = 0.5;    % roughness value
angle = 60;     % incident angle in degree
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

filename = [num2str(angle),'outputdepth_', num2str(alpha),'.txt'];
fileID = fopen(filename);
C5 = textscan(fileID,'%f');
fclose(fileID);

x = C1{1};
y = C2{1};
z = C3{1};
weight = C4{1};
depth = C5{1};

% parameter for gaussian heightfield test
observe = 6000;
total = sum(weight);

x = x/observe;
y = y/observe;
z = z/observe;

if mirror
    if sum(z<0)/length(z)>0.01
        msg = 'More than 1% rays go down, not mirror data.';
        error(msg)
    end
    % only take z>=0 data
    x = x(z>=0);
    y = y(z>=0);
    z = z(z>=0);
end

% visualization parameters
munum = 100;
phinum = 400;
mu_unit = 1/munum;
phi_unit = 2*pi/phinum;
result = zeros(phinum,munum);
epsilon = 1e-6;
for i = 1:length(x)
    if z(i) >=0
        phi = atan2(y(i),x(i));
        if phi<0
            phi = phi + 2*pi;
        end
        result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) = result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) + weight(i);
    end
end

% energy plot
figure
imagesc(result/total)
if mirror
    title(['relfect energy, alpha = ',num2str(alpha),' angle = ', num2str(angle),' mirror'])
else
    title(['relfect energy, alpha = ',num2str(alpha),' angle = ', num2str(angle)])
end
xlabel('mu_o')
ylabel('phi_o')
colorbar

% brdf*cos plot
figure
imagesc(result/total*(2*pi/(phinum*munum)))
if mirror
    title(['reflect brdf*cos, alpha = ',num2str(alpha),' angle = ', num2str(angle),' mirror'])
else
    title(['reflect brdf*cos, alpha = ',num2str(alpha),' angle = ', num2str(angle)])
end
xlabel('mu_o')
ylabel('phi_o')
colorbar

if ~mirror
    % transmission
    result2 = zeros(phinum,munum);
    for i = 1:length(x)
        if z(i)<0
            if abs(x(i))<1e-6 && y(i)>=0
                theta = pi/2;
            elseif (abs(x(i))<1e-6 && y(i)<0)
                theta = 3*pi/2;
            else
                theta = atan(y(i)/x(i));
                if x(i)< 0
                    theta = theta + pi;
                else
                    if y(i) < 0
                        theta = theta + 2*pi;
                    end
                end
            end
            result2(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) = result2(ceil(theta/phi_unit),abs(floor(z(i)/mu_unit))) + weight(i);
        end
    end
    % energy plot
    figure
    imagesc(result2/total)
    title(['transmit energy, alpha = ',num2str(alpha),' angle = ', num2str(angle)])
    xlabel('mu_o')
    ylabel('phi_o')
    colorbar
    
    % brdf*cos plot
    figure
    imagesc(result2/total*(2*pi/(phinum*munum)))
    title(['transmit brdf*cos, alpha = ',num2str(alpha),' angle = ', num2str(angle)])
    xlabel('mu_o')
    ylabel('phi_o')
    colorbar
end

% hitogram of energy percentage
depth_r = depth(z>0);
yt = 0:0.1:1;
figure
histogram(depth,'Normalization','cdf');
title(['accumulated energy percentage r+t, alpha=', num2str(alpha)])
xlabel('depth')
ylabel('energy ratio')
figure
histogram(depth_r,'Normalization','cdf');
title(['accumulated energy percentage from reflectance, alpha =',num2str(alpha)])
xlabel('depth')
ylabel('energy ratio')

if ~mirror
    depth_t = depth(z<0);
    figure
    histogram(depth_t,'Normalization','cdf');
    title(['accumulated energy percentage from transmission, alpha =',num2str(alpha)])
    xlabel('depth')
    ylabel('energy ratio')
end
