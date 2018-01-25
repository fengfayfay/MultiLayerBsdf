%
%   Visualize varying incident angle heightfield data in an angle range
%

close all
clear
clc

mirror = true;

if mirror
    datadir = '/Users/mandy/Github/MultiLayerBsdf/build_clang';
%     datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
else
    datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting';

cd(datadir)

trainnum = 1e6;
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
epsilon = 1e-6;

if mirror
    if sum(z<0)/length(z)>0.01
        msg = 'More than 1% rays go down, not mirror data.';
        error(msg)
    end
    % only take z>=0 data
    angle = angle(z>=0);
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
angle1 = 0;
angle2 = 1;
for i = 1:length(x)
    if angle(i)<angle2*pi/2/90 && angle(i)>angle1*pi/2/90
        phi = atan2(y(i),x(i));
        if phi<0
            phi = phi + 2*pi;
        end
        result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) = result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) + 1;
    end
end

% energy plot
figure
imagesc(result/length(x))
title('reflection lobe')
xlabel('mu_o')
ylabel('phi_o')
colorbar