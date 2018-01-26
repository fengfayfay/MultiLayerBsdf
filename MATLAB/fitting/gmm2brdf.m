%
% convert from slope domain probability to brdf*cos value
%

clear
% close all
clc

dim = 3;
numg = 50;
alpha = 0.5;
angle = 0;
theta = angle*pi/180;
wi = [sin(theta), 0, cos(theta)];
munum = 100;
phinum = 400;
mu = (linspace(0,munum,munum+1)+0.5)/munum;
mu = mu(1:munum);
phi = 2*pi*(linspace(0,phinum,phinum+1)+0.5)/phinum;
phi = phi(1:phinum);
[MU,PHI] = meshgrid(mu,phi);
brdfcos = zeros(phinum,munum);

if dim == 2
    % 2d gm data
    dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_2d/';
    filename = [dir,'half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G',num2str(numg),'.mat'];
else
    % 3d gm data
    dir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
    filename = [dir,'3dhalf_projected_z1_alpha_',num2str(alpha), '_#G',num2str(numg),'.mat'];
end
load(filename,'obj')

for i = 1:phinum
    for j = 1:munum
        sintheta = sqrt(1 - MU(i,j)^2);
        wo = [sintheta*cos(PHI(i,j)), sintheta*sin(PHI(i,j)), MU(i,j)];
        h1 = (wi + wo)/2;
        h = h1/norm(h1);
        
        if dim==2
            p = pdf(obj,[h(1)/h(3),h(2)/h(3)]);
        else
            p = pdf(obj,[h(1)/h(3),h(2)/h(3),theta]);
            % conditioned on this incident angle
            p = p/(1/(pi/2));
        end

        % Jacobian        
        detJ = (MU(i,j) * wi(3) + sintheta*sin(PHI(i,j))*wi(2) + sintheta*cos(PHI(i,j))*wi(1) + 1)/(wi(3)+MU(i,j))^3;
        brdfcos(i,j)= p*detJ;
    end
end

figure
imagesc(2*brdfcos)
colorbar()
xlabel('mu')
ylabel('phi')
title(['brdf*cos angle=',num2str(angle), ' alpha=', num2str(alpha)])