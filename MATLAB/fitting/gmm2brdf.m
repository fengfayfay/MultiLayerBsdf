% convert from slope domain probability to brdf*cos value
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

% filename = ['half_projected_z1',num2str(angle),'_alpha_',num2str(alpha), '_#G5.mat'];
% load(filename,'obj')

for i = 1:phinum
    for j = 1:munum
        sintheta = sqrt(1 - MU(i,j)^2);
        wo = [sintheta*cos(PHI(i,j)), sintheta*sin(PHI(i,j)), MU(i,j)];
        h1 = (wi + wo)/2;
        h = h1/norm(h1);
%         % 2d gm
%         p = pdf(obj,[h(1)/h(3),h(2)/h(3)]);
        % 3d gm
        p = pdf(obj,[h(1)/h(3),h(2)/h(3),theta]);
        
        % conditioned on this incident angle
        p = p/(1/(pi/2));

        % Jacobian
        dxdphi = -sintheta*sin(PHI(i,j))/((wi(3)+MU(i,j)));
        dxdmu = (-cos(PHI(i,j))*MU(i,j)/sintheta*(wi(3)+MU(i,j)) - (wi(1)+sintheta*cos(PHI(i,j))))/((wi(3)+MU(i,j))^2);
        dydphi = sintheta*cos(PHI(i,j))/(wi(3) + MU(i,j));
        dydmu = (-sin(PHI(i,j))*MU(i,j)/sintheta*(wi(3)+MU(i,j)) - (wi(2)+sintheta*sin(PHI(i,j))))/((wi(3)+MU(i,j))^2);
        J = [dxdphi, dxdmu; dydphi, dydmu];
        
        brdfcos(i,j)= p*det(J);
    end
end

figure
imagesc(brdfcos)
colorbar()
xlabel('mu')
ylabel('phi')
title(['brdf*cos 3d gm reflect angle=0 plane, angle=', num2str(angle),' alpha=', num2str(alpha)])