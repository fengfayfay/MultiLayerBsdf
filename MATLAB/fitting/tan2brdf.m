%
% convert from tan probability to brdf*cos value
%



function tan2brdf(pd, angle)
    munum = 100;
    phinum = 400;

    scale = 2.0 * pi/(munum * phinum);
    %scale = 1;

    mu = (linspace(0,munum,munum+1)+0.5)/munum;
    mu = mu(1:munum);
    phi = 2*pi*(linspace(0,phinum,phinum+1)+0.5)/phinum;
    phi = phi(1:phinum);
    [MU,PHI] = meshgrid(mu,phi);
    theta = angle;
    
    wi = [sin(theta), 0, cos(theta)];
    brdfcosfit = zeros(phinum,munum);
    brdfcosw = zeros(phinum,munum);
    psq = pi;
    
    for i = 1:phinum
        for j = 1:munum
            sintheta = sqrt(1 - MU(i,j)^2);
            wo = [sintheta*cos(PHI(i,j)), sintheta*sin(PHI(i,j)), MU(i,j)];
            h = (wi + wo);
            h = h/norm(h);
            hxz = h(1)/h(3);
            hyz = h(2)/h(3);
            tantest2 = (hxz * hxz + hyz * hyz);
            tantest = sqrt(tantest2);
            cos2ThetaH = h(3) * h(3);
            cos4ThetaH = cos2ThetaH * cos2ThetaH;
            
            r= tantest;
            p_rphi = pdf(pd, r)/(2*pi);
            %jacobian drphi/dxy
            detJ1 = 1.0/r;
            p_xy = p_rphi*detJ1;

            %jacobian dxy/dwi
            detJ = (wo(3) * wi(3) + wo(2)*wi(2) + wo(1)*wi(1) + 1)/(wi(3)+wo(3))^3;
            brdfcosfit(i,j)= p_xy * detJ * scale;

            brdfcosw(i, j) = distribution(tantest2, cos4ThetaH, .5)/(4*wi(3)) * scale;
        end
    end
    
    figure
    imagesc(brdfcosfit)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['brdffit*cos angle=',num2str(angle * 180/pi)])
    
    figure
    imagesc(brdfcosw)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['brdfwat*cos angle=',num2str(angle * 180/pi)])
    
    fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcosfit(:)));
end %function

function [dist] = distribution(tantest2, cos4ThetaH,  alpha)
    alpha = alpha * alpha;
    dist = exp(-tantest2/alpha)/(pi * alpha * cos4ThetaH);
end
