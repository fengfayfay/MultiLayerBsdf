%
% convert from tan probability to brdf*cos value
%

function tan2brdf(pd, angle)


    munum = 100;
    phinum = 400;

    scale = 2.0 * pi/(munum * phinum);

    mu = (linspace(0,munum,munum+1)+0.5)/munum;
    mu = mu(1:munum);
    phi = 2*pi*(linspace(0,phinum,phinum+1)+0.5)/phinum;
    phi = phi(1:phinum);
    [MU,PHI] = meshgrid(mu,phi);
    theta = angle;
    
    wi = [sin(theta), 0, cos(theta)];
    brdfcos = zeros(phinum,munum);
    psq = pi;
    
    for i = 1:phinum
        for j = 1:munum
            sintheta = sqrt(1 - MU(i,j)^2);
            wo = [sintheta*cos(PHI(i,j)), sintheta*sin(PHI(i,j)), MU(i,j)];
            h1 = (wi + wo);
            h = h1/norm(h1);
            hxz = h(1)/h(3);
            hyz = h(2)/h(3);
            tantest = sqrt(hxz * hxz + hyz *hyz);
            cosThetaH = h(3) * h(3);
            cosThetaH = cosThetaH * cosThetaH;
            p = pdf(pd, tantest) / (psq *cosThetaH);
            
            detJ = 1/(4 * wi(3));
            brdfcos(i,j)= p*detJ * scale;
        end
    end
    
    figure
    imagesc(brdfcos)
    colorbar()
    xlabel('mu')
    ylabel('phi')
    title(['brdf*cos angle=',num2str(angle)])
    
    fprintf('at angle %2d, %4.4f energy is perserved\n',angle,sum(brdfcos(:)));
end %function
