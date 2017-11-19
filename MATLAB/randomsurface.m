function [f,x,y] = randomsurface(N,rL,h,lx,ly)

format long;

x = linspace(-rL/2,rL/2,N); y = linspace(-rL/2,rL/2,N);
[X,Y] = meshgrid(x,y); 

Z = h.*randn(N,N); % uncorrelated Gaussian random rough surface distribution
                   % with mean 0 and standard deviation h

% isotropic surface
if nargin == 4
    
    % Gaussian filter
    F = exp(-((X.^2+Y.^2)/(lx^2/2)));

    % correlation of surface including convolution, inverse
    % Fourier transform and normalizing prefactors
    f = 2/sqrt(pi)*rL/N/lx*ifft2(fft2(Z).*fft2(F));
    
% non-isotropic surface
elseif nargin == 5
    
    % Gaussian filter
    F = exp(-(X.^2/(lx^2/2)+Y.^2/(ly^2/2)));
    f = 2/sqrt(pi)*rL/N/sqrt(lx)/sqrt(ly)*ifft2(fft2(Z).*fft2(F));
    
end
figure
surf(f)
% cd('/Users/mandy/Github/MultiLayerBsdf/build');
% filename = ['pz', num2str(h/lx), '.txt'];
% pz = fopen(filename,'w');
% fprintf(pz,'%5f %5f %5f %5f %5f %5f %5f %5f\n',f);
% fclose(pz);