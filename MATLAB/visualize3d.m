% One Layer Gaussian fitting
close all
clear

mirror = true;

if mirror
    % mirror dir
%     datadir = '/Users/mandy/Github/MultiLayerBsdf/build_clang';
        datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SinglelayerMirror_3d/';
else
    % glass dir
    % datadir = '/Users/mandy/Github/pixar/ritest/GaussianHeightField/SingleLayer/pi:3/output';
end
fundir = '/Users/mandy/Github/MultiLayerBsdf/MATLAB/fitting';

cd(datadir)

trainnum = 1e6;
alphavec = [0.1, 0.2, 0.4, 0.5, 0.7, 0.9];
alpharange = 1:length(alphavec);
plotrange = 0;
for k = 4
    alpha = alphavec(k);
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
    
%     % check if angle is consistent
%     diff = inputangle(index+1)-angle;
%     assert(sum(diff)<1e-10);
    
    if mirror
        % for mirror
        if sum(z<0)/length(z)>0.01
            msg = 'More than 1% rays go down, not mirror data.';
            error(msg)
        end
        % only take z>=0 data
        x = x(z>=0);
        y = y(z>=0);
        z = z(z>=0);
        angle = angle(z>=0);
    end
    
    % visualization parameters
    munum = 100;
    phinum = 400;
    mu_unit = 1/munum;
    phi_unit = 2*pi/phinum;
    result = zeros(phinum,munum);
    angleslice = 60;
    for i = 1:length(x)
        if angle(i)<angleslice*pi/2/90 && angle(i)>(angleslice-1)*pi/2/90
            phi = atan2(y(i),x(i));
            if phi<0
                phi = phi + 2*pi;
            end
            result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) = result(ceil(phi/phi_unit),ceil(z(i)/mu_unit)) + 1;
        end
    end
    figure
    % energy
    % imagesc(result/total)
    % brdf*cos
    imagesc(result/sum(result(:))*(2*pi/(phinum*munum)))
    title('reflection lobe')
    xlabel('mu_o')
    ylabel('phi_o')
    colorbar
end