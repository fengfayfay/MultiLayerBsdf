close all
anglevec = [0,10,20,30,40,50,60,70,80,89]; 
W = [];
M = [];
R = [];
% fit current slice
for i = 1:length(anglevec)
    angle = anglevec(i);
    [obj_prev,W,M,R] = fit2D(angle,i>1,W,M,R);
    objcell{i} = obj_prev;
end

result = {};
numgaussian = 5;
% calculate the similarity between this slice and the previous slice
for i = 1:length(objcell)-1
    count = 1;
    simarray = zeros(numgaussian^2,3);
    for j = 1:numgaussian
        for k = 1:numgaussian
            sim = gmpairsimilarity(2, objcell{i}.mu(j,:),objcell{i+1}.mu(k,:),objcell{i}.Sigma(:,:,j),objcell{j}.Sigma(:,:,k));
            simarray(count,1) = j;
            simarray(count,2) = k;
            simarray(count,3) = sim;
            count = count+1;
        end
    end
    simarray = sortrows(simarray,-3);
    result{i} = simarray;
end