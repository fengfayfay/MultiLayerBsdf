close all
numgaussian = 5;
% objstruct = {obj_0,obj_10,obj_20,obj_30,obj_40,obj_50,obj_60,obj_70,obj_80,obj_89};
objstruct = objcell;
anglevec = [0,10,20,30,40,50,60,70,80,89];
maxx = -inf;
minx = inf;
maxy = -inf;
miny = inf;

for i = 1:length(anglevec)
    maxx = max(maxx, max(objstruct{i}.mu(:,1)));
    minx = min(minx, min(objstruct{i}.mu(:,1)));
    maxy = max(maxy, max(objstruct{i}.mu(:,2)));
    miny = min(miny, min(objstruct{i}.mu(:,2)));
end

loops = length(anglevec);
F(loops) = struct('cdata',[],'colormap',[]);
figure
xlim([minx maxx])
ylim([miny maxy])
hold on
for i = 1:length(anglevec)
    scatter(objstruct{i}.mu(:,1),objstruct{i}.mu(:,2),'filled');
    
    for j = 1:numgaussian
        txt = ['\leftarrow g',num2str(j)];
        text(objstruct{i}.mu(j,1),objstruct{i}.mu(j,2),txt)
    end
    
    legendInfo{i} = num2str(anglevec(i));
    pause
    F(i) = getframe;
end
legend(legendInfo)
hold off

sim = zeros(numgaussian^2,3);
count = 1;
for i = 1:numgaussian
    for j = 1:numgaussian
        sim(count,1) = i;
        sim(count,2) = j;
        sim(count,3) = gmpairsimilarity(2, objstruct{1}.mu(i,:),objstruct{2}.mu(j,:),objstruct{1}.Sigma(:,:,i),objstruct{2}.Sigma(:,:,j));
        count = count + 1;
    end
end