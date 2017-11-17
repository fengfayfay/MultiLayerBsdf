N = 4800;
rL =  40;
lx = 0.1;
hvec = [0.02];


for i = 1:length(hvec)
    h = hvec(i);
    [f,x,y] = randomsurface(4800,rL,h,lx);
end