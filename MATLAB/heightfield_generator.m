N = 2000;
rL =  10;
hvec = [0.05];
lx = 0.1;

for i = 1:length(hvec)
    h = hvec(i);
    [f,x,y] = randomsurface(N,rL,h,lx);
end