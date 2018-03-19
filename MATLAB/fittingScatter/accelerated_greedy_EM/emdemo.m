% demo Chunky Greedy EM and non-Chunky Greedy EM algorithm
%
% Jan Nunnink, 2003

% graphics stuff
figure(1);clf;

% create dataset and testset
n = 100000; k = 5; d = 2; c = 4; e = 10;
fprintf('------------\nCreating data set of %.0f points...\n',n);
[X,T,ff,genll] = mixgen(n,n,k,d,c,e);
fprintf('Log-likelihood of test set: %.4f\n', genll);

% build tree
fprintf('------------\nBuilding kd-tree...\n');
t=cputime;
tree = buildtree(X, 0, 0, 3, round(n/(50*k)));
fprintf('Runtime: %.2f seconds\n', cputime-t);

% perform greedy chunky em
fprintf('------------\nRunning accelerated greedy EM');
t=cputime;
[W,M,R,ff,Ws,Ms,Rs] = em(X,[],k,4,1,tree);
fprintf('\nRuntime: %.2f seconds\n', cputime-t);

logl=[];
for i=1:length(Ws)
   Ft = em_gauss(T,Ms{i},Rs{i}) * Ws{i};
   Ft(find(Ft < eps)) = eps;
   logl = [logl mean(log(Ft))];
end
logl
[L, I] = max(logl);
fprintf('Best log-likelihood of test set: %.4f, using %i components.\n', L, I);
fprintf('------------\n');


