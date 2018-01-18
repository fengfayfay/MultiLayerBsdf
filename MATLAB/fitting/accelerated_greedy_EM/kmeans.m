function M = kmeans(X,k)

% kmeans - clustering with k-means algorithm
%
% M = kmeans(X,kmax)
%
% X    - (n x d) d-dimensional input data
% kmax - (maximal) number of means
%
% returns
%
% M  - (k x d) matrix of cluster centers; k is computed dynamically
%

[n,d]     = size(X);

% initialize 
% use random subset of data as means
tmp    = randperm(n);
M      = X(tmp(1:k),:); 

% squared Euclidean distances to means; Dist (k x n)
Dist = sqdist(M',X');  

% Voronoi partitioning
[Dwin,Iwin] = min2(Dist',2);

% update VQ's
for i=1:size(M,1)
    I = find(Iwin==i);
    if size(I,1)>d
        M(i,:) = mean(X(I,:));
    end
end
