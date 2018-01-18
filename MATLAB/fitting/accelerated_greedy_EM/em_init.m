function [W,M,R,sigma] = em_init(X,k)

%  em_init - initialization of EM for Gaussian mixtures 
%
%  [W,M,R,sigma] = em_init(X,k)
%
%  X - (n x d) matrix of input data 
%  k - initial number of Gaussian components
%
%  returns
%
%  W - (k x 1) vector of mixing weights
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of components covariances in vector reshaped form
%  sigma - measure of variance
%
%  Jan Nunnink, 2003

[n,d] = size(X);

% means
M = kmeans(X,k);

[D,I]=min2(sqdist(M',X'),1);

% mixing weights
W=zeros(k,1);
for i=1:k
  W(i) = length(find(I==i))/n;  
end

% covariance matrices 
R = zeros(k,d^2);
if k > 1
  for j = 1:k
    J=find(I==j);
    if length(J)>2*d;Sj = cov(X(J,:));else Sj=cov(X);end
    R(j,:) = Sj(:)';
  end
else
  R = cov(X);
  R = R(:)';
end

sigma = 0.5 * (4/(d+2)/n)^(1/(d+4)) * sqrt(norm(cov(X)));
