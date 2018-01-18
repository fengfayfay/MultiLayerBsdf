function [W,M,R,LL,part,time] = em_step_partition(X,W,M,R,plo,part,comp,time)

% em_step_partition - EM learning step for multivariate Gaussian mixtures
%
% [W,M,R,LL,part,time] = em_step_partition(X,W,M,R,plo,part,comp,time)
%
%  X - (n x d) matrix of input data
%  W - (k x 1) vector of mixing weights
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of components covariances in vector reshaped format.
%  plo - if 1 then plot ellipses for 2-d data
%  part - partition containing input data
%  comp - vector of component indices to be updated (while keeping rest fixed)
%  time - time vector
%
% returns
%
%  W - (k x 1) matrix of components priors
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of components covariances in vector reshaped format
%  LL - loglikelihood of input data under _old_ parameters
%  part - updated partition
%  time - updated time
%
%  Jan Nunnink, 2003

[n,d] = size(X);
nr = length(comp);
k = length(W);

% increase time index of the components to be updated
time(comp)=time(comp)+1;

% compute posteriori's
SW=zeros(k,1);
SWX=zeros(k,d);
SWXX=zeros(k,d^2);
LL=0;
for i=1:length(part)
   [SWb, SWXb, SWXXb, LLb, new_t] = make_stats(part{i}, W, M, R, time);
   part{i}=new_t;
   SW=SW+SWb;
   SWX=SWX+SWXb;
   SWXX=SWXX+SWXXb;
   LL=LL+LLb;
end

LL=LL/n;		% mean loglikelihood of data under mixture

% update weights of ~comp so they sum to 1
alpha = SW(comp)/n;
ff = W(comp);
W(comp)=0;
if sum(W)>eps
   W = W*((1-sum(alpha))/sum(W));
end
W(comp)=ff;

% update covariance matrix, weights and mean of comp
for i=1:nr
   if SW(comp(i))>eps
      W(comp(i)) = SW(comp(i))/n;
      
      M(comp(i),:) = SWX(comp(i),:)./(SW(comp(i))*ones(1,d));
      
      c = reshape(SWXX(comp(i),:), d, d);
      Sj = c/SW(comp(i)) - M(comp(i),:)'*M(comp(i),:);
      
      % check cov matrices for singularities
      [U,L,V] = svd(Sj); 
      l = diag(L);
      if (min(l) > eps) & (max(l)/min(l) < 1e4)
         R(comp(i),:) = Sj(:)';
      end
   end
end

if plo
   plotmix(X, W, M, R);
   drawnow;
end
