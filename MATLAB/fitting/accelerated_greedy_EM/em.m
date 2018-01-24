function [W,M,R,Tlogl,mixW,mixM,mixR] = em(X,T,kmax,nr_of_cand,plo,tree)

% em - EM algorithm for adaptive multivariate Gaussian mixtures
%
% [W, M, R, Tlogl, mixW, mixM, mixR] = em(X, T, kmax, cands, plo, tree)
%
%  X     - (n x d) d-dimensional zero-mean unit-variance data
%  T     - (m x d) test data (optional, set [] if none)
%  kmax  - maximum number of components allowed
%  cands - number of candidates per component, zero gives non-greedy EM
%  plo   - if 1 then plot ellipses for 2-d data
%  tree  - multi-resolution kd-tree containing the data
%
%  returns:
%
%  W - (k x 1) vector of mixing weights
%  M - (k x d) matrix of components means
%  R - (k x d^2) matrix of components covariances in vector reshaped format.
%  Tlogl -  average log-likelihood of test data
%  mixW, mixM, mixR - array with parameters of mixture sequence
%
%  Jan Nunnink, 2003

[n,d] = size(X);

if d~=2 plo=0;
end

if isempty(T) test=0; else test=1; Tlogl=[]; end

if nr_of_cand
   k=1;
else
   k=kmax;
end

% initialize em using kmeans.
[W,M,R,sigma] = em_init(X,k);
sigma=sigma^2;

% initialize basic partition
basic_part{1} = tree;

THRESHOLD = 1e-5;        % convergence threshold
oldlogl = -realmax;

% time index for recording in which update step each component is.
time(1,1:k) = 0;

% count number of iterations
iter = 0;
while 1
   t1 = tic;
   veryoldlogl = -realmax; oldlogl = -realmax;
   part = make_partition(basic_part, 1, 5);
   while 1
      while 1
         [W,M,R,logl,part,time] = em_step_partition(X,W,M,R,0,part,1:k,time);
         
         % Print running information (Mandy)
         iter = iter+1;
         if mod(iter,40)==0 && iter>0
             fprintf('\nNumber of iterations: %4d \n', iter)
             fprintf('\nLast 40 iterations used time: %4.4f \n', toc(t1))
             F1(iter) = abs(logl/oldlogl-1);
             t1 = tic;
         end
         
         if abs(logl/oldlogl-1)<THRESHOLD break; end
         oldlogl=logl;
      end
      F2(iter) = abs(logl/oldlogl-1);
      if abs(logl/veryoldlogl-1)<THRESHOLD break; end
      veryoldlogl=logl;
      part = make_partition(part, 1, 1);
   end
   if test
      Ft = em_gauss(T,M,R) * W;
      Ft(find(Ft < eps)) = eps;
      Tlogl(length(Tlogl)+1) = mean(log(Ft));
   else
      Tlogl=0;
   end
   
   if plo
      plotmix(X,W,M,R);
      drawnow;
   end
      
   mixW{k}=W;
   mixM{k}=M;
   mixR{k}=R;
   
   fprintf('.');
   
   if k == kmax;    return;  end
   
   % find best new component
   [Mnew,Rnew,alpha,par] = find_component(part,W,sigma,nr_of_cand); 
   
   if alpha==0		% no candidates
      if test
         Tlogl = [Tlogl; mean(log(Ft))];
      else
         Tlogl=0;
      end  
      return;
   end
   
   % add new component
   M=[M;Mnew];
   R=[R;Rnew];
   W=[(1-alpha)*W; alpha];
   k=k+1;
   
   % add new component to time
   time(:,length(time)+1)=0;
   
   % sparse EM updating only new component and parent
   veryoldlogl = -realmax; oldlogl = -realmax;
   part = make_partition(basic_part, 1, 5);
   while 1
      while 1				
         [W,M,R,logl,part,time] = em_step_partition(X,W,M,R,0,part,[par k],time);
         
         if abs(logl/oldlogl-1)<THRESHOLD break; end
         oldlogl=logl;
      end
      if abs(logl/veryoldlogl-1)<THRESHOLD break; end
      veryoldlogl=logl;
      part = make_partition(part, 1, 1);
   end
end

