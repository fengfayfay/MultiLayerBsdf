function [Mus, Covs, Ws, parent] = find_component(part,W,sigma,nr_of_cand)

% [Mus, Covs, Ws, parent] = find_component(part, W, sigma, nr_of_cand)
%
% finds the best new component to insert in the mixture
%
% -- gets:
% part - partition
% W - weights
% sigma - measure of data variance
% nr_of_cand - number of candidates form each component
%
% -- returns:
% Ws, Mus, Covs - new component parameters
% parent - index of parent component
%
% Jan Nunnink, 2003
%

k       = length(W);
d       = length(part{1}.centroid);
THRESHOLD = 1e-2;      % threshold in relative loglikelihood improvement for convergence in sparse EM 
B       = length(part);

% select for each box the 'owner' component
% NOTE: this only works if part is exactly the last partition
%       used in any global EM steps. (else posteriors arent up-to-date)
q = zeros(B, k);
lb = zeros(B, 1);
nb = zeros(B, 1);
for b=1:B
   [q(b,:),lb(b),nb(b)]= get_posterior(part{b}, W);
end
[tmp, I] = max2(q,2);
N = sum(nb);

par=[];
Mus=[];Covs=[];Ws=[];

for i=1:k
   
   AI        = find(I==i);
   nonAI     = find(I~=i);

   start     = size(Mus,1);    
   j=0;
   nA = length(AI);
   
   if nA > 2*d+1               % generate candidates for this parent only if it owns enough boxes
      while j < nr_of_cand    % number of candidates per parent component
         r = randperm(nA); r = r(1:2);
         
         % split AI into 2 parts by calculating the euclidian distance
         % between chosen boxes' centroids and other centroids
         
         p{1}=[];p{2}=[];m{1}=[];m{2}=[];np{1}=[];np{2}=[];
         for c=1:nA
            d1 = sum((part{AI(c)}.centroid - part{AI(r(1))}.centroid).^2);
            d2 = sum((part{AI(c)}.centroid - part{AI(r(2))}.centroid).^2);
            if d1>d2
               p{1} = [p{1} c];
               m{1} = [m{1}; part{AI(c)}.centroid];
               np{1} = [np{1}; part{AI(c)}.numpoints]; 
            else
               p{2} = [p{2} c];
               m{2} = [m{2}; part{AI(c)}.centroid];
               np{2} = [np{2}; part{AI(c)}.numpoints];
            end
         end
         
         % compute mean, covariance and weight of new candidates
         % weight = parent weight / 2
         
         for guy = 1:2
            lp = length(p{guy});
            if lp > d+1 & j<nr_of_cand
               j = j + 1;
               par = [par i];			% parent
               Mus = [Mus; sum(m{guy}.*(np{guy}*ones(1,d)))/sum(np{guy})];
               Rloc = zeros(d,d);
               for h=1:lp
                  Rloc = Rloc + np{guy}(h)*(part{AI(p{guy}(h))}.cov - ...
                     Mus(end,:)'*Mus(end,:));
               end
               Rloc = Rloc/sum(np{guy});
               Covs = [Covs; Rloc(:)'];
               Ws = [Ws W(i)/2];
            end
         end
      end
   end
   
   % do sparse incremental EM
   
   for j=start+1:length(Ws)
      iter=0;oldlogl=-realmax;done=0;
      while ~done
         % E step
         qbn=0;qbm=0;qbc=0;logl=0;
         for c=1:nA
            a_e_logphi = Ws(j)*exp(opt_approx(part{AI(c)}, Mus(j,:), Covs(j,:))); 
            % average likelihood of complete mixture of box:
            ll = a_e_logphi + (1-Ws(j))*(exp(part{AI(c)}.ll(1,:))*W);
            ll(find(ll==0))=realmin;
            logl = logl + log(ll)*part{AI(c)}.numpoints;
            qb(c) = a_e_logphi / ll; 
            qbn = qbn + qb(c)*part{AI(c)}.numpoints;
            qbm = qbm + qb(c)*part{AI(c)}.numpoints*part{AI(c)}.centroid;
            qbc = qbc + qb(c)*part{AI(c)}.numpoints*part{AI(c)}.cov;
         end
         % M step
         Ws(j) = qbn/N;
         Mus(j,:) = qbm/qbn;
         Rj = qbc/qbn - Mus(j,:)'*Mus(j,:);
         Covs(j,:) = Rj(:)';
         % convergence check
         done = abs(logl/oldlogl-1)<THRESHOLD;
         oldlogl=logl;
         if iter>30; done=1; end; iter=iter+1;
      end
      fll(j) = logl + sum(log((1-Ws(j))*lb(nonAI)).*nb(nonAI));
   end
   %fprintf('I');
end

% select best candidate

I=[];
for i=1:length(Ws) % remove some candidates that are unwanted
   S=reshape(Covs(i,:),d,d);S=min(eig(S));
   if (S<sigma/400 | Ws(i)<2*d/N  | Ws(i)>.99) I=[I i]; end
end

par(I)=[];
Ws(I)=[];
fll(I)=[];
Mus(I,:)=[];
Covs(I,:)=[];
parent=[];

if isempty(Ws)
   Ws=0;
else
   [logl sup]=max(fll);sup=sup(1);
   parent=par(sup);
   Mus=Mus(sup,:); Covs=Covs(sup,:);Ws=Ws(sup);
end


