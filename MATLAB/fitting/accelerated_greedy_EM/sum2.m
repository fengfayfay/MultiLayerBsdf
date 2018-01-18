function S = sum2(A,dim);

% calculate sum of A over dimension dim
%
% Jan Nunnink, 2002

[d1, d2] = size(A);

if dim==1
   if d1==1
      S=A;
   else
      S=sum(A);
   end
else
   if d2==1
      S=A;
   else
      S=sum(A')';
   end
end
