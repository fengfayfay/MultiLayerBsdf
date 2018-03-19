function [M,I] = max2(A,dim);

% calculate maximum of A over dimension dim
%
% Jan Nunnink, 2002

[d1, d2] = size(A);

if dim==1
   if d1==1
      M=A; I=ones(1,d2);
   else
      [M,I]=max(A);
   end
else
   if d2==1
      M=A; I=ones(d1,1);
   else
      [M,I]=max(A');
      M=M';
      I=I';
   end
   
end
