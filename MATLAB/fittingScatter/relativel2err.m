function [err,maxl1err,minl1err,varl1err] = relativel2err(A,B)

diff = abs(A - B);
diff2 = diff.*diff;
A2 = A.*A;
err = sqrt(sum(diff2(:))/sum(A2(:)));
maxl1err = max(diff(:));
minl1err = min(diff(:));
varl1err = var(diff(:));

end