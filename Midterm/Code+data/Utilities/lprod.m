function [DA]=lprod(D,A)
% function [DA]=lprod(D,A)
% k by n matrix A, k by 1 vector D
% result is diag(D)*A using faster element-wise multiplication
 [k,n]=size(A); 
 if (length(D)~=k)
   error('incompatible dimensions'); 
 end;
 DA=repmat(reshape(D,k,1),1,n).*A;

