function [AD]=rprod(A,D)
% function [AD]=rprod(A,D)
% k by n matrix A, n by 1 vector D
% result is A*diag(D) using faster element-wise multiplication
 [k,n]=size(A); 
 if (length(D)~=n)
   error('incompatible dimensions'); 
 end;
 AD=A.*repmat(reshape(D,1,n),k,1);

