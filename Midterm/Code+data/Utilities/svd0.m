
function [A,D,F]=svd0(x)
% produce the standard reduced form SVD of p-by-n matrix x
% as x=A*diag(D)*F
%   column j of A is the loadings on the j-th factor
%   row j of F is the j-th factor
%   ordered in decreasing importance as measured by elements in D
%
 [p,n]=size(x);
 [A,D,F]=svd(x,0);
 D=diag(D)';
 if p<n
  D=D(1:p);
  F=F(:,1:p);
 end;
 F=F';

  


