function [y]=subtract_mean_rows(x)
% subtract sample mean by rows of x
%
 [m,n]=size(x);
 y=x-repmat(mean(x,2),1,n);

%

