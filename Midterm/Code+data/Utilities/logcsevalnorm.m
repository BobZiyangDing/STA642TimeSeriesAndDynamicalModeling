function logpdf=logcsevalnorm(xy,m,s)
% evaluate the log of the multivariate normal in d-dims at n points in 
% input dxn matrix xy
% m is dx1 mean vector, s is dxd variance (spd) matrix
%
[d,n]=size(xy); 
xy=xy-repmat(m,1,n);
q=diag(xy'*inv(s)*xy)'; 
logpdf=-(d*log(2*pi)+log(det(s))+q)/2;
% 

