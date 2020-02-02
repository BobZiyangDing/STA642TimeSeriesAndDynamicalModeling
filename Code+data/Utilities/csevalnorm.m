function pdf=csevalnorm(xy,m,s)
% evaluate the multivariate normal in d-dims at n points in 
% input dxn matrix xy
% m is dx1 mean vector, s is dxd variance (spd) matrix
%
[d,n]=size(xy); 
xy=xy-repmat(m,1,n);
q=diag(xy'*inv(s)*xy)'; 
a=(2*pi)^(d/2)*sqrt(det(s));
pdf=exp(-q/2)/a;
% 


