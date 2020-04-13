function [logp] = lpostAR1(theta,x)
%
% log posterior density for theta=(phi,v) in the AR(1) model for x=x_{1:n}
%  evaluated at the m points in the mx2 array theta of (phi,v) values
% Assumes stationarity - |phi|<1, v>0

n=length(x); x=reshape(x,n,1);
m=size(theta,1); phi=theta(:,1); v=theta(:,2);
logp=zeros(m,1); 

% compute only for stationary cases, all others deliver 0 so set logp=NaN
% for those: 
i=find(abs(phi)<1);  j=1:m; j(i)=[]; im=length(i);  % finds stationary cases

a=1-phi(i).*phi(i);  e=repmat(x(2:n)',im,1)-outer(phi(i),x(1:(n-1)),'*');
Q=a*x(1)^2+sum(e.*e,2);
logp(i)=-Q./(2*v(i))+log(a)/2-log(v(i))*(n/2+1); 

% the above all works fine, but will have generated complex numbers for any
% nonstationary phi values: set them to zero ..
logp(j)=NaN; 