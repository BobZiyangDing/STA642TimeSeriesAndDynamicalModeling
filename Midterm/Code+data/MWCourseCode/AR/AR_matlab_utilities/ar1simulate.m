function [x]=ar1simulate(phi,s,n)
% simulate a simple AR(1) process
% n observations from stationary AR(1|(phi,v)) - note input includes the marginal variance s
% rather than the innovations variance v
%
rv=sqrt((1-phi^2)*s); 
x=randn(n,1);                 % Normal innovations
%x=trnd(1.1,n,1);  

x(1)=x(1)*sqrt(s); x(2:end)=x(2:end)*rv;

for t=2:n, x(t) = phi*x(t-1)+x(t); end

plot(x); xlabel('Time t'); ylabel('x_t')


% the above is efficient - try simulating longer series the "obvious" and inefficient
%   way that calls rv generation each step ... e.g., compare cputimes the two ways ...

%a=cputime; x=randn(10000,1); for t=2:10000, x(t)=0.9*x(t-1)+x(t); end; cputime-a
%a=cputime; x=randn(1,1);     for t=2:10000, x=[x;0.9*x(end)+randn(1,1)]; end; cputime-a
 

