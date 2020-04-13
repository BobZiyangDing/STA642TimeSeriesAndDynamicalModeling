function [x] = qnorm(p,m,v)
% inverse cdf at p of N(m,v^2)
% 
  r=p; eps=1e-12; r(r<=eps)=eps; r(r>1-eps)=1-eps; 
  x = m+sqrt(2)*v.*erfinv(2*r-1);
