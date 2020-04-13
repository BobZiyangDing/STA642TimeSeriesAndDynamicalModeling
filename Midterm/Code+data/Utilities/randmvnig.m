function [theta,v] = randmvnig(p,m,C,n,s,I)
%  Generates random sample of size I for a multivariate normal-inverse
%  gamma distribution for theta, v
%  Specifically: 
%           1/v        ~ Ga(n/2, ns/2),  or v = n.s/chi_n^2        
%           theta | v  ~ N(m, C v/s)     [ with margin theta ~ T_n(m,C) ] 
%
%  Inputs: 
%   p  --   dimension of normal
%   m  --   mean p-vector of normal
%   C  --   pxp variance matrix of conditional normal
%   n  --   inv gamma dof
%   s  --   inv gamma estimate (harmonic mean of v) 
%  Outputs: 
%   theta  -- pxI array with sampled vectors as columns
%   v      -- 1xI vector of samples v values paired with cols of theta

v     = 1./gamrnd(n/2,2/(n*s),1,I);  
e     = randn(p,I).*repmat(v/s,p,1);     % scaled normals underlying p(theta|v)
theta = repmat(m,1,I) + chol(C)'*e;  
%
