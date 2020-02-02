function [Oms] = Wishartrnd(s,S,mc);
% MW cover function for generating mc draws from Wishart posteriors  
%   s:  degree of freedom s>0 
%   S:  sum of squares parameter matrix 
% Distribution is Omega ~ W_p(d,A) with d=s+p-1, A=inv(S)
% pdf:  p(Omega) = cons. |Omega|^(s/2-1) exp(-trace(Omega S)/2) 
%                = cons. |Omega|^((d-p-1)/2) exp(-trace(Omega A^{-1})/2) 
%   and Sigma=inv(Omega) ~ IW_p(s,S)
%
% Returns array of sampled matrices Oms of dimension (mc,p,p) 
%
p=size(S,1); Oms=zeros(mc,p,p);
A=inv(S);  A=(A+A')/2; % symmetrize to deal with rounding errors on inverse
[c,C]=wishrnd(A,s+p-1); 
Oms(1,:,:)=c; for i=2:mc, Oms(i,:,:)=wishrnd(A,s+p-1,C); end
% OR: 
% [c,C]=iwishrnd(S,s+p-1); 
% Oms(1,:,:)=inv(c); for i=2:mc, Oms(i,:,:)=inv(iwishrnd(S,s+p-1,C)); end
Oms=squeeze(Oms); 
