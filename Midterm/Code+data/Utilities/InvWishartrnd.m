function [Sigs] = InvWishartrnd(s,S,mc);
% MW cover function for generating mc draws from Inverse Wishart posteriors  
%   s:  degree of freedom s>0 
%   S:  sum of squares parameter matrix 
% Distribution is Sigma ~ IW_p(s,S) with s>0, S=sum of squares param
% pdf:  p(Sigma) = cons. |Sigma|^{-(p+s/2)} exp(-trace(Sigma^{-1} S)/2) 
%
%   and Omega=inv(Sigma) ~ W_p(d=s+p-1,A=inv(S))
%       p(Omega) = cons. |Omega|^((d-p-1)/2) exp(-trace(Omega A^{-1})/2) 
% Returns array of sampled matrices Oms of dimension (mc,p,p) 
%
p=size(S,1); Sigs=zeros(mc,p,p);
[c,C]=iwishrnd(S,s+p-1); 
Sigs(1,:,:)=c; for i=2:mc, Sigs(i,:,:)=iwishrnd(S,s+p-1,C); end
Sigs=squeeze(Sigs); 
