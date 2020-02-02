function [sd]=mcmcsd(x)
% compute sample acf of MCMC sequence x and the corresponding
%  asymptotic estimate of MC sd. 

y=acovf(x); m=length(x);
sd=y(1)+2*sum(y(2:end).*(1-(1:(m-1))/m)); 
sd=sqrt(sd/m); 
