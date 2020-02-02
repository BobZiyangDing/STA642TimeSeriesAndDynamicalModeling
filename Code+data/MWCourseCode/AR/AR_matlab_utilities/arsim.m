function [xpred,bsamp,vsamp] = arsim(x,p,b,B,nu,v,nmc,npred)
% simulates nmc realizations from the conditional reference posterior in the AR(p) model
% (results from ar.m) and generates sample predictions into the future
%
% first sample posterior for innovations variance
% 
vsamp=v*nu./chi2rnd(nu,1,nmc);  v2samp=sqrt(vsamp);
%
% then, given these values, sample b 
% 
bsamp=repmat(b,1,nmc)+chol(inv(B))'*randn(p,nmc).*repmat(v2samp,p,1);
% 
% then sequence through npred time steps .. doing predictions in parallel
%  for each sampled parameter
%

X=repmat(x(end:-1:(end-p+1)),1,nmc);  
xpred=randn(npred,nmc);                  % Normal innovations
%k=1; xpred=trnd(k,npred,nmc);            % T_k innovations

for t=1:npred,
    xpred(t,:) = sum(X.*bsamp,1)+v2samp.*xpred(t,:); 
    X=[xpred(t,:);X]; X(end,:)=[];      
end

