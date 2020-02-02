function [thetasamp,vsamp] = tvarFFBS(x,p,del,m0,C0,s0,n0,nmc)
% FFBS for TVAR(p) model of univariate series x
% Discount factors del(1) for state and del(2) for obs var
% Inputs:
%   x --  series of length T
%   p --  model order
%  m0 --  px1 vector prior mean for state
%  C0 --  pxp prior var matrix
%  n0 --  prior df
%  s0 --  prior estimate of obs var
%  nmc--  Monte Carlo sample size
% Outputs: 
%  thetasamp --  pxTxnmc   posterior samples of TVAR state vectors
%  vsamp     --   Txnmc       and observation variances
% 
% organise data ...
 d=del(1); b=del(2);
 arx=x-mean(x);
 T=length(arx);
 arx=reshape(arx,T,1); 
 m=repmat(m0,[1 T]);  C=repmat(C0,[1 1 T]);  
 s=repmat(s0,[1 T]);  n=repmat(n0,[1 T]);
 mt=m0; Ct=C0; st=s0; nt=n0;

% forward filtering
 for t=(p+1):T
     F=arx((t-1):-1:(t-p));
     A=Ct*F/d; q=F'*A+st; A=A/q; e=arx(t)-F'*mt; 
     mt=mt+A*e; m(:,t)=mt; 
     r=b*nt+e*e/q; nt=b*nt+1; r=r/nt; st=st*r;
     n(t)=nt; s(t)=st;   
     Ct=r*(Ct/d-A*A'*q); Ct=(Ct+Ct')/2;
     C(:,:,t)=Ct;
 end;

% backward sampling 
 thetasamp=zeros(p,T,nmc); vsamp=zeros(T,nmc); 

 % first, sample at end time T: NOTE SAMPLING ALL NMC AT ONCE - EFFICIENT
 vt     = 1./gamrnd(n(T)/2,2/(n(T)*s(T)),1,nmc);  
 e      = randn(p,nmc).*repmat(vt/s(T),p,1);     % scaled normals underlying p(theta|v)
 thetat = repmat(m(:,T),1,nmc) + chol(C(:,:,T))'*e;  
 vsamp(T,:)=vt; thetasamp(:,T,:)=thetat;         % save 

  % then, recurse back sampling each time point conditional on the previous values: 
 for t=(T-1):-1:(p+1)
     if (b<1)
         vt = b./vt + gamrnd((1-b)*n(t)/2, 2/(n(t)*s(t)), 1,nmc); vt=1./vt; 
     end
     ht =repmat((1-d)*m(:,t),1,nmc) + d*thetat;   % conditional mean vector, expanded nmc times 
     Ht =(1-d)*C(:,:,t);                          % conditional variance matrix 
     e  = randn(p,nmc).*repmat(vt/s(t),p,1);      % underlying scaled normals 
     thetat = ht + chol(Ht)'*e; 
     thetasamp(:,t,:)=thetat; vsamp(t,:)=vt;        % save 
 end;
 % now need to add something for ad-hoc treatment of 1st p values ... ? 
 for t=1:p
    thetasamp(:,t,:)=thetat; vsamp(t,:)=vt; 
 end
%



