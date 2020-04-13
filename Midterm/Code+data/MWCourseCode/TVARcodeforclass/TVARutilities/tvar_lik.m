function [popt,delopt,likp] = tvar_lik(x,pvals,dn,bn,m0,C0,s0,n0)
% Compute lik fn for TVAR model orders 1,...,p 
% Inputs:
%   x --  series of length T
% pvals -- [pmin,pmax] -- range for model order
% ndel -- number of discount factors to consider in range 0.95-1
%  m0 --  px1 vector prior mean for state
%  C0 --  pxp prior var matrix
%  n0 --  prior df
%  s0 --  prior estimate of obs var
% Outputs: 
% likp  --  px ndel x ndel array with log-lik fn
% popt   -- MLE of model order
% delopt -- MLEs of discounts 

% organise data ...
 ndel(1)=length(dn); ndel(2)=length(bn); % AR and v discounts 
 arx=reshape(x-mean(x),length(x),1);;  
 T=length(arx);
 pmax=pvals(2); pmin=pvals(1);
 likp=zeros(pmax-pmin+1,ndel(1),ndel(2)); 
 popt=0; dopt=1; bopt=1; maxlik=-1e300;
 for p=pmin:pmax
  for i=1:ndel(1)
   d=dn(i);
   for j=1:ndel(2)
     b=bn(j);
     mt=m0(1:p); Ct=C0(1:p,1:p); st=s0; nt=n0;
     llik=0;
     for t=(p+1):T
         F=arx((t-1):-1:(t-p));
         A=Ct*F/d; q=F'*A+st; A=A/q; f=F'*mt; e=arx(t)-f; 
         nt=b*nt;
         if t>2*pmax   % ignore first 2*pmax observations for comparison
            llik=llik+gammaln((nt+1)/2)-gammaln(nt/2)...
                 -log(nt*q)/2-(nt+1)*log(1+e*e/(q*nt))/2; 
         end;
         mt=mt+A*e; r=nt+e*e/q; nt=nt+1; r=r/nt; st=st*r; Ct=r*(Ct/d-A*A'*q); 
     end;
     likp(p-pmin+1,i,j)=llik;
     if (llik>maxlik)
        popt=p; dopt=d; bopt=b; maxlik=llik;
     end;
   end;
  end;
 end;
 delopt=[dopt,bopt];
%


