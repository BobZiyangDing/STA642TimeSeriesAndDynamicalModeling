function  [popt,likp] = armodelik(xdata,maxp,m0,M0,s0,n0)
% Compute approximate marginal likelihood for AR model orders 0,...,maxp 
% Inputs:
%  xdata --  series of length T
%  maxp  --  max model order
%  m0 --  px1 vector prior mean for state
%  M0 --  pxp prior var matrix
%  n0 --  prior df
%  s0 --  prior estimate of obs var
% Outputs: 
% likp  --  approx marginal likihood over 0:maxp
% popt   -- MLE of model order

% organise data ...
 arx=xdata-mean(xdata);
 T=length(arx);
 likp=zeros(maxp+1,1); 
 popt=0; maxlik=-1e300;
 for p=0:maxp
     mt=m0(1:p); Mt=M0(1:p,1:p); st=s0; nt=n0;
     llik=0;
     for t=(p+1):T
         x=arx((t-1):-1:(t-p)); 
         A=Mt*x; q=x'*A+1; A=A/q; e=arx(t)-x'*mt; 
         if t>2*maxp  % use only last n-2*maxp observations to compare models: for comparability
            llik=llik+gammaln((nt+1)/2)-gammaln(nt/2)...
                 -log(q*st*nt*pi)/2-(nt+1)*log(1+e*e/(q*nt*st))/2; 
         end;
         mt=mt+A*e; st=(nt*st+e*e/q)/(nt+1); nt=nt+1; Mt=Mt-A*A'*q; 
     end;
     likp(p+1)=llik;
     if (llik>maxlik)     
	    popt=p; maxlik=llik;
     end;
 end;
 subplot(1,1,1);
 bar(0:p,exp(likp-max(likp))); axis([0 p+1 0 1.1])
 xlabel('Model order p'); ylabel('Model Likelihood') 
 
% 
