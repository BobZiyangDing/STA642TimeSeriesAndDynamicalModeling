function  [logmlik,aic,bic] = arpcompare(x,maxp,z)
% Compute conditional marginal likelihood for AR model orders 0,...,maxp 
% ***under reference prior*** 
% Inputs:
%  x     --  series of length T  (assumed 0 mean 
%  maxp  --  max model order
%  z     --  if 0, include case p=0; if non-zero, start at p=1
% Outputs: 
% likp  --  conditional marginal likihood over 0:maxp
% popt   -- MLE of model order
% Prado & West, Section 2.3.4
%  nb. compares analyses based on using only the last T-maxp observations
% 
 T=length(x); n=T-maxp; 
 logmlik=zeros(maxp,1); aic=logmlik; bic=aic; 
 % fit reference analyses for AR(p)a - reuse code from ar.m function
 y=x(T:-1:maxp+1);                        % response data using only maxp+1:T time series
 X=hankel(x(T-1:-1:maxp),x(maxp:-1:1));   % AR predictors based on T-maxp data pts
 % 
 % now fit each model order p=1,2....
 %  
 for p=1:maxp  
   H=X(:,1:p);          % select just the last p values as predictors
   B=H'*H; b=B\(H'*y);  % reference posterior mean 
   r=y-H*b;             % residuals
   R=r'*r;              % residual sum of squares
   nu=n-p;              % residual dof
   logmlik(p) = ( -nu*log(pi) -log(det(B)) -nu*log(R) )/2 + gammaln(nu/2); 
                    %  log marginal lik under reference prior - to a const
                    %  includes all terms in H, p .. 
   ic=n*log(R/nu);  % error in P&W: Fig 2.6. That incorrectly used ic=n*log(R/n);
                    % try it and you will reproduce that figure exactly
   aic(p)= -(2*p+ic)/2;      % -0.5 AIC
   bic(p)= -(log(n)*p+ic)/2; % -0.5 BIC
                             % both times -1/2  to convert to log lik scale
 end
 %
 plim=maxp; plow=0; prange=1:maxp; 
 if (z==0) 
 % add case of p=0 - data are N(0,v) with ref prior on v
    R=y'*y; 
    logmlik0 = -n*log(pi*R)/2 + gammaln(n/2);
    ic0  = n*log(R/n); 
    logmlik = [ logmlik0 ; logmlik ]; 
    aic = [ -ic0/2 ; aic ];   bic = [ -ic0/2 ; bic ]; 
    plim=maxp+1; plow=-0.5; prange=0:maxp; 
 end
 %
 clf; set(gca,'linewidth',1.15,'fontsize',14)
 ml = logmlik-max(logmlik);    % to normalize to 0 at maximum
 ma = aic-max(aic);            % ditto
 mb = bic-max(bic);            % ditto  
 c=[ml;ma;mb];  d=range(c)/25; 
 axis([plow plim+0.5 min(c)-d max(c)+d ]); 
 xlabel('AR model order p'); ylabel('Log likelihood')
 text(prange, logmlik-max(logmlik), repmat('m',plim,1), 'color','b','fontsize',14)
 text(prange, aic-max(aic), repmat('a',plim,1), 'color','r','fontsize',14)
 text(prange, bic-max(bic), repmat('b',plim,1), 'color','k','fontsize',14)
 text(maxp*1/5,min(c)+range(c)*3/20,   'm = log reference marginal likelihood','color','b','fontsize',14)
 text(maxp*1/5,min(c)+range(c)*2/20,'a = -AIC/2','color','r','fontsize',14)
 text(maxp*1/5,min(c)+range(c)*1/20, 'b = -BIC/2','color','k','fontsize',14)

% 
