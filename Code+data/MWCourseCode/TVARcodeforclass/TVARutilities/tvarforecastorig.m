function [y,mu,phi,v] = tvarforecast(x,K,I,p,del,mT,CT,sT,nT)
                       
% Forecast ahead in a TVAR(p) model of univariate series x
% Discount factors del(1) for state and del(2) for obs var
% Inputs:
%   x --  series of length T, defining data DT
%   K --  number of steps ahead to forecast from current time T
%   I --  number of Monte Carlo samples to give simulate futures
%   p --  model order
%  mT --  px1 vector posterior mean for state given DT
%  CT --  pxp posterior var matrix
%  nT --  posterior df
%  sT --  posterior estimate of obs var
% Outputs: 
%  phi  --  pxKxI    sampled future state vectors 
%  mu   --  KxI      sampled means of future series ... F'phi 
%  v    --  KxI      sampled future volatilities 
%  y    --  KxI      sampled future series ... "synthetic futures"
% 

% Based on D_T, we use a fixed W evolution variance matrix
% for additive state evolution each step ahead, and a fixed beta 
% parameters for the beta innovations in the volatility model
% 
% organise data, initialize and define arrays for outputs ...
T=length(x);  dx=reshape(x,T,1); 
d=del(1); b=del(2); 

phi=zeros(p,K,I); v=zeros(K,I);    
y=zeros(K+p,I);  y(1:p,:)=repmat(dx(T-p+1:T),1,I); % prepend last p data points     

gT = b*nT/2; hT = (1-b)*nT/2;   % 1-step beta shock parameters 
W  = CT*(1/d-1);                % 1-step evolution variance matrix
L =  chol(W)';                  %    and its transposed Cholesky
A = [eye(p-1) zeros(p-1,1)];
phit = zeros(p,I); vt=zeros(1,I); 
 
% step-ahead forecasting: enforcing local stationarity on phit at each t
% by sampling step-ahead priors and only accepting those that correspond
% to locally stationary AR(p) 
%

% First, make initial draws from the 1-step ahead prior: 
uT = nT/2; zT = nT*sT/2;  aT=mT;  RT = CT+W;  LT=chol(RT)'; imc=0; 
while (imc<I)
    vtest = 1/gamrnd(uT,1/zT);                      % draw inv gamma for volatility  
    phitest = aT + sqrt(vtest/sT)*LT*randn(p,1);    % draw conditional normal prior for state
      %[imc vtest phitest']
      if all( abs(eig([ phitest' ; A ])) < 1 )        % local stationarity? 
        imc=imc+1; 
        phit(:,imc)=phitest; vt(imc)=vtest; 
      end
end
 
  
  
% Now move over future time steps 1,2,3, ..., k: 
 for j=1:K
     t=p+j; 
     for i=1:I
         vold = vt(i); phiold = phit(:,i); 
         imc=0;
         if (j==1)
             phitest = phit(:,i);   vtest = vt(:,i);  
         else
             while (imc==0)
             vtest = b*vold./betarnd(gT,hT,1,1);              % sample next volatility
             phitest = phiold + sqrt(vtest/sT)*L*randn(p,1);  % sample next state 
             %[imc vtest phitest']
             if all( abs(eig([ phitest' ; A ])) < 1 )         % local stationarity? 
                 imc=1;  
             end
             end
         end
         %i
         phi(:,j,i) = phitest;   v(j,i) = vtest;
         % simulate next future: 
         mu(j,i) = y((t-1):-1:(t-p),i)'*phitest; 
         y(t,i) = mu(j,i) +...                        % forecast means  
                     + randn*sqrt(vtest);             %  and add innovations 
     end 
 end;
 y(1:p,:)=[];
 return
 
 



