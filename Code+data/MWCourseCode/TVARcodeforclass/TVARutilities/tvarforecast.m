function [y,mu,phi,v,nsc] = tvarforecast(x,K,I,p,del,mT,CT,sT,nT)
                       
% Forecast ahead in a TVAR(p) model of univariate series x
% Discount factors del(1) for state and del(2) for obs var
% Inputs:
%   x --  series of length T, defining data DT
% % %   K --  number of steps ahead to forecast from current time 
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
% nsc   --  Kx1      per step ahead, percent nonstationary draws 
% 

% Based on D_T, we use a fixed W evolution variance matrix
% for additive state evolution each step ahead, and a fixed beta 
% parameters for the beta innovations in the volatility model
% 
% organise data, initialize and define arrays for outputs ...
T=length(x);  dx=reshape(x,T,1); 
d=del(1); b=del(2); 

phi=zeros(p,K,I); v=zeros(K,I);   mu=v;  nsc=zeros(K,1); 
y=zeros(K+p,I);  y(1:p,:)=repmat(dx(T-p+1:T),1,I); % prepend last p data points     
nsc=zeros(K,1);  % collect MC probs of nonstationary phi vectors at each step ahead 

gT = b*nT/2; hT = (1-b)*nT/2;   % 1-step beta shock parameters 
W  = CT*(1/d-1);                % 1-step evolution variance matrix
L =  chol(W)';                  %    and its transposed Cholesky
A = [eye(p-1) zeros(p-1,1)];
phit = zeros(p,I); vt=zeros(1,I); 
 
% step-ahead forecasting: enforcing local stationarity on phit at each t
% by sampling step-ahead priors and only accepting those that correspond
% to locally stationary AR(p) 
%
rho=0.99;  % max abs(eig) for TVAR state vector samples: constrain to locally stny

% First, make initial draws from the 1-step ahead prior: 
uT = nT/2; zT = nT*sT/2;  aT=mT;  RT = CT+W;  LT=chol(RT)'; imc=0; 
vt  = 1./gamrnd(uT,1/zT,1,I);    % draw inv gamma for volatility 
phit = repmat(aT,1,I) + repmat(sqrt(vt/sT),p,1).*(LT*randn(p,I));    % draw conditional normal prior for state
for i=1:I,              % visit each MC draw to "correct" to locally stationary
    [phit(:,i),ins] = localstnytvar(p,phit(:,i),rho,A);
    nsc(1)=nsc(1)+ins; 
end
 
% Now move over future time steps 1,2,3, ..., k: 
 for j=1:K
     t=p+j;      
     if (j>1)
         if (b<1)
             vt = b*vt./betarnd(gT,hT,1,I);              % sample next volatility
         end
         phit = phit + repmat(sqrt(vt/sT),p,1).*(L*randn(p,I));  % sample next state 
         for i=1:I,              % visit each MC draw to "correct" to locally stationary
             [ phit(:,i), ins ] = localstnytvar(p,phit(:,i),rho,A);
             nsc(j)=nsc(j)+ins; 
         end
     end
     phi(:,j,:) = phit;   v(j,:) = vt;
     % simulate next future: 
     mu(j,:) = sum(y((t-1):-1:(t-p),:).*phit,1); 
     y(t,:) = mu(j,:) + sqrt(vt).*randn(1,I);             %  and add innovations 
 end
 y(1:p,:)=[]; nsc=100*nsc/I; 
 
 
 



