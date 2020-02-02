function [f,V] = forecast(k,nmc,h,D,M,C,beta,delta,q,arp,p,Yt)
% FORECAST k steps ahead in TV-VAR(arp) + intercept model in q dimensions
% p = dimension of model paramater = 1+arp*q
% nmc = size of Monte Carlo sample from predictive distribution
% Y = data on current and last arp-1 times: Y(:,t:-1:t-arp+1) where t="now"
% Current posterior summaries: 
%    h = degrees of freedom of Wishart, > q-1
%    D = volatility sum of squares matrix
%    M = current estimate of parameters
%    C = current variance info on parameters
% Model info:  
%   beta = volatility discount
%   delta = parameter discount
%
% Outputs: f = mean vector, V = variance matrix of predictions

Yf=zeros(q,k,nmc);  % sampled futures 
n=h-q+1; 
Sigi = InvWishartrnd(n,D,1); M = rMatNorm(M,C,Sigi); 
Wt=C*(1/delta-1);  Ft = repmat([1; reshape(Yt,p-1,1)],1,nmc); iq=2:q+1; 

for t = 1:k
      if (n>q/beta), 
          h    = beta*h; D = beta*D; n=h-q+1;  % evolve Wishart parameters 
      end
      Sigs = InvWishartrnd(n,D,nmc);        % sample the set of MC volatility matrices at T+t
                       % state evolution and 1-step prior variance matrix
      for i=1:nmc
          Sigi      = squeeze(Sigs(i,:,:)); 
          M         = rMatNorm(M,Wt,Sigi);                   % sample the state matrix 
          Yf(:,t,i) = rMNorm(M'*Ft(:,i),Sigi,1);             % and the outcome 
          Ft(2:p,i) = [ Yf(:,t,i); Ft(2:p-q,i) ];            % predictors for next time point
      end 
      % display([t n])
end

Yf= squeeze(Yf(:,k,:));  % MC sample of returns at time k
f =  mean(Yf')'; V = cov(Yf');

