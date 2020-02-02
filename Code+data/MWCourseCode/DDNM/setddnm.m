
% specify parallel DLMs with steady state and volatility evolutions
%  and selected parental and lagged predictors in each 
 
% modelname='N&W14 Inspired Model';
pa=cell(1,q);     
pa{1}=2:3; % Parents of interest rate: Inflation, Unemployment   
           % Inflation & Unemployment have no parents
npa=[2 0 0]'; 

% fit model starting at time t0+1
% t0 is the number of lags
maxlag = 3; t0 = maxlag;

% indicator matrix of lagged predictors
% (series, predictor, lags)
pr=zeros(q,q,maxlag); 

% IntRate <-  IntRate(lag 1) + Infln(lag 3) + Unempl(lag 3) 
pr(1,1,1)=1; 
pr(1,2,3)=1;  
pr(1,3,3)=1;  

% Infln <- Infln(lag 1,2) + IntRate(lag 1,2,3)
pr(2,1,[1 2 3])=1;   
pr(2,2,[1 2])=1;               

% Unempl <- Unempl(lag 1, 2) + Infln(lag 1) 
pr(3,2,1)=1; 
pr(3,3,[1 2])=1; 

npr=zeros(q,1); 

% dimensions
for j=1:q, 
    npr(j)=sum(sum(squeeze(pr(j,:,:)))); 
end 

% full state dimensions are p+1 since we add intercepts
p =  npa + npr; 

% max lags for priors
maxlags = [3 3 2];



