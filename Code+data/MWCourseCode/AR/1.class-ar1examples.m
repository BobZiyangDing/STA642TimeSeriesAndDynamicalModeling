
startup


help ar1simulate   % and open/examine the function for tips ...

x=ar1simulate(0.9,1,100); 
x=ar1simulate(0.2,1,100); 
x=ar1simulate(-0.9,1,100);

% look at nonstationary cases, e.g., x=ar1simulate(1.1,1,100) 
%   but first need to edit function to directly specific constant innovations 
%    variance v .... try it. 

x=ar1simulate(0.8,1,1000); hist(x,20); histf(x,35); qqplot(x)

% look at two sections -- explore stationarity ...

x1=x(1:250); x2=x(751:1000); plot([x1 x2])

[ mean(x1) std(x1)
  mean(x2) std(x2)]

qqplot(x1)
qqplot(x2)

subplot(2,1,1); plot(x1); subplot(2,1,2); plot(x2)

scatter( sort(x1), sort(x2) ); line([-3 3],[-3 3],'color','k')
 
figure(1); clf; scatter( x1(1:(end-1)), x1(2:end) )
                line([-3 3],[-3 3],'color','k')
                b=regress(x1(2:end),[ones(249,1) x1(1:(end-1))])
                b(2)
                
figure(2); clf; scatter( x2(1:(end-1)), x2(2:end) )
                line([-3 3],[-3 3],'color','k') 
                b=regress(x2(2:end),[ones(249,1) x2(1:(end-1))])
                b(2)
                
              
% reversibility 

plot(x)
n=length(x)

plot( x(n:-1:1) )

plot( flipud(x) )
    

i=randsample(n,1);  plot( x([ 1:i n:-1:(i+1) ]) ); 
                pause
                line([i i],[-4 4],'color','r')


%--------------------------------------                
% Now here's some real data -- it does not conform 'exactly' to AR(1) model 
% assumptions, for good reason, but it is close and gives a start
% on thinking about more realistic models. Also, it gets us into 
% exploratory model estimation/fitting 
% 
y=load('soi.txt'); 
plot(y); xlabel('Time in months'); ylabel('SOI')
         title('Southern Oscillation Index')
         
        % Monthly values of the The Southern Oscillation Index (SOI) 
        % of interest in connection with warming trends and "El Nino" 
        % cycles in climate. SOI is calculated from the monthly or 
        % seasonal fluctuations in the air pressure difference between 
        % Tahiti and Darwin. There is MUCH structure in SOI beyond AR(1)
        % See http://www.cgd.ucar.edu/cas/catalog/climind/soi.html
           
size(y)
n=length(y)
hist(y,25)

x=y-mean(y); 
plot(x)

% stationary characteristics? 

i=randsample(n-150,1); figure(1); plot( x(i:i+150) ); [ mean(x(i:i+150)), std(x(i:i+150)) ]
 
scatter( x(1:(n-1)), x(2:n) )

corrcoef( x(1:(n-1)), x(2:n) )   % SAMPLE AUTOCORRELATION AT LAG 1 

b=regress(x(2:n),[ones(n-1,1) x(1:(n-1))])  % A DIFFERENT WAY 

acf(x,1,10);    % plot sample autocorrelations ACF 
pacf(x,1,10);   % plot sample partial autocorrelations PACF 

e=x(2:n)-0.64*x(1:(n-1)); plot(e)
                          qqplot(e) 
                          acf(e,1,10);  pacf(e,1,10);      
var(e)/(1-0.64^2)
                          
z=ar1simulate(0.64,1,n); 

plot( [z; x] )

plot( [z((n-199):n) ; x(1:200) ]) 

for g=1:100
    	i=10+randsample(n-20,1);  
    	plot( [ z(1:i) ; x(1:(n-i)) ]); axis([0 n+1 -3 3]) 
        text(20,3.2,'Simulated'); text(480,3.2,'Real')
        pause
        line([i i],[-4 4],'color','r')
pause; end



%--------------------------------------    
% Fit reference Bayesian analysis of AR(1) model to data x
% 
  p=1; % this is an AR(1) model - 
  
 [b,B,r,nu,v]=ar(x,p);
 b
 B
 sqrt(v.*diag(inv(B)))
 sqrt(v)
 nu
 
 
 % predict future of the series multiple times via simulation

 npred=100; % number of future time points
 nmc=2500;   % number of sampled "futures"
 
 [xpred,bsamp,vsamp] = arsim(x,p,b,B,nu,v,nmc,npred);
 
 mean(abs(bsamp)>1)      % Monte Carlo estimate of Pr(nonstationary)
  mean(bsamp)
 hist(bsamp(1,:),15); xlabel('\phi'); title('Posterior for \phi'); cmap('bw')


 % aspects of the bivariate posterior ..
 figure(1);
 scatter( bsamp(1,:), sqrt(vsamp(1,:)))
 ylabel('v^{1/2}','rotation',0); xlabel('\phi');  axis([0.5 0.8 0.65 0.9])
 
 figure(2); % now use hist2d and contour the points ... 
             % this uses mw custom function contourpoints ... 
 help contourpoints
  contourpoints( [bsamp(1,:);sqrt(vsamp(1,:))],8,14,[1 2],[],cellstr(char('\phi','v^{1/2}')))
 axis([0.55 0.75 0.7 0.85]); colormap('jet'); ylabel('v^{1/2}','rotation',0)


 
 % plot some synthetic futures 
 for j=1:20
        plot(1:n,x); xlabel('Time'); ylabel('x_t'); hold on; 
        plot((n+1):(n+npred),xpred(:,j),'color','r'); hold off; pause; end
 
 % eg, what's predicted for 5 steps ahead? 
  histf(xpred(5,:),15); xlabel('x_{n+5}'); title('Predictions at t=n+5')
  
 
 % go back a few years to the "highest" SOI and redo .... 
 [a,i]=max(x); savex=x; x=x(1:i); n=i;
 [xpred,bsamp,vsamp] = arsim(x,p,b,B,nu,v,nmc,npred);
 % then look again at some of the above plots .... 
   j=8; boxplot(xpred(1:j,:)'); line([0 j],[0 0])
%---------------------------------------










%%%%%  -------------------------------------------
% HMM AR(1) 
% Look at ACF to see how a signal gets buried or messed up in noise in HMMs 


phi=0.9; s=1; n=200;   % strong signal 
x=ar1simulate(phi,s,n); 
 
w=1; % as an example - SN is 50% so a lot of noise
y=x+sqrt(w)*randn(n,1); 
figure(1); clf
plot(1:n,x,1:n,y); xlabel('Time'); box off; legend('Signal x_y','Signal+Noise y_t'); legend boxoff
figure(2); clf
subplot(2,2,1); acf(x,1,20); title('ACF of signal'); ylim([-0.2 1])
subplot(2,2,3); acf(y,1,20); title('ACF of signal+noise'); ylim([-0.2 1])
subplot(2,2,2); pacf(x,1,20); title('PACF of signal'); ylim([-0.2 1])
subplot(2,2,4); pacf(y,1,20); title('PACF of signal+noise'); ylim([-0.2 1])

w=2; % repeat ...


% repeat with lower \phi -- signal gets buried in noise faster ....




  %-----------------------------------------------
  % -- SV simulation - numbers look like those in real studies - for
  % example and many references, check out 
  %  O. Aguilar and M. West (2000) Bayesian dynamic factor models and 
  %      variance matrix discounting for portfolio allocation.   
  %        Journal of Business and Economic Statistics, 18, 338-357. 
  %
  mu=0.01; phi=0.985; s=0.25; n=1000; 
  x=ar1simulate(phi,s,n);
  y=randn(n,1).*exp(mu+x); 
  
  figure(1)
  subplot(2,1,1); plot(mu+x); title('Log volatility process')
                            xlabel('Time'); ylabel('\mu+x_t')
  subplot(2,1,2); plot(y); title('Returns time series')
                            xlabel('Time'); ylabel('y_t')
                                                  
  figure(2); qqplot(y)
             hist(y,25) 
  scatter( x(1:(n-1)), x(2:n) )
  scatter( y(1:(n-1)), y(2:n) )
  
  
  % Here's some real data - n=2556 daily returns on exchange rates of
  % 12 currencies relative to the US dollar: web link is
  % 
  y=load('returns_exrates.txt'); [n,p]=size(y)
  
  plot(y(:,12))
  
  plot(y(2000:end,12))
  
  qqplot(y(:,12))
  
  
  
  
 
  

