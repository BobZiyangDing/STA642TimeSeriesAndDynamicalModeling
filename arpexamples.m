

% Explore basics of some AR model structures and analyses 
 
% simple AR(2) simulation - trick the arsim function ...
p=2; 
% An example with two real roots
phi=[ .99+.75, -.99*.75 ]'; 
% or
phi=[ .99+(-.75), -.99*(-.75) ]'; 

B=eye(2)*10^12; nu=10^12; v=1;
npred=500; nmc=1; x=[0 0]'; 


[xpred,bsamp,vsamp] = arsim(x,p,phi,B,nu,v,nmc,npred);
x=xpred(51:500); 

clf
subplot(3,1,1); plot(x); 
subplot(3,1,2); acf(x,1,30);
subplot(3,1,3); pacf(x,1,30); 

% An example with complex roots - quasi-cyclical processes
r=0.9; om=2*pi/12;  phi=[2*r*cos(om) -r*r]';

% run above simulation and plots ..

r=0.99999;  % repeat

r=0.1;      % repeat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% explore SOI data ...

y=load('soi.txt'); 
plot(y); xlabel('Time in months'); ylabel('SOI')
         title('Southern Oscillation Index')
     
size(y)
n=length(y)
x=y-mean(y); 

acf(x,1,36); 
acf(x,1,60);
pacf(x,1,36); 
%--------------------------------------    
% Fit reference Bayesian analysis of AR(p) model to data x

p=24;  % why 24 months?   Recall AR(2p+) "folklore" for AR(p) signal in noise
  
[b,B,r,nu,v]=ar(x,p);

sd=sqrt(v*diag(inv(B)));
a=[ b-2*sd  b+2*sd ]; [a(:,1) b a(:,2) ]

figure(1); clf
plot(1:p,b,'d'); box off; for j=1:p, line([j j],a(j,:)); end
    xlabel('AR coeff index j'); ylabel('95% interval for \phi_j')
    title(['Inference on \phi in AR(',int2str(p),') model of SOI'])
    line([0 p+1],[0 0])
    

 % sample posterior and sample predict futures 
 
 npred=n;    % number of future time points
 nmc=2500;   % number of sampled "futures"
 
 [xpred,bsamp,vsamp] = arsim(x,p,b,B,nu,v,nmc,npred);
 
 figure(2); clf; 
 histf(bsamp(1,:),25); xlabel('\phi_1'); title('Posterior for \phi_1')
  
 figure(2); clf
 boxplot(bsamp'); ylabel('AR(1) coefficients \phi_j'); 
                  xlabel(['Coefficient lag j=1:',int2str(p)]); box off
 line([0 p+1],[0 0],'color','k','linestyle',':')


 % plot some of the sample futures with some data -- random number of data
 % points and synthetic future points .. can you see where data ends and 
 % simulations beging?  Acid test of "model fit"
 
 
for g=1:100
    	i=30+randsample(n-60,1);  j=randsample(nmc,1); 
    	plot( [ x(1:i) ; xpred(1:(n-i),j) ], 'r' ); 
        axis([0 n+1 min(x) max(x)]); box off;  xlabel('Months'); ylabel('x_t'); 
        text(n-100,4.2,'Synthetic'); text(60,4.2,'Real')
        pause
        line([i i],[-4 4],'color','r')
        hold on; plot( x(1:i), 'b'); hold off
pause; end


 
npred=200;  xpred=xpred(1:npred,:); 
 % plot some of the sampled futures .. clearly now
 for j=1:50, plot(1:n,x); axis([0 n+npred+1 min(x) max(x)]); box off;  xlabel('Months'); ylabel('x_t'); hold on; 
              plot((n+1):(n+npred),xpred(:,j),'color','r'); hold off; 
 pause(0.5); end
hold on; plot((n+1):(n+npred),mean(xpred'),'k.'); hold off

 
 % find eigenvalues - characteristic roots .. first just use the posterior mean b
 phi=b 
 
 % directly via eigenstructure: 
  G=[ phi' ; [eye(p-1) zeros(p-1,1)] ]; lambda=eig(G)
  
  [ abs(lambda) 2*pi./angle(lambda) ]

  ang=angle(lambda); j=find(ang>0&ang<pi); [ abs(lambda(j)) 2*pi./angle(lambda(j)) ]
  
  % added following 2019(b) Matlab update .. new convention for eig ...
  [a,i]=sort(abs(lambda),'descend');  % make sure order down in modulii ...
  [ abs(lambda(i)) 2*pi./angle(lambda(i)) ]
  
  % note: suggested ~4year cycle +/- -- El Nino
  % Plus annual harmonics: subcycles of ~2yrs, 1yr, 1/2 yr,  ? crude but a start ...

  
 % stationary? 
  abssamp=[]; 
  for i=1:nmc 
      abssamp = [abssamp sort(abs(1./roots([-flipud(bsamp(:,i));1])),'descend')];
  end
  boxplot(abssamp'); ylabel('AR(24) modulii, ordeeed down');
 
  
  % MC estimate of Prob(non-stationary): 
  sum(sum(abssamp'>1))/nmc
  
 % explore decomposition of data at fitted model parameters:   components
 % ordered in terms of decreasing wavelength then modulii for real eigs
  [waves,mods,decomp]=ar_decomp(x,p,p/2);
  
  [waves mods]
  
  clf;  decomp_plot(x,decomp,4);
        decomp_plot(x,decomp,2);

 % what about uncertainty about wavelengths, components ...
 % simulate posterior for decomposition of AR under reference analysis: 
 nmc==5000;
 [swaves,smods,sdecomp]=ar_decomp_sim(x,p,nmc);
 figure(1); clf
 subplot(2,1,1); boxplot(smods'); xlabel('moduli'); box off
 subplot(2,1,2); boxplot(swaves(1:4,:)'/12); ylim([0 10]); xlabel('wavelengths (years)'); box off
 prctile(swaves(1:4,:)'/12,[5 25 50 75 95 99])      

 figure(3); clf
 decomp_plot_sim(x,sdecomp,3)

      
  % repeat:  different values of p .... 
  % hard to find longer-term structure using shorter order models
  % often need to use longer AR orders, overfitting to find structure
  % which means additional low-frequency/low modulus noise terms too
  
   
% AIC, BIC and reference marginal likelihood for model order
maxp=25;  
[logmlik,aic,bic] = arpcompare(x,maxp,0);   % the input 1 means ignore case p=0
        
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % EEG series 
  x=load('eeg400.txt'); % or 800, or more ... 
  x=x-mean(x); n=length(x); plot(x)
  acf(x,1,40);
  [waves,mods,decomp]=ar_decomp(x,8,8);
  decomp_plot(x,decomp,3);
  [waves mods]
  
  [waves,mods,decomp]=ar_decomp(x,4,2);
  decomp_plot(x,decomp,2);
  
    
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % US Macro econ data  
    fprintf('q=3 series, quarterly data, 1960/Q1 - 2014/Q4\n')
    fprintf('  - Inflation = GDP Deflator, annual percent change = (Nominal GDP)/(Real GDP)*100\n')
    fprintf('  - Interest Rate: Short Term = 3 month Treasury bill\n')
    fprintf('  - Unemployment Rate\n')
    %
    [Y,Ynames]=xlsread('USmacrodata1960-2014.xlsx');
    [T q]=size(Y); q=q-2; zq=zeros(q,1); oq=ones(q,1); Iq=eye(q);
    yrqtr=Y(:,1:2); Y(:,1:2)=[]; Y=Y'; Ynames(1:2)=[]; names=char(Ynames);
    %
    %  difference all series:
    origY=Y; Y=[zeros(q,1) diff(Y')'];
    %
    nticks=14; h=round(T/(4*nticks)); tticks=1:4*h:T; % plot nticks axis ticks & labels
    tdates=cell(size(tticks)); tdi=int2str(yrqtr(tticks,1));
    for i=1:length(tticks), tdates{i}=tdi(i,3:4); end
    xa=['set(gca,''Xtick'',tticks);set(gca,''XtickLabel'',tdates);xlabel(''Year/Q1''); xlim([0 T+1]); box off;'];
    % plot ...
    mmy=[ min(Y')' max(Y')']; 
    clf; plot(Y'); eval(xa); legend(names); legend boxoff
% ---------------

%%%% Select inflation series: 
iy = 1; x=Y(iy,:)';  name=['\Delta ',names(iy,:)]

figure(1); clf
subplot(2,1,1); plot(origY(iy,:)); eval(xa)
    title('Inflation'); ylabel('Annual %') 
subplot(2,1,2); plot(x); eval(xa) 
    title(name); ylabel(['Change %'])
 
figure(2); clf
subplot(2,1,1); acf(x,1,48); xlabel('Quarters') 
    title(['ACF ',name]); ylabel('ACF') 
subplot(2,1,2); pacf(x,1,48); xlabel('Quarters') 
    title(['PACF ',name]); ylabel('PACF') 

% AR models for "business cycles" ? 

 p=4;  % quarterly data
 [waves,mods,decomp]=ar_decomp(x,p,4);
  decomp_plot(x,decomp,4);
  [waves mods]      % waves in quarters
  [waves/4 mods]    % waves in years
  
p=8;  % how to think about choice of p? 
       % business cycles - folklore 3.5-6 years
       
 % uncertainty about wavelengths, components ...
 % simulate posterior for decomposition of AR under reference analysis: 
  nmc=10000;
 [swaves,smods,sdecomp]=ar_decomp_sim(x,p,nmc);
 figure(1); clf;
        subplot(2,1,1); boxplot(smods'); xlabel('moduli'); box off
        subplot(2,1,2); boxplot(swaves'/4); xlabel('wavelengths (years)'); box off
 figure(2); clf;
        decomp_plot_sim(x,sdecomp,3); eval(xa)
        
        
              
  
  
  