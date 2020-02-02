
% example of use of TVAR functions 
% TV-VAR(p) models for US differenced inflation data 
%  
startup
% --------------- Read in and transform data: 
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
    mmy=[ min(Y')' max(Y')']; plot(Y'); eval(xa); legend(names); legend boxoff
% ---------------



%%%% Selected inflation series: 
iy = 1; x=Y(iy,:);  name=['\Delta ',names(iy,:)]

figure(1); clf
subplot(2,1,1); plot(origY(iy,:)); eval(xa)
    title('Inflation'); ylabel('Annual %') 
subplot(2,1,2); plot(x); eval(xa) 
    title(name); ylabel(['Change %'])
 
figure(2); clf
subplot(2,1,1); acf(x,1,16); xlabel('Quarters') 
    title(['ACF ',name]); ylabel('ACF') 
subplot(2,1,2); pacf(x,1,16); xlabel('Quarters') 
    title(['PACF ',name]); ylabel('PACF') 
    
figure(3); clf
[logmlik,aic,bic] = arpcompare(x',20,1);

% ---------------------------------------------------------
% Fit TVAR(p)
p=4;
fprintf(['\n Fit TVAR(',int2str(p),') to Delta Inflation: Hit a key to continue ...\n'])
pause
m0=zeros(p,1); n0=1; s0=0.05; C0=eye(p); % initial priors 
del=[0.98 .98];  %  AR and variance discount factors
% 
[m,C,n,s,e,mf,Cf,sf,nf,ef,qf] = tvar(x,p,del,m0,C0,s0,n0);

% look at forward filtered summaries 
figure(1); clf
plot(mf'); eval(xa); title('TVAR coeffs \phi_t')
                legend( int2str((1:p)'), 'location','south','orientation','horizontal'); legend boxoff
figure(2); clf
subplot(2,1,1); plot(sqrt(sf)); eval(xa); title('S.D. v_t')
subplot(2,1,2); plot(abs(x)); eval(xa); title(['|',name,'|'])      


% look at smoothed estimates of TVAR coeffs and innovations SD over time 
figure(3); clf
plot(m'); eval(xa); title('TVAR coeffs \phi_t'); 
   ylim([-.8 .8])
   legend( int2str((1:p)'), 'location','south','orientation','horizontal'); legend boxoff
figure(4); clf
subplot(2,1,1); plot(sqrt(s)); eval(xa); title('S.D. v_t'); 
  ylim([.1 .6])
subplot(2,1,2); plot(abs(x)); eval(xa); title(['|',name,'|'])              
                
%----------------
fprintf('\n Compute estimated decomposition: Hit a key to continue...\n')
pause
% compute "plug-in" estimate of decomposition, wavelengths etc, then plot
[waves,mods,decomp,nr,nc] = tvar_decomp(x,m);   
         decomp(:,1:p)=NaN; waves(:,1:p)=NaN; mods(:,1:p)=NaN; % just to remove impact of initial priors
figure(1); clf 
decomp_plot(x,decomp,4); eval(xa) 

ic=min(3,size(waves,1))      % select 3 quasi-prediodic components if they exist
                              % cut-back to max number otherwise

figure(2); clf
 plot(waves(1:ic,:)');    % complex -> real at high freqs -> noise? 
  eval(xa); ylabel('Quarters'); title('Component wavelengths: Quarters')
  legend( int2str((1:ic)'), 'location','east'); legend boxoff
% figure(3); clf
% plot(1./waves(1:ic,:)');   %   frequency scale 
%   eval(xa); ylabel('Cycles/Quarter'); title('Component frequencies: Cycles/quarter')
%   legend( int2str((1:ic)'), 'location','east'); legend boxoff
%  
  

% ---------------- 
% FULL posterior sampling using FFBS and MC inference on components 
 % -- trajectories rather than just marginal posteriors at each past time point

I = 1000;  % MC sample size 
[thetasamp,vsamp] = tvarFFBS(x,p,del,m0,C0,s0,n0,I);
%

% plot FFBS with median, 50% and 90% intervals for TVAR parameters
iphi=[1 p]; % select just 2 TVAR params for example
figure(1); clf;
for i=1:2
    subplot(2,1,i); 
    pr = prctile(squeeze(thetasamp(iphi(i),:,:))',[5 25 50 75 95])';
    ciplot(pr(p+1:T,1),pr(p+1:T,5),p+1:T,[0.95 0.95 0.95]); hold on
    ciplot(pr(p+1:T,2),pr(p+1:T,4),p+1:T,[0.85 0.85 0.85]);
    line(p+1:T,pr(p+1:T,3),'color','k')
    axis tight; box off; eval(xa)
    hold on; plot([1 T],[0 0],'b:');hold off
    ylabel('Qtrs'); title(['TVAR parameter \phi_{',int2str(iphi(i)),'}']) 
    if (i==1), xlabel(''); end 
end


% Now compute eigendecomp and extra posterior samples on wavelengths & moduli
nc=2;   % select only the 2 largest wavelength components here 
waves=zeros(nc,T,I); mods=waves; 
for i=1:I
    [wa,mo] = tvar_decomp(x,thetasamp(:,:,i));  
    waves(:,:,i)=wa(1:nc,:); mods(:,:,i)=mo(1:nc,:); 
    i
end
 

% plot FFBS with median, 50% and 90% intervals for wavelengths
figure(2); clf;
yl =[ 24 10 ]; % upper limits in qtrs for plots y-axis 
for i=1:2
    subplot(2,1,i); 
    pr = prctile(squeeze(waves(i,:,:))',[5 25 50 75 95])';
    ciplot(pr(p+1:T,1),pr(p+1:T,5),p+1:T,[0.95 0.95 0.95]); hold on
    ciplot(pr(p+1:T,2),pr(p+1:T,4),p+1:T,[0.85 0.85 0.85]);
    line(p+1:T,pr(p+1:T,3),'color','k')
    axis tight; ylim([0 yl(i)]); box off; eval(xa)
    ylabel('Qtrs'); title(['Wavelength ',int2str(i)]) 
    if (i==1), xlabel(''); end 
end
% .. and for moduli 
figure(3); clf 
for i=1:2
    subplot(2,1,i); 
    pr = prctile(squeeze(mods(i,:,:))',[5 25 50 75 95])';
    ciplot(pr(p+1:T,1),pr(p+1:T,5),p+1:T,[0.95 0.95 0.95]); hold on
    ciplot(pr(p+1:T,2),pr(p+1:T,4),p+1:T,[0.85 0.85 0.85]);
    line(p+1:T,pr(p+1:T,3),'color','k')
    axis tight; ylim([0 1.1]); box off; eval(xa)
    hold on; plot([1 T],[1 1],'b:');hold off
    ylabel('Qtrs'); title(['Modulus ',int2str(i)]) 
    if (i==1), xlabel(''); end 
end
 
 
 
%----------------
% Prediction: forecast the next h quarters .... 
h = 12; % steps ahead
I = 2500; % MC samples
mT=m(:,T); CT=squeeze(C(:,:,T)); sT=s(T); nT=n(T); 
[ y,mu phi,v,probsns] = tvarforecast(x,h,I,p,del,mT,CT,sT,nT);
probsns'       

figure(1); clf
nx=75;    
pr = prctile(y',[10 25 50 75 90])'; 
ciplot(pr(:,1),pr(:,5),T+1:T+h,[0.95 0.95 0.95]); hold on
ciplot(pr(:,2),pr(:,4),T+1:T+h,[0.85 0.85 0.85]);
scatter(T+1:T+h,pr(:,3),10,'ro')
scatter(1:T,x,10,'+'); hold on; ylim([min(x) max(x)])
%%errorbar(T+1:T+h,mean(y,2),1.5*std(y,0,2),'m*')  
hold off; eval(xa); xlim([T-nx T+h+1])
title(['Predictions of ',name]);  ylabel(['Change %'])

% convert to the original case rather than changes scale .... 
figure(2); clf
Y = origY(iy,end)+cumsum(y);
pr = prctile(Y',[10 25 50 75 90])'; 
ciplot(pr(:,1),pr(:,5),T+1:T+h,[0.95 0.95 0.95]); hold on
ciplot(pr(:,2),pr(:,4),T+1:T+h,[0.85 0.85 0.85]);
scatter(T+1:T+h,pr(:,3),10,'ro')
scatter(1:T,origY(iy,:),10,'+'); hold on; % ylim([min() max(x)])
%%errorbar(T+1:T+h,mean(y,2),1.5*std(y,0,2),'m*')  
hold off; eval(xa); xlim([T-nx T+h+1]); line([T-nx T+h+1],[0 0])
title(['Predictions of ',names(iy,:)]);  ylabel(['Rate'])


% 
% % look at a few synthetic predicted realities 
figure(3); clf
r=8; ir=randsample(I,r);
scatter(1:T,origY(iy,:),10,'+'); hold on
plot(T+1:T+h,Y(:,ir),'+-');     
hold off; eval(xa); xlim([T-nx T+h+1])
title(['Synthetic futures of ',names(iy,:)]);  ylabel(['Rate'])

% 


% 


fprintf('\n Compare multiple models: Hit a key to continue...\n')
pause
pvals=[2 12]; p=pvals(2);  dn=0.985:.005:0.995; bn=0.8:.005:0.995;  
m0=zeros(p,1); n0=1; s0=0.05; C0=eye(p); 
[popt,delopt,likp] = tvar_lik(x,pvals,dn,bn,m0,C0,s0,n0);
 fprintf('Marginal likelihood max at (delta,beta,p)=(%4.3f,%4.3f,%i)\n',delopt,popt)  % results in:  0.9950    0.9450   12.0000 -- rerun

 
 
