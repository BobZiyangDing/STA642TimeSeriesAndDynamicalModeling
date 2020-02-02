startup 

% example of use of TVAR functions 
% Long/more recent SOI series with lots of action 

% ---------------------------------------------------------
% read in data   
[x,months]=xlsread('soi1876-2012.xlsx'); months=char(months(1,:)); years=x(:,1); x(:,1)=[]; x=reshape(x',prod(size(x)),1); 
i=find(isnan(x)); x(i)=[]; T=length(x)
x=(x-mean(x))/std(x); plot(x)

tticks=1880:10:2012;
tdates= reshape([ int2str(tticks') repmat('|',length(tticks),1) ]',1,5*length(tticks));
 r=T/length(tticks); tticks=(tticks(1)-years(1))*12:r:T; 
xa='xlim([0 T+1]); set(gca,''Xtick'',tticks);set(gca,''XtickLabel'',tdates);box off';
eval(xa)


figure(1); clf
plot(x); box off; eval(xa); ylabel('SOI')
%
acf(x,1,75)
pacf(x,1,75)
sx=x-mean(x); sx=sx/std(sx);


% ---------------------------------------------------------
% Fit TVAR(p) to standardized data 
x=sx; 
p=12; m0=zeros(p,1); n0=5; s0=0.25; C0=10*eye(p);  % initial priors 
del=[0.99,0.99];  %  AR and variance discount factors
% 
[m,C,s,n,e]= tvar(x,p,del,m0,C0,s0,n0);

% look at smoothed estimates of TVAR coeffs and innovations SD over time 
figure(1); clf
plot(m'); eval(xa); title('TVAR coeffs \phi_t')
                legend( int2str((1:min(4,p))'), 'location','southwest'); legend boxoff
figure(2); clf
subplot(2,1,1); plot(sqrt(s)); eval(xa); title('S.D. v_t')
subplot(2,1,2); plot(abs(x)); eval(xa); title('|SOI|')              
                
%----------------
fprintf('\n Compute estimated decomposition: Hit a key to continue...\n')
pause
% compute "plug-in" estimate of decomposition, wavelengths etc, then plot
[waves,mods,decomp,nr,nc] = tvar_decomp(x,m);
figure(1); clf 
decomp_plot(x,decomp,4); eval(xa) 

figure(2); clf
 plot(waves(1:4,:)');    % complex -> real at high freqs -> noise? 
  eval(xa); ylabel('Months'); title('Component wavelengths: Months')
  legend( int2str((1:min(4,p))'), 'location','northwest'); legend boxoff

  
%----------------
% sample posterior N times at selected points
%  and plot trajectories of first k (at most) frequencies ...
%  overlay the "plug-in" estimates on approximate 95% posterior intervals
fprintf('\n Sample posterior for wavelengths etc: Hit a key to continue...\n')
pause
N=500; k=4; times=75*(1:T/75);
[freqs] = tvar_sim(m,C,n,times,k,N); 
 %
 % Frequency trajectory plots plus intervals ...  select a few components 
 k=[1 2];   % components selected ...
 figure(1); clf
 sim_plot_wave(times,freqs,waves,k);
     for j=1:length(k), subplot(length(k),1,j);  eval(xa);   ylabel('Months'); title(['Wavelength comp ',int2str(k(j))]); end
%    


%return
%----------------
% explore likelihood fn for discounts and model order

fprintf('\n Compare multiple models: Hit a key to continue...\n')
pause
pvals=[8 30]; p=pvals(2);  dn=0.98:.01:1; bn=0.98:.01:1;  
m0=zeros(p,1); n0=5; s0=0.25; C0=10*eye(p); 
[popt,delopt,likp] = tvar_lik(x,pvals,dn,bn,m0,C0,s0,n0);
 fprintf('Marginal likelihood max at (delta,beta,p)=(%4.3f,%4.3f,%i)\n',delopt,popt)  % results in:  0.9950    0.9450   12.0000 -- rerun

 
 
 
 