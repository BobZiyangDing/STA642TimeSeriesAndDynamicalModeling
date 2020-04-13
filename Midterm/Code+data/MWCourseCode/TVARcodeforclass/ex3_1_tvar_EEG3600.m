startup 
addpath('../../Data/EEG-data+utilities/')
   % folder with EEG data & scripst for EEG layouts, saved .mat files, movies etc 

% example of use of TVAR functions 
% Length 3600 subsampled-by-6 EEG series from one channel 
% 
% These commented-out lines setup of initial data
% clear; set(0,'DefaultFigureWindowStyle','docked')
%  eeg  % reads in all 19 channels of EEG data etc; 256 Hz recordings
%  x=Y(:,10); Y=x; channames(10,:)  % select Cz channel
%  t0=1001; t=t0:6:(t0+6*3600-1); t1=max(t); T=length(t);x=x(t); 
%  % % selects subsample-by-6 
%  % now timing ..
%  nsecs=3600*6/256; % full time frame in seconds and scale for wavelengths to seconds
%  hscale=256/6;   % subsampled # obsns per second  
%  ti=(1:3600)/hscale;  % 3600 time points in seconds 
%  xticks=round((10:10:80)*hscale); xticklabels='10|20|30|40|50|60|70|80'; 
%  xa='set(gca,''Xtick'',xticks,''XtickLabel'',xticklabels); xlabel(''Time in seconds''); box off; xlim([0 T+1])';
%  save 'eegCz-3600.mat' x T t0 t1 Y nsecs hscale xticks xticklabels xa ti


% ---------------------------------------------------------
% read in data   
load 'eegCz-3600.mat' 

% general plots 
fprintf('\n Reading and plotting EEG data: \n Using figures 1:3- Hit key to continue...\n')
close all; figure(1); clf; figure(2); clf; figure(3); clf
pause 
figure(1); clf
plot(x); box off; eval(xa); ylabel('EEG Cz potential')
%
figure(2); clf 
sx=x-mean(x); sx=sx/std(sx);
s=[200 900 1600 2300 3000]; t=1:500; 
ix=cell(1,11); for i=1:11, ix{i}=int2str(i); end
plot(t,sx(s(1)+t),t,7+sx(s(2)+t), t,14+sx(s(3)+t), t,19+sx(s(4)+t), t,23+sx(s(5)+t)); box off
 text(250,3,[int2str(s(1)) ':' int2str(s(1)+t(end))]); text(250,10,[int2str(s(2)) ':' int2str(s(2)+t(end))])
 text(250,16.5,[int2str(s(3)) ':' int2str(s(3)+t(end))]); text(250,21,[int2str(s(4)) ':' int2str(s(4)+t(end))])
 text(250,24,[int2str(s(5)) ':' int2str(s(5)+t(end))])
 set(gca,'Xtick', hscale*(1:11),'XtickLabel',ix); xlabel('Time in seconds')
 set(gca,'ycolor',get(gcf,'color'))

 
 



% ---------------------------------------------------------
% Fit TVAR(p) to standardized data 
fprintf('\n Fit TVAR(8) to EEG: Hit a key to continue ...\n')
pause
x=sx; 
p=8; m0=zeros(p,1); n0=1; s0=0.01; C0=eye(p);  % initial priors 
del=[0.99,0.95];  %  AR and variance discount factors
% 
[m,C,n,s,e,mf,Cf,sf,nf,ef,qf] = tvar(x,p,del,m0,C0,s0,n0);

% look at smoothed estimates of TVAR coeffs and innovations SD over time 
figure(1); clf
plot(m'); eval(xa); title('TVAR coeffs \phi_t')
                legend( int2str((1:min(4,p))'), 'location','southwest'); legend boxoff
figure(2); clf
subplot(2,1,1); plot(sqrt(s)); eval(xa); title('S.D. v_t')
subplot(2,1,2); plot(abs(x)); eval(xa); title('|EEG|')              
                
%----------------
fprintf('\n Compute estimated decomposition: Hit a key to continue...\n')
pause
% compute "plug-in" estimate of decomposition, wavelengths etc, then plot
[waves,mods,decomp,nr,nc] = tvar_decomp(x,m);
figure(1); clf 
decomp_plot(x,decomp,4); eval(xa) 

figure(2); clf
 plot(waves(1:4,:)'/hscale);    % complex -> real at high freqs -> noise? 
  eval(xa); ylabel('Seconds'); title('Component wavelengths: seconds')
  legend( int2str((1:min(4,p))'), 'location','northwest'); legend boxoff
                                % brain-wave wavelengths in secs
figure(3); clf
plot(hscale./waves(1:4,:)');   % brain-wave frequency scale 
  eval(xa); ylabel('Hz'); title('Component frequencies: Cycles per second')
  legend( int2str((1:min(4,p))'), 'location','northwest'); legend boxoff
  line(repmat([0 T+1],3,1)',[[4 4];[8 8];[13 13]]','linestyle',':','color',[0.7 0.7 0.7])

  
%----------------
% sample posterior N times at selected points
%  and plot trajectories of first k (at most) frequencies ...
%  overlay the "plug-in" estimates on approximate 95% posterior intervals
fprintf('\n Sample posterior for wavelengths etc: Hit a key to continue...\n')
pause
N=2000; k=4; times=75*(1:T/75);
[freqs] = tvar_sim(m,C,n,s,times,k,N); 
 %
 % Frequency trajectory plots plus intervals ...  select a few components 
 k=[1 2];   % components selected ...
 figure(1); clf
 sim_plot_wave(times,freqs*hscale,waves/hscale,k);
     for j=1:length(k), subplot(length(k),1,j);  eval(xa);   ylabel('Seconds'); title(['Wavelength comp ',int2str(k(j))]); end
%    
 figure(2); clf
 sim_plot_freq(times,freqs*hscale,waves/hscale,k);
     for j=1:length(k), subplot(length(k),1,j);  eval(xa);   ylabel('Hz'); title(['Frequency comp ',int2str(k(j))]); end
    

%return
%----------------
% explore likelihood fn for discounts and model order

fprintf('\n Compare multiple models: Hit a key to continue...\n')
pause
pvals=[7 13]; p=pvals(2);  dn=0.985:.005:1; bn=0.94:.005:0.955;  
m0=zeros(p,1); n0=1; s0=0.01; C0=eye(p); 
[popt,delopt,likp] = tvar_lik(x,pvals,dn,bn,m0,C0,s0,n0);
 fprintf('Marginal likelihood max at (delta,beta,p)=(%4.3f,%4.3f,%i)\n',delopt,popt)  % results in:  0.9950    0.9450   12.0000 -- rerun

 
 
 
 