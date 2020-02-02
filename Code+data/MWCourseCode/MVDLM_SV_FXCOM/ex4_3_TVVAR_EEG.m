
startup

% TV-VAR(p) models for EEG data at various subsampled time scales
% 

load eeg_alldata
  % Y=Xeeg(4001:8:16000,3:21)';  % subsample every 8 - mlik likes arp=1, delta=beta=0.99
  Y=Xeeg(2001:4:18000,3:21)'; % for finer time scale mlik likes arp=2, delta=beta=0.99 
  
[q T]=size(Y); names=channames; Ynames=cellstr(names);
Y=subtract_mean_rows(Y); saveY=Y; 


%r=[4 7 8 10 19]; % select fewer channels out of the 19 for pilot examples
%Y=Y(r,:); Ynames=Ynames(r); names=names(r,:); q=length(r); saveY=Y; 

arp=2; p=arp*q; % ---- set up TV-VAR(arp) with zero mean
F=zeros(q*arp,T-arp); Xnames=cell(1,p);  r=0; 
for j=1:arp, 
    F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
    for i=1:q, r=r+1; Xnames{r}=[ char(Ynames(i)),':Lag-',int2str(j) ]; end
end
Y(:,1:arp)=[]; T=T-arp; 
 

% choose a pair of discounts and run with it: 
delta=0.99; beta=0.99; 
n0=1; h0=n0+q-1;  D0=h0*eye(q)*1000;     % priors
z = zeros(p,q);  zq=zeros(q,1); 
M0=z; r=0.99; 
        if (p>1)
              M0(1:q,:)=r*eye(q);   % sets prior mean to be zero everywhere but r on lag-1 of same series
        end
Mt = M0; C0=eye(p);  Ct=C0/10;         % initial Theta prior 
n = n0; h=h0; D = D0;  St=D/h;         % initial Sigma prior
sMt=zeros(p,q,T);  sCt=zeros(p,p,T);  sdCt=zeros(p,q,T); sSt=zeros(q,q,T); sEt=zeros(q,T); 
sft=zeros(q,T);    sQt=zeros(q,q,T); snt=zeros(1,T); sloglik=zeros(1,T); 

% forward filtering: 
        for t = 1:T
            ft = Mt'*F(:,t);            sft(:,t)=ft; 
            et = Y(:,t) - ft;
            Rt = Ct/delta; 
            h  = beta*h;  n=h-q+1;  D = beta*D;  snt(t)=n;     
            qvt = 1 + F(:,t)'*Rt*F(:,t); sEt(:,t) = et./sqrt(qvt*diag(St)); 
                                         sQt(:,:,t) = St*qvt; 
                                         sloglik(t) = ltpdf(et,zq,qvt,n,D); 
            At = Rt*F(:,t)/qvt;
            h=h+1; n=n+1; D = D+et*et'/qvt;  St=D/h; St=(St+St')/2;  
            Mt = Mt + At*et'; Ct = Rt - At*At'*qvt;   sCt(:,:,t)=Ct;
            sSt(:,:,t)=St; sMt(:,:,t) = Mt; sdCt(:,:,t) = sqrt(diag(Ct)*diag(St)'); 
        end

% reverse smoothing  -
K=inv(sSt(:,:,T)); n=snt(T); Mt = sMt(:,:,T); Ct = sCt(:,:,T); sEt(:,T) = Y(:,T) - Mt'*F(:,T);
for t=(T-1):-1:1
    n = (1-beta)*snt(t)+beta*snt(t+1);          snt(t)=n; 
    K=(1-beta)*inv(sSt(:,:,t))+beta*K;          St = inv(K); sSt(:,:,t)=St;  
    Mt = (1-delta)*sMt(:,:,t) +delta*Mt;        sMt(:,:,t) = Mt; 
    Ct = (1-delta)*sCt(:,:,t) + delta^2*Ct;     sCt(:,:,t) = Ct; sdCt(:,:,t) = sqrt(diag(Ct)*diag(St)');
    sEt(:,t) = Y(:,t) - Mt'*F(:,t);
end
  
%clf
%plot(sloglik);  title('log marginal likelihood')
       

%---------- PLOTS --------------
% some TV-VAR coeffs over time 
for j=1:q
    display(['Coeffs of ',char(Ynames(j)),'...']); clf
    for i=1:3
        subplot(3,1,i); 
        mt=squeeze(sMt(i,j,:)); sct=squeeze(sdCt(i,j,:)); 
        plot(1:T,mt,'r',1:T,mt+sct,'b:',1:T,mt-sct,'b:')
        line([0 T],[0 0],'color','k','linestyle','--');  
        title([char(Ynames(j)),' on ',char(Xnames(i)),':   \Theta_{t,',int2str(i),',',int2str(j),'}'])
        ylabel('Coeff'); box off
        pause
    end
    display('Hit a key to see trajectories for next series ... ')
end
%
%  Some volatilities over time 
%
figure(1); clf; c=char('a','b','c','d'); 
JJ=[4 10 18];
for r=1:length(JJ)
    J=JJ(r); subplot(3,1,r)
        plot(sqrt(squeeze(sSt(J,J,:))),'k')
       title(['(',c(r),') ','SD(',deblank(names(J,:)),')'])
end
%
%  Some correlations over time 
%
figure(2); clf; c=char('a','b','c','d'); 
JJ=[[1 2];[1 3];[2 3]];
for r=1:size(JJ,1)
    J=JJ(r,:); subplot(3,1,r)
        plot(squeeze(sSt(J(1),J(2),:))./sqrt(squeeze(sSt(J(1),J(1),:)).*squeeze(sSt(J(2),J(2),:))),'k')
        title(['(',c(r),') ','Corr(',deblank(names(J(1),:)),',',deblank(names(J(2),:)),')'])
end
%
% image plots of estimated TV-VAR coefficient matrices
j=10; clf
for la=1:arp, 
    imagesc(squeeze(sMt(1+(la-1)*q:la*q,j,:))); cmap('gr'); box off
    title([names(j,:),' : lag ',int2str(la),' coefficients'])
    set(gca,'Ytick',1:q);set(gca,'YtickLabel',names)
    xlabel('Time'); colorbar; pause
end 
%
%
% better - in-context movie of EEG estimated TV-VAR coefficiants on head shots 
% Show time images of countours of coefficients at lags 1-arp of all EEG
% channels on one chosen channel -- to look at the feed-forward connectivity
% of channel i on the rest
load('eeg_electrode_layout')
i=10;       % choose channel i
figure(1); clf; cmap(jet)
np=ceil(arp/2); t=0:pi/20:2*pi; a=0.08; b=sqrt(1-a*a); 
for r = 250:10:T-250 
   M=squeeze(sMt(:,:,r));    
  for j=1:arp
    subplot(np,2,j);    
    plot(sin(t),cos(t),'k'); hold on; plot([-a 0 0 a],[b 1+2*a 1+2*a b],'k'); axis off; hold on; axis square
    [x,y]=meshgrid(linspace(-1,1,100),linspace(-1,1,100));
    z=griddata(eeg_electrode_layout(:,1),eeg_electrode_layout(:,2),...
           M((j-1)*q+(1:q),i),x,y,'v4');
    z((x.*x+y.*y)>1)=NaN;  
    contour(x,y,z,12); axis square; axis tight
    text(chanlayout(:,1),chanlayout(:,2),channames); hold off
    title([char(Ynames(i)),' : lag ',int2str(j)])
    %sub_pos = get(gca,'position'); % get subplot axis position
    %set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height 
    colorbar
  end
    display(int2str(r));   pause(0.0001)
end




% ------------ model comparison ------------------------------------

% for saved results:
 load ex4_3_modelLikelihoods.mat  
 [a,i]=max(plik(end,:))
 [ char('delta  ','beta  ','p  ') char( num2str(pars(1:2,i),2), int2str(pars(3,i)) ) ]
 return
 
 
 
 
 
 
 
 
%%%%%%               
% the above summaries were saved having been computed by the following: 

% comparing different discount factors and TV-VAR model order
 Y=saveY; T=size(Y,2); % reset full data set
 %
 Delta=[0.94:0.01:1]; Beta=Delta; ndelta=length(Delta); nbeta=length(Beta); 
 Arp=0:6; narp=length(Arp);
 nm=ndelta*nbeta*narp; mlik=zeros(T,nm); pars=zeros(3,nm); im=0;
 for id=1:ndelta
     delta=Delta(id); 
     for ib=1:nbeta
         beta=Beta(ib);
         for iar=1:narp
             Y=saveY; T=size(Y,2);
             arp=Arp(iar); im=im+1; pars(:,im)=[delta beta arp]; 
             p=arp*q; % ---- set up TV-VAR(arp)with zero mean
             F=zeros(q*arp,T-arp); r=1; 
             for j=1:arp, 
                F((j-1)*q+(1:q),:)=Y(:,(arp-j+1):(T-j)); 
             end
             Y(:,1:arp)=[]; T=size(Y,2); % 
             n0=1; h0=n0+q-1;  D0=h0*eye(q)*1000;     % priors
             z = zeros(p,q);  zq=zeros(q,1); M0=z; r=0.99;
             if (p>1)
                 M0(1:q,:)=r*eye(q);   % sets prior mean to be zero everywhere but r on lag-1 of same series
             end
             Mt = M0; C0=eye(p);  Ct=C0/10;         % initial Theta prior
             n = n0; h=h0; D = D0;  St=D/h;         % initial Sigma prior
             sloglik=zeros(1,T); 
             % forward filter 
             for t = 1:T
                ft = Mt'*F(:,t);  et = Y(:,t) - ft; Rt = Ct/delta; 
                n  = beta*n;  b=(n+q-1)/h; D = b*D;  
                qvt = 1 + F(:,t)'*Rt*F(:,t);  sloglik(t) = ltpdf(et,zq,qvt,n,D); 
                At = Rt*F(:,t)/qvt;
                n=n+1; h=n+q-1; D = D+et*et'/qvt;  St=D/n; St=(St+St')/2; 
                Mt = Mt + At*et'; Ct = Rt - At*At'*qvt;  
             end
            ['Model ',int2str(im),' of ',int2str(nm)]
            mlik(arp+1:end,im) = sloglik ; % save log marginal likelihood 
         end
     end
 end
 Y=saveY; T=size(Y,2); % reset full data set
 savemlik=mlik; 
 
 ilik=250;  mlik(1:ilik,:)=[]; Tlik=T-ilik; % remove initial values of marg lik to drop impact of initial priors
 iT=ilik+1:T; plot(iT,mlik);
 plik=cumsum(mlik);  plik=exp(plik-repmat(max(plik')',1,nm)); 
                     plik=plik./repmat(max(plik')',1,nm);
 [a,i]=max(plik(end,:))
  
 save ex4_3_modelLikelihoods.mat mlik ilik iT Tlik pars savemlik Y T nm Beta Delta Arp

 
 

