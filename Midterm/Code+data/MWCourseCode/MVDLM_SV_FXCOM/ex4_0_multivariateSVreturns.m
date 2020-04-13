% FX returns data 
% Look at filtered and smoothed estimates of volatilities 
% q=13 dimensional series 2000 - 2011 of daily data on FX and commodities 

startup 

fprintf('Daily FX & commodity returns .. volatility tracking: filtering & smoothing\n\n')

for i=3:-1:1, figure(i); clf; end
fprintf('Using figures 1, 2 & 3  ...\n\n'); pause

[Y,Ynames]=xlsread('FXCommData2000-2011.xlsx');
time=Y(:,1)'; Y(:,1)=[]; Ynames(9)=[]; Y(:,9)=[]; T=size(Y,1); % deletes SGD 
for j=1:9, a=char(Ynames{j}); a(1:4)=[]; Ynames{j}=a; end; Ynames{5}='JPY'; Ynames{11}='GOL'; Ynames{12}='NSD';  % tidy up name

% returns:
Prices=Y'; 
Y=100*(Y(2:T,:)./Y(1:T-1,:)-1)'; Prices(:,1)=[]; T=T-1;  [q,T]=size(Y);  names=char(Ynames); names=names(:,1:3);
% timing axis labels ...
tticks=1:365:T; ti=int2str(time'); yr=ti(tticks,3:4); mo=ti(tticks,5:6); da=ti(tticks,7:8);
tdates=cell(size(tticks)); for i=1:length(tticks), 
    tdates{i}=[mo(i,:)  '/' yr(i,:)]; end; %'/' da(i,:) '|' ]; end; 
xa=['set(gca,''Xtick'',tticks);set(gca,''XtickLabel'',tdates);xlabel(''mo/yr''); box off;'];
% 


% % example of plot
r=[4 5 11]; figure(1); clf
for j=1:3
    subplot(3,1,j); plot(Y(r(j),:)'); ylabel('returns'); title(names(r(j),:))
    eval(xa) ; end

p=1;    F=ones(1,T);           % simple local mean model with volatility 

%%%%% - change discounts to assess, compare marginal likelihood

delta=0.99; % discount level
beta =0.98;   % discount volatility

%n0=3; h0=n0+q-1;  
h0=1/(1-beta); n0=h0-q+1; 
D0=h0*eye(q);  z = zeros(p,q);   zq=zeros(q,1); M0=z; r=0.99; % priors
Mt = M0; C0=.01*eye(p);  Ct=C0;        % initial Theta prior 
n = n0; h=h0; D = D0;  St=D/h;         % initial Sigma prior
sMt=zeros(p,q,T);  sCt=zeros(p,p,T);  sdCt=zeros(p,q,T); sSt=zeros(q,q,T);  
snt=zeros(1,T); sloglik=zeros(1,T); 


fprintf('Data and model set up. Hit a key for forward filtering ...\n')
pause


% forward filtering: 
        for t = 1:T
            ft = Mt'*F(:,t);             
            et = Y(:,t) - ft;
            Rt = Ct/delta; 
            h  = beta*h;  n=h-q+1;  D = beta*D;   
            qvt = 1 + F(:,t)'*Rt*F(:,t); sEt(:,t) = et./sqrt(qvt*diag(St)); 
                                         sloglik(t) = ltpdf(et,zq,qvt,n,D); 
            At = Rt*F(:,t)/qvt;
            h=h+1; n=n+1; D = D+et*et'/qvt;  St=D/h; St=(St+St')/2; 
            Mt = Mt + At*et'; Ct = Rt - At*At'*qvt;   sCt(:,:,t)=Ct;
            snt(t)=n; 
            sSt(:,:,t)=St; sMt(:,:,t) = Mt; sdCt(:,:,t) = sqrt(diag(Ct)*diag(St)'); 
        end
        sht=snt+q-1; 


  
%---------- PLOTS --------------
%  Some volatilities over time 
%
figure(1); clf; figure(3); clf
JJ=[4 5 11];
for r=1:3
    j=JJ(r); 
    figure(1); subplot(3,1,r)
        plot(sqrt(squeeze(sSt(j,j,:))),'b')
        eval(xa); title([ ,'SD(',names(j,:),')'])
    figure(3); subplot(3,1,r)  
    mt=squeeze(sMt(1,j,:)); sct=squeeze(sdCt(1,j,:)); 
    plot(1:T,mt,'b'); hold on
    plot(1:T,mt+sct,'color',[.6 .6 .6 ]);  plot(1:T,mt-sct,'color',[.6 .6 .6 ])
    line([0 T],[0 0],'color','k','linestyle',':');  
    title(names(j,:))
    eval(xa); ylabel('Level +/-1 s.d.')
end

%  Some correlations over time 
%
figure(2); clf
JJ=[[4 5];[4 11];[4 12]];
for r=1:3
    j=JJ(r,:); subplot(3,1,r)
        plot(squeeze(sSt(j(1),j(2),:))./sqrt(squeeze(sSt(j(1),j(1),:)).*squeeze(sSt(j(2),j(2),:))),'b')
        line([0 T],[0 0],'color','k','linestyle',':');  
        eval(xa); title([ 'Corr(',names(j(1),:),',',names(j(2),:),')'])
end


fprintf('Hit a key to perform retrospective smoothing ...\n')
pause
      

% reverse smoothing  -
K=inv(sSt(:,:,T)); n=snt(T); Mt = sMt(:,:,T); Ct = sCt(:,:,T); 
for t=(T-1):-1:1
    K=(1-beta)*inv(sSt(:,:,t))+beta*K;          St = inv(K); sSt(:,:,t)=St;  
    Mt = (1-delta)*sMt(:,:,t) +delta*Mt;        sMt(:,:,t) = Mt; 
    Ct = (1-delta)*sCt(:,:,t) + delta^2*Ct;     sCt(:,:,t) = Ct; sdCt(:,:,t) = sqrt(diag(Ct)*diag(St)');
end
   

% add to plots ...

JJ=[4 5 11];
for r=1:3
    j=JJ(r); 
    figure(1); subplot(3,1,r); hold on
        plot(sqrt(squeeze(sSt(j,j,:))),'r')
    figure(3); subplot(3,1,r); hold on  
    mt=squeeze(sMt(1,j,:)); sct=squeeze(sdCt(1,j,:)); 
    plot(1:T,mt,'r')
end

%  Some correlations over time 
%
figure(2)
JJ=[[4 5];[4 11];[4 12]];
for r=1:3
    j=JJ(r,:); subplot(3,1,r); hold on
    plot(squeeze(sSt(j(1),j(2),:))./sqrt(squeeze(sSt(j(1),j(1),:)).*squeeze(sSt(j(2),j(2),:))),'r')
end





%-------------------------------------------------------------------------------------
% Now, Backward Sampling and its use in dynamic PCA ... aka empirical dynamic factor analysis
% Uses Monte Carlo sample size for retrospective trajectories of the matrix volatility process 
nmc=100; mcS=zeros(q,q,T,nmc); 
% now start at time T: 
snt=sht-q+1; D = sht(T)*squeeze(sSt(:,:,T)); mcK=Wishartrnd(snt(T),D,nmc); 
for im=1:nmc, mcS(:,:,T,im)=inv(squeeze(mcK(im,:,:))); end
for t=(T-1):-1:1 
    t
    D = sht(t)*squeeze(sSt(:,:,t)); 
    mcK = beta*mcK + Wishartrnd((1-beta)*sht(t)-q+1,D,nmc); 
    for im=1:nmc, mcS(:,:,t,im) = inv(squeeze(mcK(im,:,:))); end
end
mS=mean(mcS,4); 
TT=T;


%---------- PLOTS --------------
% eigendecomositions  
% 
S=mS;
E=S; D=zeros(q,T); L=Y;
for t=1:T
   [e,d]=eig(S(:,:,t)); d=sqrt(diag(d)); [d,i]=sort(d,'descend'); e=e(:,i);    
   D(:,t)=d;  
   % now, fix arbirary sign of evecs .... 
   if (t>1)
       a=e-E(:,:,t-1); f=-e-E(:,:,t-1);
       [a,i]=min([sum(a.*a);sum(f.*f)]); j=find(i==2); e(:,j)=-e(:,j);
   end
   %
   E(:,:,t)= e; % rprod(e(:,i),d); 
   L(:,t)=rprod(e,1./d)'*Y(:,t);
end


% posterior samples of pcs - look at variation explained by first npca
c=q; pc = zeros(q,T,nmc); 
for j=1:nmc
   for t=1:T
       [e,d]=eig(mcS(:,:,t,j)); d=sqrt(diag(d)); [d,i]=sort(d,'descend'); 
       pc(:,t,j)=d(1:c); 
   end
end

npca=3; 
% variation: 
cd  = cumsum(pc); cd=100*pc./repmat(cd(q,:,:),q,1); cd=squeeze(mean(cd,npca));
% plots ....  
clf; axv=[0 TT+1 0 2.5];
subplot(4,1,1,'replace'); plot(1:TT,sqrt(squeeze(pc(1,1:TT,:))),'color',[0.75 0.75 0.75]); axis tight; box off
 hold on; plot(sqrt(squeeze(D(1,:))),'k'); hold off
 axis(axv); ylabel('volatility')
 set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates);
 title('(a) First volatility component') 
subplot(4,1,2,'replace'); plot(1:TT,sqrt(squeeze(pc(2,1:TT,:))),'color',[0.75 0.75 0.75]); axis tight; box off
 hold on; plot(sqrt(squeeze(D(2,:))),'k'); hold off
 axis(axv); ylabel('volatility')
 set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates);
 title('(b) Second volatility component') 
subplot(4,1,3,'replace'); plot(1:TT,sqrt(squeeze(pc(3,1:TT,:))),'color',[0.75 0.75 0.75]); axis tight; box off
 hold on; plot(sqrt(squeeze(D(3,:))),'k'); hold off
 axis(axv); ylabel('volatility')
 set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates);
 title('(c) Third volatility component') 
subplot(4,1,4,'replace'); 
 plot(1:TT,cd(1,1:TT),'k-',1:TT,cd(2,1:TT),'k:',1:TT,cd(3,1:TT),'k--'); axis([0 TT+1 0 100]); box off
 legend('comp 1','comp 2','comp 3','Orientation','horizontal'); legend boxoff; ylabel('% variation');
 set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates);xlabel('days within month/year')
 title('(d) Contributions of components')
  

% ------------ 
% multiscale components in returns - just first C and 'the rest' ... 
%  try with C=3
Dv=(D./repmat(max(D')',1,T)).^2; 
C=3; axr=[0 TT+1 -5 5]; axv=[0 TT+1 0 0.25];  % You can change the plot ranges displayer
for j=1:q
    names(j,:)
    e=squeeze(E(j,:,:)); vj=sum((e.*D).^2); 
    axr([3 4]) = [ min(Y(j,:)) max(Y(j,:))]; axv([3 4]) = [ 0 max(vj) ]; 
    rjc=D(1:C,:).*e(1:C,:).*L(1:C,:); 
    resrj=sum(D((C+1):q,:).*e((C+1):q,:).*L((C+1):q,:));
    vjc=(D(1:C,:).*e(1:C,:)).^2; resvj=sum((D((C+1):q,:).*e((C+1):q,:)).^2);
    figure(1); clf
    subplot(C+2,1,1); plot(1:TT,Y(j,1:TT),'k'); box off; 
    title(['(a) ',names(j,:),' returns']); axis(axr)
    set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates)
    subplot(C+2,1,C+2); plot(1:TT,resrj(1:TT),'k'); box off;  
    title('(d) Other returns components'); axis(axr)
    set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates)
    xlabel('days within month/year')
    figure(2); clf
    subplot(C+2,1,1); plot(1:TT,vj(1:TT),'k'); box off;  
    title(['(a) ',names(j,:),' volatility']); axis(axv)
    set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates)
    subplot(C+2,1,C+2); plot(1:TT,resvj(1:TT),'k'); box off;  
    title('(d) Other volatility components'); axis(axv)
    set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates)
    xlabel('days within month/year')
    cc=char('(b) ','(c) ','(d) ','(e) ','(f) ','(g) '); 
    for c=1:C, 
        figure(1); subplot(C+2,1,1+c); plot(1:TT,rjc(c,1:TT),'k'); box off; axis(axr)
        set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates)
        title([cc(c,:), 'Returns component ',int2str(c)]) 
        figure(2); subplot(C+2,1,1+c); plot(1:TT,vjc(c,1:TT),'k'); box off; axis(axv)
        set(gca,'Xtick',tticks);set(gca,'XtickLabel',tdates)
        title([cc(c,:) 'Volatility component ',int2str(c)])
            hold on; plot(axv(4)*Dv(c,:),'color',[0.5 0.5 0.5],'linestyle',':'); hold off
    end
    %
    input('Hit a key to see decomposition of next series\n');   
end




 