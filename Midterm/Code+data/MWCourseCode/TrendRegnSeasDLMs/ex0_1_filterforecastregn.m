
startup

figure(1) 
% load data set and plot 
load salesindex
plot(1:T,sales,1:T,index); eval(xa);    
legend('Sales','Index','location','east'); legend boxoff
    
figure(2)
scatter(index,sales); box off
xlabel('Index'); ylabel('Sales')

    pause
% components of a DLM:
% trend: 
Ftrend=[1]; Gtrend=[1];
ntrend=length(Ftrend); itrend=1:ntrend;          
% regn:
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% overall DLM:
F = [Ftrend;Fregn]; n=length(F);  G = blkdiag(Gtrend,Gregn); 
Y=sales;

% priors p(\theta_1,v|D_0): 
nu=6; S=0.15^2;   % IG(nu/2, nu.S/2) 
a=[9.5 -0.7]';
R=blkdiag(0.09,0.01);  % N(a,Rv) 
% discount factor to determine W_t
delta=0.99; 

sm=zeros(n,T); sa=zeros(n,T); sC=zeros(n,n,T); sR=zeros(n,n,T); 
sf=zeros(1,T); sQ=zeros(1,T); sS=zeros(1,T); snu=zeros(1,T); 

% forward filtering and forecasting at each time point if iforecast is set
iforecast=1; k=24; % yes to forecast, k months ahead 
ftk=zeros(1,k); Qtk=zeros(1,k); Yk=[Y ones(1,k)*NaN]; Xk=[X ones(1,k)*NaN];

for t=1:T
    
    if (t>1) a = G*m; 
             R = G*C*G'; R=R/delta; 
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    if iforecast
        W=R*(1-delta)/delta; 
        F(iregn)=X(:,t)';  atk=a; Rtk=R; Qtk(1)=Q; ftk(1)=f;
        for h=2:k
            F(iregn)=Xk(:,t+h-1)';  atk=G*atk; Rtk=G*Rtk*G'+W; Qtk(h)=F'*Rtk*F+S; ftk(h)=F'*atk;
        end    
        figure(1); clf; 
        subplot(2,1,1); plot(1:t-1,Y(:,1:t-1),'ko'); eval(xa); hold on
        sq=tinv(.95,nu);  h=sqrt(Qtk).*sq; ciplot(ftk-h,ftk+h,t:t+k-1,[.85 .85 .85])
        plot(t:t+k-1,ftk); 
        plot(t,Y(:,t),'ro',t+1:t+k-1,Yk(:,t+1:t+k-1),'b+'); hold off
        title('90% prediction intervals and step ahead forecasts'); ylabel('Sales')
    end
    y = Y(:,t);  e=y-f; r=(nu+e^2/Q)/(nu+1);
    m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;  
    % save:
    sm(:,t)=m; sC(:,:,t)=C; sa(:,t)=a; sR(:,:,t)=R; sf(t)=f; sQ(t)=Q; sS(t)=S; snu(t)=nu;
    if iforecast
    subplot(2,1,2); hold on;  axis([0 T+1 -2 0]);
        sq=tinv(.95,nu);  
        h=sqrt(squeeze(sC(iregn,iregn,1:t)))'*sq; ciplot(sm(iregn,1:t)-h,sm(iregn,1:t)+h,1:t,[.85 .85 .85]); hold on
        plot(1:t,sm(iregn,1:t),'b+',t+1:t+k,ones(1,k)*NaN); eval(xa); hold off
         title('90% on-line posterior intervals for Index regn '); ylabel('Parameter')
        display('hit a key to proceed to next time point')
        pause
    end
 
end
display('hit a key to proceed to summary learning on v')
pause
clf
%sm=zeros(n,T); sa=zeros(n,T); sC=zeros(n,n,T); sR=zeros(n,n,T); 
%sf=zeros(1,T); sQ=zeros(1,T); sS=zeros(1,T); snu=zeros(1,T); 
subplot(2,1,1); plot(1:T,sqrt(sS)); eval(xa); ylabel('S_t^{1/2}')
subplot(2,1,2); plot(1:T,snu); eval(xa); ylabel('n_t')


