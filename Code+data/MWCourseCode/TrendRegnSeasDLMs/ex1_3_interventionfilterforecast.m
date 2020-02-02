
startup


for i=1:6, figure(i); clf; end; figure(1)
fprintf('Using figure 1 only ...\n\n'); pause

% load data set and plot 
load salesindex
plot(1:T,sales,1:T,index); 
eval(xa);   legend('Sales','Index','location','east'); legend boxoff
    
    pause
% components of a DLM:
% trend: 
Ftrend=[1]; Gtrend=[1];
  % Ftrend=[1 0]'; Gtrend=[[1 1];[0 1]]; 
ntrend=length(Ftrend); itrend=1:ntrend;          
% regn:
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1 3 4]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 
Y=sales;

% priors p(\theta_1,v|D_0): 
nu=6; S=0.15^2; 
a=[9.5 -0.7 0.691 1.159 0.283 -0.050 -0.217 0.144]';
R=blkdiag(0.09,0.01,0.0067*eye(nseas)); 
deltrend=0.9; delregn=0.98; delseas=0.95; 

sm=zeros(n,T); sa=zeros(n,T); sC=zeros(n,n,T); sR=zeros(n,n,T); 
sf=zeros(1,T); sQ=zeros(1,T); sS=zeros(1,T); snu=zeros(1,T); 

% forward filtering and forecasting at each time point if iforecast is set
iforecast=1; k=2*p; ftk=zeros(1,k); Qtk=zeros(1,k); Yk=[Y ones(1,k)*NaN]; Xk=[X ones(1,k)*NaN];

for t=1:T
    
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         
    end
    
    if (t==31)
        fprintf('Intervention time! ...\n'); pause; pause
        fprintf('Intervention at time t=%2i - Decreased seasonal discount\n',t)
        R(iseas,iseas)=R(iseas,iseas)/delseas^12;
        pause
    end
    
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    if iforecast
        W=blkdiag(R(itrend,itrend)*(1-deltrend)/deltrend, R(iregn,iregn)*(1-delregn)/delregn, R(iseas,iseas)*(1-delseas)/delseas);
        F(iregn)=X(:,t)';  atk=a; Rtk=R; Qtk(1)=Q; ftk(1)=f;
        for h=2:k
            F(iregn)=Xk(:,t+h-1)';  atk=G*atk; Rtk=G*Rtk*G'+W; Qtk(h)=F'*Rtk*F+S; ftk(h)=F'*atk;
        end    
        figure(1); clf; subplot(2,1,1); plot(1:t-1,Y(:,1:t-1),'ko'); eval(xa); hold on
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
    subplot(2,1,2); hold on;  axis([0 T+1 -1.2 -0.1]);
        sq=tinv(.95,nu);  
        h=sqrt(squeeze(sC(iregn,iregn,1:t)))'*sq; ciplot(sm(iregn,1:t)-h,sm(iregn,1:t)+h,1:t,[.85 .85 .85]); hold on
        plot(1:t,sm(iregn,1:t),'b+',t+1:t+k,ones(1,k)*NaN); eval(xa); hold off
         title('90% on-line posterior intervals for Index regn '); ylabel('Parameter')
        display('hit a key to proceed to next time point')
        pause
    end
    
end


