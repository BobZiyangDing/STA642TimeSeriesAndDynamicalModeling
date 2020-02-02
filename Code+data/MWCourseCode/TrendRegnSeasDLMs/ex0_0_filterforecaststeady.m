

startup

for i=1:6, figure(i); clf; end; figure(1)
fprintf('Using figure 1 only ...\n\n')

% load data set and plot 
Y=load('uk-usa1000.txt'); Y=Y(700:900); T=length(Y); 
plot(1:T,Y); 
pause
% components of a DLM:
% trend: 
Ftrend=[1]; Gtrend=[1]; F=Ftrend; G=Gtrend; 
ntrend=length(Ftrend); itrend=1:ntrend;     n=1; 


% priors p(\theta_1,v|D_0): 
nu=6; S=0.12; 
a=0; R=blkdiag(0.09); 
deltrend=0.9;  

sm=zeros(n,T); sa=zeros(n,T); sC=zeros(n,n,T); sR=zeros(n,n,T); 
sf=zeros(1,T); sQ=zeros(1,T); sS=zeros(1,T); snu=zeros(1,T); 

% forward filtering and forecasting at each time point if iforecast is set
iforecast=1; k=20; ftk=zeros(1,k); Qtk=zeros(1,k); Yk=[Y ones(1,k)*NaN];  

for t=1:T
    
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
    end
    A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    if iforecast
        W=blkdiag(R(itrend,itrend)*(1-deltrend)/deltrend);
        atk=a; Rtk=R; Qtk(1)=Q; ftk(1)=f;
        for h=2:k
            atk=G*atk; Rtk=G*Rtk*G'+W; Qtk(h)=F'*Rtk*F+S; ftk(h)=F'*atk;
        end    
        figure(1); clf; 
        plot(1:t-1,Y(:,1:t-1),'ko'); xlim([0 T+1]); ylim([min(Y) max(Y)]);
        hold on
        sq=tinv(.95,nu);  h=sqrt(Qtk).*sq; ciplot(ftk-h,ftk+h,t:t+k-1,[.85 .85 .85])
        plot(t:t+k-1,ftk); 
        plot(t,Y(:,t),'ro',t+1:t+k-1,Yk(:,t+1:t+k-1),'b+'); hold off
        title('90% prediction intervals and step ahead forecasts'); ylabel('Y')
    end
    y = Y(:,t);  e=y-f; r=(nu+e^2/Q)/(nu+1);
    m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;  
    % save:
    sm(:,t)=m; sC(:,:,t)=C; sa(:,t)=a; sR(:,:,t)=R; sf(t)=f; sQ(t)=Q; sS(t)=S; snu(t)=nu;
    pause
 
end


