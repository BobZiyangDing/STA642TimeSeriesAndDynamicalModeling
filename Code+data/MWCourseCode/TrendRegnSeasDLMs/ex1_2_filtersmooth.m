startup

% this version does filtering and smoothing AND adds variance discounting -
% delvar

for i=1:6, figure(i); clf; end
fprintf('Using all figures 1-6 ...\n\n'); pause

% DLM analysis:  Set the following flag to 0 or 1 prior to running 
% iforecast=1;     % to stop and produce, plot step-ahead forecasts for t+1,...,t+k at each t
                   %       in the foreward filtering sequential analysis 
 iforecast=1;      % to just perform the forward updating with no forecasts at each time ...               
              
% load data set and plot 
load salesindex
xa='set(gca,''Xtick'',1:p:T,''XtickLabel'',tdates); xlabel(''Date (month/year)''); box off; xlim([0 T+1])';

    figure(1); clf
    plot(1:T,sales,1:T,index); eval(xa)
    legend('Sales','Index','location','east'); legend boxoff
    
% components of a DLM:
% trend: % e.g. locally constant:  
     Ftrend=[1]; Gtrend=[1];
%    Ftrend=[1 0]'; Gtrend=[[1 1];[0 1]];    % locally linear  
    ntrend=length(Ftrend); itrend=1:ntrend;          
% regn:
    q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); 
    nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
    p=12; rseas=[1 3 4]; pseas=length(rseas); 
    nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas;    
    Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
% DLM matrices F & G: 
    F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 
% response time series: 
    Y=sales;
% priors p(\theta_1,v|D_0) & 3 discount factors for component DLMs: 
    nu=6; S=0.15^2; 
    a=[9.5 -0.7 0.691 1.159 0.283 -0.050 -0.217 0.144]';
    R=blkdiag(0.09,0.01,0.0067*eye(nseas)); 
    deltrend=0.9; delregn=0.98; delseas=0.95; delvar=0.9; 
% storage arrays for filtering summaries: 
    sm=zeros(n,T); sa=zeros(n,T); sC=zeros(n,n,T); sR=zeros(n,n,T); 
    sf=zeros(1,T); sQ=zeros(1,T); sS=zeros(1,T); snu=zeros(1,T); 
    mlik=zeros(1,T); % to save the t densities at each point that give the model marginal likelihood
% 
display('Hit a key to start analysis ...')
pause

if iforecast
        k=p; ftk=zeros(1,k); Qtk=zeros(1,k); Yk=[Y ones(1,k)*NaN]; Xk=[X ones(1,k)*NaN]; 
        display(['Running DLM updating analysis with ',int2str(k),' step-ahead forecasts at each t ...'])
else
        display('Running DLM updating analysis - no step-ahead forecasts  ...')
end
% now analysis ...
    for t=1:T
        if (t>1) a = G*m; 
                 R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                             R(iregn,iregn)=R(iregn,iregn)/delregn; 
                             R(iseas,iseas)=R(iseas,iseas)/delseas;
                             nu=delvar*nu;% discount variance dof
        end
        Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
        if iforecast
            W=blkdiag(R(itrend,itrend)*(1-deltrend)/deltrend, R(iregn,iregn)*(1-delregn)/delregn, R(iseas,iseas)*(1-delseas)/delseas);
            F(iregn)=X(:,t)';  atk=a; Rtk=R; Qtk(1)=Q; ftk(1)=f;
            for h=2:k
                F(iregn)=Xk(:,t+h-1)';  atk=G*atk; Rtk=G*Rtk*G'+W; Qtk(h)=F'*Rtk*F+S; ftk(h)=F'*atk;
            end    
            figure(2); subplot(2,1,1); plot(1:t-1,Y(:,1:t-1),'ko'); eval(xa); hold on
            sq=tinv(.95,nu);  h=sqrt(Qtk).*sq; ciplot(ftk-h,ftk+h,t:t+k-1,[.85 .85 .85])
            plot(t:t+k-1,ftk); plot(t,Y(:,t),'ro',t+1:t+k-1,Yk(:,t+1:t+k-1),'b+')
            xlabel('Months ahead)'); hold off
            title('90% prediction intervals and step ahead forecasts'); ylabel('Sales')
        end
        y = Y(:,t);  
        % updating equations: 
            e=y-f; rQ=sqrt(Q); 
            mlik(t) = tpdf(e/rQ,nu)/rQ;
            r=(nu+e^2/Q)/(nu+1);
            m = a+A*e; C = r*(R-A*A'*Q); 
            nu=nu+1; S=r*S;  

        % save posterior at time t:
            sm(:,t)=m; sC(:,:,t)=C; sa(:,t)=a; sR(:,:,t)=R; sf(t)=f; sQ(t)=Q; sS(t)=S; snu(t)=nu;
        if iforecast
            subplot(2,1,2); eval(xa); hold on;  
            h=sqrt(squeeze(sC(iregn,iregn,1:t)))'*sq; ciplot(sm(iregn,1:t)-h,sm(iregn,1:t)+h,1:t,[.85 .85 .85]);  
            plot(1:t,sm(iregn,1:t),'b+',t+1:t+k,ones(1,k)*NaN); xlabel('Months');  hold off
             title('90% on-line posterior intervals for Index regn '); ylabel('Parameter')
            display('hit a key to proceed to next time point')
            pause
            %subplot(2,1,1); axis tight; subplot(2,1,2); axis tight
        end
    end

% -------------- Look at some summaries of the online posteriors/trajectories ..    
display('Finished updating: hit a key to plot some online trajectories  ...')
pause
% vector of quantiles of posterior t distributions at each time: 
    sq=tinv(.95,snu);
% plot 1-step forecasting and error summaries over time
    figure(1); clf   
    subplot(2,1,1)
    h=sqrt(sQ).*sq; ciplot(sf-h,sf+h,1:T,[.85 .85 .85]); hold on
    plot(1:T,sf,'b',1:T,Y,'r+'); eval(xa)
    title('90% prediction intervals and 1-step forecasts'); ylabel('Sales')
    subplot(2,1,2)
    ciplot(-sq,+sq,1:T,[.85 .85 .85]); hold on
    plot(1:T,(Y-sf)./sqrt(sQ),'r+'); plot([0 T+1],[0 0],'k:'); eval(xa)
    title('90% prediction intervals and standardized 1-step errors'); ylabel('Sales')
% plot trajectory of online estimated dynamic regn coeff and error SD: 
    figure(2); clf  
    subplot(2,1,1)
    h=sqrt(squeeze(sC(iregn,iregn,:)))'.*sq; ciplot(sm(iregn,:)-h,sm(iregn,:)+h,1:T,[.85 .85 .85]); hold on
    plot(1:T,sm(iregn,:),'b-'); eval(xa)
    title('90% on-line posterior intervals for Index regn '); ylabel('Parameter')
    subplot(2,1,2); plot(sqrt(sS)); eval(xa)
    title('On-line estimate of observation error s.d.'); ylabel('S.D.')
 % plot trajectory of online estimated harmonic seasonal components:   
    figure(3); clf
    for j=1:pseas
        i=iseas(2*j-1); 
        subplot(pseas,1,j)
        h=sqrt(squeeze(sC(i,i,:)))'.*sq; ciplot(sm(i,:)-h,sm(i,:)+h,1:T,[.85 .85 .85]); hold on
        plot(1:T,sm(i,:),'b-');  plot([0 T+1],[0 0],'k:'); eval(xa); hold off
        title(['90% on-line posterior intervals for harmonic ',int2str(rseas(j))]); ylabel('Harmonic')
    end
% plot marginal likihood and its per time poitn components: 
    figure(6); clf
    subplot(2,1,1); hold on
    plot(1:T,log(mlik)); title('Log marginal likelihood components'); ylabel('log marg lik'); eval(xa)
    subplot(2,1,2); hold on
    plot(1:T,cumsum(log(mlik))); title('Log cumulative marginal likelihood'); ylabel('log marg lik'); eval(xa)
% -------------------------------------------------------------------    
display('hit a key to perform retrospective updating ...')
pause
% backward smoothing for retrospection - n.b. overwrites saved information
% first save filtered summaries for later backward sampling below: 
    Sm=sm; SC=sC; SR=sR;  Sa=sa; SS=sS; Snu=snu; 
    % then perform backward (Viterbi style) updating - overwriting
    %  online posterior summaries at each time point with full posteriors: 
    C=sC(:,:,T); St=sS(T); nu=snu(T); 
    for t=T-1:-1:1
        B         = sC(:,:,t)*G'*inv(sR(:,:,t+1));              
        sm(:,t)   = sm(:,t)+B*(sm(:,t+1)-sa(:,t+1));        
        C         = sC(:,:,t)+B*(C-sR(:,:,t+1))*B';  
        St        = (1-delvar)/sS(t)+delvar/St; St=1/St; SS(t)=St; 
        nu        = (1-delvar)*snu(t)+delvar*nu; Snu(t)=nu; 
        sC(:,:,t) = C*St/sS(t);
    end
% posterior quantile of time T final posterior t distribution:         
    sq=tinv(.95,snu(T));
% plot smoothed/retrospective estimates and intervals for dynamic regn coeff: 
    figure(4); clf
    subplot(2,1,1)
    h=sqrt(squeeze(sC(iregn,iregn,:)))'.*sq; ciplot(sm(iregn,:)-h,sm(iregn,:)+h,1:T,[.85 .85 .85]); hold on
    plot(1:T,sm(iregn,:),'b-'); eval(xa); hold off
    title('90% smoothed posterior intervals for Index regn '); ylabel('Parameter')
% plot smoothed/retorspective estimates and intervals for harmonic seasonal components:      
    figure(5); clf
    for j=1:pseas
        i=iseas(2*j-1); 
        subplot(pseas,1,j)
        h=sqrt(squeeze(sC(i,i,:)))'*sq; ciplot(sm(i,:)-h,sm(i,:)+h,1:T,[.85 .85 .85]); hold on
        plot(1:T,sm(i,:),'b-');  plot([0 T+1],[0 0],'k:'); eval(xa); hold off
        title(['90% smoothed posterior intervals for harmonic ',int2str(rseas(j))]); ylabel('Harmonic')
    end
     
   
% -------------------------------------------------------------------    
display('hit a key to do some retrospective sampling ...')
pause
% backward sampling -- as above, no sampling full joint posterior of trajectories 
%       of state vectors over t=1:T - uses saved info on filtered
%       posteriors:  see  ****** West & Harrison section 15.2.3
%                    and  ****** Prado & West sections 4.3.5 and 4.5
    figure(4) % plot smoothed mean estimate trajectory as reference:
    subplot(2,1,2); plot(1:T,sm(iregn,:),'b-'); eval(xa); hold on
    title('Posterior samples of trajectories of Index regn '); ylabel('Parameter')
    col='rcmbg';   % colours for each of the 5 random samples of trajectories: 
    rm=sm; 
    for isam=1:5     % simulate and add trajectory to the plot: 
        rv = sqrt(1/gamrnd(Snu(T)/2,2/(Snu(T)*SS(T)),1,1)); 
        rm(:,T) = Sm(:,T)+rv*chol(SC(:,:,T)/SS(T))'*randn(n,1); 
        for t=T-1:-1:1
            B         = SC(:,:,t)*G'*inv(SR(:,:,t+1)); 
            C         = SC(:,:,t)-B*SR(:,:,t+1)*B'; 
            rm(:,t)   = Sm(:,t)+B*(rm(:,t+1)-Sa(:,t+1)) + rv*chol(C/SS(t))'*randn(n,1); 
        end
        plot(1:T,rm(iregn,:),col(isam)); 
    end
    hold off

% -------------------------------------------------------------------    
display('hit a key for detrending & deseasonalization ...')
pause
    
% plot smoothed/retrospective estimates and intervals for trend, regn, trend+regn 
%  and seasonal separately
    clf
    subplot(5,1,1);  
    plot(1:T,Y,'b+-'); eval(xa); title('Sales');  hold off; axis off
    subplot(5,1,3)
    h=sqrt(squeeze(sC(1,1,:)))'.*sq; g=sm(1,:);
    ciplot(g-h,g+h,1:T,[.85 .85 .85]); box off; ylim(mean(g)+[5 12]-mean([5 12])); hold on
    plot(1:T,g,'b-'); eval(xa); title('Local level');  hold off; axis off
    subplot(5,1,4)
    h=abs(X).*sqrt(squeeze(sC(iregn,iregn,:)))'.*sq; g=X.*sm(iregn,:);
    ciplot(g-h,g+h,1:T,[.85 .85 .85]); box off; ylim([5 12]-mean([5 12])); hold on
    plot(1:T,g,'b-'); eval(xa); title('Index regn');  hold off; axis off
    subplot(5,1,2)
    g=zeros(1,T); h=g; i=[itrend iregn];
    for t=1:T, 
        F= [ Ftrend; X(:,t) ]; g(t)=F'*sm(i,t); 
        h(t)=sqrt(F'*squeeze(sC(i,i,t))*F)*sq; 
    end; 
    ciplot(g-h,g+h,1:T,[.85 .85 .85]); box off; ylim([5 12]); hold on
    plot(1:T,g,'b-'); eval(xa); title('Deseasonalized'); hold off; axis off
    subplot(5,1,5)
    g=Fseas'*sm(iseas,:);
    for t=1:T, 
        h(t)=sqrt(Fseas'*squeeze(sC(iseas,iseas,t))*Fseas)*sq; 
    end; 
    ciplot(g-h,g+h,1:T,[.85 .85 .85]); box off; ylim(mean(g)+[5 12]-mean([5 12])); hold on
    plot(1:T,g,'b-'); eval(xa); title('Seasonal'); hold off; axis off
    
    
    
         