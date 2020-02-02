% FX data 
% TV-VAR model for actual FX and commodity prices rather than returns
% Portfolio examples
% q=13 dimensional series 2000 - 2011 of daily data on FX and commodities 

% This file runs summaries of analysis based on the compute-intensive 
% full analysis in ex4_2_fullanalysisForecasting+Portfolios.m

startup

fprintf('\n Daily FX & commodity prices .. TVVAR analysis & porfolios \n\n')

for i=1:3, figure(i); clf; end; figure(1)
fprintf('Using figures 1, 2 & 3 ...\n\n')

load ex4SAVEportfoliosARP=3k=5beta=98.mat
fprintf('Analysis results loaded:\n')
fprintf(['   TVAR(',int2str(arp),') model with MV volatility\n'])
fprintf(['   Discount rates:  delta(TVVAR)=',num2str(delta),', beta(SV)=',num2str(beta,2),'\n\n'])

%load ex4SAVEportfoliosARP=3k=5beta=1.mat
%fprintf('Analysis results loaded:\n')
%fprintf(['   TVAR(',int2str(arp),') model with MV volatility\n'])
%fprintf(['   Discount rates:  delta(TVVAR)=',num2str(delta),', beta(SV)=',num2str(beta,2),'\n\n'])

% % example of plot
r=[4 5 11]; figure(1); clf
for j=1:3
    subplot(3,1,j); plot(Y(r(j),:)'); ylabel('Log prices'); title(names(r(j),:))
    eval(xa) ; end

fprintf('Hit return to continue ...\n\n'); pause
fprintf(['Hit return to see some of the ',int2str(k),'-step ahead forecasts ...\n\n'])

% look at some k-step ahead forecasts: 
 for j=1:q
        figure(1); clf; 
        plot(1:T,Y(j,:),tf+1+k:T+k,ffore(j,:)); title(char(Ynames(j)));
        legend('Data',[int2str(k),'-step ahead forecasts'],'location','southeast'); legend boxoff
        eval(xa); 
        figure(2); clf; 
        plot(1:T,Y(j,:),tf+1+k:T+k,ffore(j,:)); title(char(Ynames(j)));
        legend('Data',[int2str(k),'-step ahead forecasts'],'location','southeast'); legend boxoff
        eval(xa);   xlim([T-250 T+k])      
        pause; end

    
 
fprintf(['Hit return to see examples of portfolios at the end time ...\n\n']); pause
% target mean r and minvar portfolios ..
% We are working here with data y that are LOG PRICES so need to transform back to p = exp(y) for decisions
% If y ~ N(f,V) then p=exp(y) is multivariate lognormal with mean vector m
%             and variance matrix M where:   m = exp(f+diag(V)/2) and M = (m*m').*(exp(V)-1).
% As we need to use (m,M) as the mean and variance matrix to predict actual prices
% Then,  do move to portfolio optimizations focussed on returns,  this converts to exact calculations for returns as follows:

% e.g.:  choose portfolios k steps ahead from final time point:

f=ffore(:,end); V=Vfore(:,:,end);       V=(V+V')/2;                     % these are forecast mean and variance matrix for log prices
m = exp(f+diag(V)/2); M=(m*m').*(exp(V)-1);                             % mean & var mx of actual prices
% now find mean vector and variance matrix of *forecast returns* ....
pr = exp(Y(:,end));                                      % current actual prices
mr = 100*(m./pr-1);                                      % exact mean vector of forecast returns on percent scale
Mr = 10^4*lprod(1./pr,rprod(M,1./pr));   Mr=(Mr+Mr')/2;  % exact variance matrix of forecast returns on percent scale
fprintf('Final %i-step ahead forecast mean and sds for percent returns:\n',k)
[mr sqrt(diag(Mr))]
fprintf('\n\n Plotting correlations ...\n\n')
    figure(1); clf
    s=1./sqrt(diag(Mr)); C=lprod(s,rprod(Mr,s)); imagesc(C);  colormap hot; colorbar  
    xticks(1:p); xticklabels(names); xtickangle(90); yticks(1:p); yticklabels(names);
% Now portfolios ...  
r=min(.5,0.8*max(mr)); %  r is the target percent return 
fprintf('Hit return to see optimal portfolios (target return=%3.2f):\n',r); pause
[ wtarget wtargetSD wtargcon wtargconSD wminvar wminvarSD] = portfolios(mr,Mr,r,char(Ynames),1);

 

fprintf(['Hit return for results of sequential cumulative portfolios ...\n\n']); pause
% target mean r and minvar portfolios ..
% reallocate each k time points and see how much we make ... 
%
ftime = T-size(Yreturns,2)+1:T;  tticks=1:365:length(ftime); 
ti=int2str(time(ftime)'); yr=ti(tticks,3:4); mo=ti(tticks,5:6); da=ti(tticks,7:8);
tdates=cell(size(tticks)); for i=1:length(tticks), 
    tdates{i}=[mo(i,:)  '/' yr(i,:)]; end; %'/' da(i,:) '|' ]; end; 
xa=['set(gca,''Xtick'',tticks);set(gca,''XtickLabel'',tdates);xlabel(''mo/yr''); box off;']; 

figure(1); clf
%subplot(3,1,1); plot(wealth'); eval(xa)
subplot(2,1,1); plot(100*cumprod(1+wealth(1:3,:)'/100)-1); eval(xa)
    title(['Portfolios on ',int2str(k),'-day returns: TV-VAR(',int2str(arp),'),\beta=',num2str(beta,2)])
        legend('targcon','target','minvar','location','southeast'); legend boxoff; ylabel('%','rotation',0) 
           text (50,150,'Cumulative returns'); %ylim([0 160])
subplot(2,1,2); plot(wsdall'); eval(xa); text(50,9,'Risk'); ylabel('%','rotation',0); ylim([0 30])

figure(2); clf
 subplot(3,1,1); plot(squeeze(wall(4,1,:))'); eval(xa); title('GBP: targcon'); line([0 1+length(ftime)],[0 0],'color','k','linestyle',':')
 subplot(3,1,2); plot(squeeze(wall(4,2,:))'); eval(xa); title('target'); line([0 1+length(ftime)],[0 0],'color','k','linestyle',':')
 subplot(3,1,3); plot(squeeze(wall(4,3,:))'); eval(xa); title('minvar'); line([0 1+length(ftime)],[0 0],'color','k','linestyle',':')
 
figure(3); clf
 subplot(3,1,1); plot(squeeze(wall(11,1,:))'); eval(xa); title('GOLD: targcon'); line([0 1+length(ftime)],[0 0],'color','k','linestyle',':')
 subplot(3,1,2); plot(squeeze(wall(11,2,:))'); eval(xa); title('target'); line([0 1+length(ftime)],[0 0],'color','k','linestyle',':')
 subplot(3,1,3); plot(squeeze(wall(11,3,:))'); eval(xa); title('minvar'); line([0 1+length(ftime)],[0 0],'color','k','linestyle',':')

  
