
% Dynamic dependency network models: Example with US macro data
startup
readddnmdata     % read in data and plot series over time 
setddnm          % define parental sets and lagged predictors for the series in a chosen order
setpriors        % specify priors and discount factors

names
icontrol=1; % index of T-bill/Fed funds target interest rate- potential control variable to change
 
% Now perform FF over an initial time period, then look at k-step ahead forecasts at
% a specific time point tfore 
%  
 
 
k =16;                  % generate k-step ahead forecasts at t=tfore
%tfore = 156; 
tfore = randsample(100:T-k-1,1);   % stop FF and forecast at this time; your choice 
lastt = max(1,tfore+1-5*k);   % past data to show on plots 
I = 1000;               % Monte Carlo simulation sample size for all forecasts  

for j=1:q             % learning phase 1:tfore
    p0.m=prior{1,j}; p0.C=prior{2,j}*prior{4,j}; p0.n=prior{3,j}; p0.s=prior{4,j}; 
    p0.delta=prior{5,j}; p0.beta=prior{6,j}; 
    F = ones(1,tfore); if (npa(j)>0), F=[F; Y(pa{j},1:tfore)]; end
    if (npr(j)>0), 
        for h=1:t0, 
            i=find(squeeze(pr(j,:,h))); 
            if (length(i)>0), F = [ F; [zeros(length(i),h) Y(i,(h+1:tfore)-h) ] ]; end
        end
    end    
    pq1tfore(j) = ff(Y(j,1:tfore),F,tfore,t0,p0);   
end

% now we are standing at time t=tfore looking ahead  
%
t = tfore; 
tpred = tfore+(1:k);             % next k time points 
spred = zeros(1,k,I);            % saved synthetic futures 
ypred = zeros(q,7,k);            % saved forecast summaries  
                                 %  5-fig forecast summaries: [ mean, median, SD, 2.5%, 97.5% ]                                 
                                 
rr = Y(icontrol,tfore);          % selected value of T Bill r_{t+1}=rr - may choose any value
kstepforecast                % perform 1:k step predictive simulations at current t = tfore   

 
ifore=find(yrmoq==yrmoq(tfore)); moqtfore=find(ifore==tfore);
for j=1:q      % plot 1:k-step ahead forecast summaries ...
    fmean=squeeze(ypred(j,1,:))';fmedian=squeeze(ypred(j,2,:))';
    figure(j); clf   
    ciplot(squeeze(ypred(j,4,:)),squeeze(ypred(j,7,:)),tpred,[.95 .95 .95]);
    hold on
    ciplot(squeeze(ypred(j,5,:)),squeeze(ypred(j,6,:)),tpred,[.85 .85 .85]);       
    plot((lastt+1):tfore,Y(j,(lastt+1):tfore),'k+',tpred,fmean,'r-',tpred,fmedian,'b--')
    eval(xa); ylabel(deblank(names(j,:)))
    title({[deblank(names(j,:)),': ','1:',int2str(k),'-moq ahead predictions'],
           ['from baseline t=moq ',int2str(moqtfore),'/',int2str(yrmoq(tfore))],
           [deblank(names(icontrol,:)),'_t\leftarrow',num2str(rr,2),'%']}) 
    xlim([lastt-1 tfore+k+1]); %ylim(mmy(j,:))
    hold on; scatter(tpred,Y(j,tpred),20,'r+'); hold off
end
figure(icontrol); hold on; 
  scatter(tfore,rr,20,'ro'); scatter(tfore,rr,20,'k+');
  hold off
   

    
    
  