 
% read in and plot the detrended O2 isotope time series - which is in 3,000 years reverse time

y=load('o2plusdetrended.txt'); ti=y(:,1); x=y(:,3); T=length(ti); 
xa='box off; axis tight; xlabel(''Reverse time (in kyears)''); ylabel('''')';
figure(1); clf 
subplot(2,1,1); plot(ti,x); 
title('Detrended O2 Ratio Series: t=1:866'); eval(xa)


% now make forward time as we are interested in forecasting from now ...
x=flipud(x); ti=-flipud(ti);
xa='box off; axis tight; xlabel(''Time in kyears''); ylabel('''')';
subplot(2,1,2); plot(ti,x); eval(xa)



 