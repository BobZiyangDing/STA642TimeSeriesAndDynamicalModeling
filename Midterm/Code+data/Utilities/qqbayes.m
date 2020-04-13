
% Example of a Normal quantile plot of data x to provide a visual
% assessment of its conformity with a normal (data is standardised first)
%
% The ordered data values are posterior point estimates of the underlying
% quantile function. So, if you plot the ordered data values (y-axis) against
% the exact theoretical quantiles (x-axis), you get a scatter that should be close
% to a straight line if the data look like a random sample from the theoretical
% distribution. This function chooses the normal as the theory, to provide a
% graphical/visual assessment of how normal the data appear.

% Look at the S-Plus functions qqplot() and qqnorm()
% To help with assessing the relevance of sampling variability
% on just "how close" to the normal the data appears, we add (very) approximate
% posterior 95% intervals for the uncertain quantile function at each point
% (Based on approximate theory)

% qqbayes below is an matlab function to do this
  
% You can use the function with any data. e.g., to see how a random sample of size 250
% from a Beta(5,3) distribution compares with a normal (a normal with the same mean and
% variance as the Beta sample), just enter
%                                 qqbayes( betarnd(5,3,1,250) )
% Another eg:
%                                 qqbayes( trnd(4,1,100) )
% Use
%                                 qqbayes( normrnd(0,1,1,100) )
% several times for a specific sample size n -- to get a "feel" for how much
% sampling variability there is in a sample of size n that DOES come from a normal.
% Repeat with different values of n, e.g., n=30, n=100, n=500, etc%
%

function qqbayes(mydata)

   n=length(mydata);
   p=(1:n)/(n+1);
   %% x=mydata; 
    x=(mydata-mean(mydata))/sqrt(var(mydata));
   x=sort(x);
   z=qnorm(p,0,1);

   plot(z,x,'r+'); axis square
   hold on	
   xlabel('Standard Normal Quantiles','fontsize',12);
   ylabel('Ordered Residuals','fontsize',12);

   s=2.5758*sqrt(p.*(1-p)/n);
   pl=p-s; i=pl<1&pl>0;
   plot(z(i),prctile(x,100.*pl(i)),'b--');
   pl=p+s; i=pl<1&pl>0;
   plot(z(i),prctile(x,100.*pl(i)),'b--');

   plot([-3,3], [-3,3],'k--');
   plot(z,x,'r+');

   hold off; box off
   
   legend('residuals',' ~ 99% intervals','location','northwest')
   legend boxoff
 %  title('Normal quantile plot with approx 95% posterior intervals');


