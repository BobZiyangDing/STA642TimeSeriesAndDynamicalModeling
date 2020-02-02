function y=acf(x,plott,k)
% y=acf(x,plot,k) gives autocorrelation function 
% if third argument is given, plot only up to lag k
% if second argument is given, plot is drawn and appr 95% limits

y=acovf(x); y=y/y(1); 

if nargin>1
   n=length(y); 
   if nargin==2
       k=n; 
   end
   plot([0:(k-1);0:(k-1)],[zeros(1,k);y(1:k)],'k'); line([0 k],[0 0],'color','k')
   hold on;
   plot([0 k],[1.96/sqrt(n) 1.96/sqrt(n)],':b',[0 k],-ones(1,2)*1.96/sqrt(n),':b');
   hold off;
   xlabel('Lag'); ylabel('Autocorrelation')
end

 