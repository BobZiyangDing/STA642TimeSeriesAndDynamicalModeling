function pr=pacf(x,plfg,m)
%
% pr=pacf(x,m,plfg)
%
% computes the pacf of the series in x
%
% x = input series
% m = maximum lag to compute
% pa = pacf(j) j=1:m

% first compute all the sample correlations needed
nx=length(x);
%rho=acf(x,3);  % divide by sample var * no samples in kth product
rho=acf(x);  % divide by sample var * tot samples
% r=xcov(x,'unbiased'); % using xcov has 2nx-1 el
% sigsq=r(nx);% the middle one is the variance
% rho=r(nx:nx+m)/sigsq;% lag 0 to end
pr(1)=rho(2); % lag 1 rho
%pr=zeros(1,m);
for k=2:m
   pmat=toeplitz(rho(1:k));% lag 0 to k-1
   rhovec=rho(2:k+1);		% lag 1 to k
   phi=inv(pmat)*rhovec';
   pr(k)=phi(k);
end
if plfg
   ALPHA=.05; 
   plot([1:m;1:m],[zeros(1,m);pr(1:m)],'k'); line([0 m],[0 0],'color','k')
   hold on;
   thr=norminv(1-ALPHA/2,0,1/sqrt(nx));% later replace with rinv
   plot([0:m],thr*ones(1+size(pr)),':b',[0:m],-thr*ones(1+size(pr)),':b');
   %title(['PACF, nx=',int2str(nx),' alpha = ',num2str(ALPHA)]);
   hold off; xlabel('Lag'); ylabel('Partial autocorrelation'); box off
end

      
   