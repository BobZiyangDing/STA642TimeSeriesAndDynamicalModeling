function [waves,mods,decomp] = ar_decomp(x,p,ncomp)
% Fit a reference posterior AR(p) model to the time series x
% and produce the decomposition into latent components 
% Outputs: 
%   waves  -- the estimated wavelengths of ARMA(2,1) components
%   mods   -- the estimated moduli of all components, the first 
%   		ones corresponding to complex roots in order in waves	
%   decomp  -- the matrix of estimated components, in same order from top town, 
%               starting at t=p and running to t=length(x) 
%
% organise data ...
 arx=x-mean(x); mi=min(arx); ma=max(arx); ndat=length(arx);
 y=arx(ndat:-1:p+1); n=length(y); X=hankel(arx(ndat-1:-1:p),arx(p:-1:1));
% compute AR fit ...
 xtx=X'*X; b=xtx\(X'*y); r=y-X*b; nu=n-p; s=r'*r/nu;
% now decomposition ... 
 G=[b,[diag(ones(p-1,1)),zeros(p-1,1)]']';
 [e,l]=eig(G); l=diag(l); arg=angle(l); mods=abs(l); H=diag(e(1,:))*inv(e);
%
% reorder in terms of increasing arguments and remove neg args
% nc is # complex,  nk is total, num real is nk-nc
 arg(arg==0)=pi;
 [arg,i]=sort(arg);   mods=mods(i); e=e(i,:); H=H(i,:);
 nc=sum(arg<0); 
 i=arg>0; arg=arg(i); mods=mods(i); e=e(i,:); H=H(i,:);
 nk=sum(i);  
 H(1:nc,:)=2*real(H(1:nc,:)); i=(nc+1):nk; H(i,:)=real(H(i,:));
% 
 decomp=[zeros(p-1,nk)        
        fliplr([x(ndat-(0:(p-1)))'
                X]*H')];
 decomp=fliplr(H*[arx(ndat-(0:(p-1))) X' zeros(p,p-1)]);
 waves =2*pi./arg;
%

