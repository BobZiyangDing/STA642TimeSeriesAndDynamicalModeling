function [swaves,smods,sdecomp] = ar_decomp_sim(x,p,nmc)
% Fit a reference posterior AR(p) model to the time series x
% and produce posterior samples of the decomposition into latent components 
% Inputs:
% x     -- ndat dimensional column vector time series, centered
% p     -- AR model order
% nmc   -- number of Monte Carlo samples 
%
% Outputs: 
%   swaves  -- the sampled wavelengths of quasi-periodic ARMA(2,1) components
%   smods   --  sampled moduli of all components, the first 
%      	        	ones corresponding to complex roots in order in swaves	
%   sdecomp  -- p.nmc matrix of estimated components, in same order from top town, 
%               starting at t=p and running to t=length(x) 
%
% organise data ...
 arx=x-mean(x); mi=min(arx); ma=max(arx); ndat=length(arx);
 y=arx(ndat:-1:p+1); n=length(y); X=hankel(arx(ndat-1:-1:p),arx(p:-1:1));
% compute AR fit ...
 B=X'*X; b=B\(X'*y); nu=n-p; r=y-X*b; nu=n-p; v=r'*r/nu;
% simulate posterior ...
 [~,bsamp,~] = arsim(x,p,b,B,nu,v,nmc,0);
% now decomposition ... repreat for each sampled phi vector 
swaves=zeros(p,nmc); smods=zeros(p,nmc); sdecomp=NaN(ndat,p,nmc); 
ncompmax=0; 
for imc=1:nmc
    G=[bsamp(:,imc),[diag(ones(p-1,1)),zeros(p-1,1)]']';
    [e,l]=eig(G); l=diag(l); arg=angle(l); mods=abs(l); H=diag(e(1,:))*inv(e);
    %
    % reorder in terms of increasing arguments and remove neg args
    % nc is # complex,  nk is total, num real is nk-nc
    arg(arg==0)=pi;
    [arg,i]=sort(arg);   mods=mods(i); e=e(i,:); H=H(i,:);
    nc=sum(arg<0);
    i=arg>0; arg=arg(i); mods=mods(i); e=e(i,:); H=H(i,:);
    nk=sum(i); ncompmax = max(nk,ncompmax); 
    H(1:nc,:)=2*real(H(1:nc,:)); i=(nc+1):nk; H(i,:)=real(H(i,:));
    %
    decomp=[zeros(p-1,nk)
        fliplr([x(ndat-(0:(p-1)))'
        X]*H')];
    decomp=fliplr(H*[arx(ndat-(0:(p-1))) X' zeros(p,p-1)]);
    waves =2*pi./arg;
    swaves(1:nk,imc)=waves; smods(1:nk,imc)=mods; sdecomp(:,1:nk,imc)=decomp'; 
end
%
swaves=swaves(1:ncompmax,:); smods=smods(1:ncompmax,:); 
sdecomp=sdecomp(:,1:ncompmax,:); 




