function [waves,mods,decomp,maxnr,maxnc] = tvar_decomp(x,m)
% Decomposition of TVAR model based on post means in matrix m
%
% ** ignores issue of changing numbers of complex and real comps! **
%
% Inputs:
%   x --  series of length T
%   m --  pxT array of post means for states
% 
% Outputs: 
%   waves  -- floor(p/2)xT matrix of wavelengths of ARMA(2,1) components
%   mods   -- pxT matrix of moduli of all comps
%   decomp -- pxT trajectories of all components 
%   maxnr  -- max number of real comps
%   maxnc  -- max number of complex comps


[p,T]=size(m); arx=reshape(x,T,1); 
decomp=zeros(p,p); waves=zeros(p,p); mods=zeros(p,p);
G=[m(:,1),[diag(ones(p-1,1)),zeros(p-1,1)]']';
maxnr=0; maxnc=0;

% decomposition ... 
for t=(p+1):T
 G(1,:)=m(:,t)'; [e,l]=eig(G); l=diag(l); arg=angle(l); mod=abs(l);
 H=diag(e(1,:))*inv(e);
%
% reorder roots in terms of 
%  - decreasing wavelengths (for complex) and remove neg args
%  - followed by reals in terms of decreasing moduli
arg(arg==0)=pi; [arg,i]=sort(arg);  mod=mod(i); H=H(i,:);
i=find(arg<0); nc=length(i); nr=p-2*nc; nk=nc+nr;
% now remove negative mods
    arg(i)=[]; mod(i)=[]; H(i,:)=[]; 
 H=real(H);  H(1:nc,:)=2*H(1:nc,:);
%
 decomp(1:nk,t)=H(1:nk,:)*arx((t-1):-1:(t-p));
 waves(1:nc,t) =2*pi./arg(1:nc);
 mods(1:nk,t)  =mod(1:nk);
 maxnr=max(nr,maxnr);  maxnc=max(nc,maxnc); 
end;
 % now ad-hoc treatment of 1st p values 
 waves(:,1:p)=repmat(waves(:,p+1),1,p);
 mods(:,1:p) =repmat(mods(:,p+1),1,p);
 % and trim down the output arrays
 waves=waves(1:maxnc,:);
 mods=mods(1:min(p,(maxnc+maxnr)),:);
 decomp=decomp(1:min(p,(maxnc+maxnr)),:);
%

