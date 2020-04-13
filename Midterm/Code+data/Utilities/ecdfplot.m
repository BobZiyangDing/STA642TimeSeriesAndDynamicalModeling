
function [pz,py]=ecdfplot(x,iplot,sur,lc)
%
% computes ECDF of a set of k samples in columns of matrix X 
% and delivers both x values evaluated (pz), and ecdf values (py), at these 
% Plots on same frame if iplot=1
% Plots survival functions instead of CDFs if sur=1
% Uses colors if lc  is an input; if not, uses standard cols

if (nargin==2) 
   sur=0; lc='krbgmcykrgbmc'; 
end
if (nargin==3)
   lc='krbgmcykrgbmc';
end

[n,r]=size(x); d=1/(n+1); py=[]; pz=[];

if (iplot), clf; end
for i=1:r
    z=sort(x(:,i)); 
    z=reshape(repmat(z,1,2)',1,2*n)'; 
    y=d*reshape(repmat((0:n)',1,2)',1,2*(n+1))'; y([1 end])=[]; 
	py=[py y]; pz=[pz z]; 
    if (iplot) 
        plot(z,(1-sur)*y+sur*(1-y),'Color',lc(i)); 
        axis([min(z) max(z) 0 1])
        hold on
    end	
end
if (iplot), hold off; end
py=unique(py,'rows'); py(1,:)=[]; pz=unique(pz,'rows'); 
% 


