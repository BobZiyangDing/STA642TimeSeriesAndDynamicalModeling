function [freqs,phis] = tvar_sim(m,C,n,s,times,k,N)
% Sample posteriors in TVAR model
% 
% Inputs:
%  m  --  pxT array of post means for states
%  C  --  pxpxT post var matrices
%  n  --  T post dfs
%  s  --  T post obs var estimate
%  times  --  vector of times to sample, length nt 
%  k  --  max number of freqs to sample
%  N  --  MC sample size
%
% Output: 
%   phis   -- p x nt x N array of post sampled \phi_t
%   freqs  -- k x nt x N array of post sampled TVAR comp frequencies
%   
%   
nt=length(times); [p,T]=size(m);
k=min(k,floor(p/2)); freqs =zeros(k,nt,N); phis=zeros(p,nt,N);

% sample frequencies ...
for it=1:nt
 t=times(it);
 theta=repmat(m(:,t),1,N)+(chol(C(:,:,t))'*randn(p,N))...
         ./repmat(sqrt(gamrnd(n(t)/2,2/n(t),1,N)), p,1);
 phis(:,it,:)=theta; 
 for mc=1:N
   arg=angle(roots([1 -theta(:,mc)']));
   arg(arg==0)=pi; [arg,i]=sort(arg); nc=sum(arg<0); i=arg>0; arg=arg(i); 
   freqs(:,it,mc) =arg(1:k);
 end;
end;
%
