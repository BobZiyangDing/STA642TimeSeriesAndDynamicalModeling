function [m,C,n,s,e,mf,Cf,sf,nf,ef,qf] = tvar(x,p,del,m0,C0,s0,n0)
% Fit TVAR(p) model to univariate series x
% Discount factors del(1) for state and del(2) for obs var
% Inputs:
%   x --  series of length T
%   p --  model order
%  m0 --  px1 vector prior mean for state
%  C0 --  pxp prior var matrix
%  n0 --  prior df
%  s0 --  prior estimate of obs var
% Outputs: 
%  forward filtered: 
%  mf  --  pxT   post mean vectors
%  Cf  --  pxpxT post var matrices
%  nf  --  T post dfs
%  sf  --  T post obs var estimates
%  ef  --  1-step forecast errors (zero to t=p) 
%  qf  --  1-step forecast variance factors
%  retrospectively updated/smoothed: 
%  m   --  pxT   post mean vectors
%  C   --  pxpxT post var matrices
%  n   --  T post dfs
%  s   --  T post obs var estimates
%  e   --  estimated innovations (zero to t=p) 
% 
% organise data ...
 d=del(1); b=del(2);
 arx=x-mean(x);
 T=length(arx);
 arx=reshape(arx,T,1); 
 m=repmat(m0,[1 p]);  C=repmat(C0,[1 1 p]);  
 s=repmat(s0,[1 p]);  n=repmat(n0,[1 p]); 
 e=zeros(T,1);  q=zeros(T,1);  
 mt=m0; Ct=C0; st=s0; nt=n0;

% forward filtering
 for t=(p+1):T
     F=arx((t-1):-1:(t-p));
     A=Ct*F/d; qt=F'*A+st; A=A/qt; et=arx(t)-F'*mt; e(t)=et; q(t)=qt; 
     mt=mt+A*et; m(:,t)=mt; 
     r=b*nt+et*et/qt; nt=b*nt+1; r=r/nt; st=st*r;
     n(t)=nt; s(t)=st;   
     Ct=r*(Ct/d-A*A'*qt); Ct=(Ct+Ct')/2;
     C(:,:,t)=Ct;
 end;
% save filtered values ...
mf=m; Cf=C; sf=s; nf=n; ef=e; qf=q; 

% backward smoothing
 for t=(T-1):-1:1
     m(:,t)=(1-d)*m(:,t)+d*m(:,t+1);
     if t>p
        e(t)=arx(t)-m(:,t)'*arx((t-1):-1:(t-p));
     end;
     n(t)  =(1-b)*n(t)+b*n(t+1);  
     st=s(t); s(t)=1/((1-b)/st+b/s(t+1)); 
     C(:,:,t)=s(t)*((1-d)*C(:,:,t)/st + d*d*C(:,:,t+1)/s(t+1)); 
 end;
 % now ad-hoc treatment of 1st p values 
 m(:,1:p)=repmat(m(:,p+1),1,p);  C(:,:,1:p)=repmat(C(:,:,p+1),[1 1 p]);
 n(1:p)=repmat(n(p+1),1,p); s(1:p)=repmat(s(p+1),1,p);
%



