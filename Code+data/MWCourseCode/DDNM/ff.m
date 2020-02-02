function p1T = ff(y,F,T,t0,p0)
%   Forward filter in univariate dynamic regression 
%    y = 1xT series, F=pxT predictors
%   p0 = structure with prior info  m0, C0, n0, d0, delta, beta 
%      = initial prior at time t0+1
%  p1T = cell array with forward updated m1:T, C1:T, n1:T, d1:T
%          and log of 1-step predictive pdfs for marg lik, lp1:T
%  Perform filtering starting at time t0+1 
p=size(F,1);
p1T.a=zeros(p,T); p1T.R=zeros(p,p,T); 
p1T.m=zeros(p,T); p1T.C=zeros(p,p,T); p1T.n=zeros(1,T); p1T.d=zeros(1,T); 
p1T.f=zeros(1,T); p1T.q=zeros(1,T); p1T.lp=zeros(1,T); 

mt = p0.m; Ct=p0.C; nt=p0.n; st=p0.s; dt=nt*st; delta=p0.delta; beta=p0.beta;         


% forward filtering: z
        for t = t0+1:T
            Ft=F(:,t); ft = Ft'*mt; et = y(t) - ft; stm1=st; 
            Rt = Ct/delta; At = Rt*Ft; qt = st+Ft'*At; 
            p1T.a(:,t)=mt; p1T.R(:,:,t)=Rt; nt=beta*nt;
            p1T.f(t)=ft; p1T.q(t)=qt; p1T.lp(t) = lsttpdf(et,nt,0,qt); 
            At=At/qt; 
            nt = nt+1; dt = beta*dt+stm1*et*et/qt;  st=dt/nt;    
            mt = mt + At*et; Ct = (Rt - At*At'*qt)*(st/stm1);  Ct=(Ct+Ct')/2; 
            p1T.m(:,t)=mt; p1T.C(:,:,t)=Ct; p1T.n(t)=nt; p1T.s(t)=st; 
        end

 % initial values frozen: 
        for t = 1:t0
            p1T.a(:,t)=p1T.a(:,t0+1); p1T.R(:,:,t)=p1T.R(:,:,t0+1); 
            p1T.f(t)=y(t); p1T.q(t)=p1T.q(t0+1); p1T.lp(t) = 0; 
            p1T.m(:,t)=p1T.m(:,t0+1); p1T.C(:,:,t)=p1T.C(:,:,t0+1); 
            p1T.n(t)=p1T.n(t0+1); p1T.s(t)=p1T.s(t0+1); 
        end
        