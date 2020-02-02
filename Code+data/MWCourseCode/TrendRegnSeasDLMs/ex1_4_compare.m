
startup

for i=1:6, figure(i); clf; end; figure(1)
fprintf('Using figure 1 only ...\n\n'); pause


load salesindex
Y=sales;  deltrend=0.9; delregn=0.98; delseas=0.95; 
 
% compare nm models using marginal likelihood ... 
nm=6; mlik=zeros(nm,T); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  M1: linear term, 3 seasonal components 
im=1; display('model 1: locally linear term, 3 seasonal components')
Ftrend=[1 0]'; Gtrend=[[1 1];[0 1]];  ntrend=length(Ftrend); itrend=1:ntrend;          
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1 3 4]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 

% priors 
a=[9.5 0 -0.7 0.691 1.159 0.283 -0.050 -0.217 0.144]';
R=100*blkdiag(0.09,0.01,0.01,0.0067*eye(nseas)); nu=6; S=0.15^2; 

% forward filtering 
for t=1:T
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         R(iseas,iseas)=R(iseas,iseas)/delseas;
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    y = Y(:,t);  e=y-f; 
    rQ=sqrt(Q); mlik(im,t) = tpdf(e/rQ,nu)/rQ;
    r=(nu+e^2/Q)/(nu+1);   m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;   
end
figure(1); clf; 
    subplot(2,1,1); hold on
    plot(1:T,mlik); title('Marginal likelihood'); ylabel('marg lik'); eval(xa)
    
pause    
display('hit a key to see model 2: linear term, 2 seasonal components')
pause
  
im=2; 
Ftrend=[1 0]'; Gtrend=[[1 1];[0 1]];  ntrend=length(Ftrend); itrend=1:ntrend;          
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1 3]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 

% priors 
a=[9.5 0 -0.7 0.691 1.159 0.283 -0.050]';
R=100*blkdiag(0.09,0.01,0.01,0.0067*eye(nseas)); nu=6; S=0.15^2; 

% forward filtering 
for t=1:T
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         R(iseas,iseas)=R(iseas,iseas)/delseas;
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    y = Y(:,t);  e=y-f; 
    rQ=sqrt(Q); mlik(im,t) = tpdf(e/rQ,nu)/rQ;
    r=(nu+e^2/Q)/(nu+1);   m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;   
end
figure(1); clf; 
    subplot(2,1,1); 
    plot(1:T,mlik); title('Marginal likelihood'); ylabel('marg lik'); eval(xa)
 
pause
display('hit a key to see model 3: linear term, 1 seasonal component')
pause
  


im=3; 
Ftrend=[1 0]'; Gtrend=[[1 1];[0 1]];  ntrend=length(Ftrend); itrend=1:ntrend;          
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 

% priors 
a=[9.5 0 -0.7 0.691 1.159]';
R=100*blkdiag(0.09,0.01,0.01,0.0067*eye(nseas)); nu=6; S=0.15^2; 

% forward filtering 
for t=1:T
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         R(iseas,iseas)=R(iseas,iseas)/delseas;
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    y = Y(:,t);  e=y-f; 
    rQ=sqrt(Q); mlik(im,t) = tpdf(e/rQ,nu)/rQ;
    r=(nu+e^2/Q)/(nu+1);   m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;   
end
figure(1); clf; 
    subplot(2,1,1); 
    plot(1:T,mlik); title('Marginal likelihood'); ylabel('marg lik'); eval(xa)
    
pause
display('hit a key to see model 4: constant term, 3 seasonal components')
pause
  
im=4; 
Ftrend=[1]'; Gtrend=[1];  ntrend=length(Ftrend); itrend=1:ntrend;          
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1 3 4]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 

% priors 
a=[9.5 -0.7 0.691 1.159 0.283 -0.050 -0.217 0.144]';
R=100*blkdiag(0.09,0.01,0.0067*eye(nseas)); nu=6; S=0.15^2; 

% forward filtering 
for t=1:T
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         R(iseas,iseas)=R(iseas,iseas)/delseas;
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    y = Y(:,t);  e=y-f; 
    rQ=sqrt(Q); mlik(im,t) = tpdf(e/rQ,nu)/rQ;
    r=(nu+e^2/Q)/(nu+1);   m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;   
end
figure(1); clf; 
    subplot(2,1,1); hold on
    plot(1:T,mlik); title('Marginal likelihood'); ylabel('marg lik'); eval(xa)
    
pause
display('hit a key to see model 5: constant term, 2 seasonal components')
pause
  

im=5; 
Ftrend=[1]'; Gtrend=[1];  ntrend=length(Ftrend); itrend=1:ntrend;          
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1 3]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 

% priors 
a=[9.5 -0.7 0.691 1.159 0.283 -0.050]';
R=100*blkdiag(0.09,0.01,0.0067*eye(nseas)); nu=6; S=0.15^2; 

% forward filtering 
for t=1:T
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         R(iseas,iseas)=R(iseas,iseas)/delseas;
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    y = Y(:,t);  e=y-f; 
    rQ=sqrt(Q); mlik(im,t) = tpdf(e/rQ,nu)/rQ;
    r=(nu+e^2/Q)/(nu+1);   m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;   
end
figure(1); clf; 
    subplot(2,1,1); hold on
    plot(1:T,mlik); title('Marginal likelihood'); ylabel('marg lik'); eval(xa)
    

pause
display('hit a key to see model 6: constant term, 1 seasonal component')
pause
  

im=6; 
Ftrend=[1]'; Gtrend=[1];  ntrend=length(Ftrend); itrend=1:ntrend;          
q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier components of seasonal:
p=12; rseas=[1]; pseas=length(rseas); nseas=2*pseas; iseas=ntrend+nregn+1:ntrend+nregn+nseas; 
Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 

% priors 
a=[9.5 -0.7 0.691 1.159]';
R=100*blkdiag(0.09,0.01,0.0067*eye(nseas)); nu=6; S=0.15^2; 

% forward filtering 
for t=1:T
    if (t>1) a = G*m; 
             R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                         R(iregn,iregn)=R(iregn,iregn)/delregn; 
                         R(iseas,iseas)=R(iseas,iseas)/delseas;
    end
    Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
    y = Y(:,t);  e=y-f; 
    rQ=sqrt(Q); mlik(im,t) = tpdf(e/rQ,nu)/rQ;
    r=(nu+e^2/Q)/(nu+1);   m = a+A*e; C = r*(R-A*A'*Q); nu=nu+1; S=r*S;   
end
figure(1); clf; 
    subplot(2,1,1); hold on
    plot(1:T,mlik); title('Marginal likelihood'); ylabel('marg lik'); eval(xa)
    

    
%%%%    
display('hit a key to compare models')
pause

llik=log(mlik); subplot(2,1,2);
    plot(1:T,llik); title('Log marginal likelihood'); ylabel('log marg lik'); eval(xa)
    
    pause
    
    llik=cumsum(llik')'; 
    hold off
    plot(1:T,llik); title('Cumulative log marginal likelihood'); ylabel('cum log marg lik')
    eval(xa)
    
    pause
    
    prob=exp(llik-repmat(max(llik),nm,1)); prob=prob./repmat(sum(prob),nm,1);
    hold off
    plot(1:T,prob); title('Cumulative model probabilities'); ylabel('probability'); eval(xa)
    legend([repmat('M',nm,1),int2str((1:nm)')],'location','east'); legend boxoff
   

