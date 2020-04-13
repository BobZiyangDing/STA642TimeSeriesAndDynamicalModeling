startup

for i=1:2, figure(i); clf; end; figure(1); figure(2) 

fprintf('Using figures 1 and 2 ...\n\n'); pause


%  uncertain period of a single Fourier component DLM ...            
%   in sales/index example               
% load data set and plot 
load salesindex
    figure(1); clf; subplot(3,1,1) 
    plot(1:T,sales,1:T,index); eval(xa)
    xlabel('Date (month/year)'); box off; xlim([0 T+1]); ylim([-3 12]);
    legend('Sales','Index','location','east'); legend boxoff
	fprintf('Hit return to compute marginal likelihood ...\n') 
    pause
    
% components of a DLM:
% trend: % e.g. locally constant:  
    Ftrend=1; Gtrend=1; ntrend=length(Ftrend); itrend=1:ntrend;          
% regn:
    q=1; X=index; Fregn=zeros(1,q); Gregn = eye(q); 
    nregn=q; iregn=ntrend+1:ntrend+q; 
% Fourier component:
    pseas=1; nseas=2; iseas=ntrend+nregn+1:ntrend+nregn+nseas;    
    Fseas=[1 0]'; Gseas = zeros(2,2); 
% DLM matrices F & G: 
    F = [Ftrend;Fregn;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gregn,Gseas); 
% response time series: 
    Y=sales;
% 3 discount factors for component DLMs: 
    deltrend=0.9; delregn=0.98; delseas=0.95; 

% range of wavelengths to explore: 
lam = 8:.1:16; nlam=length(lam);  logmlik=zeros(nlam,1); 
%
for ilam=1:nlam
    om = 2*pi/lam(ilam); G(iseas,iseas)=[[ cos(om) sin(om) ]; [ -sin(om) cos(om)]]; 
    % priors 
    nu=6; S=0.15^2; a=[9.5 -0.7 0 0]'; R=blkdiag(0.09,0.01,0.01*eye(nseas)); 
    for t=1:T
        if (t>1) a = G*m; 
                 R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                             R(iregn,iregn)=R(iregn,iregn)/delregn; 
                             R(iseas,iseas)=R(iseas,iseas)/delseas;
        end
        Fregn=X(:,t)'; F(iregn)=Fregn; A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
        y = Y(:,t);  e=y-f; rQ=sqrt(Q); 
        logmlik(ilam) = logmlik(ilam) + log(tpdf(e/rQ,nu)/rQ);
        r=(nu+e^2/Q)/(nu+1); m = a+A*e; C = r*(R-A*A'*Q);  nu=nu+1; S=r*S;    
        C=(C+C')/2;    % numerical check
    end
end


   subplot(3,1,2); plot(lam,logmlik); text(9,-64,'Log marginal likelihood for period \lambda'); box off
        xlabel('Period \lambda=2\pi/\omega'); ylabel('Log likelihood'); ylim([-100 -60])
   subplot(3,1,3); plot(lam,exp(logmlik-max(logmlik))); 
   text(9,1.1,'Marginal likelihood for period \lambda'); box off
        xlabel('Period \lambda=2\pi/\omega'); ylabel('Likelihood'); ylim([0 1.2])
   [a,i]=max(logmlik); lam(i) 
   
   
       
figure(2); clf
display('Next comes O2 example ...')
pause
clear
Ti=1:866; T=length(Ti); x=load('o2.dat'); Y=x(Ti,2); 
   subplot(3,1,1); plot(3*(1:T),Y); xlabel('kyears'); xlim([0 3*(T+1)]); box off
   text(840,5.7,'Oxygen isotope ratio series') 
  	fprintf('Hit return to compute marginal likelihood ...\n') 
    pause

   
% trend: % e.g. locally constant:  
    Ftrend=1; Gtrend=1; ntrend=length(Ftrend); itrend=1:ntrend;          
% Fourier component:
    pseas=1; nseas=2; iseas=ntrend+1:ntrend+nseas;    
    Fseas=[1 0]'; Gseas = zeros(2,2); 
% DLM matrices F & G: 
    F = [Ftrend;Fseas]; n=length(F);  G = blkdiag(Gtrend,Gseas); 
% 2 discount factors for component DLMs: 
    deltrend=1.0; delseas=1.00; 

% range of wavelengths to explore: 
lam = (10:1:150)/3; nlam=length(lam);  logmlik=zeros(nlam,1); 
%
for ilam=1:nlam
    om = 2*pi/lam(ilam); G(iseas,iseas)=[[ cos(om) sin(om) ]; [ -sin(om) cos(om)]]; 
    nu=2; S=1; a=[4.5 0 0]'; R=1000*eye(n); 
    for t=1:T
        if (t>1) a = G*m; 
                 R = G*C*G'; R(itrend,itrend)=R(itrend,itrend)/deltrend;  
                             R(iseas,iseas)=R(iseas,iseas)/delseas;
        end
        
        if any(diag(R)<0), display('uhoh'); pause; end
        
        A = R*F; Q=F'*A+S; A=A/Q; f=F'*a;
        y = Y(t);  e=y-f;   rQ=sqrt(Q); 
        logmlik(ilam) = logmlik(ilam) + log(tpdf(e/rQ,nu)/rQ);
        r=(nu+e^2/Q)/(nu+1); m = a+A*e; C = r*(R-A*A'*Q);  nu=nu+1; S=r*S;  
        C=(C+C')/2;    % numerical check
    end
end

   subplot(3,1,2); plot(3*lam,logmlik); text(50,-120,'Log marginal likelihood for period \lambda'); box off
        xlabel('Period \lambda=2\pi/\omega in kyears'); ylabel('Log likelihood')
   subplot(3,1,3); plot(3*lam,exp(logmlik-max(logmlik))); 
    text(50,0.88,'Marginal likelihood for period \lambda'); box off
        xlabel('Period \lambda=2\pi/\omega in kyears'); ylabel('Likelihood')
   [a,i]=max(logmlik); 3*lam(i) 
   
   
    

   
        