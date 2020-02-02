
% seasonal factor estimates and uncertainties:
startup


for i=1:6, figure(i); clf; end; figure(1)
fprintf('Using figure 1 only ...\n\n'); pause

q=12; 
 mseas = [ 0.8 1.3 1 1.1 1.2 0 -1.2 -1 -1 -1.5 -0.5 0 ]';   
 Cseas = 0.04*eye(q); 
 
  plot(mseas); hold on;  errorbar(1:q,mseas,2*sqrt(diag(Cseas)),'r*'); hold off; box off
        xlabel('Months (1=Jan)'); title('Unconstrained seasonals')

% constrained prior - must sum to zero as they are seaonal deviations 
%   from nonseasonal trend+regression:
 z=ones(q,1); a=Cseas*z; r=z'*a; a=a/r;  
 mseas = mseas-a*sum(mseas); Cseas = Cseas-r*a*a'; fmseas=mseas;
  plot(mseas); hold on;  errorbar(1:q,mseas,2*sqrt(diag(Cseas)),'r*'); hold off; box off
        xlabel('Months j=1:12 (1=Jan)'); title('Constrained seasonals \phi_{1,j}, j=1:12'); 
        ylabel('\phi_{1,j}','rotation',0)
        
% convert to prior on Fourier components 
 p=12; rseas=[1 3 4]; pseas=length(rseas); nseas=2*pseas;  
 Fseas=repmat([1 0],1,pseas)'; Gseas = zeros(nseas,nseas); 
    for j=1:pseas
        c=cos(2*pi*rseas(j)/p); s=sin(2*pi*rseas(j)/p); i=2*j-1:2*j; 
        Gseas(i,i)=[[c s];[-s c]]; end
 L = zeros(q,nseas);  L(1,:)=Fseas'; 
 for t=2:q, L(t,:)=L(t-1,:)*Gseas; end
 H=(L'*L)\L'; 
 
 % resulting prior parameters: 
 mseas = H*mseas;  Cseas = H*Cseas*H'; 
 fprintf('Initial prior mean and variance matrix of seasonal states:\n')
 [mseas Cseas] 
 
  
 T=1:q; f=zeros(1,q); a= mseas; 
 for t=1:q
        if (t>1) a = Gseas*a; end
        f(t)=Fseas'*a;
 end

 
 