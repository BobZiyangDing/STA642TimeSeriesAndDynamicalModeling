 

%  Initial priors for DDNM of 3 series macro example, monthly data 
% Tbill:  
j=1;    m0 = zeros(p(j)+1,1); m0(1+npa(j)+1)=0; m0(4)=0.9; % 0.9 is prior est of AR(1) coeff
        C0 = eye(p(j)+1);   C0(1,1)=0.1; % constrained prior on intercept
        s0 = .1; n0=2;  delta=0.99; beta=0.99;          
    prior{1,j}=m0; prior{2,j}=C0; prior{3,j}=n0; prior{4,j}=s0; 
    prior{5,j}=delta; prior{6,j}=beta;       
o1T=ones(1,T);  lmlik=zeros(q,T); 
%Infln:  
j=2;    m0 = zeros(p(j)+1,1); m0(1+npa(j)+1)=0; m0(3)=0.9; % 0.9 is prior est of AR(1) coeff
        C0 = eye(p(j)+1);   C0(1,1)=0.1; % constrained prior on intercept
        s0 = .1; n0=2;  delta=0.99; beta=0.99;  
    prior{1,j}=m0; prior{2,j}=C0; prior{3,j}=n0; prior{4,j}=s0; 
    prior{5,j}=delta; prior{6,j}=beta;     
% Unemp:  
j=3;    m0 = zeros(p(j)+1,1); m0(1+npa(j)+1)=0; m0(3)=0.9; % 0.9 is prior est of AR(1) coeff
        C0 = eye(p(j)+1);  C0(1,1)=0.1; % constrained prior on intercept
        s0 = .1; n0=2;  delta=0.99; beta=0.99;  
    prior{1,j}=m0; prior{2,j}=C0; prior{3,j}=n0; prior{4,j}=s0; 
    prior{5,j}=delta; prior{6,j}=beta;            
