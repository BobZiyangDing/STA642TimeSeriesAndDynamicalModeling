function [K] = wishart_InvA_rndInvA(df,S,mc);
% Generates mc draws from a Wishart - allows for singular Wishart too, in
%    cases that the distn parameter S is rank deficient of rank r
% 
% K is W_p(df,A) with sum-of-squares parameter S = A^{-1} and d.o.f. df
% Dimension p is implicit
% Usual nonsingular case: r=p<df; now allow for r<p  with r integral
% Note that  E(K)=df.S^{-1}  in this notation 
%
% Nonsing  pdf is p(K) = cons. |K|^((df-p-1)/2) exp(-trace(K S)/2) 
% Sing case usual modification 
%
% Returns matrix W of dimension (mc,p,p) of rank r=rank(S)
%
%  EXAMPLE:  reference posterior in a normal model N(0,Sigma) with precision 
%    matrix Omega = Sigma^{-1} 
%    Random sample of size n has sample var matrix V=S/n with S=\sum_i x_ix_i'
%    Ref posterior for precision mx Omega is W_p(n,A) with A=S^{-1}
%    e.g.,  K=wishart_InvA_rnd(n,n*V,1000); 
%  Useful for looking at posteriors for correlations and eigenvalues
%  and also for posterior on graph - looking for elements of Omega near 0
%
%  To deal with numerical problems when p is larger, or high collinearity
%  and/or rank deficiency, use this: 
%
p=size(S,1); 
[P,D]=svd0(S); i=find(D>max(D)*1e-9); r=length(i);
               P=rprod(P(:,i), 1./sqrt(D(i))); 
%OR: P=chol(inv(S))'; r=p;
% This means that PP'=A=inv(S)
%
% Now generate W_p(df,I) using Bartlett decomposition: 
for imc=1:mc
    U=zeros(p,r);
    for i=1:r
        U(i,i:r)=[sqrt(gamrnd(((df-i+1)/2),2,1,1)) randn(1,r-i)];
    end;
    U=U*P'; 
    K(imc,:,:)=U'*U;  
end
%
K=squeeze(K); 
