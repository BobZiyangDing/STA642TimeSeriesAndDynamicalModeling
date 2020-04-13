function pdf = mvnormixpdf(x,w,mu,Sigma)
% x  = p.n array of n values of p-dim MV normal mixture 
% w = column k vector probs
% mu = p.k matrix of component means 
% Sigma = p.p.k array of variance matrice 
% pdf = n vector of pdf values
%
[p n]=size(x); k=length(w); pdf=zeros(1,n);
for j=1:k
    C=chol(Sigma(:,:,j)); e=inv(C)'*(x-repmat(mu(:,j),1,n)); 
    if (p==1), q=e.*e; else, q=sum(e.*e); end
    pdf = pdf + w(j)*exp(-q/2)/( prod(diag(C))*(2*pi)^(p/2) );
end
