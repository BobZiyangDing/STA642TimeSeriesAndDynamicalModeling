function X=rMatNorm(M,U,V);
% rMatNorm(M,U,V)  generate a r.q random matrix normal X 
% M = r.q matrix mean
% U = r.r variance matrix of each column 
% V = q.q variance matrix of each row 

[r q]=size(M);
C=chol(U); D=chol(V);  
X = M + C'*randn(r,q)*D; 

% X=reshape(rMNorm(reshape(M,r*q,1),kron(V,U),1),r,q);  % direct, bad way!