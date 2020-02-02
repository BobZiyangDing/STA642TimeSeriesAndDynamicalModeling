function f = Eforecast(k,M,q,arp,p,Yt)
% Expectation-only  FORECAST k steps ahead in TV-VAR(arp) + intercept model in q dimensions
% p = dimension of model paramater = 1+arp*q
% Y = data on current and last arp-1 times: Y(:,t:-1:t-arp+1) where t="now"
% Current posterior summaries: 
%    M = current estimate of parameters
% Outputs: f = mean vector
f=zeros(q,1); Ft = [1; reshape(Yt,p-1,1)]; 

for t = 1:k
          f = M'*Ft;             
          Ft(2:p) = [ f ; Ft(2:p-q) ];         
      end 
end


