function [x]=outer(a,b,op)

% general outer produce operator for vectors a,b that
% are coerced to column vectors
% e.g., when op ='*', outer(a,b) is a*b'
%        when op ='==', outer(a,b) has ij element ai==bj
% etc

% needs error checks on inputs adding to code ...


na=length(a); nb=length(b);
x=repmat(reshape(a,na,1),1,nb); y=repmat(reshape(b,nb,1),1,na)'; 

if (op=='*'|op=='/'|op=='^') 
   x=eval(['x.',op,'y']);
else
   x=eval(['x',op,'y']);   
end;

