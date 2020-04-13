function [f,x]=histf(y,m)
%HISTF 
% histogram with frequency scale rather than counts
% Inputs and outputs are as in hist
if (nargin>1), [n,x]=hist(y,m); 
else, [n,x]=hist(y); m=length(n); end;
f=n/sum(n*(x(2)-x(1))); bar(x,f,'histc'); 
box off
