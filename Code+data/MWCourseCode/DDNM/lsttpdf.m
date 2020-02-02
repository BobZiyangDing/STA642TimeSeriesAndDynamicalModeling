function lp = lsttpdf(x,n,m,v)
% log Student T_n(x|m,v) density  
lp = gammaln((n+1)/2)-gammaln(n/2)-log(n*v*pi)/2-log(1+(x-m).^2/(n*v))*(n+1)/2;

