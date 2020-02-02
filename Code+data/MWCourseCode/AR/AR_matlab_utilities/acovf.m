function y=acovf(x)
% acvf(Z) sample autocovariances - uses the filter function in matlab
z=x(:)'-mean(x);
n = length(z);
y = filter(z(n:-1:1),1,z);
y = y(n:-1:1)/n; 
