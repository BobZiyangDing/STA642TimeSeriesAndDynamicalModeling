% setup for 7 component normal mixture approximation to log chi-square for
%  SV models

J=7;
q=[  0.00730  0.00002  0.10556   0.25750  0.34001 0.24566 0.04395 ]';
b=[ -5.7002   -4.9186   -2.6216   -1.1793   -0.3255    0.2624    0.7537 ]';
w=[  1.4490    1.2949    0.6534    0.3157    0.1600    0.0851    0.0418 ]';

x=-6+8*((0:500)/501)';
fx=zeros(size(x));; 
for j=1:J, fx=fx+q(j)*normpdf(x,b(j),sqrt(w(j))); end % *exp(-(x-b(j)).^2/(2*w(j)))/sqrt(2*pi*w(j)); end
px=sqrt(2/pi)*exp(x-exp(2*x)/2); 
plot(x, fx); line(x,px,'color','r'); 