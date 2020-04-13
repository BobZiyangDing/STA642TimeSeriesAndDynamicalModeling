function decomp_plot(x,decomp,k)
% Plots x plus the first k (at most) series in decomp
%
  T=length(x); arx=reshape(x,T,1); 
  nk=min(k,size(decomp,1));
  mi=min(arx);  ma=max(arx);  d=1/(ma-mi);  h=0.02*T;
  cen=mean(d*(arx-mi));
  plot(d*(arx-mi)-cen,'k');  text(-6*h,0,'Data');
  xlabel('Time'); 
  box off;
  axis([-h T+h -1 nk+cen]);
  set(gca,'Xgrid','on');
  hold on; 
  for i=1:nk
      plot(-cen+i+d*(decomp(i,:)-mi));
      text(-5.5*h,i,'comp');
  end;
  hold off;
%
