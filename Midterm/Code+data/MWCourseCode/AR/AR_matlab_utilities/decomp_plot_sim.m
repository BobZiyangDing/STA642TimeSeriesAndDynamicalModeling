function decomp_plot_sim(x,sdecomp,k)
% Plots x plus the samples of first k (at most) component series in sdecomp
%
  T=length(x); nk=min(k,size(sdecomp,2));
  mi=min(x);  ma=max(x);  d=1/(ma-mi);  h=0.02*T;
  cen=mean(d*(x-mi));
  plot(d*(x-mi)-cen,'k');  text(-6*h,0,'Data');
  xlabel('Time'); 
  box off;
  axis([-h T+h -1 nk+cen]);
  set(gca,'Xgrid','on');
  hold on; 
  for i=1:nk
      decomp=squeeze(sdecomp(:,i,:))'; 
      pr= -cen+i+d*(prctile(decomp,[5 25 50 75 95])'-mi);
      ciplot(pr(:,1),pr(:,5),1:T,[0.8 0.8 0.8]); 
      ciplot(pr(:,2),pr(:,4),1:T,[0.5 0.65 0.65]);
      scatter(1:T,pr(:,3),10,'.');
      text(-5.5*h,i,'comp');
  end;
  hold off;
%