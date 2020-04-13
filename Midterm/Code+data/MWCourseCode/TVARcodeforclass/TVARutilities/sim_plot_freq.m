function sim_plot_freq(times,freqs,waves,k)
% Trajectory plots plus intervals ...
%
nk=length(k);
for i=1:nk
  subplot(nk,1,i);
  plot(1./waves(k(i),:)'); hold on;
  errorbar(times,mean(freqs(k(i),:,:),3)/(2*pi),1.5*std(freqs(k(i),:,:)/(2*pi),0,3),'r*');
  hold off;
end;

