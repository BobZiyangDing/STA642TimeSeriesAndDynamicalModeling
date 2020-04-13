function sim_plot_wave(times,freqs,waves,k)
% Trajectory plots plus intervals ...
%
nk=length(k); w=2*pi./freqs;
for i=1:nk
  subplot(nk,1,i);
  plot(waves(k(i),:)'); hold on;
  errorbar(times,mean(w(k(i),:,:),3),1.5*std(w(k(i),:,:),0,3),'r*')
end;

