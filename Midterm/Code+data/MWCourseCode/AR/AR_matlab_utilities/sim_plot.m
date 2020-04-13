function sim_plot(times,freqs,waves,k,fw)
% Posterior estimates - i.e., trajectories over time 
%       - of frequencies (fw=0) 
%       - or wavelengths (fw~=0) 
%  in TVAR(p) model of series, with approximate 95% posterior intervals ...
%  inputs come from call to tvar_sim function 
nk=min(k,size(freqs,1)); w=2*pi/freqs;
%
for i=1:nk
  subplot(nk,1,i);
  if (fw==0) 
     plot(2*pi./waves(i,:)'); hold on;
     errorbar(times,mean(freqs(i,:,:),3),2*std(freqs(i,:,:),0,3),'r*');
  else
     plot(waves(i,:)'); hold on;
     errorbar(times,mean(w(i,:,:),3),2*std(w(i,:,:),0,3),'r*');
  end
  hold off;
end;
    
 
