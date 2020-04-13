
startup

load eeg_alldata.mat
load eeg_electrode_layout

Y=Xeeg(:,2+ilayout); Y_1_2=Xeeg(:,1:2); clear Xeeg 

figure(1); clf
t=0:pi/20:2*pi; 
plot(sin(t),cos(t),'k'); hold on; a=0.08; b=sqrt(1-a*a); plot([-a 0 0 a],[b 1+2*a 1+2*a b],'k');
axis off; hold on
text(chanlayout(:,1),chanlayout(:,2),channames)
hold off




