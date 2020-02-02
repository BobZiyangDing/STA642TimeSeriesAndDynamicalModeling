startup

np=ceil(arp/2); t=0:pi/20:2*pi; a=0.08; b=sqrt(1-a*a); 

for i=1:q
clf
for j=1:arp
    subplot(np,2,j);    
    plot(sin(t),cos(t),'k'); hold on; plot([-a 0 0 a],[b 1+2*a 1+2*a b],'k'); axis off; hold on; axis square
    [x,y]=meshgrid(linspace(-1,1,100),linspace(-1,1,100));
    z=griddata(eeg_electrode_layout(:,1),eeg_electrode_layout(:,2),...
           M((j-1)*q+(1:q),i),x,y,'v4');
    z((x.*x+y.*y)>1)=NaN;  
    contour(x,y,z,12); text(chanlayout(:,1),chanlayout(:,2),channames); hold off
    title([char(Ynames(i)),' : lag ',int2str(j)])
end
display( [ char(Ynames(i)),'...' ])
pause; end


i=10;
clf
for r = 100:10:T 
   M=squeeze(sMt(:,:,r));    
  for j=1:arp
    subplot(np,2,j);    
    plot(sin(t),cos(t),'k'); hold on; plot([-a 0 0 a],[b 1+2*a 1+2*a b],'k'); axis off; hold on; axis square
    [x,y]=meshgrid(linspace(-1,1,100),linspace(-1,1,100));
    z=griddata(eeg_electrode_layout(:,1),eeg_electrode_layout(:,2),...
           M((j-1)*q+(1:q),i),x,y,'v4');
    z((x.*x+y.*y)>1)=NaN;  
    contour(x,y,z,12); text(chanlayout(:,1),chanlayout(:,2),channames); hold off
    title([char(Ynames(i)),' : lag ',int2str(j)])
  end
    display(int2str(r));   pause(0.0001)
end


