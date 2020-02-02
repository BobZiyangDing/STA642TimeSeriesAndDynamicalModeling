startup

load eeg_electrode_layout
load eeg_alldata.mat
Xeeg=Xeeg(1000:8:end-1000,3:end)'; [q,T]=size(Xeeg);  Ti=1:T; 
u=0:pi/20:2*pi; 
[x,y]=meshgrid(linspace(-1,1,100),linspace(-1,1,100));
a=0.08; b=sqrt(1-a*a);  
    z=griddata(eeg_electrode_layout(:,1),eeg_electrode_layout(:,2),...
          Xeeg(:,Ti(1)),x,y,'v4');
    z((x.*x+y.*y)>1)=NaN;
    f=figure(1); clf; plot(sin(u),cos(u),'k'); hold on; box off; axis off
    surfc(x,y,z); colormap(hot); shading interp
    text(chanlayout(:,1),chanlayout(:,2),'+')
    %title(['Time=',int2str(Ti(1))]); hold off
    rect = get(gcf,'Position'); rect(1:2) = [0 0];
    vidObj = VideoWriter('Movie.avi'); 
    vidObj.FrameRate=2; vidObj.Quality=100; open(vidObj);
    set(gcf,'Renderer','zbuffer'); 
for t=1:20:T
    z=griddata(eeg_electrode_layout(:,1),eeg_electrode_layout(:,2),...
          Xeeg(:,Ti(t)),x,y,'v4');
    z((x.*x+y.*y)>1)=NaN;
    %figure(1); 
    clf
    plot(sin(u),cos(u),'k'); hold on; box off; axis off
    a=0.08; b=sqrt(1-a*a); plot([-a 0 0 a],[b 1+2*a 1+2*a b],'color','k');
    surfc(x,y,z); colormap(hot); shading interp
    text(chanlayout(:,1),chanlayout(:,2),'+')
    title(['Time=',int2str(Ti(t))]); hold off
    %drawnow
    writeVideo(vidObj,getframe(f,rect)); pause(0.1)
end
close(vidObj);

