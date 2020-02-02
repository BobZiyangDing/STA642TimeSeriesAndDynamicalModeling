startup

load('eeg_electrode_layout')
load('eeg_channels_data')
clf;
t=0:pi/20:2*pi;
plot(sin(t),cos(t))
%axis equal
%set(gca,'xtick',((0:20)-10)/10)
%set(gca,'ytick',((0:20)-10)/10)
%grid on
hold on
axis off

[x,y]=meshgrid(linspace(-1,1,100),linspace(-1,1,100));

clf;
for i=1:4
    subplot(2,2,i);
    plot(sin(t),cos(t),'k');
    hold on;
    a=0.08; b=sqrt(1-a*a); plot([-a 0 0 a],[b 1+2*a 1+2*a b]);
    axis off;
    z=griddata(eeg_electrode_layout(:,1),eeg_electrode_layout(:,2),...
          eeg_channels_data(:,i),x,y,'v4');
    z((x.*x+y.*y)>1)=NaN;
%    contour(x,y,z,12);
%    contourf(x,y,z,12);
    surf(x,y,z); colormap hot; shading interp 
%    mesh(x,y,z);
end;    

