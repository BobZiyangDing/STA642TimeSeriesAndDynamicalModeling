function contourpoints(x,nc,ng,ind,xscat,labs)
% Contours sample points: for one or several pairs of columns of x 
% 
% x     =     p.m array of m samples in p dimensions
% nc    =     number of contours chosen (at fractions of maximum) 
% ng    =     number of grid points in each dimensions 
% ind   =     dimensions to be contoured - all pairs specified in ind are contoured
% xscat =     index vector (if any) of points to scatter on contours
% labs  =     p dim cell array of character string variable labels (if present)
%
    if (nargin<4), 
        error('Not enough arguments'); 
    end
    if (nargin<6),
        labs=cell(0,1); 
        for i=1:p, labs(i,:)=cellstr(['x(',int2str(i),')']); end
        if (nargin<5),
            xscat=[];  
        end
    end
ixscat=(length(xscat)>0); 
clf
%
[p m]=size(x); 
n=length(ind); 
pic=1;
for i=1:(n-1)
  for j=(i+1):n
    subplot(n-1,n-1,pic) 
    ix=ind(i); iy=ind(j);   
    ax=[min(x(ix,:)) max(x(ix,:))]+[-1 1]*range(x(ix,:))/8;
    ay=[min(x(iy,:)) max(x(iy,:))]+[-1 1]*range(x(iy,:))/8;
    vX = linspace(ax(1),ax(2),ng);
    vY = linspace(ay(1),ay(2),ng);
    h=hist2d(x([iy ix],:)',vY,vX);
    vX=(vX(2:end)+vX(1:(end-1)))/2; vY=(vY(2:end)+vY(1:(end-1)))/2;
    contourf(vX,vY,h,nc); cmap('summer');
    hold on; 
    if (ixscat)
        scatter(x(ix,xscat),x(iy,xscat),'b+')
    end
    %axis([ax ay]);
    xlabel(char(labs(i,:))); ylabel(char(labs(j,:)));      
    hold off;
    pic=pic+1;
  end;
  pic=pic+i;
end;
% 
hold off;
% 

