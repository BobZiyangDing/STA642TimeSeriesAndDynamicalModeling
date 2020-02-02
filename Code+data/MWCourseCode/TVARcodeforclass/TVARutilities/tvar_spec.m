function [times,omega,spec]=tvar_spec(m,s,times,freqs)
% plots the time-varying spectrum of AR(p) model
% Inputs:
%   m     -- pxT   post mean vectors of model
%   s     -- T vector of innov vars
%   times -- nt vector of time indices to select 
%   omega -- freq range to select
% Output:
%   times -- nt vector of time indices to select 
%   omega -- freq range to select
%   spec  -- 200 x nt array of spectra
%
nt=length(times); [p,T]=size(m);
omega=min(freqs)+(1:200).*(max(freqs)-min(freqs))/201; spec=zeros(200,nt);
eom=exp(-j*omega);
clf;
if (nt>1)
    for it=1:nt
        t=times(it);
        sp=abs(polyval([-fliplr(m(:,t)') 1],eom));
        spec(:,it)=s(t)./(2*pi*(sp.*sp))';
    end;
    % meshc(times,omega,spec);   view(75,-30);
    surfc(times,omega,spec);   view(-110,-40);
    colormap colorcube; shading interp
    title(['TVAR(',int2str(p),') spectral density']);
    ylabel('Frequency'); xlabel('Times'); zlabel('Density')
    %
else
    sp=abs(polyval([-fliplr(m') 1],eom));
    spec=s./(2*pi*(sp.*sp))';
    plot(omega,spec);  title(['AR(',int2str(p),') spectral density']);
    xlabel('Frequency'); ylabel('Density');
end
