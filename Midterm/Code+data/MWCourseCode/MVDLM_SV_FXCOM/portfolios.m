function [ wtarget wtargetSD wtargcon wtargconSD wminvar wminvarSD] = portfolios(f,V,r,names,iprint)
%PORTFOLIOS 
% r = target return for target portfolio
% f,V are forecast mean vector and variance matrix
% names is character array of variable names
% iprint = 1 to print summary results, 0 (or anything else) otherwise
q=length(f); f=reshape(f,q,1); 
K = inv(V); l = ones(q,1); h = K*l; e = (l'*h)*(f'*K*f)-(f'*h)^2; g = (l*r-f)/e; z = -f'*K*g;  u = h'*g; 
wtarget = K*(u*f+z*l); wminvar=h/(l'*h);

%%% constrained portfolio with non-negative weights ...    
%%%  fmincon is a general optimization fn, so can be used for other things:
% options = optimset('Algorithm','sqp','Display','notify'); 
% [wtargcon,val] = fmincon(@nestedfun,ones(q,1)/q,[],[],[ones(1,q);f'],[1 r]',zeros(q,1),[],[],options);
%    function y = nestedfun(w); y=w'*V*w; end

% this is direct constrained quadratic programming: simpler, faster ...
% options = optimset('Algorithm','active-set','Display','final','LargeScale','off'); 
options = optimset('Algorithm','interior-point-convex','Display','off','LargeScale','off'); 
[wtargcon, val] = quadprog(2*V,zeros(q,1),[],[],[ones(1,q);f'],[1 r]',zeros(q,1),[],ones(q,1)/q,options); 

wtargcon(wtargcon<0)=0; % just to make -0.000... into 0
%  risks:    
wtargetSD = sqrt(wtarget'*V*wtarget); wminvarSD = sqrt(wminvar'*V*wminvar); wtargconSD = sqrt(wtargcon'*V*wtargcon);

if iprint
    display(' ')
    a= char([ num2str(.001*floor(1000*[wtargcon wtarget wminvar f sqrt(diag(V))]),2) repmat('     ',q,1) names ]);
    display(char('wtargcon  wtarget  wminvar    f       sd      var',a))
    display('  ')
    a=[char('tarcon :  ','target :  ','minvar :  ') num2str([[wtargcon'*l wtargcon'*f  wtargconSD];[wtarget'*l wtarget'*f  wtargetSD];[wminvar'*l wminvar'*f wminvarSD]],4)];
    display(char('w       wsum      mean        sd',a))
end
%
end

