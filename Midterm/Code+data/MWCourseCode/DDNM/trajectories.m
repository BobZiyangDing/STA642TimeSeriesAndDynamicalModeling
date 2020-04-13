
% Look at FF trajectories of coefficients & volatilities
% This follows running of FF in the DDNMexample.m code
% 
  
xatraj=[xa(1:71),'xlim([0 tfore+1]);box off'];
i=1; 
for j=1:q                          % look at series j
    display(['Hit a key for ',deblank(names(j,:)),' model state vector'])
    pause    
    for ifig=1:i, figure(ifig); clf; end
    yj = Y(j,1:tfore)';            % data over this time period 1:tfore
    fj = pq1tfore(j).f';           % one-step forecast means
    rqj = sqrt(pq1tfore(j).q)';    % one-step forecast SD scale factors (T distns)
    nj = pq1tfore(j).n';           % filtered dofs
    rsj = sqrt(pq1tfore(j).s)';    % filtered estimates of volatility SD
    ej = yj-fj;                    % one step forecast errors
    mj = pq1tfore(j).m;            % filtered means of coefficients
    Cj = pq1tfore(j).C;            % filtered covariance of coefficents, up to T dofs...
    tq95 = tinv(0.95,nj); tq75 = tinv(0.75,nj); % T %iles for intervals    
    % intercept:
    i=1; coeffij='Intercept'; plottrajectory
    % parental coeffs: nothing happens if there are no parents, i.e. if npa(j) is 0
    for ipa=1:npa(j) 
        i=1+ipa; coeffij=[deblank(names(pa{j}(ipa),:)),' lag-0']; plottrajectory
    end
    % lagged coeffs: nothing happens if there are no lagged predictors, i.e. if npr(j) is 0
    % move through lags h=1:t0 where t0=max lag allowed 
    for h=1:t0, 
        ijh = find(squeeze(pr(j,:,h)));  % finds lag h coefficients, if any
        for r=1:length(ijh)        
            i=i+1; coeffij=[deblank(names(ijh(r),:)),' lag-',int2str(h)]; plottrajectory
        end
    end % end loop over lags h
end  % end loop over series j


return

figure(1); clf; subplot(2,1,1); plot(abs(ej./rqj),'h-'); eval(xatraj)
                    subplot(2,1,2); plot(rsj); eval(xatraj)
             
  