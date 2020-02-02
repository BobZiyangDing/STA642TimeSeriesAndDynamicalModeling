
%  Intervene to adjust current realized Int rate r_t to control value r
%  Then simulate evolution from t=tfore to generate forecasts for times t+1, ..., t+k 
% 
%  
X = Y(:,t-t0+1:t);                 % recent lagged Ys (real or simulated)
X(icontrol,end) = rr;              % adjust current end-qtr t value to control
X = [ X NaN(q,1) ];                % add time t+1 column  
X = repmat(X,1,1,I);               % replicate to create array for simulated futures 

sv = zeros(q,I);             % working array of sampled volatilies 
stheta = cell(1,q);          % working cell array of sampled states 
g = zeros(q,2); W=cell(1,q); L=W; 

ik=1; % one-step ahead first: 
for j=q:-1:1             % sample 1-step ahead for series j
       
    p1 = pq1tfore(j);    % extract saved posterior at last time point t-1
    m = p1.m(:,end); C=p1.C(:,:,end); n=p1.n(end); s=p1.s(end); d=n*s;
    delta=prior{5,j}; beta=prior{6,j};
    
    C=C/s;                            % scale free condikional post var at tfore
    g(j,1) = beta*n/2; g(j,2) = (1-beta)*n/2;   % 1-step beta shock parameters
    Wj = C*(1/delta-1);                 % 1-step evolution variance matrix, scale free
    Lj =  chol(Wj)';                    %    and iks transposed Cholesky
    W{j}=Wj; L{j}=Lj; W{j}=Wj; L{j}=Lj;
    
    %  stheta;            % [ 1+p(j) I ] - cols are MC sampled states
    %  sv;                % [ 1 I ] - MC sampled volatilities
    
    svj = 1./gamrnd(n/2,2/d,1,I);                  % draw inv gamma for volatiliy
    svj = beta*svj./betarnd(g(j,1),g(j,2),1,I);    % evolve ik to 1-step ahead
    a=m; R=C+Wj;                                   % evolve and draw 1-step state
    sthetaj = repmat(a,1,I) + repmat(sqrt(svj),1+p(j),1).*( chol(R)'*randn(1+p(j),I) );    % draw 1-step ahead state
    sv(j,:)=svj; stheta{j}=sthetaj;

    % now draw the predictive sample at time t+1
    % first, set up predictor vectors including past sampled predictors
    sft=ones(1+p(j),I);                % to hold step-ahead regression vectors
    isf=1;                             % index to add predictors to sft
    if (npa(j)>0),
        sft(isf+(1:npa(j)),:)= squeeze(X(pa{j},end,:));
        isf=isf+npa(j);
    end
    if (npr(j)>0),  % lagged Y value predictors
        for h=1:t0
            i=find(squeeze(pr(j,:,h))); ni=length(i);
            if (ni>0)
                sft(isf+(1:ni),:) = squeeze(X(i,1+t0-h,:));
                isf=isf+ni;
            end
        end
    end
    % now sample the forecast distribution:
    syj = sum(sft.*sthetaj) + sqrt(svj).*randn(1,I);
   
    % and save ....      
    ypred(j,:,ik) = [ mean(syj) median(syj) std(syj) prctile(syj,[5 25 75 95]) ];
    X(j,end,:) = syj;
    spred(j,ik,:) = syj; 
end % now have samples 1-step ahead made at t = tfore for time t+1 for all series j=1:q


% now sequence over steps ahead up to k .... 
for ik=2:k           % move forecast one step at a time
    % first, update "current" sampled values in working array X, to update sft
    X(:,1:end-1,:)=X(:,2:end,:); X(:,end,:)=NaN; % shift time step up by 1
    
    for j=q:-1:1             % sample ik-step ahead for series j
        
        % evolve one step and predict ...
        svj = prior{6,j}*sv(j,:)./betarnd(g(j,1),g(j,2),1,I);
        sthetaj = stheta{j} + repmat(sqrt(svj),1+p(j),1).*(L{j}*randn(1+p(j),I));
        
        %  set up predictor vectors including past sampled predictors
        sft=ones(1+p(j),I);                % to hold step-ahead regression vectors
        isf=1;                             % index to add predictors to sft
        if (npa(j)>0),
            sft(isf+(1:npa(j)),:)= squeeze(X(pa{j},end,:));
            isf=isf+npa(j);
        end
        if (npr(j)>0),
            for h=1:t0
                i=find(squeeze(pr(j,:,h))); ni=length(i);
                if (ni>0)
                    sft(isf+(1:ni),:) = squeeze(X(i,1+t0-h,:));
                    isf=isf+ni;
                end
            end
        end
        
        syj = sum(sft.*sthetaj) + sqrt(svj).*randn(1,I);
        
        % and save ...
        ypred(j,:,ik) = [ mean(syj) median(syj) std(syj) prctile(syj,[5 25 75 95]) ];
        X(j,end,:) = syj;
        spred(j,ik,:) = syj; 
    
    end % end loop over series j
end  % end loop over steps ahead

clear  stheta sthetaj sv svj X % clean up a bit


