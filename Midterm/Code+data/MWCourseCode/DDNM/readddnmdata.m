 
 

[Y,Ynames]=xlsread('USMacroData1965_2016updated.xlsx');
yr1=1965; yrT=2016; 
fprintf('q=9 series, monthly data, 1/1965 - 9/2016 via FRED\n\n')

month=Ynames(:,1); month(1,:)=[]; Ynames=Ynames(1,2:end);  names=char(Ynames); 
 
% choose just 3 series for pilot example 
j=[6 1 3]; Y=Y(:,j); Ynames=Ynames(j); names=names(j,:)  
% % then difference so model is for monthly changes
% Y=diff(Y')'; month=month(2:end); 

[T q]=size(Y); zq=zeros(q,1); oq=ones(q,1); Iq=eye(q); 
yr=yr1:yrT; nyr=yrT-yr1+1; yrmoq=reshape(repmat(yr1:yrT,12,1),nyr*12,1);
    % whole number of years in data set, counts years
nticks=14; h=round(T/(12*nticks)); tticks=1:12*h:T; % plot nticks axis ticks & labels
tdates= int2str(yrmoq(tticks,1)); 
xa=['set(gca,''Xtick'',tticks);set(gca,''XtickLabel'',tdates);xlabel(''year''); xlim([0 T+1]); box off;'];
 
mmy=[ min(Y')' max(Y')'];
figure(1); plot(Y); eval(xa); legend(names); legend boxoff 

 
% 
%-----------------------------------------------------------
%%%%   Below, select quarterly data -- comment out if want monthly
% Qrtly data starts at 1/1965 and ends at 4/2015 ... use cistart to choose starting qtr
fprintf('q=3 selected series, converting to quarterly data\n\n')
%yr1=1982; 
cistart=['1/1/',int2str(yr1)];  mo=char([month]); 
istart = find(sum(mo(:,1:8)==repmat(cistart,T,1),2)==8);
month{[istart]} ;% check
saveallY = Y; 
Y=Y(istart:T,:); T=size(Y,1); nyr=yrT-yr1+1; 
% reduce to quarterly data:    
Y = Y(1:3:end,:); T=size(Y,1);    
yrmoq=reshape(repmat(yr1:yrT,4,1),nyr*4,1);  yrmoq=yrmoq(1:T,:); 
nticks=14; h=round(T/(4*nticks)); tticks=1:4*h:T; % plot nticks axis ticks & labels
tdates= int2str(yrmoq(tticks,1)); 
xa=['set(gca,''Xtick'',tticks);set(gca,''XtickLabel'',tdates);xlabel(''Q1/year''); xlim([0 T+1]); box off;'];
figure(2); plot(Y); eval(xa); legend(names); legend boxoff 
Y=Y'; 



 