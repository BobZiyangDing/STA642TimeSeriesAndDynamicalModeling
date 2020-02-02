% This scripts is called from trajectories.m     
% plots trajectory of selected coefficient i for series j, coeff name in coeffij

    mji=mj(i,:)'; rCji=sqrt(squeeze(Cj(i,i,:))); 
    figure(i); clf
    ciplot(mji-tq95.*rCji,mji+tq95.*rCji,1:tfore,[.9 .9 .9]);
    eval(xatraj); hold on
    ciplot(mji-tq75.*rCji,mji+tq75.*rCji,1:tfore,[.75 .75 .75]);
    line([0 tfore+1],[0 0],'linestyle','--','color','k'); eval(xa); xlim([0 tfore+1])
    plot(1:tfore,mji,'k-')
    title(['\rm ',deblank(names(j,:)),' \leftarrow ',coeffij, ' trajectory'])
    hold off 
    
    
    