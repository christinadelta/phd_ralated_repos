function fp = plotInd(subject, subj, thisvol)
% FUNCTION that plots data for one subject

% subject number:
sID = subject;



% recode choices to 0 & 1 (1 = correct, 0 = incorrect)
choices = 2-subj.choice(:)';
smoothingkernel = 6; % average across n trials

fp = figure; box off; hold on;

set(fp,'position', [10 60 700 400],'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w')
title(sprintf('%s volatility, subject %02d', thisvol, sID));
plot(subj.prep.feedbackprob,'k-','linewidth',2);    % plot the thick black line (fore feedback probability)
plot(choices,'k*','markersize',10);                 % plot choices on every trial as astersisks
plot(plotSmooth(choices,smoothingkernel,[],'backward'),'r:','linewidth',2)
ylabel('p(reward|blue)');
ylim([-0.1 1.1]); % distance 
xlabel('trial');
legend({'p(reward|blue)','choice','smoothed choice'},'location','northeastoutside')
legend boxoff

end