function subfig = plot_onesub(allQs, allPs, choices, feedbackprob, subfig, params)

% plot simulated data/model output for one dataset

subfig = figure; box off; hold on;
set(subfig,'position', [10 60 700 400],'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w');
plot(feedbackprob,'k','linewidth',2)
plot(allPs(:,1),'b-','linewidth',2) % plot only choice probabilities for vertical gabor p(vertical)
plot(allQs(:,1),':','color',[.4 0.4 1],'linewidth',2) % plot Q vals for vertical gabor
plot(allQs(:,2),':','color',[1 .35 0.1],'linewidth',2) % plot Q vals for horizontal gabor
plot(choices, 'k*')

legend({'p(reward|verical)','p(choose vertical)','value(vertical)','value(horizontal)',...
    'choice'},'location','northeastoutside');

legend boxoff
ylabel('probability');
ylim([-0.1 1.1]);
xlabel('trial');
title(sprintf('simulated data, alpha = %0.02f, beta = %02.01f',...
    params(1),params(2)));

end