function hprior = plot_modelBeliefBL_v1(est_q,q_space,prior_p_q,y)
% function hprior = plot_modelBeliefBL_v1(est_q,q_space,prior_p_q,y)

% plot beliefs results:
% outcomes (y)
% estimated q
% q space 
% prior (trial-by-trial)

%% ---------

% plot the result
figure; hold on;

for c = 1:15

    hprior = subplot(3,5,c); hold on;

    t = 100:114; % 30:44; % representative trials
    plot([est_q(t(c)) est_q(t(c))],[0 1],'r-.','MarkerSize',15);
    plot(q_space,prior_p_q(t(c),:));
    
    if y(t(c)) == 1
        title(['t' int2str(t(c)) ': vertical']);
    else
        title(['t' int2str(t(c)) ': horizontal']);
    end

    xlabel('value of q'); ylabel('p(value of q)'); set(gcf,'color','w'); set(gca,'XLim',[0 1],'YLim',[0 0.05]);

end % end of trials c

end % end of function