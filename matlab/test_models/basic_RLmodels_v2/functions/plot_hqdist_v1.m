function hdist = plot_hqdist_v1(H_candidates,q_candidates,prior_p_qH)
% function hdist = plot_hqdist_v1(H_candidates,q_candidates,prior_p_qH)

% plot beliefs results:
% H_candidates
% q candidates
% prior (trial-by-trial)

%% ---------
% plot the result
figure; hold on;

for c = 1:15

    hdist = subplot(3,5,c); hold on;

    t = 100:114; % 30:44; % representative trials
    imagesc(H_candidates,q_candidates,prior_p_qH(:,:,t(c))); 
    ylabel('value of q'); 
    xlabel('value of H'); 
    title(['trial ' int2str(t(c))]);

end % end of trials c

end % end of function