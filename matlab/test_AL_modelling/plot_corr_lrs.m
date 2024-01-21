function lr_timeseries = plot_corr_lrs(fitLrs,simlr)

% correlate simulated lrs vs fitted learning rates and plot as a time
% series: 1) across trials, 2) across repetitions 

% how many reps?
reps = size(simlr,2);

%% first re-arrange data and average over repetitions 

% loop over simulations and extract simulated and fitted learning rates
for i = 1:reps

    this_rep_lrs_sim    = simlr{1,i};
    this_rep_lrs_fit    = fitLrs{1,i};
   
    % re-arrange simulated lrs
    all_sim_lrs(:,1,i)  = this_rep_lrs_sim{1,1}(:,1); % small stc/stable lrs
    all_sim_lrs(:,2,i)  = this_rep_lrs_sim{1,2}(:,1); % small stc/volatile lrs
    all_sim_lrs(:,3,i)  = this_rep_lrs_sim{2,1}(:,1); % medium stc/stable lrs
    all_sim_lrs(:,4,i)  = this_rep_lrs_sim{2,2}(:,1); % medium stc/volatile lrs
    all_sim_lrs(:,5,i)  = this_rep_lrs_sim{3,1}(:,1); % large stc/stable lrs
    all_sim_lrs(:,6,i)  = this_rep_lrs_sim{3,2}(:,1); % large stc/volatile lrs

    % re-arrange fitted lrs
    all_fit_lrs(:,1,i)  = this_rep_lrs_fit{1,1}(:,1); % small stc/stable lrs
    all_fit_lrs(:,2,i)  = this_rep_lrs_fit{1,2}(:,1); % small stc/volatile lrs
    all_fit_lrs(:,3,i)  = this_rep_lrs_fit{2,1}(:,1); % medium stc/stable lrs
    all_fit_lrs(:,4,i)  = this_rep_lrs_fit{2,2}(:,1); % medium stc/volatile lrs
    all_fit_lrs(:,5,i)  = this_rep_lrs_fit{3,1}(:,1); % large stc/stable lrs
    all_fit_lrs(:,6,i)  = this_rep_lrs_fit{3,2}(:,1); % large stc/volatile lrs

end % end of reps loop

% reshape 
for j = 1:reps

    tmp_sim                 = all_sim_lrs(:,:,j);
    tmp_fit                 = all_fit_lrs(:,:,j);
    reshaped_sim_lrs(:,j)   = tmp_sim(:);
    reshaped_fit_lrs(:,j)   = tmp_fit(:);
end

%% correlate and plot 

% correlate across trials
for k = 1:size(reshaped_fit_lrs,1)

    this_lrs_sim        = reshaped_sim_lrs(k,:)';
    this_lrs_fit        = reshaped_fit_lrs(k,:)';
    [r,p]               = corrcoef(this_lrs_fit, this_lrs_sim);
    pvals(k)            = p(1,2);
    lr_timeseries(k)    = r(1,2);

end

end % end of function 