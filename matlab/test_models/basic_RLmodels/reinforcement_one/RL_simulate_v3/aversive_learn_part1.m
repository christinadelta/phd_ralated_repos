%  PART 1 OF THE AVERSIVE - LEARNING VERSION. 

% Date created: 29/12/2022
% modified 1: 15/1/2023
% modified 2: 29/1/2023 (added trial-by-trial plots, and more simulated
% participants)
% modified 3: 1/2/2023 (added heatmap plots to teack performance in
% different parameter combinations)

% In part 1 I only simulate data for the stable condition and let's start
% with simulating data one participant only.
% The steps are:
% 1. simulate dataset(s)
% 2. model simulated dataset(s)
% 3. visualise simulation

% Information about the task: Aversive learning task Version 1:
% A simple 2-armed bandit choice. More details in the simulation function:
% avlearn_simulate_m1.m

clear all  
clc

%% initialise variables

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);

% initialise variables 
subjects        = 1;
params          = [.25 4];                % alpha and beta values 
condition       = 2;                      % only stable condition for now (if 2 = stable and volatile)
task            = 1;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];          % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20]; % probabilities of the stable + volatile condition
end

labels          = {'alpha', 'beta'};      % for ploting
condstring      = {'stable', 'volatile'}; % for ploting 
trials          = 100;                    % per volatility condition

%% simulate dataset(s)

for sub = 1:subjects
    % simulate dataset
    data                = avlearn_simulate_v1(condition, probs, trials, outpath, task); % data is a structure containaining the simulation output
    allsub_data{1,sub}  = data; % if many subjects add them to cell

end
    
%% model the dataset(s)

for sub = 1:subjects

    this_data               = allsub_data{1,sub};
    [modelout]              = modelRW_v1(params, this_data, outpath);
    allsub_modelout{1,sub}  = modelout;

end

%% visualise one simulated dataset/subject

% scoring will be used to visualise the average nb times that subjects
% chose the low aversive outcome 
score           = nan(2,subjects);

% unpack the info we need for visualisation 
allQs           = modelout.Qvals;
allPs           = modelout.allPs;
choices         = 2 - modelout.a; % convert to 0/1 for ploting 
feedbackprob    = data.feedbackprob;

subfig          = nan(2,subjects); % figure handle for ploting 1 simulated dataset

% plot one simulated dataset
subfig          = plot_onesub(allQs, allPs, choices, feedbackprob, subfig, params);

% store the plot
filename = fullfile(figpath, 'stable_noswitch_plot.fig');
saveas(subfig, filename)

%% simulate, model and visualise the 2 stable conditions  

% first simulate the two datasets 
for task = 1:2
    % simulate dataset
    data                = avlearn_simulate_v1(condition, probs, trials, outpath, task); % data is a structure containaining the simulation output
    cond_data{1,task}   = data; % if many subjects add them to cell

end

% model the two datasets 
for task = 1:2

    thiscond_data               = cond_data{1,task};
    [modelout]                  = modelRW_v1(params, thiscond_data, outpath);
    allcond_modelout{1,task}    = modelout;

end

% visualise the two datasets 
subfig          = nan(2,subjects); % figure handle for ploting 1 simulated dataset

for task = 1:2
    
    % extract structures
    this_data       = cond_data{1,task};
    this_modelout   = allcond_modelout{1,task};

    % extract data from structures
    allQs           = this_modelout.Qvals;
    allPs           = this_modelout.allPs;
    choices         = 2 - this_modelout.a; % convert to 0/1 for ploting 
    feedbackprob    = this_data.feedbackprob;
    subfig(task)    = plot_onesub(allQs, allPs, choices, feedbackprob, subfig(task), params); % plot data

    % save figures 
    filename = fullfile(figpath, sprintf('stable_%d_plot.fig', task));
    saveas(subfig(task),filename)

end

%% explore parameter values

% run combparams.m function to simulate data and model with different
% parameter values of alpha and beta

% parameter values to be used:
alphas  = [0.20 0.5 0.75 1];
betas   = [3 5 9 15];
[trlbytrl_choices] = combparams(condition,probs,trials, outpath,task); % run simulations
% [trlbytrl_choices, correctchoicese, correctchoicesl] = combparams(condition,probs,trials, outpath,task); % run simulations

%% if stable condition look at trial-by-trial choices-vertical 

% loop over alphas and betas to extract reps and average them 
for i = 1:length(alphas)

    for j = 1:length(betas)

        avtrl_choices{1,i}(j,:) = mean(trlbytrl_choices{1,1}{1,i}{1,j});
    end
end

% for every alpha parameter seperately, plot trial-by-trial choices for
% each beta parameter value
trl_plts = trlplots(avtrl_choices,alphas,betas);

% store the figure 
filename = fullfile(figpath, 'trialbytrial_params_plot.fig');
saveas(trl_plts, filename)


%% if stable condition look at choice-vertical (early-late) trials

% average choices over beta values for early and late trials
earlystable = nanmean(correctchoicese{1,1},3);
latestable = nanmean(correctchoicesl{1,1},3);

plts = figure; box off; hold on;
% plot the early and late trials 
% for early and late trials
subplot(2,1,1)
l1 = plot(alphas, earlystable);
xlabel('learning rate, \alpha')
ylabel('choice p(vertical)')
title('parameter combinations stable-early', 'FontWeight','normal')

subplot(2,1,2)
l2 = plot(alphas, latestable);
xlabel('learning rate, \alpha')
ylabel('choice p(vertical)')
title('parameter combinations stable-late', 'FontWeight','normal')

% make legend 
for i = 1:length(betas)

    leg{i} = ['\beta = ' num2str(betas(i))];

end

leg1 = legend(l1(end:-1:1), {leg{end:-1:1}});

set(leg1, 'fontsize', 12)
set(leg1, 'position', [0.6267    0.6453    0.1500   0.2000]);

% store the figure 
filename = fullfile(figpath, 'earlyVslatetrls.fig');
saveas(plts, filename)

%% visualise many subjects

% for each subject, extract Q values, choice probabilities, actions for ploting 
for sub = 1:subjects

    allQs(sub,:,:)  = allsub_modelout{1,sub}.Qvals;
    allPs(sub,:)    = allsub_modelout{1,sub}.allPs(:,1); % only extract choice probabilities for vertical 
    choices(sub,:)  = 2 - allsub_modelout{1,sub}.a; % convert to 0-1 for plotting

end % end of subjects loop

fprob               = data.feedbackprob; % needed for plotting 
volatility          = data.volatility;

% plot the averaged sub model
figh = plot_manysubs(allQs, allPs, choices, fprob, params, condition, volatility);

% save the figure
filename = fullfile(figpath, 'plot_manysubs_stable.fig');
saveas(figh, filename)

%% add the volatility component 

% loop over subject and over condition to simulate datasets
for sub = 1:subjects

    for cond = 1:condition

        % simulate data for one sub and condition
        data                        = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        allsub_data{1,sub}{1,cond}  = data; % if many subjects add them to cell
    end % end of condition loop

end % end of subjects loop

% model the stable and volatile conditions 
for sub = 1:subjects
    for cond = 1:condition

        this_data                       = allsub_data{1,sub}{1,cond};               % extract condition dataset
        [modelout]                      = modelRW_v1(params, this_data, outpath);   % run model
        allsub_modelout{1,sub}{1,cond}  = modelout;                                 % store model output
    end
end

%% check the model tables (optional)

% concatenate the tables (both for conditions and simulated data and model
% output)
for sub = 1:subjects

    for cond = 1:condition

        sub_sim             = allsub_data{1,sub}{1,cond}.datatable;
        sub_mod             = allsub_modelout{1,sub}{1,cond}.modeltable;

        condtable{1,cond}   = [sub_sim, sub_mod];

    end % end of condition loop

    tmptbl = [condtable{1,1}; condtable{1,2}];
    subtables{1,sub} = [condtable{1,1}; condtable{1,2}];


end % end of subjects loop

% store sample table for visualisation
filename            = 'sample_table.xlsx'; % save table as xlsx file

writetable(tmptbl,filename, 'Sheet', 1) 
movefile('*.xlsx', outpath) % move file to output dir 

%% visualise model results for one dataset (2 conditions)

% visualise the two conditions
subfig          = nan(2,subjects); % figure handle for ploting simulated dataset for both conditions

for cond = 1:condition

    % extract structures
    this_data       = allsub_data{1,1}{1,cond};
    this_modelout   = allsub_modelout{1,1}{1,cond};

    % extract data from structures
    allQs           = this_modelout.Qvals;
    allPs           = this_modelout.allPs;
    choices         = 2 - this_modelout.a; % convert to 0/1 for ploting 
    feedbackprob    = this_data.feedbackprob;
    subfig(cond)    = plot_onesub(allQs, allPs, choices, feedbackprob, subfig(cond), params); % plot data

    % save figures 
    filename = fullfile(figpath, sprintf('plot_onesub_%s.fig', this_data.volatility));
    saveas(subfig(cond),filename)

end % end of condition loop

%% visualie the results for many subjects (averaged - 2 conditions)

% now that the volatility component is added 
% loop over conditions 
for cond = 1:condition

    % loop over subjects 
    for sub = 1:subjects

        cond_allQs{1,cond}(sub,:,:) = allsub_modelout{1,sub}{1,cond}.Qvals;
        cond_allPs{1,cond}(sub,:)   = allsub_modelout{1,sub}{1,cond}.allPs(:,1);
        cond_choices{1,cond}(sub,:) = 2 - allsub_modelout{1,sub}{1,cond}.a;
        cond_fprob{1,cond}          = allsub_data{1,1}{1,cond}.feedbackprob;
        cond_volatility{1,cond}      = allsub_data{1,1}{1,cond}.volatility;
    end

end 

% plot the averaged sub model
figh = plot_manysubs(cond_allQs, cond_allPs, cond_choices, cond_fprob, params, condition, cond_volatility);

% save figures 
filename = fullfile(figpath, 'plot_manysubs_allconds.fig');
saveas(figh,filename)

%% make trial-wise plots with the two conditions (stable-volatile)

% run many simulations first 
% parameter values to be used:
alphas  = [0.25 0.5 0.75 1];
betas   = [3 5 9 15];
[trlbytrl_choices] = combparams(condition,probs,trials, outpath,task); % run simulations

%% visualise stable and volatile condition look at trial-by-trial choices-vertical 

% loop over alphas and betas to extract reps and average them 
for cond = 1:condition
    for i = 1:length(alphas)
    
        for j = 1:length(betas)
    
            allavtrl_choices{cond}{1,i}(j,:) = mean(trlbytrl_choices{1,cond}{1,i}{1,j});
        end
    end
end

% for every alpha parameter seperately, plot trial-by-trial choices for
% each beta parameter value
for cond = 1:condition

    avtrl_choices = allavtrl_choices{1,cond};

    trl_plts = trlplots(avtrl_choices,alphas,betas);
    
    % store the figure 
    filename = fullfile(figpath, sprintf('trialbytrial_params_plot_%d.fig',cond));
    saveas(trl_plts, filename)

end

%% plot heatmap with combinations of alpha and beta parameter values 


% plot heatmap with averaged choice-probabilities and performance of all possible combinations of beta and alpha parameter
% values
% parameter values to be used:
alphas  = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
betas   = [1 3 5 7 9 11 13 15];

for rep = 1:1000 

    for alpha = 1:length(alphas)

        for beta = 1:length(betas)
            for cond = 1:condition

                tmp_params                              = [alphas(alpha) betas(beta)]; % combine parameter values
                performance                             = heatmap_combparams(tmp_params,cond,probs,trials,outpath,task); % simulate data, run model
                allperformance{1,cond}(alpha,beta,rep)  = performance;

            end % end of condtion loop
        end % end of betas loop
    end % end of alpha loop

end % end of reps loop

% for each condition, average performance over the repetitions
for i = 1:condition

    avhm{1,i} = nanmean(allperformance{1,i},3);

end

% plot the heatmap 
for cond = 1:condition

    cond_hm = avhm{1,cond}; cond_hm = cond_hm';
    t = plotheatmap(cond_hm);
    % thiscond = condstring{1,cond}

    % store the figures
    % filename = fullfile(figpath, sprintf('plot_performance_%s_heatmap.fig',thiscond));
    % saveas(t, filename)

end % end of condition loop

