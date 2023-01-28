%  PART 1 OF THE AVERSIVE - LEARNING VERSION. 

% Date created: 29/12/2022
% modeified 1: 15/1/2023

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
condition       = 1;                      % only stable condition for now (if 2 = stable and volatile)
task            = 1;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];          % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20]; % probabilities of the stable + volatile condition
end

labels          = {'alpha', 'beta'};      % for ploting
trials          = 100;                    % per volatility condition

%% simulate one dataset

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
    filename = fullfile(figpath, sprintf('stable_%d.fig', task));
    saveas(subfig(task),filename)

end

%% explore parameter values

% use different parameter values and re-run the model 
alphas  = [0.25 0.5 0.75 1]; % use different alpha values and re-run the model 
betas   = [3 5 9 15];

for rep = 1:500
    for alpha = 1:length(alphas)
    
        for beta = 1:length(betas)
    
            for cond = 1:condition
    
                temp_params = [alphas(alpha) betas(beta)]; 
                
                % simulate data 
                data        = avlearn_simulate_v1(condition, probs, trials, outpath, task); % data is a structure containaining the simulation output
                
                % run model
                modelout    = modelRW_v1(temp_params, data, outpath);
                
                % store mean probability of choosing the vertical gabor (high
                % prob gabor) for all combinations of alphas and betas
                % prob_vertical{cond}(alpha,beta,rep)     = nanmean(modelout.allPs(:,1)); % average choice prob for vertical over trials
                

                % extract actions/choices for the high probability
                if cond == 1
                    [~, imax] = max(probs(1,:));
                else
                    [~, imax] = max(probs(2,:)); % probabilities are different for volatile condition 
                end

                % store trial-by-trial vertical-gabor choices for ploting 
                trlbytrl_choices{alpha}{beta}(rep,:) = modelout.a == imax;
                
                
                % store correct choices for the vertical gabor (high
                % prob gabor) for all combinations of alphas and betas
                % correctchoices{cond}(alpha,beta,rep) = nanmean(modelout.a == imax);
                correctchoicese{cond}(alpha,beta,rep) = nanmean(modelout.a(1:10) == imax); % only get the first 10 trials
                correctchoicesl{cond}(alpha,beta,rep) = nanmean(modelout.a(end-9:end) == imax); % only get the last 10 trials

            end % end of condition loop
    
        end % end of betas

    end % end of alphas loop
end % end of repetition loop

%% if stable condition look at trial-by-trial choices-vertical 

% loop over alphas and betas to extract reps and average them 
for i = 1:length(alphas)

    for j = 1:length(betas)

        avtrl_choices{1,i}(j,:) = mean(trlbytrl_choices{1,i}{1,j});
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

% visualise model results
% visualise the two datasets 
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
    filename = fullfile(figpath, sprintf('plot_%s.fig', this_data.volatility));
    saveas(subfig(cond),filename)

end % end of condition loop

