%  VERSION 2 OF THE 2-armed bandit (aversive-learning like task) MODELlING  

% Date created: 29/12/2022 -- Version 1
% modified 1: 15/1/2023
% modified 2: 29/1/2023 (added trial-by-trial plots, and more simulated
% participants)
% modified 3: 1/2/2023 (added heatmap plots to teack performance in
% different parameter combinations)
% modified 4: 20/2/2023 -- Version 2

% In version 2 I simulate data for both conditions and plot them together (combined).
% I keep track of choice probabilities (computed using the softmax
% function) and learning rates
% The steps are:
% 1. simulate dataset(s)
% 2. model simulated dataset(s)
% 3. visualise simulation:


% Information about the task: Aversive learning task Version 2:
% A simple 2-armed bandit choice task setting is used. More details in the simulation function:
% avlearn_simulate_v2.m

clear all  
clc

%%  initialise variables

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
alltrials       = trials*condition;       % total trials

%% simulate dataset(s) for both conditions

for sub = 1:subjects
    
    for cond = 1:condition
        % simulate dataset
        data{1,cond}                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        
    end

    allsub_data{1,sub}  = data; % if many subjects add them to cell

    % for plotting we we need outcomes and underlying probabilities both conditions, so, concatinate the
    % trials of the two conditions
    feedback{sub}           = cat(1,data{1,1}.feedback, data{1,2}.feedback);
    feedbackprob{sub}       = cat(1,data{1,1}.feedbackprob, data{1,2}.feedbackprob); % for ploting

end
    
%% model the dataset(s)

for sub = 1:subjects
    
    % extract this subject dataset
    this_data               = allsub_data{1,sub};
    for cond = 1:condition
        [modelout{1,cond}]              = modelRW_v1(params, this_data{1,cond}, outpath);
        
    end
    
    % add all subjects model output to one cell
    allsub_modelout{1,sub}  = modelout;

    % for plotting we we need Qvalues, choice probabilities and choices both conditions, so, concatinate the
    % trials of the two conditions
    allQs{sub}      = cat(1, modelout{1,1}.Qvals, modelout{1,2}.Qvals);
    allPs{sub}      = cat(1, modelout{1,1}.allPs(:,1), modelout{1,2}.allPs(:,1));
    choices         = cat(2, modelout{1,1}.a, modelout{1,2}.a);
    allcoices{sub}  = 2 - choices'; % convert to be 1 and 0;

end

%% plot model results for one sublect/dataset

% prepare ploting inputs (for one dataset)
probability = feedbackprob{1};
actions     = allcoices{1};
choiceProb  = allPs{1};
Vvals       = allQs{1};

% plot
h = plotRW_V2(probability,choiceProb,Vvals,actions);

% store the plot
filename = fullfile(figpath, 'RW_singlesubPlot.fig');
saveas(h, filename)






