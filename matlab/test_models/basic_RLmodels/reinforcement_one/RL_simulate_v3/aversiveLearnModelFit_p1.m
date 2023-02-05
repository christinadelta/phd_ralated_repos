% Model fitting for the simple RW model (Part of aversive-learning
% modelling) -- PART 1

% Model fitting using fmincon for the computation of LOG-LIKELIHOOD (LL) 

% created: 04/02/2023

% STEPS:
% simulate dataset(s)
% run model 
% using actions and rewards, fit model to obtain -ll


% -----------------------------------------

clear all  
clc

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);
addpath(fullfile(pwd,'modelfit'))

%% define parameters and init variables

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
condstring      = {'stable', 'volatile'}; % for ploting 
trials          = 100;                    % per volatility condition

%% simulate dataset(s)

for sub = 1:subjects

    for cond = 1:condition
    % simulate dataset
        data                        = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        allsub_data{1,sub}{1,cond}  = data; % if many subjects add them to cell
    end

end % end of subjects loop

%% model the dataset(s)

for sub = 1:subjects
    
    for cond = 1:condition
        this_data                       = allsub_data{1,sub}{1,cond};
        [modelout]                      = modelRW_v1(params, this_data, outpath);
        allsub_modelout{1,sub}{1,cond}  = modelout;
    end

end % end of subjects loop

%% fit model (for each (simulated) subject and condition)

% loop over subjects and conditions to fit the model
for sub = 1:subjects

    for cond = 1:condition

        this_data   = allsub_modelout{1,sub}{1,cond};

        % extract actions and rewards 
        actions     = this_data.a;
        rewards     = this_data.r;
        [xfit ll]   = modelfitRW_v1(params, actions, rewards);

    end
end % end of subjects loop





