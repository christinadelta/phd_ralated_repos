%  AVERSIVE LEARNING MODELING USING THE VOLATILE KALMAN FILTER -- VERSION 1. 

% Date created: 16/02/2023
% modified 1: 

% The script simulates data the aversive learning task, runs vkf model and
% plots outpus. 
% All functions are called using this script.
% Outputs are also stored in tables for easy access and reading
% figures are saved in .svg and .fig formats
% 

% Information about the task: Aversive learning task Version 1:
% A simple 2-armed bandit choice. More details in the simulation function:
% aversivelearn_sim_v1.m

clear all  
clc

%% set paths & initialise variables

% set paths
addpath(fullfile(pwd,'functions'))
addpath(fullfile(pwd,'main_scripts'))
outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
hgfpath         = fullfile(pwd, 'tapas_HGF');   addpath(hgfpath);
cbmpath         = fullfile(pwd, 'cbm_local');   addpath(cbmpath);

% initialise variables 
subjects        = 1;
condition       = 2;                        % only stable condition for now (if 2 = stable and volatile)
task            = 1;                        % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];            % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20];   % probabilities of the stable + volatile condition

end

model_params    = {};                       % init struct to store all model parameters 

% init parameters fo the VKF model 
model_params.lambda             = 0.1;
model_params.init_vol           = 0.1;
model_params.variance           = 0.1;
model_params.omega              = 0.1;

labels                          = {'lambda', 'initial volatility'...
    'variance', 'omega'};                                   % for ploting
condstring                      = {'stable', 'volatile'};   % for ploting 
trials                          = 100;                      % per volatility condition

% ----------------- FROM NOW ONE DIFFERENT MODELING PARTS ARE PRERFOMED ----------------- %

%% simulate dataset 

for sub = 1:subjects

    for cond = 1:condition
    % simulate dataset
        data                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        cond_data{1,cond}   = data; % if many subjects add them to cell
    end

    % for modelling with VKF we only need the outcomes for both conditions, so, concatinate the
    % trials of the two conditions
    feedback{sub}           = cat(1,cond_data{1,1}.feedback, cond_data{1,2}.feedback);
    feedbackprob{sub}       = cat(1,cond_data{1,1}.feedbackprob, cond_data{1,2}.feedbackprob); % for ploting
%     feedback{sub}           = cat(1,cond_data{1,2}.feedback, cond_data{1,1}.feedback);
%     feedbackprob{sub}       = cat(1,cond_data{1,2}.feedbackprob, cond_data{1,1}.feedbackprob); % for ploting
end
 
%% model the dataset(s)

for sub = 1:subjects
    
    % run vkf model with init parameters
    [predictions, signals]  = vkf_v1(feedback{sub},model_params);
    allsubs_signals{1,sub}  = signals; % all model output

end % end of subjects loop

% look at model output in table
probability                 = feedbackprob{1};
outc                        = feedback{1};
predicted_state             = signals.predictions;
volatility                  = signals.volatility;
PE                          = signals.prediction_error;
learning_rate               = signals.learning_rate;
m1                          = 1./(1+exp(-predictions(:,1))); % bound predicted state between 0-1

signals_table               = table(probability,outc, predicted_state,volatility,learning_rate,PE); % 

% store table in .xlsx format
% filename = 'VKFbin_modelout.xlsx';
% writetable(signals_table,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

%% plot VKF bin results 

h           = plotVFK_bin(volatility(:,1),learning_rate(:,1),m1,probability,outc(:,1));

% save figure
filename    = fullfile(figpath, 'averslearn_VFKbin_plot.fig');
saveas(h, filename)

%% test VKF with different combinations of parameter values  

% simulate dataset for stable & voaltile conditions and extract outcome and true underlying probabilities
for sub = 1:subjects

    for cond = 1:condition
    % simulate dataset
        data                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        cond_data{1,cond}   = data; % if many subjects add them to cell
    end

    % for modelling with VKF we only need the outcomes for both conditions, so, concatinate the
    % trials of the two conditions
    feedback{sub}           = cat(1,cond_data{1,1}.feedback, cond_data{1,2}.feedback);
    feedbackprob{sub}       = cat(1,cond_data{1,1}.feedbackprob, cond_data{1,2}.feedbackprob); % for ploting
%     feedback{sub}           = cat(1,cond_data{1,2}.feedback, cond_data{1,1}.feedback);
%     feedbackprob{sub}       = cat(1,cond_data{1,2}.feedbackprob, cond_data{1,1}.feedbackprob); % for ploting
end

% define 4 different combinations of parameter values 
pv = [.1 .5 .1; % lambda, initial volatility, omega 
    .5 .5 .1;
    .1 .2 .1;
    .1 .5 .2];

% run the model for the different combinations of parameter values 
for sub = 1:subjects

    % loop over the 4 combs
    for i = 1:size(pv,1)

        tmp_params.lambda       = pv(i,1);
        tmp_params.init_vol     = pv(i,2);
        tmp_params.omega        = pv(i,3);
        tmp_params.variance     = 0.1;
    
        % run vkf model with init parameters
        [predictions, signals]  = vkf_v1(feedback{sub},tmp_params);
        % allsubs_signals{1,sub}  = signals; % all model output

        predicted_state(:,i)    = signals.predictions(:,1); % extract only predicted-state only for probability of vertical 
        volatility(:,i)         = signals.volatility(:,1);
        learning_rate(:,i)      = signals.learning_rate(:,1);

        clear tmp_params

    end

end % end of subjects loop

%% plot VKF bin results for different combinations of param values 

% prepare true probability and outcomes for plotting 
probability                 = feedbackprob{1};
outc                        = feedback{1}(:,1); % plot outcomes for probability of vertical shapes only 

fh                          = plotVFK_PVbin(volatility,learning_rate,predicted_state,probability,outc,pv);

% save figure
% filename                    = fullfile(figpath, 'averslearn_VFKbin_PVplot.fig');
% saveas(fh, filename)

%% compare VKF & HGF 

% simulate datasets for fitting and comparison
simulations = 500;

for sim = 1:simulations

    for cond = 1:condition
    % simulate dataset
        data                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        cond_data{1,cond}   = data; % if many subjects add them to cell
    end

    % for modelling with VKF we only need the outcomes for both conditions, so, concatinate the
    % trials of the two conditions
    f           = cat(1,cond_data{1,1}.feedback, cond_data{1,2}.feedback); % outcomes
    p           = cat(1,cond_data{1,1}.feedbackprob, cond_data{1,2}.feedbackprob); % true probabilities

    % only keep outcomes for vertical shape
    ys(:,sim) = f(:,1);
    xs(:,sim) = p;

end


 




