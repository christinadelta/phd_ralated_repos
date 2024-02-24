% AL Simple RL model sensitivity analysis
% Created February 2024

%%%% Sensitivity Analysis with Parameter Recovery:   %%%%%%

% In sensitivity analysis with parameter recovery, we simulate data with 
% only varying one parameter at a time but also fit the model to the 
% simulated data generated with that varied parameter values. This involves 
% "freeing" one parameter (letting it vary across a range of values), 
% generating data based on those values, and then trying to recover those 
% values by fitting the model to the generated data. The other parameters 
% can be set to typical or baseline values during the data generation phase 
% but are also estimated during the fitting process to see how well the 
% model can recover the known (simulated) parameter values.

%%%%% Objective and Interpretation:  %%%%%

% The objective of varying one parameter at a time in sensitivity analysis 
% is to understand the influence of each parameter on the model's behavior 
% or performance. In parameter recovery, the goal is to assess the model's 
% ability to accurately estimate parameters from data, particularly how 
% the estimation of one parameter might be influenced by its true value 
% and the values of other parameters.

% If the model consistently recovers the correct values for the varied 
% parameter across its range, this suggests good parameter identifiability 
% and robustness. If recovery is poor, this might indicate issues such as 
% parameter correlations, insufficient sensitivity of the model outputs to 
% the parameter, or other model specification issues.

% ----------

%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
% datadir         = fullfile(startpath,'data');
% subs            = dir(fullfile(datadir, '*sub-*'));
% nsubs           = length(subs);
% 
% figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
total_trials    = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
nCues           = 2;
% params          = [0.3 1.5 0.1];

data            = ALsimdata_v2(probabilities, total_trials,condtrials);
x               = data.x; % extract true reward rate for plotting

%% start sensitivity analysis: ALPHA FREE

% how many repetitions?
repetitions = 100;

% define ranges for alpha the values of the fixed parameters  
alpha_range     = [0 1]; % range for alpha
beta_fix        = 1.5;
decay_fix       = 0.1;

for i       = 1:reps

    % Generate random true values for alpha and beta within the specified ranges
    true_alpha              = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand();
    allsim_alpha(i)         = true_alpha;

    % simulate data
    d                       = ALsimdata_v2(probabilities, total_trials,condtrials);
    out                     = modelRW_v2([true_alpha beta_fix decay_fix], d); % simulate model
    allsim_output{i}        = out;
    sim_actions             = out.a;

    

    













end % end of reps loop



