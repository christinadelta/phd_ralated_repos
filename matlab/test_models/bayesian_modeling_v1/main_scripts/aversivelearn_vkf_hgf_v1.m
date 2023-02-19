%  AVERSIVE LEARNING MODELING USING THE VOLATILE KALMAN FILTER & HIERARCHICAL GAUSSIAN FILTERING -- VERSION 1. 

% Date created: 19/02/2023 --- adapted from Piray & Daw (2020). 

% STEPS PERFOMED IN THIS SCRIPT:
% 1. simulation of 500 datasets 
% 2. run models: hgf & vkf 
% 3. fit models using the cbm toolbox 
% 4. plot results 

clear all
clc

% ----------------
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

labels                          = {'lambda', 'initial volatility'...
    'variance', 'omega'};                                   % for ploting
condstring                      = {'stable', 'volatile'};   % for ploting 
trials                          = 100;                      % per volatility condition


%% simulate one dataset

for cond = 1:condition
% simulate dataset
    data                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
    cond_data{1,cond}   = data; % if many subjects add them to cell
end

% for modelling with VKF we only need the outcomes for both conditions, so, concatinate the
% trials of the two conditions
feedback          = cat(1,cond_data{1,1}.feedback, cond_data{1,2}.feedback);
feedbackprob      = cat(1,cond_data{1,1}.feedbackprob, cond_data{1,2}.feedbackprob); % for ploting



