% Testing the HGF model with the initial parameters for the aversive
% learning task 
% Created: 7/3/2023

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


% extract feedback/outcome for vertical gabor (for simulated sub 1) for
% modelling:
u                           = feedback{1,1}(:,1);

%% run initial model just to take a look at volatility, state predictions and learning rates

% to run the hierarchical learning model we need:
% 1. a perceptual model (e.g., hgf_binary)
% 2. a response model (e.g., hgf_softmax)
% 3. and an optimisation algorithm 
outp = tapas_fitModel([],...
                         u,...
                         'tapas_hgf_binary_config',...
                         'tapas_bayes_optimal_binary_config',...
                         'tapas_quasinewton_optim_config');

% quickly plot to check the outputs
tapas_hgf_binary_plotTraj(outp)

% outp.mat is a structure that contains all the info we need. For plotting
% and comparing with VKF we need volatility -- mu(:,3), state predictions --
% mu(:,2) and dynamic learning rate estimates -- wt(:,1). 

%% 

