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

%% simulate datasets for fitting and comparison
simulations      = 500;

for sim = 1:simulations

    for cond = 1:condition
    % simulate dataset
        data                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        cond_data{1,cond}   = data; % if many subjects add them to cell
    end

    % for modelling with VKF we only need the outcomes for both conditions, so, concatinate the
    % trials of the two conditions
    f                       = cat(1,cond_data{1,1}.feedback, cond_data{1,2}.feedback); % outcomes
    p                       = cat(1,cond_data{1,1}.feedbackprob, cond_data{1,2}.feedbackprob); % true probabilities

    % only keep outcomes for vertical shape
    simy(:,sim)             = f(:,1);
    simx(:,sim)             = p;

end

%% run the two models 

% simulate learning rate and volatility values for the two models 
[tx_vkf, tx_hgf] = sim_params(simy,simx,feedback(:,1));
%[mc nc tx]  = mainf(xs,ys);

% RUN VKF MODEL
% get params for model
lambda              = tx_vkf(1);
v0                  = tx_vkf(2);
omega               = tx_vkf(3);

% run vkf model
[m, k, v]   = vkf_bin(feedback(:,1),lambda,v0,omega);

val1        = m;
vol1        = v;
lr1         = k;

% RUN HGF MODEL
% get hgf parameters
nu                  = tx_hgf(1);
kappa               = tx_hgf(2);
omega               = tx_hgf(3);

% run hgf model
[~, ~, mu2, mu3, sigma2] = hgf_bin(feedback(:,1),nu,kappa,omega);  

m                   = mu2(1:end-1);
v                   = (mu3(1:end-1));
sigma2              = sigma2(2:end);    
val2                = m;
vol2                = v;
lr2                 = sigma2;

%% plot VKF and HGF bin results 

% get all params in a struct
allsignals.m1       = val1;
allsignals.m2       = val2;
allsignals.v1       = vol1;
allsignals.v2       = vol2;
allsignals.lr1      = lr1;
allsignals.lr2      = lr2;

% prepare true probability and outcomes for plotting 
probability         = feedbackprob;
outc                = feedback(:,1); % plot outcomes for probability of vertical shapes only 

% plot results
plotVFK_HGF(probability,outc,val1,vol1,lr1,val2,vol2,lr2);

%% make table with the results
results_table = table(predicted_state, predicted_state2, volatility, volatility2, learning_rate, learning_rate2);

% store table in .xlsx format
filename = 'VKF_HGFbin_modelout.xlsx';
writetable(results_table,filename, 'Sheet', 1)
movefile('*.xlsx', outpath) % move file to output dir 
