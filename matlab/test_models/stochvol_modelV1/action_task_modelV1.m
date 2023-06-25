% first version of the aversive learning task with the
% stochasticity-volatility model 

% STEPS:

% created May 2023 

clear all
clc

%% set paths and init variables

% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'volstoch');    addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

% initialise variables 
subjects        = 1;
condition       = 4;                        % stable & volatile / small & large stochasticity
% task            = 1;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.88 .10;                 % small stochasticity probabilities
    .65 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 200;                      % total trials
condtrials      = [100 25];                 % 100 per stochasticity condition in stable env and 25 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes

% define initial learning rate inverse temperature parameters (don't think
% that this will be used)
params          = [.2 .6;
    3 6]; %alpha and beta parameteres 

% models to run:
mtorun = 3; % (first test the PF using the exact same params as Piray (2021), then run healthy model and stochasticity model

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
    alldata{sub,1}  = data;


end % end of subjects loop


%% run stoch-vol model with the Piray params 

% this will be used to define parameter values etc..
model = 1;

% define parameters and config for the particle filter 
nsim        = 1000;
x           = data.state;
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);

params      = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);

% RUN PF (firts the healthy and then the lesioned) with the parameters
% defined above
stochVolSim = runModel(data, params, config, model, condition,...
    probabilities, trials,condtrials, outpath, outtype);

% plot the results 

fsiz        = [0 0 .45 1];

figure; 

nr          = 2;
nc          = 2;
subplots    = 1:4;

% plot estimated reward, learning rates, estimated volatility and estimated
% stochasticity 
h = plotPF(nr,nc,subplots,stochVolSim,x);



%% run healthy model only - with different parameter values

%% run the stochasticity lesion model with different parameters 



