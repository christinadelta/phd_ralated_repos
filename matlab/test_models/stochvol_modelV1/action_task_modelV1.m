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
probabilities   = [.88 .12;                 % small stochasticity probabilities
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 200;                      % total trials
condtrials      = [100 25];                 % 100 per stochasticity condition in stable env and 25 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes

% define initial learning rate inverse temperature parameters (don't think
% that this will be used)
params          = [.2 .6;
    3 6]; %alpha and beta parameteres 

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
    alldata{sub,1}  = data;


end % end of subjects loop

%% run particle filter 

% unpack data struct 
o           = data.outcome; 
x           = data.state;
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);

% define parameters and config for the particle filter 
params      = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',2,'model_parameters',params);

rng(config.rng_id); 
nsim = config.nsim;

% init params to be estimated 
N           = length(o);
outcome     = nan(N,nsim);
vol         = nan(N,nsim);
stc         = nan(N,nsim);
lr          = nan(N,nsim);
val         = nan(N,nsim);  

% init cells to store the above parameters
outcomes    = cell(1,2);
vols        = cell(1,2); % stable, volatile
stcs        = cell(1,2); % small, large
lrs         = cell(1,2); % 
vals        = cell(1,2); % vals --> action values? 
v_example   = nan(N,2);
glabels     = {'Control','Anxious'};
lnames      = {'Healthy', sprintf('%s lesion',def_actions('stc'))};      




%% generate responses

%% plot results
