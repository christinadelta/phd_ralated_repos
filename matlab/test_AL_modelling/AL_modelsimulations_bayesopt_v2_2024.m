% RBPF MODEL SIMULATIONS AND FITTING USING BAYESOPT 

% CREATE JANUARY 2024 FOR TESTING PURPOSES ONLY 

%% 

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
trials          = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
nCues           = 2;
beta            = 1;

data            = ALsimdata_v2(probabilities, trials,condtrials);

%% run inerence model 

% run pf_core function which includes both the kalman and particle filters 
output      = infModel_v2(data, probabilities, trials,condtrials,beta);

%% plot learning rates 

% extract averaged learning rates 
lrs         = output.mdata(1).mlr; % for option A 
[h, g, f]   = plotLRs(lrs);

%% plot model performance 

o                       = output.binary_o;
actions                 = output.actions;

% compute correct responses
[stable_corr, vol_corr] = getCorrect(o,actions);
[h, g, f]               = plotPerf(stable_corr,vol_corr);

%% prepare data for model fititng using bayesopt

state(:,1)  = data.x;
state(:,2)  = 1 - state(:,1);
ss          = data.stcind;
vv          = data.t;
test_o      = output.test_o;
test_a      = output.test_a;
test_a(find(test_a==0)) = 2; % convert action 2 to 0

config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);




