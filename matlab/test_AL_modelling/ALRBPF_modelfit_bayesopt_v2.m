% AL RBPF model fit version 2 -- using BayesOpt toolbox 


% Optimisation Algorithm: BayesOpt routine: 
% https://www.mathworks.com/help/stats/bayesopt.html
% basic optimisation algorithm is used for now. 

% free parameters: 
% Full model: lambda_s, lambda_v, beta, v0 and s0
% Lambdas + beta model: lambda_s, lambda_v, beta
% Lambdas model: lambda_s, lambda_v

% changes introduced:


%%

clc 
clear all

%% dock figures first

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% get subject data
% 
[all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL, allsub_mrtsStc, allsub_mrtsVol,allsub_mrtsAL] = clean_ALdata(datadir,subs,nsubs);

%% extract participant choice data and rts

for sub = 1:nsubs 
    allchoices{1,sub}       = all_ALTdata{1,sub}(:,8); % get responses 
    alloutcomes{1,sub}      = all_ALTdata{1,sub}(:,5); % get outcomes (p(loss|blue)) 
    allrts{1,sub}           = all_ALTdata{1,sub}(:,9);

    alldata{sub,1}.actions  = allchoices{1,sub};
    alldata{sub,1}.outcome  = alloutcomes{1,sub};
    alldata{sub,1}.rts      = allrts{1,sub};

end 

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
stdvals         = [0.1 0.15 0.2];       % define vector of stDev values for each stc level
nCues           = 2;
beta            = 1;

data            = ALsimdata_v3(probabilities, trials,condtrials);


%%