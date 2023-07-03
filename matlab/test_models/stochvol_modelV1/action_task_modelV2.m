% second version of the aversive learning task with the
% stochasticity-volatility model 

% STEPS:

% created May July 2023
% changes introduced:
% changed the environment of the experiment
% changed the number of trials

clear all
clc

%% TESTING STCS:
% For large stochasticity will test contigences: [68% 32%] & [60% 40%]
% Stable condition will include 40 trials. The environment looks like this:
% Low stc trials: [40 16 16 16 16 16 40] -- first 40 & last 40 constitute stable env
% High stch: [40 16 16 16 16 16 40]

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
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.80 .20;                 % small stochasticity probabilities
    .68 .32];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 160;                      % total trials
condtrials      = [40 20];                  % 80 per stochasticity condition in stable env and 16 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;

% define initial learning rate inverse temperature parameters (don't think
% that this will be used)
% params          = [.2 .6;
%     3 6]; %alpha and beta parameteres 

% models to run:
mtorun = 3; % (first test the PF using the exact same params as Piray (2021), then run healthy model and stochasticity model

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV2(condition, probabilities, trials,condtrials, outpath, outtype, task);
    alldata{sub,1}  = data;


end % end of subjects loop

%% run stoch-vol model with the Piray params 

% this will be used to define parameter values etc..
model = 1;
trls = 161:320; % either: 160, 320 or all

% define parameters and config for the particle filter 
nsim        = 100;

tvolatile   = data.t(trls,2);
tstable     = data.t(trls,1);

for cue = 1:nCues

    x               = data.state(trls,1); % 
    o               = data.outcome(trls,cue); % outcomes computed using outcome variance

    params          = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
    config          = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);

    % RUN PF (firts the healthy and then the lesioned) with the parameters
    % defined above
    [stochVolSim, vals] = runModel(cue, config, model, condition,...
        probabilities, trials,condtrials, outpath, outtype,trls);
    
    % store estimated values for each option/cue
    vvals{1,cue}        = vals;         % values for both options are required to simulate responses 
    allSVsims{1,cue}    = stochVolSim;  % only estimated values for option A are needed for now

end


%% simulate responses and compute performance

% compute responses using vals
ngroups     = 2;

xstate      = data.state(trls,:);
o           = data.feedback(trls,:);
xx          = nan(nsim, ngroups); % I guess I wont include volatility for now 
mR          = nan(nsim, ngroups); 
mxx         = nan(ngroups, 2); % mean correct 
exx         = nan(ngroups, 2); % se correct 
beta        = 5;

% use softmax function to simulate responses for binary choices
% loop over (volatility) conditions 
for j = 1:ngroups
    
    % extract vals for group j and each option/cue
    valA = vvals{1,1}{1,j};
    valB = vvals{1,2}{1,j};

    % loop over simulations
    for i = 1:nsim

        val(:,1) = valA(:,i);
        val(:,2) = valB(:,i);

        % generate responses for control and clinical simulated groups using the softmax function 
        % [xx(i,:),mR(i,:),~,~] = responseModel_v1(xstate, val,tvolatile, tstable, beta);
        [xx(i,:), mR(i,:),~,~] = responseModel_v1(xstate, val,tvolatile, tstable,o, beta);

    end % end of simulations loop

%     mxx(j,:) = median(xx);
%     exx(j,:) = se_median(xx);
    mRR(j,:) = median(mR);
    eRR(j,:) = se_median(mR);
    
    clear xx mR
end % end of volatility condition

%% plot particle filter results 
fsiz        = [0 0 .45 1];

figure; 

nr          = 2;
nc          = 2;
subplots    = 1:4;

% plot estimated reward, learning rates, estimated volatility and estimated
% stochasticity 
h = plotPF(nr,nc,subplots,allSVsims{1,1},x);

%% plot performance 

% plot performance 
col = def_actions('col');
fsy = def_actions('fsy');
subplots = 1;

labels = {'Stable','Volatile'};


h = plot_bar(1,1,subplots(1),{mRR},{eRR},{'control','ASD'},{'Reward','Reward'},'',col);
set(h,'ylim',[0 1]);
legend(h,labels,'fontsize',fsy,'location','north','box','off');
title(h,'Model');

%%


