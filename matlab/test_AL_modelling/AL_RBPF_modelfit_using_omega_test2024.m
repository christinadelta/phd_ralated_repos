% AL model fitting test using the same method is in Piray and Daw (2021)
% simulations.

% Generation of linear outcomes using some omega value and volatility 

% Optimisation Algorithm: BADS https://github.com/acerbilab/bads
% basic optimisation algorithm is used for now. See bads_examples.m 

% free parameters: lambda_s, lambda_v

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

%% test the objective function with one subject 

testsub     = alldata{4,1};

% exatact what we need from data struct and prepare data for model fitting 
state(:,1)  = data.x;
state(:,2)  = 1 - state(:,1);
ss          = data.stcind;
vv          = data.t;
omega       = .01;
N           = length(state);

o(:,1)          = state(:,1) + sqrt(omega)*randn(N,1); % outcomes for p(loss|option A)
o(:,2)          = state(:,2) + sqrt(omega)*randn(N,1); % outcomes for p(loss|option B)

testsub.o       = o;                                    % add linear outcomes to the data struct
parameters      = [0.3 0.3 1 0.8 1 100 1 1 0.001];      % parameters for modelling 
x               = [0.3 0.3];                            % Starting points of free parameters (lambda_s and lambda_v)

% parameters: 
% 1st: lambda_s, -- free (not used here, see x(1))
% 2nd: lambda_v, -- free (not used here, see x(2))
% 3rd: beta,
% 4rd: s0,
% 5th: v0,
% 6th: nparticles
% 7th: x0_unc
% 8th: model type (1=core, 2=stc lesioned, 3=vol lesioned)
% 9th: prior for lesion models

config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
o               = testsub.o;

% I also need to reverse the actions. The output for Gorilla, saves as
% correct action the choice for the low probability cue ( because this is
% the good/no-loss option) we need to swap this and test the model 
a               = testsub.actions;
previousTwos    = a == 2;
a(a == 1)       = 2;
a(previousTwos) = 1;

% which model should we run?
NLL = rbpf_lambdas(x,o,a,parameters,config);

%% test model: fit the RBPF using the BADS toolbox and outcomes generated using omega and state 

% fit lambda RBPF model for each subject seperately 
% exatact what we need from data struct and prepare data for model fitting 
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

% loop over subjects 
for sub = 1:nsubs 

    dat = alldata{sub,1};
    omega       = .01;
    N           = length(state);

    o(:,1)          = state(:,1) + sqrt(omega)*randn(N,1); % outcomes for p(loss|option A)
    o(:,2)          = state(:,2) + sqrt(omega)*randn(N,1); % outcomes for p(loss|option B)

    % define parameters 
    parameters      = [0.4 0.2 1 1.2 1 100 1 1 0.001];  % parameters for modelling 
    config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);

    x0  = [parameters(1) parameters(2) parameters(3) parameters(4) parameters(5)];  % Starting point
    lb  = [0.1 0.1 0.1 0.1 0.1];                                                    % Lower bounds
    ub  = [0.9 0.9 10 2 2];                                                             % Upper bounds
    plb = [0.2 0.3 0.5 0.2 0.2];                                                    % Plausible lower bounds
    pub = [0.9 0.8 5 2 2];                                                          % Plausible upper bounds

    dat.o           = o;                                    % add linear outcomes to the data struct
    a               = dat.actions;
    previousTwos    = a == 2;
    a(a == 1)       = 2;
    a(previousTwos) = 1;
   

    % Set BADS options
    options                         = bads('defaults');
    options.MaxObjectiveEvaluations = 50;
    options.MaxTime                 = 3600;


    % Run BADS, which returns the minimum X and the NLL.
    bFunc       = @(x) rbpf_core_full(x,o,a,parameters,config);
    [x,fval]    = bads(bFunc,x0,lb,ub,plb,pub,[],options);
   
    % now run the model with the minimised x values to get subject learning
    % rates
    [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(x,o,parameters,config);

end % end of subjects loop


%% plot learning rates 

% split learning rates in stc and vol conditions
for j = 1:3

    stc_lr = lr(ss(:,j),1); % extract only lrs for p(loss|option A)

    for k = 1:2

        cond_lrs{j,k} = stc_lr(vv(:,k),:);

    end 
end % end of stc 

[h, g, f] = plotLRs(cond_lrs);


%%
