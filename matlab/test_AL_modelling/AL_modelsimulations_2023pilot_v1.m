% MODEL SIMULATION USING GAUSSIAN NOISE AUGMENTATION FOR BINARY OUTCOMES 

% Created: December 2023
% @christinadelta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Method Name:          %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gaussian Noise Augmentation for Binary Outcomes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Purpose:               %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To transform binary outcome data into a continuous format suitable for 
% statistical analyses or modeling techniques that require or benefit 
% from continuous data inputs.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Method Description:     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method involves adding Gaussian-distributed noise to binary 
% outcomes. The binary outcomes typically represent discrete events or 
% states (e.g., success/failure, win/loss, presence/absence). By 
% introducing Gaussian noise, each binary outcome is transformed into a 
% continuous value, while still reflecting its original state. This 
% transformation is particularly useful for applying statistical models or 
% analytical techniques that assume continuous data distributions, such as 
% regression models, Kalman filters, or other Gaussian-based methods.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%    Mathematical Notation:  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's denote:

% B = {b1,b2,...,bn} as the set of binary outcomes, where $b(i) ∈{0,1} for i=1,2,...,n.
% σ(j) as the standard deviation of the Gaussian noise corresponding to the j(th) stochasticity level.
% G = {g_1,g_2,...,g_n} as the set of Gaussian noise values, where g(i) is drawn from N(0,σ(j)^2), 
% the normal distribution with mean 0 and variance σ(j)^2
% C = {c_1,c_2,...,c_n} as the set of continuous outcomes resulting from this transformation.
% The transformation from binary to continuous outcomes can be mathematically represented as: 
% c(i) = b(i) + g(i) for i=1,2,...,n.

% STEPS:
% 1. simulate data and add Gaussian noise to the outcomes
% 2. run inference model
% 3. run response model to simulate responses?
% 4. Plot results

% Additional stuff to do:
% use different approach to convert binary to continuous outcomes
% use binary outcomes as they are to run the inference model 


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

data            = ALsimdata_v3(probabilities, trials,condtrials);
% data            = ALsimdata_v2(probabilities, trials,condtrials);

%%  test the distribution of the newly simulated outcomes 

std_o = data.std_o;
std_oR = data.std_oR;


% plot outcomes with a histogram 
% Plot a histogram
figure;
histogram(std_o, 'Normalization', 'pdf','NumBins', 120);
title('Histogram of Continuous Outcomes');
xlabel('Value');
ylabel('Probability Density');

%% Q-Q plot
figure;
qqplot(std_o);
title('Q-Q Plot of Continuous Outcomes');

%% Shapiro-Wilk test (if the Statistics and Machine Learning Toolbox is installed)
[h, p] = kstest(zscore(std_o)) % zscore standardizes the data

%% run inerence model 

% run pf_core function which includes both the kalman and particle filters 
output      = infModel(data, probabilities, trials,condtrials,beta);

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


%% prepare data for model fititng 
state(:,1)  = data.x;
state(:,2)  = 1 - state(:,1);
ss          = data.stcind;
vv          = data.t;
test_o      = output.test_o;
test_a      = output.test_a;
test_a(find(test_a==0)) = 2; % convert action 2 to 0

config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
%% fit model to the simulated outcomes and actions 

parameters      = [0.7 0.3 1 1 0.6 100 1 1 0.001];  % parameters for modelling 

%---                       ---%
%--- fit simple BADS model ---%

% We specify the bounds and plausible bounds that (hopefully) 
% contain the solution. Plausible bounds represent your best guess at 
% bounding the region where the solution might lie.

x0  = [parameters(1) parameters(2) parameters(3) parameters(4) parameters(5)];  % Starting point
lb  = [0.1 0.1 0.1 0.2 0.2];                                                    % Lower bounds
ub  = [0.95 0.95 10 2 2];                                                       % Upper bounds
plb = [0.3 0.3 0.5 0.7 0.5];                                                    % Plausible lower bounds
pub = [0.8 0.8 3 1.8 1.8];                                                      % Plausible upper bounds

% Set BADS options
options                         = bads('defaults');
options.MaxObjectiveEvaluations = 50;
options.MaxTime                 = 3600;

% Run BADS, which returns the minimum X and the NLL.
bFunc                       = @(x) rbpf_core_full(x,test_o,test_a,parameters,config);
[x,fval,result]             = bads(bFunc,x0,lb,ub,plb,pub,[],options);

% now run the model with the minimised x values to get subject learning
% rates
[val,vol,unp,lr,unc,choice] = rbpf_coreb_full(x,test_o,parameters,config);

%% plot learning rates 

% split learning rates in stc and vol conditions
for j = 1:3

    stc_lr = lr(ss(:,j),1); % extract only lrs for p(loss|option A)

    for k = 1:2

        cond_lrs{j,k} = stc_lr(vv(:,k),:);

    end 
end % end of stc 

[h, g, f] = plotLRs(cond_lrs);

%% model 2: lambdas_beta RBPF model

parameters      = [0.7 0.3 0.2 1 1 100 1 1 0.001];  % parameters for modelling 

%---                       ---%
%--- fit simple BADS model ---%

% We specify the bounds and plausible bounds that (hopefully) 
% contain the solution. Plausible bounds represent your best guess at 
% bounding the region where the solution might lie.

x0  = [parameters(1) parameters(2) parameters(3)];  % Starting point
lb  = [0.1 0.1 0.1];                                                    % Lower bounds
ub  = [0.95 0.95 10];                                                   % Upper bounds
plb = [0.3 0.3 0.5];                                                    % Plausible lower bounds
pub = [0.8 0.8 3];                                                      % Plausible upper bounds

% Set BADS options
options                         = bads('defaults');
options.MaxObjectiveEvaluations = 50;
options.MaxTime                 = 3600;

% Run BADS, which returns the minimum X and the NLL.
bFunc                           = @(x) rbpf_lambdaBeta_full(x,test_o,test_a,parameters,config);
[x,fval,result]                 = bads(bFunc,x0,lb,ub,plb,pub,[],options);

% now run the model with the minimised x values to get subject learning
% rates
[val,vol,unp,lr,unc,choice]     = rbpf_lambdaBetab_full(x,test_o,parameters,config);

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



