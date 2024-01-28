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

% loop over stc and volatility levels
for i = 1:3

    % extract this stc outcomes and actions:
    stc_o = test_o(ss(:,i),:);
    stc_a = test_a(ss(:,i),:);

    for j = 1:2

        % extract this volatility outcomes and actions
        vol_o           = stc_o(vv(:,j),:);
        vol_a           = stc_a(vv(:,j),:);

        % re-generate parameter values for recovery
        % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
        lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
        lambda_s    = lambda(1);
        lambda_v    = lambda(2);
        
        % 2. softmax temperatyre beta -- sould be between 0.1 and 10
        beta        = 0.1 + (5-0.1).* rand(1,1); 
        
        % 3. s0 and v0 -- should be between 0.1 and 2
        init_vals   = 0.1 + (2-0.1).* rand(2,1); 
        s0          = init_vals(1);
        v0          = init_vals(2);

        % Define the initial point as an array
        initialValues = [lambda_s lambda_v beta s0 v0];
        otherParams   = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed

        % Evaluate your objective function at the initial point
        initialObjective = rbpf_core_full(initialValues,vol_o,vol_a,otherParams,config); 

        % Define the parameter space
        parameters = [
            optimizableVariable('lambda_s',[0.1, 0.9],'Type','real');
            optimizableVariable('lambda_v',[0.1, 0.9],'Type','real');
            optimizableVariable('beta',[0.1, 5],'Type','real');
            optimizableVariable('s0',[0.1, 2],'Type','real');
            optimizableVariable('v0',[0.1, 2],'Type','real');
            % ... define all your parameters similarly
        ];

        % Setup the optimization function
        optimFcn = @(T) rbpf_core_full([T.lambda_s, T.lambda_v, T.beta, T.s0, T.v0],vol_o,vol_a,otherParams,config);

        % Define the initial point as a row in a table with 5 columns
        initialValuesTable = table(initialValues(1), initialValues(2), initialValues(3), initialValues(4), initialValues(5), 'VariableNames', {'lambda_s', 'lambda_v', 'beta', 's0', 'v0'});


        % Run Bayesian optimization with initial points
        results = bayesopt(optimFcn, parameters, ...
                       'MaxObjectiveEvaluations', 30, ...
                       'AcquisitionFunctionName', 'expected-improvement-plus', ...
                       'InitialX', initialValuesTable);

        % View the results
        bestParams = results.XAtMinObjective;
        minNLL = results.MinObjective;

        fitParams(1) = bestParams.lambda_s(1);
        fitParams(2) = bestParams.lambda_v(1);
        fitParams(3) = bestParams.beta(1);
        fitParams(4) = bestParams.s0(1);
        fitParams(5) = bestParams.v0(1);

         % now run the model with the optimised x values to get subject learning
         % rates
        [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(fitParams,vol_o,otherParams,config);

        % store optimised params 
        allparams(i,j,1) = fitParams(1);
        allparams(i,j,2) = fitParams(2);
        allparams(i,j,3) = fitParams(3);
        allparams(i,j,4) = fitParams(4);
        allparams(i,j,5) = fitParams(5);

        all_lr{i,j} = lr;

    end % end of stc levels
end % end of volatility levels 

%% plot learning rates 

% plot learning rates for each sub seperately
[h, g, f] = plot_fitLRs(all_lr,'simsub');

%%

