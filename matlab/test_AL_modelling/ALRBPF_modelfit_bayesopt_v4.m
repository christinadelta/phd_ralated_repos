% AL RBPF model fit version 3 -- using BayesOpt toolbox 

% Created: January 2024
% Last update: February 10 2024

% Optimisation Algorithm: BayesOpt routine: 
% https://www.mathworks.com/help/stats/bayesopt.html
% I am using more complex settings where I also add initial parameter
% values (normally BayesOpt only needs lower and upper bounds of parameters but I
% modified this part to also require initial parameter values). 

% free parameters: 
%%%% Full model: lambda_s, lambda_v, beta, v0 and s0 ( and lambda parameter
%%%% for regularisation)
% Lambdas + beta model: lambda_s, lambda_v, beta
% Lambdas model: lambda_s, lambda_v

% changes introduced:
% Added regularisation of the objective function via:
% 1: 
%%%% Ridge Regularisation %%%%:
% a. Shrinkage: Ridge is suitable if you expect that most parameters contribute 
% to the outcome but prefer to shrink their effect sizes to prevent overfitting. 
% It tends to give small, non-zero coefficients.

% b. Multicollinearity: If your model suffers from multicollinearity (high correlations 
% between parameters), Ridge can help mitigate this issue by allowing the inclusion of 
% all features but penalizing their coefficients.

% c. Bias-Variance Trade-Off: Ridge regularization introduces bias into the parameter 
% estimates but reduces their variance, which can lead to better prediction 
% accuracy on unseen data.

% 2: 
%%%% Lasso Regularisation %%%%:

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
    allchoices{1,sub}       = all_ALTdata{1,sub}(:,9); % get responses 
    alloutcomes{1,sub}(:,1) = all_ALTdata{1,sub}(:,5); % get outcomes (low prob -- good option outcomes) 
    alloutcomes{1,sub}(:,2) = all_ALTdata{1,sub}(:,6); % get outcomes (high prob -- bad option outcomes) 
    allrts{1,sub}           = all_ALTdata{1,sub}(:,10);

    alldata{sub,1}.actions  = allchoices{1,sub};
    alldata{sub,1}.outcome  = alloutcomes{1,sub};
    alldata{sub,1}.rts      = allrts{1,sub};

end 

%% extract participant choice data and rts

for sub = 1:nsubs 
    allchoices{1,sub}       = all_ALTdata{1,sub}(:,9); % get responses 
    alloutcomes{1,sub}(:,1) = all_ALTdata{1,sub}(:,5); % get outcomes (low prob -- good option outcomes) 
    alloutcomes{1,sub}(:,2) = all_ALTdata{1,sub}(:,6); % get outcomes (high prob -- bad option outcomes) 
    allrts{1,sub}           = all_ALTdata{1,sub}(:,10);

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
% stdvals         = [0.1 0.15 0.2];       % define vector of stDev values for each stc level
% stdvals         = [0.50 0.75 1.0];       % define vector of stDev values for each stc level
stdvals         = [1 1.5 2];      % define vector of stDev values for each stc level
nCues           = 2;
beta            = 1;

true_stc        = [1 2 3];
true_vol        = [1 2];

data            = ALsimdata_v2(probabilities, trials,condtrials);

%% test the BayesOpt objective function with one subject 

testsub     = alldata{1,1};
randreps    = 10;

% initialize an array to store results from each run
allResults  = cell(1, randreps);


% exatact what we need from data struct and prepare data for model fitting 
state(:,1)  = data.x;
state(:,2)  = 1 - state(:,1);
ss          = data.stcind;
vv          = data.t;
binoutcome  = testsub.outcome; binoutcome(find(binoutcome==2)) = 0; % convert outcome 2 to 0

% apply GNA to binary outcomes 
for i = 1:3

    this_stc        = binoutcome(ss(:,i),:);
    oT(:,i)         = this_stc(:,1) + stdvals(i) * randn(size(this_stc(:,1)));
    oR(:,i)         = this_stc(:,2) + stdvals(i) * randn(size(this_stc(:,2)));

end % end of stc loop

% concatinate outcomes
o(:,1) = oT(:);
o(:,2) = oR(:);

% extract outcomes
a      = testsub.actions;

% loop over repetitions
for rep = 1:randreps

    % set the random seed for reproducibility
    rng(rep); % Each iteration uses a different seed
   
    % loop over stc and vol phases 
    for i = 1:3
    
        % extract this block data:
        stc_o = o(ss(:,i),:);
        stc_a = a(ss(:,i),:);
    
        for j = 1:2
    
            % extract this volatility outcomes and actions
            vol_o           = stc_o(vv(:,j),:);
            vol_a           = stc_a(vv(:,j),:);
            config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);

            % set up the BayesOpt configuration
            % Define the parameter space
            params = [
                optimizableVariable('lambda_s',[0.1, 0.9],'Type','real');
                optimizableVariable('lambda_v',[0.1, 0.9],'Type','real');
                optimizableVariable('beta',[0.1, 5],'Type','real');
                optimizableVariable('s0',[1, 3],'Type','real');
                optimizableVariable('v0',[1, 3],'Type','real');
                optimizableVariable('lambda', [1e-4, 1e4], 'Transform', 'log'); 
                % 
            ];
    
            otherParams   = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed
    
            % Define the objective function as an anonymous function if needed
            objectiveWithInputs = @(params) bayesOptObjective(params, vol_o, vol_a, otherParams, config);
    
            % Example initial values for the parameters
            initialLambdaS      = 0.1;      % Example initial value for lambda_s
            initialLambdaV      = 0.1;      % Example initial value for lambda_v
            initialBeta         = 0.5;      % Example initial value for beta
            initialS0           = 1;        % Example initial value for s0
            initialV0           = 1;        % Example initial value for v0
            initialLambda       = 0.5;
            
            % Creating the table with initial guesses
            initialGuessTable = table([initialLambdaS], [initialLambdaV], [initialBeta], [initialS0], [initialV0], [initialLambda], ...
                'VariableNames', {'lambda_s', 'lambda_v', 'beta', 's0', 'v0', 'lambda'});
    
    
            % Then use this anonymous function in the BayesOpt call
            results = bayesopt(objectiveWithInputs, params, ...
                'IsObjectiveDeterministic', false, ...
                'MaxObjectiveEvaluations', 50, ...
                'InitialX', initialGuessTable, ...
                'PlotFcn', []); % This suppresses the display of plots
    
            % View the results
            bestParams          = results.XAtMinObjective;
            minNLL              = results.MinObjective;
    
            fitParams(1)        = bestParams.lambda_s(1);
            fitParams(2)        = bestParams.lambda_v(1);
            fitParams(3)        = bestParams.beta(1);
            fitParams(4)        = bestParams.s0(1);
            fitParams(5)        = bestParams.v0(1);
            bestLambdaReg(1)    = bestParams.lambda(1);
    
             % now run the model with the optimised x values to get subject learning
             % rates
            [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(fitParams,vol_o,otherParams,config);

            % Store the optimization results
            allResults{rep}{i,j} = results;
    
            % store output for each subejct, block
            allX{1,rep}{i,j}           = fitParams;
            allLambdaReg{1,rep}{i,j}   = bestLambdaReg;
            allNLL{1,rep}(i,j)         = minNLL;
            all_lr{1,rep}{i,j}         = lr;
            all_choice{1,rep}{i,j}     = choice;
    
        end % end of volatility phases 
    end % end of stc levels
end % end of repetitions loop

%% average NLL and LR

[meanNLL,allrepsNLL]    = getMeanNLL(allNLL);
[meanLR,allrepsLR]      = getMeanLR(all_lr);

%% plot learning rates

% plot learning rates for test sub
[h, g, f] = plot_fitLRs(meanLR,'test');

clear allX allNLL all_lr all_choice fitParams allLambdaReg allResults

%% Fit the rbpf full model using the BayesOpt toolbox and Ridge regularisation

% exatact what we need from data struct and prepare data for model fitting 
randreps        = 10;
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

for sub = 1:nsubs

    dat         = alldata{sub,1};
    binoutcome  = dat.outcome; binoutcome(find(binoutcome==2)) = 0; % convert outcome 2 to 0

    % apply GNA to binary outcomes 
    for i = 1:3
    
        this_stc    = binoutcome(ss(:,i),:);
        oT(:,i)     = this_stc(:,1) + stdvals(i) * randn(size(this_stc(:,1)));
        oR(:,i)     = this_stc(:,2) + stdvals(i) * randn(size(this_stc(:,2)));
    
    end % end of stc loop

    % concatinate outcomes 
    o(:,1)          = oT(:);
    o(:,2)          = oR(:);

    dat.o           = o;                                    % add linear outcomes to the data struct
    a               = dat.actions;

    % loop over repetitions
    for rep = 1:reps

        % set the random seed for reproducibility
        rng(rep); % Each iteration uses a different seed

        for i = 1:3
    
            % extract this stc outcomes and actions:
            stc_o = o(ss(:,i),:);
            stc_a = a(ss(:,i),:);
    
            for j = 1:2
    
                % extract this volatility outcomes and actions
                vol_o           = stc_o(vv(:,j),:);
                vol_a           = stc_a(vv(:,j),:);
                config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);

                % set up the BayesOpt configuration
                % Define the parameter space
                params = [
                    optimizableVariable('lambda_s',[0.1, 0.9],'Type','real');
                    optimizableVariable('lambda_v',[0.1, 0.9],'Type','real');
                    optimizableVariable('beta',[0.1, 5],'Type','real');
                    optimizableVariable('s0',[1, 3],'Type','real');
                    optimizableVariable('v0',[1, 3],'Type','real');
                    optimizableVariable('lambda', [1e-4, 1e4], 'Transform', 'log'); 
                    % 
                ];
        
                otherParams   = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed
        
                % Define the objective function as an anonymous function if needed
                objectiveWithInputs = @(params) bayesOptObjective(params, vol_o, vol_a, otherParams, config);
        
                % Example initial values for the parameters
                initialLambdaS      = 0.1;      % Example initial value for lambda_s
                initialLambdaV      = 0.1;      % Example initial value for lambda_v
                initialBeta         = 0.5;      % Example initial value for beta
                initialS0           = 1;        % Example initial value for s0
                initialV0           = 1;        % Example initial value for v0
                initialLambda       = 0.5;
                
                % Creating the table with initial guesses
                initialGuessTable = table([initialLambdaS], [initialLambdaV], [initialBeta], [initialS0], [initialV0], [initialLambda], ...
                    'VariableNames', {'lambda_s', 'lambda_v', 'beta', 's0', 'v0', 'lambda'});
        
        
                % Then use this anonymous function in the BayesOpt call
                results = bayesopt(objectiveWithInputs, params, ...
                    'IsObjectiveDeterministic', false, ...
                    'MaxObjectiveEvaluations', 50, ...
                    'InitialX', initialGuessTable, ...
                    'PlotFcn', []); % This suppresses the display of plots
        
                % View the results
                bestParams          = results.XAtMinObjective;
                minNLL              = results.MinObjective;
        
                fitParams(1)        = bestParams.lambda_s(1);
                fitParams(2)        = bestParams.lambda_v(1);
                fitParams(3)        = bestParams.beta(1);
                fitParams(4)        = bestParams.s0(1);
                fitParams(5)        = bestParams.v0(1);
                bestLambdaReg(1)    = bestParams.lambda(1);
        
                 % now run the model with the optimised x values to get subject learning
                 % rates
                [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(fitParams,vol_o,otherParams,config);
    
                % Store the optimization results
                allResults{sub}{rep}{i,j} = results;

                % store output for each subejct, block
                allX{1,sub}{1,rep}{i,j}            = fitParams;
                allLambdaReg{1,sub}{1,rep}(i,j)    = bestLambdaReg;
                allNLL{1,sub}{1,rep}(i,j)          = minNLL;
                all_lr{1,sub}{1,rep}{i,j}          = lr;
                all_choice{1,sub}{1,rep}{i,j}      = choice;
                all_vol{1,sub}{1,rep}{i,j}         = vol;
                all_val{1,sub}{1,rep}{i,j}         = val;
                all_stc{1,sub}{1,rep}{i,j}         = unp;
                all_outcome{1,sub}{1,rep}{i,j}     = binary_outcome; % store binary outcome here? 
    
                clear vol_a vol_o fitParams minNLL lr choice val vol unp bestLambdaReg initialGuessTable results
            end % end of volatility loop
        end % end of stochasticity loop
    end % end of repetitions loop
end% end of subjects loop

%% average NLL and LR
for sub = 1:nsubs

    subNLLs                 = allNLL{1,sub};
    subLRs                  = all_lr{1,sub};
    [meanNLL,allrepsNLL]    = getMeanNLL(subNLLs);
    [meanLR,allrepsLR]      = getMeanLR(subLRs);

    subs_lr_ridge{1,sub}    = meanLR;
    subs_NLL_ridge{1,sub}   = meanNLL;

    clear meanNLL allrepsLR allrepsNLL meanLR

end

%% plot learning rates

% loop over participants to plot learning rates
for sub = 1:nsubs
    sub_lr = subs_lr_ridge{1,sub};
    subnum = sprintf('%d',sub);

    % plot learning rates for each sub seperately
    [h, g, f] = plot_fitLRs(sub_lr,subnum);
end 

%% Fit the rbpf full model using the BayesOpt toolbox and Lasso regularisation

% exatact what we need from data struct and prepare data for model fitting 
randreps        = 10;
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

for sub = 1:nsubs

    dat         = alldata{sub,1};
    binoutcome  = dat.outcome; binoutcome(find(binoutcome==2)) = 0; % convert outcome 2 to 0

    % apply GNA to binary outcomes 
    for i = 1:3
    
        this_stc    = binoutcome(ss(:,i),:);
        oT(:,i)     = this_stc(:,1) + stdvals(i) * randn(size(this_stc(:,1)));
        oR(:,i)     = this_stc(:,2) + stdvals(i) * randn(size(this_stc(:,2)));
    
    end % end of stc loop

    % concatinate outcomes 
    o(:,1)          = oT(:);
    o(:,2)          = oR(:);

    dat.o           = o;                                    % add linear outcomes to the data struct
    a               = dat.actions;

    % loop over repetitions
    for rep = 1:reps

        % set the random seed for reproducibility
        rng(rep); % Each iteration uses a different seed

        for i = 1:3
    
            % extract this stc outcomes and actions:
            stc_o = o(ss(:,i),:);
            stc_a = a(ss(:,i),:);
    
            for j = 1:2
    
                % extract this volatility outcomes and actions
                vol_o           = stc_o(vv(:,j),:);
                vol_a           = stc_a(vv(:,j),:);
                config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);

                % set up the BayesOpt configuration
                % Define the parameter space
                params = [
                    optimizableVariable('lambda_s',[0.1, 0.9],'Type','real');
                    optimizableVariable('lambda_v',[0.1, 0.9],'Type','real');
                    optimizableVariable('beta',[0.1, 5],'Type','real');
                    optimizableVariable('s0',[1, 3],'Type','real');
                    optimizableVariable('v0',[1, 3],'Type','real');
                    optimizableVariable('lambda', [1e-4, 1e4], 'Transform', 'log'); 
                    % 
                ];
        
                otherParams   = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed
        
                % Define the objective function as an anonymous function if needed
                objectiveWithInputs = @(params) bayesOptObjectiveLasso(params, vol_o, vol_a, otherParams, config);
        
                % Example initial values for the parameters
                initialLambdaS      = 0.1;      % Example initial value for lambda_s
                initialLambdaV      = 0.1;      % Example initial value for lambda_v
                initialBeta         = 0.5;      % Example initial value for beta
                initialS0           = 1;        % Example initial value for s0
                initialV0           = 1;        % Example initial value for v0
                initialLambda       = 0.5;
                
                % Creating the table with initial guesses
                initialGuessTable = table([initialLambdaS], [initialLambdaV], [initialBeta], [initialS0], [initialV0], [initialLambda], ...
                    'VariableNames', {'lambda_s', 'lambda_v', 'beta', 's0', 'v0', 'lambda'});
        
        
                % Then use this anonymous function in the BayesOpt call
                results = bayesopt(objectiveWithInputs, params, ...
                    'IsObjectiveDeterministic', false, ...
                    'MaxObjectiveEvaluations', 50, ...
                    'InitialX', initialGuessTable, ...
                    'PlotFcn', []); % This suppresses the display of plots
        
                % View the results
                bestParams          = results.XAtMinObjective;
                minNLL              = results.MinObjective;
        
                fitParams(1)        = bestParams.lambda_s(1);
                fitParams(2)        = bestParams.lambda_v(1);
                fitParams(3)        = bestParams.beta(1);
                fitParams(4)        = bestParams.s0(1);
                fitParams(5)        = bestParams.v0(1);
                bestLambdaReg(1)    = bestParams.lambda(1);
        
                 % now run the model with the optimised x values to get subject learning
                 % rates
                [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(fitParams,vol_o,otherParams,config);
    
                % Store the optimization results
                allResults{sub}{rep}{i,j} = results;

                % store output for each subejct, block
                lasso_allX{1,sub}{1,rep}{i,j}            = fitParams;
                lasso_allLambdaReg{1,sub}{1,rep}(i,j)    = bestLambdaReg;
                lasso_allNLL{1,sub}{1,rep}(i,j)          = minNLL;
                lasso_all_lr{1,sub}{1,rep}{i,j}          = lr;
                lasso_all_choice{1,sub}{1,rep}{i,j}      = choice;
                lasso_all_vol{1,sub}{1,rep}{i,j}         = vol;
                lasso_all_val{1,sub}{1,rep}{i,j}         = val;
                lasso_all_stc{1,sub}{1,rep}{i,j}         = unp;
    
                clear vol_a vol_o fitParams minNLL lr choice val vol unp bestLambdaReg results initialGuessTable
            end % end of volatility loop
        end % end of stochasticity loop
    end % end of repetitions loop
end% end of subjects loop

%% average NLL and LR
for sub = 1:nsubs

    subNLLs                 = lasso_allNLL{1,sub};
    subLRs                  = lasso_all_lr{1,sub};
    [meanNLL,allrepsNLL]    = getMeanNLL(subNLLs);
    [meanLR,allrepsLR]      = getMeanLR(subLRs);

    subs_lr_lasso{1,sub}    = meanLR;
    subs_NLL_lasso{1,sub}   = meanNLL;

    clear meanNLL allrepsLR allrepsNLL meanLR

end

%% plot learning rates

% loop over participants to plot learning rates
for sub = 1:nsubs
    sub_lr = subs_lr_lasso{1,sub};
    subnum = sprintf('%d',sub);

    % plot learning rates for each sub seperately
    [h, g, f] = plot_fitLRs(sub_lr,subnum);
end 
