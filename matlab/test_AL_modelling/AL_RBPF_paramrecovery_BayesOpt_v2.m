% AL PARAMETER RECOVERY USING THE BAYESOPT TOOLBOX
% Created February 2024

% Parameter recovery is done using BayesOpt with:
% 1. Ridge Regularisation for better fitting and perfomance of the model
% with the assumption that all parameters contribute equally to the model's
% predictions 
% 2. Lasso Regularisation for better fitting and perfomance of the model
% with the assumption that there are subsets of paraeteres that do a better
% job at predicting the data

% Hassian matrix is also computed for all parameter recovery
% regularisations 

%% clear all
clc
clear all

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

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
true_stc        = [0.1 1 2;
                    1 2 3];
true_vol        = [0.1 1;
                    1 2];

data            = ALsimdata_v2(probabilities, trials,condtrials);

%% run parameter recovery (full model) - with ridge regularisation (hessian matrix included)

% how many reps?
reps = 50; 

for rep = 1:reps

    % set the random seed for reproducibility
    rng(rep); % Each iteration uses a different seed

    for i = 1:3
        for j = 1:2

            % PART 1 SIMULATING
            % define model parameters  
            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            sim_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            sim_lambda_s    = sim_lambda(1);
            sim_lambda_v    = sim_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            sim_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            init_vals   = 0.1 + (3-0.1).* rand(2,1); 
            sim_s0          = init_vals(1);
            sim_v0          = init_vals(2);
            parameters      = struct('nparticles',500,'x0_unc',1,'lambda_v',sim_lambda_v,'lambda_s',sim_lambda_s,'v0',sim_v0,'s0',sim_s0, 'beta',sim_beta);
    
            % run pf_core function which includes both the kalman and particle filters 
            output          = simModel_v2(data, probabilities, trials,condtrials,parameters,i,j);

            % store simulated X 
            ridge_simX{rep}{i,j}(1)         = sim_lambda_s;
            ridge_simX{rep}{i,j}(2)         = sim_lambda_v;
            ridge_simX{rep}{i,j}(3)         = sim_beta;
            ridge_simX{rep}{i,j}(4)         = sim_s0;
            ridge_simX{rep}{i,j}(5)         = sim_v0;
            ridge_sim_lr{1,rep}{i,j}        = output.lr;
            ridge_sim_action{1,rep}{i,j}    = output.a;
            ridge_sim_bo{1,rep}{i,j}        = output.binary_o;
            ridge_sim_o{1,rep}{i,j}         = output.o;
            
            % PART 2 FITTING
            % extract what we need and prepare data for model fititng 
            state(:,1)              = data.x;
            state(:,2)              = 1 - state(:,1);
            ss                      = data.stcind;
            vv                      = data.t;
            test_o                  = output.o;
            test_a                  = output.a;
            test_a(find(test_a==0)) = 2; % convert action 0 to 2
            config                  = struct('tvol',vv,'tstc',ss,'state',state,'nsim',1);

            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            fit_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            fit_lambda_s    = fit_lambda(1);
            fit_lambda_v    = fit_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            fit_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            init_vals       = 0.1 + (3-0.1).* rand(2,1); 
            fit_s0          = init_vals(1);
            fit_v0          = init_vals(2);

            lambda          = 1e-4 + (1e4 - 1e-4).* rand(1,1); 
            otherParams     = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed

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

            % Define the objective function as an anonymous function if needed
            objectiveWithInputs = @(params) bayesOptObjective(params, test_o,test_a, otherParams, config);

            % Creating the table with initial guesses
            initialGuessTable = table([fit_lambda_s], [fit_lambda_v], [fit_beta], [fit_s0], [fit_v0], [lambda], ...
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
            fitParams(6)        = bestParams.lambda(1);
            bestLambdaReg(1)    = bestParams.lambda(1);

            % now run the model with the optimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(fitParams,test_o,otherParams,config);

            % store fitted parameter values
            ridge_fitX{rep}{i,j}(1)       = fitParams(1);
            ridge_fitX{rep}{i,j}(2)       = fitParams(2);
            ridge_fitX{rep}{i,j}(3)       = fitParams(3);
            ridge_fitX{rep}{i,j}(4)       = fitParams(4);
            ridge_fitX{rep}{i,j}(5)       = fitParams(5);
            ridge_fitLambda{1,rep}(i,j)   = bestLambdaReg;
            ridge_fitNLL{1,rep}(i,j)      = minNLL;
            ridge_fit_lr{1,rep}{i,j}      = lr;
            ridge_fit_choice{1,rep}{i,j}  = choice;
            ridge_fit_vol{1,rep}{i,j}     = vol;
            ridge_fit_val{1,rep}{i,j}     = val;
            ridge_fit_stc{1,rep}{i,j}     = unp;
            % ridge_allHessian{i,j}(:,:,rep)= H;

            clear test_a test_o fitParams minNLL lr choice val vol unp bestLambdaReg H output state
            clear fit_lambda_s fit_lambda_v fit_beta fit_s0 fit_v0 sim_v0 sim_s0 sim_lambda_v sim_lambda_s sim_beta fit_lambda
            clear initialGuessTable params results fitParams bestParams best

        end % end of volatility loop
    end % end of stochasticity 
end % end of reps loop

%% run parameter recovery (full model) - with lasso regularisation (hessian matrix included)

% how many reps?
reps = 50; 

for rep = 1:reps

    % set the random seed for reproducibility
    rng(rep); % Each iteration uses a different seed

    for i = 1:3
        for j = 1:2

            % PART 1 SIMULATING
            % define model parameters  
            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            sim_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            sim_lambda_s    = sim_lambda(1);
            sim_lambda_v    = sim_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            sim_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            init_vals   = 0.1 + (3-0.1).* rand(2,1); 
            sim_s0          = init_vals(1);
            sim_v0          = init_vals(2);
            parameters      = struct('nparticles',500,'x0_unc',1,'lambda_v',sim_lambda_v,'lambda_s',sim_lambda_s,'v0',sim_v0,'s0',sim_s0, 'beta',sim_beta);
    
            % run pf_core function which includes both the kalman and particle filters 
            output          = simModel_v2(data, probabilities, trials,condtrials,parameters,i,j);

            % store simulated X 
            lasso_simX{rep}{i,j}(1)         = sim_lambda_s;
            lasso_simX{rep}{i,j}(2)         = sim_lambda_v;
            lasso_simX{rep}{i,j}(3)         = sim_beta;
            lasso_simX{rep}{i,j}(4)         = sim_s0;
            lasso_simX{rep}{i,j}(5)         = sim_v0;
            lasso_sim_lr{1,rep}{i,j}        = output.lr;
            lasso_sim_action{1,rep}{i,j}    = output.a;
            lasso_sim_bo{1,rep}{i,j}        = output.binary_o;
            lasso_sim_o{1,rep}{i,j}         = output.o;
            
            % PART 2 FITTING
            % extract what we need and prepare data for model fititng 
            state(:,1)              = output.x;
            state(:,2)              = 1 - state(:,1);
            ss                      = data.stcind;
            vv                      = data.t;
            test_o                  = output.o;
            test_a                  = output.a;
            test_a(find(test_a==0)) = 2; % convert action 0 to 2
            config                  = struct('tvol',vv,'tstc',ss,'state',state,'nsim',1);

            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            fit_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            fit_lambda_s    = fit_lambda(1);
            fit_lambda_v    = fit_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            fit_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            init_vals       = 0.1 + (3-0.1).* rand(2,1); 
            fit_s0          = init_vals(1);
            fit_v0          = init_vals(2);

            lambda          = 1e-4 + (1e4 - 1e-4).* rand(1,1); 
            otherParams     = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed

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

            % Define the objective function as an anonymous function if needed
            objectiveWithInputs = @(params) bayesOptObjectiveLasso(params, test_o,test_a, otherParams, config);

            % Creating the table with initial guesses
            initialGuessTable = table([fit_lambda_s], [fit_lambda_v], [fit_beta], [fit_s0], [fit_v0], [lambda], ...
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
            fitParams(6)        = bestParams.lambda(1);
            bestLambdaReg(1)    = bestParams.lambda(1);

            % now run the model with the optimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(fitParams,test_o,otherParams,config);

            % store fitted parameter values
            lasso_fitX{rep}{i,j}(1)       = fitParams(1);
            lasso_fitX{rep}{i,j}(2)       = fitParams(2);
            lasso_fitX{rep}{i,j}(3)       = fitParams(3);
            lasso_fitX{rep}{i,j}(4)       = fitParams(4);
            lasso_fitX{rep}{i,j}(5)       = fitParams(5);
            lasso_fitLambda{1,rep}(i,j)   = bestLambdaReg;
            lasso_fitNLL{1,rep}(i,j)      = minNLL;
            lasso_fit_lr{1,rep}{i,j}      = lr;
            lasso_fit_choice{1,rep}{i,j}  = choice;
            lasso_fit_vol{1,rep}{i,j}     = vol;
            lasso_fit_val{1,rep}{i,j}     = val;
            lasso_fit_stc{1,rep}{i,j}     = unp;
            % lasso_allHessian{i,j}(:,:,rep)= H;

            clear test_a test_o fitParams minNLL lr choice val vol unp bestLambdaReg H output state
            clear fit_lambda_s fit_lambda_v fit_beta fit_s0 fit_v0 sim_v0 sim_s0 sim_lambda_v sim_lambda_s sim_beta fit_lambda
            clear initialGuessTable params results fitParams bestParams best
        end % end of volatility loop
    end % end of stochasticity 
end % end of reps loop

%% plot parameter recovery (full model)

% first re-arrange the parameters 
[Xfit, fParams, sParams] = rearrangeX_v2(lasso_simX, lasso_fitX,reps);

%% plot simulated vs fitted parameters for each parameter seperately 
% plot lambdas_s
y           = fParams.lambda_s;
x           = sParams.lambda_s;
ylm         = [0 1];
colours     = 'magenta'; % let's start with magenta?
param_title = 'lambda_s';
a           = plotScatter(y, x, ylm, colours,param_title);

%% plot lambda_v 

y           = fParams.lambda_v;
x           = sParams.lambda_v;
ylm         = [0 1];
colours     = "red"; % let's start with magenta?
param_title = 'lambda_v';
b           = plotScatter(y, x, ylm, colours,param_title);

%% plot beta 

y           = fParams.betas;
x           = sParams.betas;
ylm         = [0 10];
colours     = "blue"; % let's start with magenta?
param_title = 'beta';
c           = plotScatter(y, x, ylm, colours,param_title);

%% plot s0 

y           = fParams.s0;
x           = sParams.s0;
ylm         = [0 2];
colours     = "cyan"; % let's start with magenta?
param_title = 's_0';
d           = plotScatter(y, x, ylm, colours,param_title);

%% plot v0 

y           = fParams.v0;
x           = sParams.v0;
ylm         = [0 2];
colours     = "green"; % let's start with magenta?
param_title = 'v_0';
e           = plotScatter(y, x, ylm, colours,param_title);


%% %% plot correlations of simulated vs recovered learning rates

[lr_timeseries,pvals,f] = plot_corr_lrs(lasso_fit_lr, lasso_sim_lr);

%% run parameter recovery (lambdas_beta model) - with ridge regularisation (hessian matrix included)

% how many reps?
reps = 50; 

for rep = 1:reps

    % set the random seed for reproducibility
    rng(rep); % Each iteration uses a different seed

    for i = 1:3
        for j = 1:2

            % PART 1 SIMULATING
            % define model parameters  
            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            sim_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            sim_lambda_s    = sim_lambda(1);
            sim_lambda_v    = sim_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            sim_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            sim_s0          = 1;
            sim_v0          = 1;
            parameters      = struct('nparticles',500,'x0_unc',1,'lambda_v',sim_lambda_v,'lambda_s',sim_lambda_s,'v0',sim_v0,'s0',sim_s0, 'beta',sim_beta);
    
            % run pf_core function which includes both the kalman and particle filters 
            output          = simModel_v2(data, probabilities, trials,condtrials,parameters,i,j);

            % store simulated X 
            ridge_simX{rep}{i,j}(1)         = sim_lambda_s;
            ridge_simX{rep}{i,j}(2)         = sim_lambda_v;
            ridge_simX{rep}{i,j}(3)         = sim_beta;
            ridge_simX{rep}{i,j}(4)         = sim_s0;
            ridge_simX{rep}{i,j}(5)         = sim_v0;
            ridge_sim_lr{1,rep}{i,j}        = output.lr;
            ridge_sim_action{1,rep}{i,j}    = output.a;
            ridge_sim_bo{1,rep}{i,j}        = output.binary_o;
            ridge_sim_o{1,rep}{i,j}         = output.o;
            
            % PART 2 FITTING
            % extract what we need and prepare data for model fititng 
            state(:,1)              = output.x;
            state(:,2)              = 1 - state(:,1);
            ss                      = data.stcind;
            vv                      = data.t;
            test_o                  = output.o;
            test_a                  = output.a;
            test_a(find(test_a==0)) = 2; % convert action 0 to 2
            config                  = struct('tvol',vv,'tstc',ss,'state',state,'nsim',1);

            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            fit_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            fit_lambda_s    = fit_lambda(1);
            fit_lambda_v    = fit_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            fit_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            fit_s0          = 1;
            fit_v0          = 1;

            lambda          = 1e-4 + (1e4 - 1e-4).* rand(1,1); 
            otherParams     = [0.3 0.6 0.5 1 1 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed

            % set up the BayesOpt configuration
            % Define the parameter space
            params = [
                optimizableVariable('lambda_s',[0.1, 0.9],'Type','real');
                optimizableVariable('lambda_v',[0.1, 0.9],'Type','real');
                optimizableVariable('beta',[0.1, 5],'Type','real');
                optimizableVariable('lambda', [1e-4, 1e4], 'Transform', 'log'); 
                % 
            ];

            % Define the objective function as an anonymous function if needed
            objectiveWithInputs = @(params) bayesOptObjective(params, test_o,test_a, otherParams, config);

            % Creating the table with initial guesses
            initialGuessTable = table([fit_lambda_s], [fit_lambda_v], [fit_beta], [lambda], ...
                'VariableNames', {'lambda_s', 'lambda_v', 'beta','lambda'});
    
    
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
            fitParams(6)        = bestParams.lambda(1);
            bestLambdaReg(1)    = bestParams.lambda(1);

            % now run the model with the optimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_lambdaBetab_full(fitParams,test_o,otherParams,config);

            % store fitted parameter values
            ridge_fitX{rep}{i,j}(1)       = fitParams(1);
            ridge_fitX{rep}{i,j}(2)       = fitParams(2);
            ridge_fitX{rep}{i,j}(3)       = fitParams(3);
            ridge_fitLambda{1,rep}(i,j)   = bestLambdaReg;
            ridge_fitNLL{1,rep}(i,j)      = minNLL;
            ridge_fit_lr{1,rep}{i,j}      = lr;
            ridge_fit_choice{1,rep}{i,j}  = choice;
            ridge_fit_vol{1,rep}{i,j}     = vol;
            ridge_fit_val{1,rep}{i,j}     = val;
            ridge_fit_stc{1,rep}{i,j}     = unp;
            % ridge_allHessian{i,j}(:,:,rep)= H;

            clear test_a test_o fitParams minNLL lr choice val vol unp bestLambdaReg H output state
            clear fit_lambda_s fit_lambda_v fit_beta fit_s0 fit_v0 sim_v0 sim_s0 sim_lambda_v sim_lambda_s sim_beta fit_lambda
            clear initialGuessTable params results fitParams bestParams best

        end % end of volatility loop
    end % end of stochasticity 
end % end of reps loop

%% run parameter recovery (lambda_beta model) - with lasso regularisation (hessian matrix included)

% how many reps?
reps = 50; 

for rep = 1:reps

    % set the random seed for reproducibility
    rng(rep); % Each iteration uses a different seed

    for i = 1:3
        for j = 1:2

            % PART 1 SIMULATING
            % define model parameters  
            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            sim_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            sim_lambda_s    = sim_lambda(1);
            sim_lambda_v    = sim_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            sim_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            sim_s0          = 1;
            sim_v0          = 1;
            parameters      = struct('nparticles',500,'x0_unc',1,'lambda_v',sim_lambda_v,'lambda_s',sim_lambda_s,'v0',sim_v0,'s0',sim_s0, 'beta',sim_beta);
    
            % run pf_core function which includes both the kalman and particle filters 
            output          = simModel_v2(data, probabilities, trials,condtrials,parameters,i,j);

            % store simulated X 
            lasso_simX{rep}{i,j}(1)         = sim_lambda_s;
            lasso_simX{rep}{i,j}(2)         = sim_lambda_v;
            lasso_simX{rep}{i,j}(3)         = sim_beta;
            lasso_simX{rep}{i,j}(4)         = sim_s0;
            lasso_simX{rep}{i,j}(5)         = sim_v0;
            lasso_sim_lr{1,rep}{i,j}        = output.lr;
            lasso_sim_action{1,rep}{i,j}    = output.a;
            lasso_sim_bo{1,rep}{i,j}        = output.binary_o;
            lasso_sim_o{1,rep}{i,j}         = output.o;
            
            % PART 2 FITTING
            % extract what we need and prepare data for model fititng 
            state(:,1)              = output.x;
            state(:,2)              = 1 - state(:,1);
            ss                      = data.stcind;
            vv                      = data.t;
            test_o                  = output.o;
            test_a                  = output.a;
            test_a(find(test_a==0)) = 2; % convert action 0 to 2
            config                  = struct('tvol',vv,'tstc',ss,'state',state,'nsim',1);

            % randomly generate values for the free parameterers of the full model:
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            fit_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            fit_lambda_s    = fit_lambda(1);
            fit_lambda_v    = fit_lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            fit_beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 1-3 and 1-2
            fit_s0          = 1;
            fit_v0          = 1;

            lambda          = 1e-4 + (1e4 - 1e-4).* rand(1,1); 
            otherParams     = [0.3 0.6 0.5 1 1 100 1 1 0.001];  % parameters 6-9 are also used for modelling, are fixed

            % set up the BayesOpt configuration
            % Define the parameter space
            params = [
                optimizableVariable('lambda_s',[0.1, 0.9],'Type','real');
                optimizableVariable('lambda_v',[0.1, 0.9],'Type','real');
                optimizableVariable('beta',[0.1, 5],'Type','real');
                optimizableVariable('lambda', [1e-4, 1e4], 'Transform', 'log'); 
                % 
            ];

            % Define the objective function as an anonymous function if needed
            objectiveWithInputs = @(params) bayesOptObjectiveLasso(params, test_o,test_a, otherParams, config);

            % Creating the table with initial guesses
            initialGuessTable = table([fit_lambda_s], [fit_lambda_v], [fit_beta], [lambda], ...
                'VariableNames', {'lambda_s', 'lambda_v', 'beta', 'lambda'});
    
    
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
            fitParams(4)        = fit_s0;
            fitParams(5)        = fit_v0;
            fitParams(6)        = bestParams.lambda(1);
            bestLambdaReg(1)    = bestParams.lambda(1);

            % now run the model with the optimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_lambdaBetab_full(fitParams,test_o,otherParams,config);

            % store fitted parameter values
            lasso_fitX{rep}{i,j}(1)       = fitParams(1);
            lasso_fitX{rep}{i,j}(2)       = fitParams(2);
            lasso_fitX{rep}{i,j}(3)       = fitParams(3);
            lasso_fitLambda{1,rep}(i,j)   = bestLambdaReg;
            lasso_fitNLL{1,rep}(i,j)      = minNLL;
            lasso_fit_lr{1,rep}{i,j}      = lr;
            lasso_fit_choice{1,rep}{i,j}  = choice;
            lasso_fit_vol{1,rep}{i,j}     = vol;
            lasso_fit_val{1,rep}{i,j}     = val;
            lasso_fit_stc{1,rep}{i,j}     = unp;
            % lasso_allHessian{i,j}(:,:,rep)= H;

            clear test_a test_o fitParams minNLL lr choice val vol unp bestLambdaReg H output state
            clear fit_lambda_s fit_lambda_v fit_beta fit_s0 fit_v0 sim_v0 sim_s0 sim_lambda_v sim_lambda_s sim_beta fit_lambda
            clear initialGuessTable params results fitParams bestParams best
        end % end of volatility loop
    end % end of stochasticity 
end % end of reps loop

%% plot parameter recovery (full model)

% first re-arrange the parameters 
[Xfit, fParams, sParams] = rearrangeX_v2(lasso_simX, lasso_fitX,reps);

%% plot simulated vs fitted parameters for each parameter seperately 
% plot lambdas_s
y           = fParams.lambda_s;
x           = sParams.lambda_s;
ylm         = [0 1];
colours     = 'magenta'; % let's start with magenta?
param_title = 'lambda_s';
a           = plotScatter(y, x, ylm, colours,param_title);

%% plot lambda_v 

y           = fParams.lambda_v;
x           = sParams.lambda_v;
ylm         = [0 1];
colours     = "red"; % let's start with magenta?
param_title = 'lambda_v';
b           = plotScatter(y, x, ylm, colours,param_title);

%% plot beta 

y           = fParams.betas;
x           = sParams.betas;
ylm         = [0 10];
colours     = "blue"; % let's start with magenta?
param_title = 'beta';
c           = plotScatter(y, x, ylm, colours,param_title);


%% %% plot correlations of simulated vs recovered learning rates

[lr_timeseries,pvals,f] = plot_corr_lrs(ridge_fit_lr, ridge_sim_lr);
