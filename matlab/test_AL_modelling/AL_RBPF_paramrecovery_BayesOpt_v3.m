% AL PARAMETER RECOVERY USING THE BAYESOPT TOOLBOX
% Created February 2024 -- version 3


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
    output                      = simModel(data, probabilities, trials,condtrials,parameters);

     % store simulated X 
    ridge_simX(1,rep)           = sim_lambda_s;
    ridge_simX(2,rep)           = sim_lambda_v;
    ridge_simX(3,rep)           = sim_beta;
    ridge_simX(4,rep)           = sim_s0;
    ridge_simX(5,rep)           = sim_v0;
    ridge_sim_lr(:,:,rep)       = output.lrs;
    ridge_sim_action(:,rep)     = output.test_a;
    ridge_sim_bo(:,:,rep)       = output.binary_o;
    ridge_sim_o(:,:,rep)        = output.test_o;

    % prepare output from simulations that are needed for fitting 
    state                   = data.x;
    state(:,2)              = 1 - state(:,1);
    ss                      = data.stcind;
    vv                      = data.t;
    test_o                  = output.test_o;
    test_a                  = output.test_a;
    test_a(find(test_a==0)) = 2; % convert action 0 to 2
    config                  = struct('tvol',vv,'tstc',ss,'state',state,'nsim',1);

    % randomly generate values for the free parameterers of the full model:
    % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
    fit_lambda              = 0.1 + (0.9-0.1).* rand(2,1); 
    fit_lambda_s            = fit_lambda(1);
    fit_lambda_v            = fit_lambda(2);
    
    % 2. softmax temperatyre beta -- sould be between 0.1 and 10
    fit_beta                = 0.1 + (5-0.1).* rand(1,1); 
    
    % 3. s0 and v0 -- should be between 1-3 and 1-2
    init_vals               = 0.1 + (3-0.1).* rand(2,1); 
    fit_s0                  = init_vals(1);
    fit_v0                  = init_vals(2);

    lambda                  = 1e-4 + (1e4 - 1e-4).* rand(1,1); 
    otherParams             = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001]; 

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
    ridge_fitX(1,rep)       = fitParams(1);
    ridge_fitX(2,rep)       = fitParams(2);
    ridge_fitX(3,rep)       = fitParams(3);
    ridge_fitX(4,rep)       = fitParams(4);
    ridge_fitX(5,rep)       = fitParams(5);
    ridge_fitLambda(1,rep)  = bestLambdaReg;
    ridge_fitNLL(1,rep)     = minNLL;
    ridge_fit_lr(:,:,rep)   = lr;
    ridge_fit_choice(:,rep) = choice;
    ridge_fit_vol(:,:,rep)  = vol;
    ridge_fit_val(:,:,rep)  = val;
    ridge_fit_stc(:,:,rep)  = unp;

    clear test_a test_o fitParams minNLL lr choice val vol unp bestLambdaReg H output state
    clear fit_lambda_s fit_lambda_v fit_beta fit_s0 fit_v0 sim_v0 sim_s0 sim_lambda_v sim_lambda_s sim_beta fit_lambda
    clear initialGuessTable params results fitParams bestParams best

end % end of repetitions loop


%% plot ridge fit output

% Number of parameters to plot
numParams = 5;

% Create a figure for subplots
figure;

% Titles for each subplot
paramTitles = {'lambda_s', 'lambda_v', 'beta', 's0', 'v0'};

% Loop through each parameter
for i = 1:numParams
    % Select subplot position
    subplot(2, 3, i);
    
    % Extract simulated and fitted values for the current parameter
    sim_values = ridge_simX(i, :);
    fit_values = ridge_fitX(i, :);
    
    % Scatter plot for the current parameter
    scatter(sim_values, fit_values, 'filled');
    hold on; % Keep the scatter plot visible when adding the line of best fit
    
    % Calculate and plot line of best fit
    coeffs = polyfit(sim_values, fit_values, 1); % Linear fit
    fittedX = linspace(min(sim_values), max(sim_values), 200); % Generating X values for the fit line
    fittedY = polyval(coeffs, fittedX); % Y values for the line
    plot(fittedX, fittedY, 'r-', 'LineWidth', 1); % Plot line of best fit in red
    
    % Formatting the subplot
    xlabel(['Simulated ' paramTitles{i}]);
    ylabel(['Fitted ' paramTitles{i}]);
    title(['Scatter plot of Simulated vs. Fitted ' paramTitles{i}]);
    
    % Calculate and display correlation coefficient
    [R, P] = corrcoef(sim_values, fit_values);
    corrText = sprintf('R = %.2f, p = %.3f', R(2), P(2));
    legend('Data', ['Fit: ' corrText], 'Location', 'best');
    
    hold off; % Release the plot for the next subplot
end

% Adjust layout to make everything fit nicely
set(gcf, 'Position', [100, 100, 1200, 600]); % Resize figure to make it wider


%%




