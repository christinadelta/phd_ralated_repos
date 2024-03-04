% AL Simple RL model sensitivity analysis
% Created February 2024

%%%% Sensitivity Analysis with Parameter Recovery:   %%%%%%

% In sensitivity analysis with parameter recovery, we simulate data with 
% only varying one parameter at a time but also fit the model to the 
% simulated data generated with that varied parameter values. This involves 
% "freeing" one parameter (letting it vary across a range of values), 
% generating data based on those values, and then trying to recover those 
% values by fitting the model to the generated data. The other parameters 
% can be set to typical or baseline values during the data generation phase 
% but are also estimated during the fitting process to see how well the 
% model can recover the known (simulated) parameter values.

%%%%% Objective and Interpretation:  %%%%%

% The objective of varying one parameter at a time in sensitivity analysis 
% is to understand the influence of each parameter on the model's behavior 
% or performance. In parameter recovery, the goal is to assess the model's 
% ability to accurately estimate parameters from data, particularly how 
% the estimation of one parameter might be influenced by its true value 
% and the values of other parameters.

% If the model consistently recovers the correct values for the varied 
% parameter across its range, this suggests good parameter identifiability 
% and robustness. If recovery is poor, this might indicate issues such as 
% parameter correlations, insufficient sensitivity of the model outputs to 
% the parameter, or other model specification issues.

% ----------

%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
% datadir         = fullfile(startpath,'data');
% subs            = dir(fullfile(datadir, '*sub-*'));
% nsubs           = length(subs);
% 
% figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
total_trials    = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
nCues           = 2;
% params          = [0.3 1.5 0.1];

data            = ALsimdata_v2(probabilities, total_trials,condtrials);
x               = data.x; % extract true reward rate for plotting

%% start sensitivity analysis: ALPHA FREE

% how many repetitions?
repetitions = 100;

% define ranges for alpha the values of the fixed parameters  
alpha_range     = [0 1]; % range for alpha
beta_fix        = 1.5;
decay_fix       = 0.1;

for i       = 1:repetitions

    % Generate random true values for alpha and beta within the specified ranges
    true_alpha              = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand();
    allsim_alpha(i)         = true_alpha;

    % simulate data
    d                       = ALsimdata_v2(probabilities, total_trials,condtrials);
    out                     = modelRW_v2([true_alpha beta_fix decay_fix], d); % simulate model
    allsim_output{i}        = out;
    sim_actions             = out.simulated_actions;
    sim_reward              = out.reward;
    sim_outcome             = d.o;
    sim_pp(:,:,i)           = out.allPs;
    sim_Qs(:,:,i)           = out.Qvals;

    % generate starting points using a random generator for odd iterations
    % assuming known plausible ranges for alpha 
    alphaRange              = [0 1]; 
    starting_alpha          = alphaRange(1) + (alphaRange(2) - alphaRange(1)) * rand();
    fixedParams             = [beta_fix decay_fix];

    % define the objective function
    obFunc                  = @(x) lik_modelRW_v2_sensAnalysis(x(1),fixedParams,out,1);

    X0                      = starting_alpha;
    LB                      = 0;
    UB                      = 1;
    [Xfit, NegLL]           = fmincon(obFunc, X0, [], [], [], [], LB, UB);

    allfit_alpha(i)         = Xfit(1);
    allNLL(i)               = NegLL;

    % prepare simulated output and fitted values to fit the model to data
    dat.simulated_actions   = sim_actions;
    dat.outcome             = sim_outcome;
    params                  = [Xfit(1) beta_fix decay_fix];
    [fit_model_output, ~]   = fit_modelRW_v2(params, dat);
    fit_pp(:,:,i)           = fit_model_output.allPs;
    fit_Qs(:,:,i)           = fit_model_output.Qvals;
    
end % end of reps loop

%% plot correlations with alpha free

% plot simulated vs estimate parameter values
% calculate correlation coefficient for Alpha
alphaCorr       = corrcoef(allsim_alpha, allfit_alpha);
alphaCorrValue  = alphaCorr(1,2); % extract the correlation value

% average choice probabilties over 3rd dimension
avg_sim_pp      = mean(sim_pp, 3);
avg_fit_pp      = mean(fit_pp, 3);
ppCorr          = corrcoef(avg_sim_pp(:,1), avg_fit_pp(:,1));
ppCorrValue     = ppCorr(1,2);

% For Alpha values
figure; 
subplot(1,2,1); % alpha values
scatter(allsim_alpha, allfit_alpha, 'filled');
hold on; % keep the scatter plot

p   = polyfit(allsim_alpha, allfit_alpha, 1); % fit a linear polynomial (degree 1)
x1  = linspace(min(allsim_alpha), max(allsim_alpha), 100); % generate x values
y1  = polyval(p, x1); % calculate y values based on the fitted polynomial
plot(x1, y1, '-r'); % plot the line
legend(['Correlation: ', num2str(alphaCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Alpha'); ylabel('Fitted Alpha');
title('Simulated and Fitted Alpha Values');
set(gca, 'FontSize', 14); % additional plot formatting
grid on;

% For choice probabiltiies
subplot(1,2,2); % Beta values
scatter(avg_sim_pp(:,1), avg_fit_pp(:,1), 'filled');
hold on; % keep the scatter plot

p   = polyfit(avg_sim_pp(:,1), avg_fit_pp(:,1), 1); % fit a linear polynomial (degree 1)
x2  = linspace(min(avg_sim_pp(:,1)), max(avg_sim_pp(:,1)), 100); % generate x values
y2  = polyval(p, x2); % calculate y values based on the fitted polynomial
plot(x2, y2, '-r'); % plot the line
legend(['Correlation: ', num2str(ppCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Choice Probabilites'); ylabel('Fitted Choice Probabilities');
title('Simulated and Fitted choice Probs');
set(gca, 'FontSize', 14); % additional plot formatting
grid on; 

%% 2. run sensitivity analysis: BETA FREE

% how many repetitions?
repetitions     = 100;

% define ranges for alpha the values of the fixed parameters  
alpha_fix       = 0.3; % range for alpha
beta_range      = [0 5];
decay_fix       = 0.1;

for i       = 1:repetitions

    % Generate random true values for alpha and beta within the specified ranges
    true_beta               = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
    allsim_beta_m2(i)       = true_beta;

    % simulate data
    d                       = ALsimdata_v2(probabilities, total_trials,condtrials);
    out                     = modelRW_v2([alpha_fix true_beta decay_fix], d); % simulate model
    allsim_output_m2{i}     = out;
    sim_actions             = out.simulated_actions;
    sim_outcome             = out.outcome;
    sim_pp_m2(:,:,i)        = out.allPs;
    sim_Qs_m2(:,:,i)        = out.Qvals;

    % generate starting points using a random generator for odd iterations
    % assuming known plausible ranges for alpha 
    betaRange               = [0 5]; 
    startingBeta            = betaRange(1) + (betaRange(2) - betaRange(1)) * rand();
    fixedParams             = [alpha_fix decay_fix];

    % define the objective function
    obFunc                  = @(x) lik_modelRW_v2_sensAnalysis(x(1),fixedParams,out,2);

    X0                      = startingBeta;
    LB                      = 0;
    UB                      = 5;
    [Xfit, NegLL]           = fmincon(obFunc, X0, [], [], [], [], LB, UB);

    allfit_beta_m2(i)      = Xfit(1);
    allNLL_m2(i)            = NegLL;

    % prepare simulated output and fitted values to fit the model to data
    dat.simulated_actions   = sim_actions;
    dat.outcome             = sim_outcome;
    params                  = [alpha_fix Xfit(1) decay_fix];
    [fit_model_output, ~]   = fit_modelRW_v2(params, dat);
    fit_pp_m2(:,:,i)        = fit_model_output.allPs;
    fit_Qs_m2(:,:,i)        = fit_model_output.Qvals;
    
end % end of reps loop


%% plot correlations with beta free

% plot simulated vs estimate parameter values
% calculate correlation coefficient for Alpha
betaCorr       = corrcoef(allsim_beta_m2, allfit_beta_m2);
betaCorrValue  = betaCorr(1,2); % extract the correlation value

% average choice probabilties over 3rd dimension
avg_sim_pp_m2      = mean(sim_pp_m2, 3);
avg_fit_pp_m2      = mean(fit_pp_m2, 3);
ppCorr_m2          = corrcoef(avg_sim_pp_m2(:,1), avg_fit_pp_m2(:,1));
ppCorrValue_m2     = ppCorr_m2(1,2);

% For Alpha values
figure; 
subplot(1,2,1); % alpha values
scatter(allsim_beta_m2, allfit_beta_m2, 'filled');
hold on; % keep the scatter plot

p   = polyfit(allsim_beta_m2, allfit_beta_m2, 1); % fit a linear polynomial (degree 1)
x1  = linspace(min(allsim_beta_m2), max(allsim_beta_m2), 100); % generate x values
y1  = polyval(p, x1); % calculate y values based on the fitted polynomial
plot(x1, y1, '-r'); % plot the line
legend(['Correlation: ', num2str(betaCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Beta'); ylabel('Fitted Beta');
title('Simulated and Fitted Beta Values');
set(gca, 'FontSize', 14); % additional plot formatting
grid on;

% For choice probabiltiies
subplot(1,2,2); % Beta values
scatter(avg_sim_pp_m2(:,1), avg_fit_pp_m2(:,1), 'filled');
hold on; % keep the scatter plot

p   = polyfit(avg_sim_pp_m2(:,1), avg_fit_pp_m2(:,1), 1); % fit a linear polynomial (degree 1)
x2  = linspace(min(avg_sim_pp_m2(:,1)), max(avg_sim_pp_m2(:,1)), 100); % generate x values
y2  = polyval(p, x2); % calculate y values based on the fitted polynomial
plot(x2, y2, '-r'); % plot the line
legend(['Correlation: ', num2str(ppCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Choice Probabilites'); ylabel('Fitted Choice Probabilities');
title('Simulated and Fitted choice Probs');
set(gca, 'FontSize', 14); % additional plot formatting
grid on; 

%% 3. run sensitivity analysis: DECAY FREE

% how many repetitions?
repetitions     = 100;

% define ranges for alpha the values of the fixed parameters  
alpha_fix       = 0.3; % range for alpha
beta_fix        = 1.5;
decay_range     = [0 1];

for i       = 1:repetitions

    % Generate random true values for alpha and beta within the specified ranges
    true_decay              = decay_range(1) + (decay_range(2) - decay_range(1)) * rand();
    allsim_decay_m3(i)      = true_decay;

    % simulate data
    d                       = ALsimdata_v2(probabilities, total_trials,condtrials);
    out                     = modelRW_v2([alpha_fix beta_fix true_decay], d); % simulate model
    allsim_output_m3{i}     = out;
    sim_actions             = out.simulated_actions;
    sim_outcome             = out.outcome;
    sim_pp_m3(:,:,i)        = out.allPs;
    sim_Qs_m3(:,:,i)        = out.Qvals;

    % generate starting points using a random generator for odd iterations
    % assuming known plausible ranges for alpha 
    decayRange              = [0 1]; 
    starting_decay          = decayRange(1) + (decayRange(2) - decayRange(1)) * rand();
    fixedParams             = [alpha_fix beta_fix];

    % define the objective function
    obFunc                  = @(x) lik_modelRW_v2_sensAnalysis(x(1),fixedParams,out,3);

    X0                      = starting_decay;
    LB                      = 0;
    UB                      = 1;
    [Xfit, NegLL]           = fmincon(obFunc, X0, [], [], [], [], LB, UB);

    allfit_decay_m3(i)      = Xfit(1);
    allNLL_m3(i)            = NegLL;

    % prepare simulated output and fitted values to fit the model to data
    dat.simulated_actions   = sim_actions;
    dat.outcome             = sim_outcome;
    params                  = [alpha_fix beta_fix Xfit(1)];
    [fit_model_output, ~]   = fit_modelRW_v2(params, dat);
    fit_pp_m3(:,:,i)        = fit_model_output.allPs;
    fit_Qs_m3(:,:,i)        = fit_model_output.Qvals;
    
end % end of reps loop

%% plot correlations with decay free

% plot simulated vs estimate parameter values
% calculate correlation coefficient for Alpha
decayCorr       = corrcoef(allsim_decay_m3, allfit_decay_m3);
decayCorrValue  = decayCorr(1,2); % extract the correlation value

% average choice probabilties over 3rd dimension
avg_sim_pp_m3      = mean(sim_pp_m3, 3);
avg_fit_pp_m3      = mean(fit_pp_m3, 3);
ppCorr_m3          = corrcoef(avg_sim_pp_m3(:,1), avg_fit_pp_m3(:,1));
ppCorrValue_m3     = ppCorr_m3(1,2);

% For Alpha values
figure; 
subplot(1,2,1); % alpha values
scatter(allsim_decay_m3, allfit_decay_m3, 'filled');
hold on; % keep the scatter plot

p   = polyfit(allsim_decay_m3, allfit_decay_m3, 1); % fit a linear polynomial (degree 1)
x1  = linspace(min(allsim_decay_m3), max(allsim_decay_m3), 100); % generate x values
y1  = polyval(p, x1); % calculate y values based on the fitted polynomial
plot(x1, y1, '-r'); % plot the line
legend(['Correlation: ', num2str(decayCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Decay'); ylabel('Fitted Decay');
title('Simulated and Fitted Decay Values');
set(gca, 'FontSize', 14); % additional plot formatting
grid on;

% For choice probabiltiies
subplot(1,2,2); % Beta values
scatter(avg_sim_pp_m3(:,1), avg_fit_pp_m3(:,1), 'filled');
hold on; % keep the scatter plot

p   = polyfit(avg_sim_pp_m3(:,1), avg_fit_pp_m3(:,1), 1); % fit a linear polynomial (degree 1)
x2  = linspace(min(avg_sim_pp_m3(:,1)), max(avg_sim_pp_m3(:,1)), 100); % generate x values
y2  = polyval(p, x2); % calculate y values based on the fitted polynomial
plot(x2, y2, '-r'); % plot the line
legend(['Correlation: ', num2str(ppCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Choice Probabilites'); ylabel('Fitted Choice Probabilities');
title('Simulated and Fitted choice Probs');
set(gca, 'FontSize', 14); % additional plot formatting
grid on; 

%%
