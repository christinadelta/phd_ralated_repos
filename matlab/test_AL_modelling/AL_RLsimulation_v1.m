% RW MODEL SIMULATIONS AND FITTING USING BAYESOPT 

% CREATE JANUARY 2024 FOR TESTING PURPOSES ONLY 

% Here, I test the AL task with a simple RW model with two parameters:
% alpha - learning rate
% beta - inverse temperature for the softmax funcyion

% ---------
% ANALYSIS AND DIAGNOSTICS RAN:
% 1. Simulate Aversive Learning data (outcomes, probabilities)
% 2. Simulate Rescorla Wagner model using two parameters: alpha and beta
% (choice probabilities, Q values)
% 3. Perform Grid search with a range of alpha and beta combinations and 
% simulated data set to explore parameter space 
% 4. run parameter recovery with a focus on the good alpha and beta
% combinations obtained in Grid Search
% 5. Fit model participant data using Grid Search approach 

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
params          = [0.3 1.5];

data            = ALsimdata_v2(probabilities, total_trials,condtrials);

%% simulate model 

modelout = modelRW_v1(params, data);

%% extract model output

% for plotting we we need Qvalues, choice probabilities and choices both conditions, so, concatinate the
% trials of the two conditions
Qvals   = modelout.Qvals;
Ps      = modelout.allPs; %
choices = modelout.a;
a       = 2 - choices; % convert to be 1 and 0
x       = data.x;
correct = modelout.correct;

% get score
score   = mean(correct);

%% plot RL model

h               = plotRW_v1(x,Ps,Qvals,a);

%% PERFORM SEARCH GRID %%

% WHY PERFORM GRID SEARCH WITH SIMULATED DATA?
% To explore the parameter space and understand how changes in 
% parameters affect model predictions. This can help identify plausible 
% ranges for parameters and reveal how sensitive the model's outcomes are 
% to different parameter values.
% This preliminary exploration can inform you about the behavior of your 
% model under controlled conditions, where you know the "true" parameters 
% used to generate the data. It's a way to familiarize yourself with the 
% model's dynamics and to ensure that the grid search process and your 
% model implementation are working as expected.

% -------
% initial variables and parameters
init_alpha      = 0.3;
init_beta       = 1.5;
bounds          = [0 1;   % alpha range
                    0 5]; % beta range
bins            = [20 25];
nparam          = length(params);

simalpha_bin    = init_alpha/bounds(1,2)*bins(1);
simbeta_bin     = init_beta/bounds(2,2)*bins(2);
best_negLL      = inf; % initialize with a large number

% create the parameter space for the grid search
for iParam = 1:nparam
    range       = linspace(bounds(iParam,1),bounds(iParam,2),bins(iParam)+1);
    p{iParam}   = range(2:end); % stay just off the zero bounds
end

% estimate nll using different combinations of parameter values 
params          = nan(1,2);
for t = 1:bins(1)
    params(1)   = p{1}(t);

    for tt = 1:bins(2)
        params(2)           = p{2}(tt);
        [~,negLL]           = fit_modelRW_v1(params,modelout);

        % store this iteration's nll
        nll(t,tt)           = negLL;

        % update best parameters if current nll is lower
        if negLL < best_negLL
            best_negLL  = negLL;
            best_alpha  = params(1);
            best_beta   = params(2);
        end
    end
end

% display best parameter values 
fprintf('Best alpha: %f, Best beta: %f, with negative LL: %f\n', best_alpha, best_beta, best_negLL);

% for visualisation, compute the likelihood instaed of log
lik = exp(nll); 

% visualise the nll space using heatmap
figure; % create a new figure
imagesc(p{1,2}, p{1,1}, nll); % create heatmap
colorbar; % show color scale
xlabel('\beta (Softmax Inverse Temperature)'); % x-axis label
ylabel('\alpha (Learning Rate)'); % y-axis label
title('Grid Search Results: Negative Log-Likelihood'); % Title
fontsize(gcf,16,"points")
set(gca, 'YDir', 'normal'); % Ensure the y-axis starts from the bottom

% Add a marker at the best alpha-beta point
hold on; % Keep the heatmap visible
% Plot the marker using the best parameter values directly
plot(best_beta, best_alpha, 'r*', 'MarkerSize', 30); % Red star marker for best parameters
hold off;

% fit the model with the best parameter values from the search grid 
best_params                 = [best_alpha best_beta];
[grid_model_output, ~]      = fit_modelRW_v1(best_params, modelout);


%% estimate marginal likelihoods of parameters alpha and beta

% Compute marginal likelihoods because they can provide a more nuanced 
% understanding of the parameter space and how individual parameters 
% contribute to model performance. Why considering marginal likelihoods can be beneficial:

% 1. Marginal likelihoods are computed by integrating (or summing) the likelihood over all 
% possible values of the other parameters, holding one parameter of interest fixed at a time. 
% This process effectively averages the performance of the model across a range of 
% scenarios for the non-fixed parameters, providing a clearer picture of how a 
% single parameter influences model fit independent of specific values of other parameters. 

% 2. by examining the marginal likelihoods, you can identify which parameters have a 
% strong influence on model performance and which do not. A parameter with a marginal 
% likelihood that sharply peaks at certain values suggests that the model is sensitive 
% to that parameter—meaning accurate estimation is crucial for good model fit. Conversely, 
% a flat marginal likelihood indicates that the model's performance is robust to changes in 
% that parameter.

% How to get Marginal NLLs:
% Marginal NLL for α: To find the best marginal α, look for the α value 
% that corresponds to the minimum of the marginal NLL calculated across all β values

% Marginal NLL for β: Similarly, to find the best marginal β, identify the β 
% value that corresponds to the minimum of the marginal NLL calculated across all α values.

% Interpretation:
% Best Marginal α: This is the value of α that, on average, best fits the 
% data when considering all possible values of β. It represents a "generalized" 
% optimal learning rate across the range of decision noise levels modeled by β.

% Best Marginal β: Conversely, this is the value of β that best fits the data 
% on average across all considered values of α, representing an optimal level 
% of decision noise or temperature across different learning rates.


% -------
% convert NLL to likelihoods for marginalization purposes
likelihoods                 = exp(-nll);

% marginalize (sum) likelihoods over beta to get marginal likelihood for alpha
marginal_likelihood_alpha   = sum(likelihoods, 2); % Sum over columns (beta)

% marginalize (sum) likelihoods over alpha to get marginal likelihood for beta
marginal_likelihood_beta    = sum(likelihoods, 1); % Sum over rows (alpha)

% optionally, convert marginal likelihoods back to NLL for visualization or analysis
marginal_nll_alpha          = -log(marginal_likelihood_alpha);
marginal_nll_beta           = -log(marginal_likelihood_beta);

% visualize Marginal NLLs
% plot marginal NLL for alpha
figure;
plot(p{1}, marginal_nll_alpha, 'LineWidth', 2);
xlabel('\alpha (Learning Rate)');
ylabel('Marginal Negative Log-Likelihood');
title('Marginal NLL for Alpha');
fontsize(gcf,16,"points")
grid on;

% plot marginal NLL for beta
figure;
plot(p{2}, marginal_nll_beta, 'LineWidth', 2);
xlabel('\beta (Softmax Inverse Temperature)');
ylabel('Marginal Negative Log-Likelihood');
title('Marginal NLL for Beta');
fontsize(gcf,16,"points")
grid on;

% estimate best marginal parameter values
[~, idx_min_alpha] = min(marginal_nll_alpha);
[~, idx_min_beta] = min(marginal_nll_beta);

% Retrieve the best marginal values for alpha and beta from their respective ranges
best_marginal_alpha = p{1}(idx_min_alpha);
best_marginal_beta = p{2}(idx_min_beta);

fprintf('Best marginal alpha: %f\n', best_marginal_alpha);
fprintf('Best marginal beta: %f\n', best_marginal_beta);

%% estimate expected values of the parameter space 

% The expected value of a random variable is a fundamental concept in 
% probability theory, representing the long-run average value of repetitions 
% of the experiment it represents. It is calculated as the sum of possible 
% outcomes weighted by their probabilities.
% In statistical modeling, when parameters are treated probabilistically, 
% the expected value of a parameter can be considered as an average value of the parameter,\
% weighted by the probability of each value occurring. This is akin to a "mean" parameter 
% value in light of the model's uncertainty about that parameter.

% Adapting to our approach
% our goal is to compute an expected value of α and β that reflects a weighted 
% average based on their contribution to the model's fit, so, we would first need 
% to convert the NLLs into a form that can act as weights in a meaningful 
% way (typically, probabilities).

% computing expected values in this manner can provide an insightful summary of 
% the parameter space, reflecting not just the single best-fitting value but 
% incorporating the entire distribution of how well different values fit the data.

% convert marginal NLLs to probabilities (for illustrative purposes)
% this is a simplified approach...
prob_alpha                  = exp(-marginal_nll_alpha - min(marginal_nll_alpha));
prob_alpha                  = prob_alpha / sum(prob_alpha);
prob_beta                   = exp(-marginal_nll_beta - min(marginal_nll_beta));
prob_beta                   = prob_beta / sum(prob_beta);

% Compute expected values
expected_alpha              = sum(p{1}(:) .* prob_alpha(:));
expected_beta               = sum(p{2}(:) .* prob_beta(:));

%% run the RW model with marginal and expected values and compare their resulting choice probs

% fit the model with the marginal parameter values from the search grid 
best_params             = [best_marginal_alpha best_marginal_beta];
[grid_model_output, ~]  = fit_modelRW_v1(best_params, modelout);
marginal_pp             = grid_model_output.allPs; clear grid_model_output

% fit the model with the expected parameter values 
expected_params         = [expected_alpha expected_beta];
[grid_model_output, ~]  = fit_modelRW_v1(expected_params, modelout);
expected_pp             = grid_model_output.allPs; clear grid_model_output


%% plot marginal vs expected choice probabilties 

%  let's look at how the model's predictions (based on marginal and expected 
% parameter values) align with the simulated data given that simulated data
% represent idealized or controlled behavior according to the model's assumptions

figure;
[~,~,h(1)] = myScatter(Ps(:,1),marginal_pp(:,1),false,[0 0 1],'x');
[~,~,h(2)] = myScatter(Ps(:,1),expected_pp(:,1),false,[1 0 0],'o');
h(3) = plot([0 1],[0 1],'k:','linewidth',2);
legend(h(1:2),{'Maximum Likelihood','Expected Value'},'location','best');legend boxoff
xlim([0 1]); ylim([0 1]);
xlabel('p(good otpion red) - simulated'); ylabel ('p(good otpion red) - estimated');
title('simulated vs. estimated choice probability');
fontsize(gcf,16,"points")

% display correlation coefficients 
corrCoeff = corrcoef(Ps, expected_pp);
disp(['Correlation (Simulated vs. Expected): ', num2str(corrCoeff(1,2))]);

corrCoeff = corrcoef(marginal_pp, expected_pp);
disp(['Correlation (Marginal vs. Expected): ', num2str(corrCoeff(1,2))]);


%%  RUN SEARCH GRID MULTIPLE TIMES %%

% INITIALIZATION FOR MULTIPLE SEARCHES
repetitions     = 50; % Number of repetitions for the grid search
all_best_params = zeros(repetitions, 2);    % to store best alpha and beta for each repetition
all_best_negLL  = inf(repetitions, 1);      % to store best negative log-likelihood for each repetition

%% MULTIPLE GRID SEARCHES

% loop over reps
for rep = 1:repetitions

    % initial variables and parameters
    init_alpha      = 0.3;
    init_beta       = 1.5;

    % start by simulating data and model choice probabilties
    data            = ALsimdata_v2(probabilities, total_trials,condtrials);

    % simulate model 
    modelout            = modelRW_v1([init_alpha init_beta], data);
    all_sim_pp(:,rep)   = modelout.allPs(:,1); 

    bounds              = [0 1;   % alpha range
                            0 5]; % beta range
    bins                = [20 25];
    nparam              = length(params);
    
    simalpha_bin        = init_alpha/bounds(1,2)*bins(1);
    simbeta_bin         = init_beta/bounds(2,2)*bins(2);
    best_negLL          = inf; % initialize with a large number
    
    % create the parameter space for the grid search
    for iParam = 1:nparam
        range           = linspace(bounds(iParam,1),bounds(iParam,2),bins(iParam)+1);
        p{iParam}       = range(2:end); % stay just off the zero bounds
    end
    
    % estimate nll using different combinations of parameter values 
    params          = nan(1,2);
    for t = 1:bins(1)
        params(1)   = p{1}(t);
    
        for tt = 1:bins(2)
            params(2)       = p{2}(tt);
            [~,negLL]       = fit_modelRW_v1(params,modelout); % Assuming 'data' is predefined or simulated earlier
    
            % store this iteration's nll
            nll(t,tt)       = negLL;
    
            % update best parameters if current nll is lower
            if negLL < best_negLL
                best_negLL  = negLL;
                best_alpha  = params(1);
                best_beta   = params(2);
            end
        end
    end
    
    % Store the results of this repetition
    all_best_params(rep, :) = [best_alpha, best_beta];
    all_best_negLL(rep)     = best_negLL;
    
    % Display results for this repetition
    fprintf('Repetition %d - Best alpha: %f, Best beta: %f, with negative LL: %f\n', rep, best_alpha, best_beta, best_negLL);

    % estimate marginal likelihoods of parameters alpha and beta
    % convert NLL to likelihoods for marginalization purposes
    likelihoods                 = exp(-nll);
    
    % marginalize (sum) likelihoods over beta to get marginal likelihood for alpha
    marginal_likelihood_alpha   = sum(likelihoods, 2); % Sum over columns (beta)
    
    % marginalize (sum) likelihoods over alpha to get marginal likelihood for beta
    marginal_likelihood_beta    = sum(likelihoods, 1); % Sum over rows (alpha)
    
    % optionally, convert marginal likelihoods back to NLL for visualization or analysis
    marginal_nll_alpha           = -log(marginal_likelihood_alpha);
    marginal_nll_beta            = -log(marginal_likelihood_beta);

    % estimate best marginal parameter values
    [~, idx_min_alpha]          = min(marginal_nll_alpha);
    [~, idx_min_beta]           = min(marginal_nll_beta);

    % Retrieve the best marginal values for alpha and beta from their respective ranges
    all_best_marginal_alpha(rep)= p{1}(idx_min_alpha);
    all_best_marginal_beta(rep) = p{2}(idx_min_beta);

    % estimate expected values of the parameter space 
    % convert marginal NLLs to probabilities (for illustrative purposes)
    % this is a simplified approach...
    prob_alpha                  = exp(-marginal_nll_alpha - min(marginal_nll_alpha));
    prob_alpha                  = prob_alpha / sum(prob_alpha);
    prob_beta                   = exp(-marginal_nll_beta - min(marginal_nll_beta));
    prob_beta                   = prob_beta / sum(prob_beta);
    
    % Compute expected values
    all_expected_alpha(rep)     = sum(p{1}(:) .* prob_alpha(:));
    all_expected_beta(rep)      = sum(p{2}(:) .* prob_beta(:));

    % run the RW model with marginal and expected values and compare their resulting choice probs
    % fit the model with the marginal parameter values from the search grid 
    best_params                 = [all_best_marginal_alpha(rep) all_best_marginal_beta(rep)];
    [grid_model_output, ~]      = fit_modelRW_v1(best_params, modelout);
    all_marginal_pp(:,rep)      = grid_model_output.allPs(:,1); clear grid_model_output

    % fit the model with the expected parameter values 
    expected_params             = [all_expected_alpha(rep) all_expected_beta(rep)];
    [grid_model_output, ~]      = fit_modelRW_v1(expected_params, modelout);
    all_expected_pp(:,rep)      = grid_model_output.allPs(:,1); clear grid_model_output

    clear expected best_params prob_alpha prob_beta idx_min_beta idx_min_alpha likelihoods
    clear p best_negLL best_alpha best_beta

end % loop over repetitions 

%% plot simulated vs expected vs marginal probabilties 

% Calculate mean probabilities across repetitions
mean_sim_pp         = mean(all_sim_pp, 2);
mean_marginal_pp    = mean(all_marginal_pp, 2);
mean_expected_pp    = mean(all_expected_pp, 2);

% Calculate the standard error of the mean (SEM) for each set of probabilities
sem_sim_pp          = std(all_sim_pp, 0, 2) / sqrt(size(all_sim_pp, 2));
sem_marginal_pp     = std(all_marginal_pp, 0, 2) / sqrt(size(all_marginal_pp, 2));
sem_expected_pp     = std(all_expected_pp, 0, 2) / sqrt(size(all_expected_pp, 2));

% Plotting
figure; hold on;
trials              = 1:length(mean_sim_pp); % Adjust this to the actual number of trials

% Error bar plots
errorbar(trials, mean_sim_pp, sem_sim_pp, 'LineWidth', 2, 'DisplayName', 'Simulated');
errorbar(trials, mean_marginal_pp, sem_marginal_pp, 'LineWidth', 2, 'DisplayName', 'Marginal');
errorbar(trials, mean_expected_pp, sem_expected_pp, 'LineWidth', 2, 'DisplayName', 'Expected');

% Add labels and legend
xlabel('Trial');
ylabel('Choice Probability');
title('Comparison of Choice Probabilities with Error Bars');
xlim([1 420]);
legend('show');

% Additional plot formatting
set(gca, 'FontSize', 16);
grid on;

%% perform parameter recovery 

% In parameter recovery, I'm using a Balanced Strategy: 
% A balanced and comprehensive parameter recovery strategy that combines 
% random and systematic exploration of parameter space.I use the best 
% parameters from the grid search as starting points for some of the optimization runs 
% (exploitation) and randomly generated starting points for others 
% (exploration). This combination ensures that while we are refining 
% the best known solutions, we're also open to discovering new solutions 
% that the grid search might have missed.

% how many repetitions?
repetitions = 100;

% define ranges for alpha and beta for true parameter values 
alpha_range = [0 1]; % range for alpha
beta_range  = [0 5]; % range for beta

% loop over repetitions 
for rep = 1:repetitions

    % Generate random true values for alpha and beta within the specified ranges
    true_alpha              = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand();
    true_beta               = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();

    allsim_alpha(rep)       = true_alpha;
    allsim_beta(rep)        = true_beta;

    % simulate data
    d                       = ALsimdata_v2(probabilities, total_trials,condtrials);
    out                     = modelRW_v1([true_alpha true_beta], d); % simulate model
    allsim_output{rep}      = out;

    % now since we employ a balanced stratedy for parameter recovery let's
    % see whether we are going to use parameter values that are randomly
    % generated or parameter values that were identified as good in the
    % previous grid search.
    if mod(rep, 2) == 0 % Check if the iteration number is even

        % select starting points randomly from the best grid search values
        all_best_alphas(1,:)    = all_best_params(:,1);
        all_best_betas(1,:)     = all_best_params(:,2);
        index                   = randi([1, length(all_best_alphas)]); % Generate a random index
        starting_alpha           = all_best_alphas(index);
        startingBeta            = all_best_betas(index);
    else
        % generate starting points using a random generator for odd iterations
        % assuming known plausible ranges for alpha and beta
        alphaRange      = [0 1]; % Adjust according to your parameter space
        betaRange       = [0 5]; % Adjust according to your parameter space
        starting_alpha  = alphaRange(1) + (alphaRange(2) - alphaRange(1)) * rand();
        startingBeta    = betaRange(1) + (betaRange(2) - betaRange(1)) * rand();
    end

    % define the objective function
    obFunc              = @(x) lik_modelRW_v1([x(1), x(2)],out);

    X0                  = [starting_alpha startingBeta];
    LB                  = [0 0];
    UB                  = [1 5];
    [Xfit, NegLL]       = fmincon(obFunc, X0, [], [], [], [], LB, UB);

    allfit_alpha(rep)   = Xfit(1);
    allfit_beta(rep)    = Xfit(2);
    allNLL(rep)         = NegLL;

end % end of reps loop

%% plot simulated and estimated parameter values

% calculate correlation coefficient for Alpha
alphaCorr       = corrcoef(allsim_alpha, allfit_alpha);
alphaCorrValue  = alphaCorr(1,2); % extract the correlation value

% calculate correlation coefficient for Beta
betaCorr        = corrcoef(allsim_beta, allfit_beta);
betaCorrValue   = betaCorr(1,2); % Extract the correlation value

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
title('Correlation between Simulated and Fitted Alpha Values');
set(gca, 'FontSize', 15); % additional plot formatting
grid on;

% for Beta values
subplot(1,2,2); % Beta values
scatter(allsim_beta, allfit_beta, 'filled');
hold on; % keep the scatter plot

p   = polyfit(allsim_beta, allfit_beta, 1); % fit a linear polynomial (degree 1)
x2  = linspace(min(allsim_beta), max(allsim_beta), 100); % generate x values
y2  = polyval(p, x2); % calculate y values based on the fitted polynomial
plot(x2, y2, '-r'); % plot the line
legend(['Correlation: ', num2str(betaCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Beta'); ylabel('Fitted Beta');
title('Correlation between Simulated and Fitted Beta Values');
set(gca, 'FontSize', 15); % additional plot formatting
grid on; 


%%




