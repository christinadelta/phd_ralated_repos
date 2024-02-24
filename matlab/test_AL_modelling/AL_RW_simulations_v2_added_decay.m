% Action Learning RW simulations version 2

% created February 2024

% added a decay/forgeting parameter to model the natural forgetting process
% or decrease in sakience of unchosen options over time 

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
params          = [0.3 1.5 0.1];

data            = ALsimdata_v2(probabilities, total_trials,condtrials);

%% simulate model 

modelout = modelRW_v2(params, data);

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
init_decay      = 0.1;
bounds          = [0 1;   % alpha range
                    0 5;  % beta range
                    0 1]; % decay range
bins            = [20 25 20];
nparam          = length(params);

simalpha_bin    = init_alpha/bounds(1,2)*bins(1);
simbeta_bin     = init_beta/bounds(2,2)*bins(2);
simdecay_bin    = init_decay/bounds(3,2)*bins(3);
best_negLL      = inf; % initialize with a large number

% create the parameter space for the grid search
for iParam = 1:nparam
    range       = linspace(bounds(iParam,1),bounds(iParam,2),bins(iParam)+1);
    p{iParam}   = range(2:end); % stay just off the zero bounds
end

% estimate nll using different combinations of parameter values 
params          = nan(1,3);
for t = 1:bins(1)
    params(1)   = p{1}(t);

    for tt = 1:bins(2)
        params(2)           = p{2}(tt);
        for ttt = 1:bins(3)

            params(3)           = p{3}(ttt);
            
            [~,negLL]           = fit_modelRW_v2(params,modelout);
    
            % store this iteration's nll
            nll(t,tt,ttt)           = negLL;
    
            % update best parameters if current nll is lower
            if negLL < best_negLL
                best_negLL  = negLL;
                best_alpha  = params(1);
                best_beta   = params(2);
                best_decay  = params(3);
            end
        end % end of decay
    end % end of beta
end % end of alpha

% display best parameter values 
fprintf('Best alpha: %f, Best beta: %f, Best decay: %f, with negative LL: %f\n', best_alpha, best_beta, best_decay, best_negLL);

% % for visualisation, compute the likelihood instaed of log
% lik = exp(nll); 

% PLOT THE NLLs AS AN IMAGESEC HIGHLIGHTING THE POSITION OF THE NLL WITH
% THE BEST PARAMETER VALUES

% specify the number of decay slices to visualize
num_decay_slices_to_visualize   = 2;

% calculate indices for evenly spaced decay slices
decay_indices                   = round(linspace(1, length(p{3}), num_decay_slices_to_visualize));

for i = 1:length(decay_indices)
    decay_index = decay_indices(i); % Get the current decay index from the selected indices
    nll_slice   = squeeze(nll(:,:,decay_index)); % Extract the NLL slice for the current decay value
    lik_slice   = exp(-nll_slice); % Convert to likelihood for visualization
    
    % Create a new figure for each selected decay slice
    figure;
    imagesc(p{1}, p{2}, lik_slice); % Adjusted for alpha (x-axis) and beta (y-axis)
    colorbar;
    xlabel('\alpha (Learning Rate)');
    ylabel('\beta (Softmax Inverse Temperature)');
    title(sprintf('Grid Search Results: Likelihood for Decay = %.2f', p{3}(decay_index)));
    fontsize(gcf,16,"points");
    set(gca, 'YDir', 'normal');
    
    % Optionally, add a marker if this slice contains the best parameters
    % This assumes 'simdecay_bin' is a way to identify the best decay index, which might need adjustment
    if any(decay_index == decay_indices) % Check if the current slice contains the best decay
        hold on;
        plot(best_alpha, best_beta, 'r*', 'MarkerSize', 30); % Note: This assumes best_alpha/beta are defined elsewhere
        hold off;
    end
end

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

% Assuming 'likelihoods' is a 3D matrix of likelihoods corresponding to alpha, beta, and decay
likelihoods = exp(-nll);

% Marginalize (sum) likelihoods over beta and decay to get marginal likelihood for alpha
marginal_likelihood_alpha = squeeze(sum(sum(likelihoods, 2), 3)); % Sum over columns (beta) and 3rd dimension (decay)

% Marginalize (sum) likelihoods over alpha and decay to get marginal likelihood for beta
marginal_likelihood_beta = squeeze(sum(sum(likelihoods, 1), 3)); % Sum over rows (alpha) and 3rd dimension (decay)

% Marginalize (sum) likelihoods over alpha and beta to get marginal likelihood for decay
marginal_likelihood_decay = squeeze(sum(sum(likelihoods, 1), 2)); % Sum over rows (alpha) and columns (beta)

% Optionally, convert marginal likelihoods back to NLL for visualization or analysis
marginal_nll_alpha = -log(marginal_likelihood_alpha);
marginal_nll_beta = -log(marginal_likelihood_beta);
marginal_nll_decay = -log(marginal_likelihood_decay);

% Visualize Marginal NLLs for all three parameters


% Plot marginal NLL for alpha
figure;
plot(p{1}, marginal_nll_alpha, 'LineWidth', 2);
xlabel('\alpha (Learning Rate)');
ylabel('Marginal Negative Log-Likelihood');
title('Marginal NLL for Alpha');
fontsize(gcf,16,"points");
grid on;

% Plot marginal NLL for beta
figure;
plot(p{2}, marginal_nll_beta, 'LineWidth', 2);
xlabel('\beta (Softmax Inverse Temperature)');
ylabel('Marginal Negative Log-Likelihood');
title('Marginal NLL for Beta');
fontsize(gcf,16,"points");
grid on;

% Plot marginal NLL for decay
figure;
plot(p{3}, marginal_nll_decay, 'LineWidth', 2);
xlabel('Decay');
ylabel('Marginal Negative Log-Likelihood');
title('Marginal NLL for Decay');
fontsize(gcf,16,"points");
grid on;

% Estimate best marginal parameter values
[~, idx_min_alpha] = min(marginal_nll_alpha);
[~, idx_min_beta] = min(marginal_nll_beta);
[~, idx_min_decay] = min(marginal_nll_decay);

% Retrieve the best marginal values for alpha, beta, and decay from their respective ranges
best_marginal_alpha = p{1}(idx_min_alpha);
best_marginal_beta = p{2}(idx_min_beta);
best_marginal_decay = p{3}(idx_min_decay);

fprintf('Best marginal alpha: %f\n', best_marginal_alpha);
fprintf('Best marginal beta: %f\n', best_marginal_beta);
fprintf('Best marginal decay: %f\n', best_marginal_decay);


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
% our goal is to compute an expected value of α and β and decay that reflects a weighted 
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
prob_decay                  = exp(-marginal_nll_decay - min(marginal_nll_decay));
prob_decay                  = prob_decay / sum(prob_decay);

% Compute expected values
expected_alpha              = sum(p{1}(:) .* prob_alpha(:));
expected_beta               = sum(p{2}(:) .* prob_beta(:));
expected_decay              = sum(p{3}(:) .* prob_decay(:));

% Now you have expected values for alpha, beta, and decay
fprintf('Expected alpha: %f\n', expected_alpha);
fprintf('Expected beta: %f\n', expected_beta);
fprintf('Expected decay: %f\n', expected_decay);

%% run the RW model with marginal and expected values and compare their resulting choice probs

% fit the model with the marginal parameter values from the search grid 
best_params             = [best_marginal_alpha best_marginal_beta, best_marginal_decay];
[grid_model_output, ~]  = fit_modelRW_v2(best_params, modelout);
marginal_pp             = grid_model_output.allPs; clear grid_model_output

% fit the model with the expected parameter values 
expected_params         = [expected_alpha expected_beta, expected_decay];
[grid_model_output, ~]  = fit_modelRW_v2(expected_params, modelout);
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
decay_range = [0 1]; % range for alpha

% loop over repetitions 
for rep = 1:repetitions

    % Generate random true values for alpha and beta within the specified ranges
    true_alpha              = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand();
    true_beta               = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
    true_decay              = decay_range(1) + (decay_range(2) - decay_range(1)) * rand();

    allsim_alpha(rep)       = true_alpha;
    allsim_beta(rep)        = true_beta;
    allsim_decay(rep)       = true_decay;

    % simulate data
    d                       = ALsimdata_v2(probabilities, total_trials,condtrials);
    out                     = modelRW_v2([true_alpha true_beta true_decay], d); % simulate model
    allsim_output{rep}      = out;

    % generate starting points using a random generator for odd iterations
    % assuming known plausible ranges for alpha and beta
    alphaRange      = [0 1]; 
    betaRange       = [0 5];
    decayRange      = [0 1]; 

    starting_alpha  = alphaRange(1) + (alphaRange(2) - alphaRange(1)) * rand();
    startingBeta    = betaRange(1) + (betaRange(2) - betaRange(1)) * rand();
    starting_decay  = decayRange(1) + (decayRange(2) - decayRange(1)) * rand();

    % define the objective function
    obFunc              = @(x) lik_modelRW_v2([x(1), x(2), x(3)],out);

    X0                  = [starting_alpha startingBeta starting_decay];
    LB                  = [0 0 0];
    UB                  = [1 5 1];
    [Xfit, NegLL]       = fmincon(obFunc, X0, [], [], [], [], LB, UB);

    allfit_alpha(rep)   = Xfit(1);
    allfit_beta(rep)    = Xfit(2);
    allfit_decay(rep)   = Xfit(3);
    allNLL(rep)         = NegLL;

 
end % end of reps loop


%% plot simulated and estimated parameter values

% calculate correlation coefficient for Alpha
alphaCorr       = corrcoef(allsim_alpha, allfit_alpha);
alphaCorrValue  = alphaCorr(1,2); % extract the correlation value

% calculate correlation coefficient for Beta
betaCorr        = corrcoef(allsim_beta, allfit_beta);
betaCorrValue   = betaCorr(1,2); % extract the correlation value

decayCorr = corrcoef(allsim_decay, allfit_decay);
decayCorrValue = decayCorr(1,2);

% For Alpha values
figure; 
subplot(1,3,1); % alpha values
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

% for Beta values
subplot(1,3,2); % Beta values
scatter(allsim_beta, allfit_beta, 'filled');
hold on; % keep the scatter plot

p   = polyfit(allsim_beta, allfit_beta, 1); % fit a linear polynomial (degree 1)
x2  = linspace(min(allsim_beta), max(allsim_beta), 100); % generate x values
y2  = polyval(p, x2); % calculate y values based on the fitted polynomial
plot(x2, y2, '-r'); % plot the line
legend(['Correlation: ', num2str(betaCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Beta'); ylabel('Fitted Beta');
title('Simulated and Fitted Beta Values');
set(gca, 'FontSize', 14); % additional plot formatting
grid on; 

% for Decay values
subplot(1,3,3); % Beta values
scatter(allsim_decay, allfit_decay, 'filled');
hold on; % keep the scatter plot

p   = polyfit(allsim_decay, allfit_decay, 1);               % fit a linear polynomial (degree 1)
x2  = linspace(min(allsim_decay), max(allsim_decay), 100);  % generate x values
y2  = polyval(p, x2);                                       % calculate y values based on the fitted polynomial
plot(x2, y2, '-r');                                         % plot the line
legend(['Correlation: ', num2str(decayCorrValue, '%.2f')]);

hold off;
xlabel('Simulated Decay'); ylabel('Fitted Decay');
title('Simulated and Fitted Decay Values');
set(gca, 'FontSize', 14); % additional plot formatting
grid on; 

%% sensitivity analysis with free parameter: alpha 



