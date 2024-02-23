% RW MODEL SIMULATIONS AND FITTING USING BAYESOPT 

% CREATE JANUARY 2024 FOR TESTING PURPOSES ONLY 

% Here, I test the AL task with a simple RW model with two parameters:
% alpha - learning rate
% beta - inverse temperature for the softmax funcyion
% for each block seperately 

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
total_trials    = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
nCues           = 2;
params_alpha    = [0.4 0.35 0.30;
                    0.7 0.6 0.5];

params_beta     = [3 2.5 2;
                    2.5 2 1.5];

data            = ALsimdata_v2(probabilities, total_trials,condtrials);

%% simulate model for each block seperately 

% extract data needed
v = data.t;
s = data.stcind;
o = data.oR;
c = data.cues;
x = data.x;

for ss = 1:3

    stc_o = o(s(:,ss),:);
    stc_x = x(s(:,ss),:);

    for vv = 1:2

        dat.oR      = stc_o(v(:,vv),:);
        vol_x       = stc_x(v(:,vv),:);
        params      = [params_alpha(vv,ss) params_beta(vv,ss)];
        dat.cues    = c;

        % simulate model
        modelout    = modelRW_v1(params, dat);

        % for plotting we we need Qvalues, choice probabilities and choices both conditions, so, concatinate the
        % trials of the two conditions
        Qvals{ss,vv}   = modelout.Qvals;
        Ps{ss,vv}      = modelout.allPs; %
        choices        = modelout.a;
        a{ss,vv}       = 2 - choices; % convert to be 1 and 0
        correct{ss,vv} = modelout.correct;
        
        % get score
        score(ss,vv)        = mean(correct{ss,vv});
        all_x{ss,vv}        = vol_x;
        all_modout{ss,vv}   = modelout;

    end % end of volatility loop
end  % end of stc loop

%% plot for each condition seperatley 

h = plotRW_v2_modified(Qvals, Ps, a, all_x);

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
% Initialize storage for best parameters and NLLs
best_params = nan(3, 2, 2); % 3 levels of stochasticity, 2 levels of volatility, 2 parameters (alpha, beta)
best_negLL = inf(3, 2); % Initialize with a large number for each condition

% Modify bins as per your requirement
bins = [20 25]; % Example bin sizes for alpha and beta

bounds = [0 1;   % alpha range
          0 5]; % beta range


% Parameter space creation
p = cell(1,2);
for iParam = 1:2 % Assuming 2 parameters: alpha, beta
    range = linspace(bounds(iParam,1), bounds(iParam,2), bins(iParam)+1);
    p{iParam} = range(2:end); % Avoiding zero bounds
end

% Grid search
nll = cell(3,2); % For storing NLLs for each condition

for stc = 1:3 % Stochasticity levels
    for vol = 1:2 % Volatility levels

        % Initialize NLL storage for the current condition
        nll{stc,vol} = inf(bins(1), bins(2));

        for t = 1:bins(1) % Alpha loop
            for tt = 1:bins(2) % Beta loop

                params = [p{1}(t), p{2}(tt)];
                mout = all_modout{stc, vol};
                [~, negLL] = fit_modelRW_v1(params, mout);

                % Store NLL
                nll{stc,vol}(t,tt) = negLL;

                % Check for best NLL and update parameters
                if negLL < best_negLL(stc, vol)
                    best_negLL(stc, vol) = negLL;
                    best_params(stc, vol, :) = params;
                end

            end % Beta loop
        end % Alpha loop

    end % Volatility loop
end % Stochasticity loop

% At this point, best_params contains the best alpha and beta for each condition,
% and best_negLL contains the corresponding NLL values.

% PLOT NLLs and mark best parameter valeus
% Define descriptive labels for stochasticity and volatility
stochasticityLabels = {'Small', 'Medium', 'Large'};
volatilityLabels = {'Stable', 'Volatile'};

% Create a figure to hold all subplots
figure;

% Loop over all conditions
for stc_idx = 1:length(stochasticityLabels)
    for vol_idx = 1:length(volatilityLabels)
        
        % Calculate subplot index
        subplotIdx = (stc_idx - 1) * length(volatilityLabels) + vol_idx;
        
        % Select subplot
        subplot(length(stochasticityLabels), length(volatilityLabels), subplotIdx);
        
        % Extract the NLL matrix for the current condition
        current_nll = nll{stc_idx,vol_idx};
        
        % Plotting
        imagesc(p{1}, p{2}, current_nll');
        colorbar;
        xlabel('Alpha');
        ylabel('Beta');
        fontsize(gcf,16,"points")
        title(sprintf('NEGLL: %s Stochasticity, %s', stochasticityLabels{stc_idx}, volatilityLabels{vol_idx}));
        
        % Highlight the best parameters
        hold on;
        best_alpha = best_params(stc_idx, vol_idx, 1);
        best_beta = best_params(stc_idx, vol_idx, 2);
        plot(best_alpha, best_beta, 'rp', 'MarkerSize', 16, 'MarkerFaceColor', 'r');
        hold off;
        
    end
end

% Adjust subplot spacing if needed
sgtitle('NEGLL for All Conditions'); % Super title for the whole figure
fontsize(gcf,16,"points")

%% estimate marginal likelihoods of parameter values

% Assuming there are 3 conditions and 2 volatilities, leading to 6 condition/volatility pairs
nConditions = 3;
nVolatilities = 2;
totalPairs = nConditions * nVolatilities;

% Assuming bins for beta are 25, initialize the matrices to store marginal NLLs
aggregated_marginal_nll_alpha = zeros(20, totalPairs); % Adjust 20 if alpha bins differ
aggregated_marginal_nll_beta = zeros(25, totalPairs);

% Counter for condition/volatility pairs
pairIndex = 0;

% Loop over conditions and volatilities to calculate marginal likelihoods
for stc = 1:nConditions
    for vol = 1:nVolatilities
        pairIndex = pairIndex + 1; % Increment pair index for each condition/volatility pair
        
        % Convert NLL to likelihoods for marginalization
        likelihoods = exp(-nll{stc, vol});
        
        % Marginalize likelihoods over parameters
        marginal_likelihood_alpha = sum(likelihoods, 2); % Sum over beta
        marginal_likelihood_beta = sum(likelihoods, 1); % Sum over alpha
        
        % Convert marginal likelihoods to NLL
        marginal_nll_alpha = -log(marginal_likelihood_alpha);
        marginal_nll_beta = -log(marginal_likelihood_beta);
        
        % Correctly aggregate the marginal NLLs
        aggregated_marginal_nll_alpha(:, pairIndex) = marginal_nll_alpha;
        aggregated_marginal_nll_beta(:, pairIndex) = marginal_nll_beta(:); % Ensure column vector
        
    end
end

% Assuming aggregated_marginal_nll_alpha and aggregated_marginal_nll_beta have the correct dimensions
% and there are 3 conditions (rows) and 2 volatilities (columns) leading to 6 plots in total

nConditions = 3; % Number of conditions
nVolatilities = 2; % Number of volatilities

% Plot Marginal NLL for Alpha across all conditions and volatilities
figure('Name', 'Marginal NLL for Alpha across Conditions and Volatilities', 'NumberTitle', 'off');
for stc = 1:nConditions
    for vol = 1:nVolatilities
        subplotIndex = (stc - 1) * nVolatilities + vol;
        subplot(nConditions, nVolatilities, subplotIndex);
        
        % Plot for the current condition/volatility pair
        plot(p{1}, aggregated_marginal_nll_alpha(:, subplotIndex), 'LineWidth', 2);
        xlabel('\alpha (Learning Rate)');
        ylabel('Marginal Negative LL');
        title(sprintf('%s Stochasticity, %s',stochasticityLabels{stc}, volatilityLabels{vol}));
        fontsize(gcf,14,"points")
        grid on;
    end
end

% Plot Marginal NLL for Beta across all conditions and volatilities
figure('Name', 'Marginal NLL for Beta across Conditions and Volatilities', 'NumberTitle', 'off');
for stc = 1:nConditions
    for vol = 1:nVolatilities
        subplotIndex = (stc - 1) * nVolatilities + vol;
        subplot(nConditions, nVolatilities, subplotIndex);
        
        % Plot for the current condition/volatility pair
        plot(p{2}, aggregated_marginal_nll_beta(:, subplotIndex), 'LineWidth', 2);
        xlabel('\beta (Softmax Inverse Temperature)');
        ylabel('Marginal Negative LL');
        title(sprintf('%s Stochasticity, %s',stochasticityLabels{stc}, volatilityLabels{vol}));
        fontsize(gcf,14,"points")
        grid on;
    end
end



%% estimate expected values of the parameter space 


% Initialize variables for expected values
expected_alpha = zeros(1, totalPairs);
expected_beta = zeros(1, totalPairs);

% Loop over all condition/volatility pairs
for pairIndex = 1:totalPairs
    
    % Convert aggregated marginal NLLs to probabilities for alpha
    prob_alpha = exp(-aggregated_marginal_nll_alpha(:, pairIndex) - min(aggregated_marginal_nll_alpha(:, pairIndex)));
    prob_alpha = prob_alpha / sum(prob_alpha);
    
    % Convert aggregated marginal NLLs to probabilities for beta
    prob_beta = exp(-aggregated_marginal_nll_beta(:, pairIndex) - min(aggregated_marginal_nll_beta(:, pairIndex)));
    prob_beta = prob_beta / sum(prob_beta);
    
    % Compute expected values for alpha and beta for the current pair
    expected_alpha(pairIndex) = sum(p{1}(:) .* prob_alpha(:));
    expected_beta(pairIndex) = sum(p{2}(:) .* prob_beta(:));
end

% 'expected_alpha' and 'expected_beta' now contain the expected values of alpha and beta
% for each condition/volatility pair

%% fit model with the best marginal and expected values from grid search

% Initialize storage for best marginal parameters
best_marginal_alpha = zeros(1, totalPairs);
best_marginal_beta = zeros(1, totalPairs);

% Loop over all condition/volatility pairs to find best marginal parameters
for pairIndex = 1:totalPairs
    
    % Find index of the best marginal alpha and beta (minimum NLL)
    [~, idx_alpha] = min(aggregated_marginal_nll_alpha(:, pairIndex));
    [~, idx_beta] = min(aggregated_marginal_nll_beta(:, pairIndex));
    
    % Store the best marginal alpha and beta values
    best_marginal_alpha(pairIndex) = p{1}(idx_alpha);
    best_marginal_beta(pairIndex) = p{2}(idx_beta);
    
    % Prepare model input for the current condition/volatility pair
    modelout = all_modout{ceil(pairIndex/nVolatilities), mod(pairIndex-1, nVolatilities)+1};
    
    % Fit the model with the best marginal parameter values for the current pair
    best_params = [best_marginal_alpha(pairIndex) best_marginal_beta(pairIndex)];
    [grid_model_output, ~] = fit_modelRW_v1(best_params, modelout);
    marginal_pp{pairIndex} = grid_model_output.allPs;
    clear grid_model_output;
    
    % Fit the model with the expected parameter values for the current pair
    expected_params = [expected_alpha(pairIndex) expected_beta(pairIndex)];
    [grid_model_output, ~] = fit_modelRW_v1(expected_params, modelout);
    expected_pp{pairIndex} = grid_model_output.allPs;
    clear grid_model_output;
end

%% plot marginal vs expected choice probabilties 

nConditions = 3; % Number of stochasticity levels
nVolatilities = 2; % Number of volatility levels
figure;

% Define descriptive labels for stochasticity and volatility
stochasticityLabels = {'Small', 'Medium', 'Large'};
volatilityLabels = {'Stable', 'Volatile'};

for stc = 1:nConditions
    for vol = 1:nVolatilities
        subplotIdx = (stc - 1) * nVolatilities + vol;
        subplot(nConditions, nVolatilities, subplotIdx);
        hold on;
        
        % Plot for marginal probabilities
        % Adjust indexing for Ps to match the 3x2 structure
        [~,~,h(1)] = myScatter(Ps{stc, vol}(:,1), marginal_pp{(stc-1)*nVolatilities+vol}(:,1), false, [0 0 1], 'x');
        
        % Plot for expected probabilities
        [~,~,h(2)] = myScatter(Ps{stc, vol}(:,1), expected_pp{(stc-1)*nVolatilities+vol}(:,1), false, [1 0 0], 'o');
        
        % Add unity line
        h(3) = plot([0 1], [0 1], 'k:', 'linewidth', 2);
        legend(h(1:2), {'Maximum Likelihood', 'Expected Value'}, 'location', 'best'); legend boxoff;
        xlim([0 1]); ylim([0 1]);
        xlabel('p(good option red) - simulated');
        ylabel('p(good option red) - estimated');
        
        % Title for each subplot with condition and volatility
        title(sprintf('%s Stochasticity, %s',stochasticityLabels{stc}, volatilityLabels{vol}));
        fontsize(gcf, 15, "points");
        
        % Display correlation coefficients for each subplot
        corrCoeffSimExpected = corrcoef(Ps{stc, vol}(:,1), expected_pp{(stc-1)*nVolatilities+vol}(:,1));
        disp(['Condition ', num2str(stc), ', Volatility ', num2str(vol), ' - Correlation (Simulated vs. Expected): ', num2str(corrCoeffSimExpected(1,2))]);
        
        corrCoeffMarginalExpected = corrcoef(marginal_pp{(stc-1)*nVolatilities+vol}(:,1), expected_pp{(stc-1)*nVolatilities+vol}(:,1));
        disp(['Condition ', num2str(stc), ', Volatility ', num2str(vol), ' - Correlation (Marginal vs. Expected): ', num2str(corrCoeffMarginalExpected(1,2))]);
        
        hold off;
    end
end

%% Initialize parameters for repetitions
nRepetitions = 100; % Number of times to repeat the grid search for each condition

% Adjust the storage structures to accommodate repetitions
best_params_repeated = nan(3, 2, 2, nRepetitions); % 3x2 conditions, 2 parameters, for each repetition
best_negLL_repeated = inf(3, 2, nRepetitions); % Best negative log-likelihood for each condition and repetition

% Modify bins as per your requirement
bins = [20 25]; % Example bin sizes for alpha and beta

bounds = [0 1;   % alpha range
          0 5]; % beta range


% Parameter space creation
p = cell(1,2);

for iParam = 1:2 % Assuming 2 parameters: alpha, beta
    range = linspace(bounds(iParam,1), bounds(iParam,2), bins(iParam)+1);
    p{iParam} = range(2:end); % Avoiding zero bounds
end

%% Perform grid search with repetitions %%

for rep = 1:nRepetitions

    fprintf('Running repetition %d...\n', rep);
    
    % Temporary storage for this repetition
    best_params = nan(3, 2, 2); % 3 levels of stochasticity, 2 levels of volatility, 2 parameters (alpha, beta)
    best_negLL = inf(3, 2); % Initialize with a large number for each condition

    % simulate data 
    data            = ALsimdata_v2(probabilities, total_trials,condtrials);

    % extract data needed
    v = data.t;
    s = data.stcind;
    o = data.oR;
    c = data.cues;
    x = data.x;

    % Your existing grid search code, slightly modified for repetition
    for stc = 1:3 % Stochasticity levels

        % extract outcomes and state 
        stc_o = o(s(:,stc),:);
        stc_x = x(s(:,stc),:);

        for vol = 1:2 % Volatility levels
            nll = inf(bins(1), bins(2)); % Initialize NLL storage for the current condition

            for t = 1:bins(1) % Alpha loop
                for tt = 1:bins(2) % Beta loop
                    params = [p{1}(t), p{2}(tt)];

                    dat.oR      = stc_o(v(:,vol),:);
                    vol_x       = stc_x(v(:,vol),:);
                    init_params = [0.3 1.5];
                    dat.cues    = c;
            
                    % simulate model and use the output for fitting 
                    modelout        = modelRW_v1(init_params, dat);
                    sim_Ps(:,:,rep) = modelout;
                    [~, negLL]      = fit_modelRW_v1(params, modelout);

                    % Store NLL
                    if negLL < nll(t,tt)
                        nll(t,tt) = negLL;
                    end

                    % Check for best NLL and update parameters
                    if negLL < best_negLL(stc, vol)
                        best_negLL(stc, vol) = negLL;
                        best_params(stc, vol, :) = params;
                    end
                end % Beta loop
            end % Alpha loop
        end % Volatility loop
    end % Stochasticity loop
    
    % Store the best parameters and NLLs for this repetition
    best_params_repeated(:,:,:,rep) = best_params;
    best_negLL_repeated(:,:,rep) = best_negLL;
    
    % Compute marginal likelihoods and expected values here if they should be updated every repetition
    % Otherwise, compute them outside of this loop if they're meant to aggregate across repetitions
end


%% reshape the best parameters matrices

% Dimension 1: Stochasticity Levels
% Dimension 2: Volatility Levels
% Dimension 3: Parameters (1 for alpha, 2 for beta)
% Dimension 4: Repetitions

nRepetitions = size(best_params_repeated, 4);

% Extract alpha values
alpha_params = best_params_repeated(:,:,1,:); % This selects all alphas
% Reshape to get rid of the singleton 3rd dimension
alpha_params_reshaped = reshape(alpha_params, [3, 2, nRepetitions]);

% Extract beta values
beta_params = best_params_repeated(:,:,2,:); % This selects all betas
% Reshape to get rid of the singleton 3rd dimension
beta_params_reshaped = reshape(beta_params, [3, 2, nRepetitions]);

% Initialize cell arrays to store the reshaped matrices for each stochasticity level
alpha_reshaped_by_stc = cell(1, 3); % Since you have 3 stochasticity levels
beta_reshaped_by_stc = cell(1, 3);

for stc = 1:3 % Loop through each stochasticity level
    % Extract the data for the current stochasticity level
    % This results in a [1, 2 volatility levels, 50 repetitions] matrix for alpha and beta
    alpha_current_stc = squeeze(alpha_params_reshaped(stc, :, :)); % Removing the singleton dimension
    beta_current_stc = squeeze(beta_params_reshaped(stc, :, :));
    
    % Since we want [50 repetitions, 2 volatility levels] and now have [2, 50],
    % we transpose the result
    alpha_reshaped_by_stc{stc} = alpha_current_stc.'; % Transpose to get [50, 2]
    beta_reshaped_by_stc{stc} = beta_current_stc.'; % Transpose to get [50, 2]
end

%% perform parameter recovery 

% In parameter recovery, this time I am randomly generating starting values

% how many repetitions?
repetitions = 100;

% define ranges for alpha and beta for true parameter values 
alpha_range = [0 1]; % range for alpha
beta_range  = [0 5]; % range for beta

% Define total number of stochasticity levels and volatility levels
nStochasticityLevels = 3;
nVolatilityLevels = 2;

for rep = 1:repetitions

    % Generate the entire task's data for this repetition
    data = ALsimdata_v2(probabilities, total_trials, condtrials);
    
    % Extract global data needed for condition-specific extraction
    v = data.t;
    s = data.stcind;
    o = data.oR;
    c = data.cues;
    x = data.x;

    % Loop over conditions within this repetition
    for ss = 1:3 % Stochasticity levels
        stc_o = o(s(:,ss),:);
        stc_x = x(s(:,ss),:);

        % extract best alpha and beta parameter values from grid search
        this_stc_best_alphas    = alpha_reshaped_by_stc{1,ss};
        this_stc_best_betas     = beta_reshaped_by_stc{1,ss};

        for vv = 1:2 % Volatility levels
            dat.oR      = stc_o(v(:,vv),:);
            dat.cues    = c;

            % Generate random true values for alpha and beta within the specified ranges
            true_alpha              = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand();
            true_beta               = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
        
            allsim_alpha{ss,vv}(rep)    = true_alpha;
            allsim_beta{ss,vv}(rep)     = true_beta;

            % Simulate model for the current condition
            modelout                        = modelRW_v1([true_alpha true_beta], dat);
            allsim_output{ss,vv}{rep}       = modelout;
            allsim_pp{ss,vv}(:,rep)         = modelout.allPs(:,1);
            allsim_qval_r{ss,vv}(:,rep)    = modelout.Qvals(:,1);
            allsim_qval_b{ss,vv}(:,rep)    = modelout.Qvals(:,2);

            % Here, I perform parameter recovery for the current condition

            % generate starying points for fititng :). Let's employ a
            % balanced startegy again!
            if mod(rep, 2) == 0 % Check if the iteration number is even

                % select starting points randomly from the best grid search values
                all_best_alphas(1,:)    = this_stc_best_alphas(:,vv);
                all_best_betas(1,:)     = this_stc_best_betas(:,vv);
                index                   = randi([1, length(all_best_alphas)]); % Generate a random index
                starting_alpha          = all_best_alphas(index);
                starting_beta           = all_best_betas(index);
            else
                % generate starting points using a random generator for odd iterations
                % assuming known plausible ranges for alpha and beta
                starting_alpha  = alpha_range(1) + (alpha_range(2) - alpha_range(1)) * rand();
                starting_beta    = beta_range(1) + (beta_range(2) - beta_range(1)) * rand();
            end
                           
            % define the objective function
            obFunc                      = @(x) lik_modelRW_v1([x(1), x(2)],modelout);
        
            X0                          = [starting_alpha starting_beta];
            LB                          = [0 0];
            UB                          = [1 5];
            [Xfit, NegLL]               = fmincon(obFunc, X0, [], [], [], [], LB, UB);

            allfit_alpha{ss,vv}(rep)    = Xfit(1);
            allfit_beta{ss,vv}(rep)     = Xfit(2);
            allNLL{ss,vv}(rep)          = NegLL;

            % get model's choice probabilities 
            [fit_model_output, ~]      = fit_modelRW_v1(Xfit, modelout);
            all_fit_pp{ss,vv}(:,rep)   = fit_model_output.allPs(:,1); % 
            allfit_qval_r{ss,vv}(:,rep) = fit_model_output.Qvals(:,1);
            allfit_qval_b{ss,vv}(:,rep) = fit_model_output.Qvals(:,2);

        end
    end
end

%% %% 

% Assuming the structure of allsim_alpha, allsim_beta, allfit_alpha, allfit_beta
% is {stochasticityLevel, volatilityLevel} and contains data for each condition

nStochasticityLevels = 3;
nVolatilityLevels = 2;

% Create a figure for alpha values
figure; 
sgtitle('Correlation between Simulated and Fitted Alpha Values per Condition');

for ss = 1:nStochasticityLevels
    for vv = 1:nVolatilityLevels
        subplotIdx = (ss - 1) * nVolatilityLevels + vv;
        subplot(nStochasticityLevels, nVolatilityLevels, subplotIdx);
        
        sim_alpha = allsim_alpha{ss,vv};
        fit_alpha = allfit_alpha{ss,vv};
        
        % Calculate and plot correlation for Alpha
        scatter(sim_alpha, fit_alpha, 'filled');
        hold on;
        
        p = polyfit(sim_alpha, fit_alpha, 1);
        x1 = linspace(min(sim_alpha), max(sim_alpha), 100);
        y1 = polyval(p, x1);
        plot(x1, y1, '-r');
        
        alphaCorr = corrcoef(sim_alpha, fit_alpha);
        alphaCorrValue = alphaCorr(1,2);
        
        legend(['Correlation: ', num2str(alphaCorrValue, '%.2f')], 'Location', 'Best');
        xlabel('Simulated Alpha'); ylabel('Fitted Alpha');
        title(sprintf('%s Stochasticity, %s',stochasticityLabels{ss}, volatilityLabels{vv}));
        set(gca, 'FontSize', 14);
        grid on;
        
        hold off;
    end
end

% Create a figure for beta values
figure; 
sgtitle('Correlation between Simulated and Fitted Beta Values per Condition');

for ss = 1:nStochasticityLevels
    for vv = 1:nVolatilityLevels
        subplotIdx = (ss - 1) * nVolatilityLevels + vv;
        subplot(nStochasticityLevels, nVolatilityLevels, subplotIdx);
        
        sim_beta = allsim_beta{ss,vv};
        fit_beta = allfit_beta{ss,vv};
        
        % Calculate and plot correlation for Beta
        scatter(sim_beta, fit_beta, 'filled');
        hold on;
        
        p = polyfit(sim_beta, fit_beta, 1);
        x2 = linspace(min(sim_beta), max(sim_beta), 100);
        y2 = polyval(p, x2);
        plot(x2, y2, '-r');
        
        betaCorr = corrcoef(sim_beta, fit_beta);
        betaCorrValue = betaCorr(1,2);
        
        legend(['Correlation: ', num2str(betaCorrValue, '%.2f')], 'Location', 'Best');
        xlabel('Simulated Beta'); ylabel('Fitted Beta');
        title(sprintf('%s Stochasticity, %s',stochasticityLabels{ss}, volatilityLabels{vv}));
        set(gca, 'FontSize', 14);
        grid on;
        
        hold off;
    end
end

%% plot probabilties and qvalues 

% Initialize arrays to hold concatenated data
concat_sim_pp       = [];
concat_fit_pp       = [];
concat_sim_val_r    = [];
concat_sim_val_b    = [];
concat_fit_val_r    = [];
concat_fit_val_b    = [];


% Loop through each condition to concatenate data
for ss = 1:3
    for vv = 1:2

        % Concatenate simulated choice probabilities
        concat_sim_pp = [concat_sim_pp; mean(allsim_pp{ss,vv}, 2)];
        
        % Concatenate fitted choice probabilities
        concat_fit_pp = [concat_fit_pp; mean(all_fit_pp{ss,vv}, 2)];

        % concatinate simulated Q values
        concat_sim_val_r = [concat_sim_val_r; mean(allsim_qval_r{ss,vv},2)];
        concat_sim_val_b = [concat_sim_val_b; mean(allsim_qval_b{ss,vv},2)];

        % concatinate fitted Q values
        concat_fit_val_r = [concat_fit_val_r; mean(allfit_qval_r{ss,vv},2)];
        concat_fit_val_b = [concat_fit_val_b; mean(allfit_qval_b{ss,vv},2)];

        
        
    end
end

% Start plotting
figure;

% Subplot 1: Choice Probabilities
subplot(2, 1, 1); % This specifies a 2-row, 1-column grid of subplots, and selects the 1st subplot
hold on;

% Plot concatenated true reward rates
plot(x, 'k-', 'LineWidth', 1, 'DisplayName', 'True Reward Rate');

% Plot concatenated simulated choice probabilities
plot(concat_sim_pp, 'b--', 'LineWidth', 2, 'DisplayName', 'Simulated Choice Probabilities');

% Plot concatenated fitted choice probabilities
plot(concat_fit_pp, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Fitted Choice Probabilities');

% Formatting the subplot
title('Simulated and Estimated Choice Probabilities Over Trials');
xlabel('Trial');
ylabel('Probability');
legend('show');
ylim([0 1]);
xlim([1 420]);
fontsize(gcf,16,"points")
hold off;

% Subplot 2: Q Values
subplot(2, 1, 2); % Selects the 2nd subplot in a 2-row, 1-column grid
hold on;

% plot lines to specify where change happens 
dx          = [0; diff(x)~=0];    
ix          = (find(dx))-1;       % trials before each switch of probabilities
tt          = 1:length(x);        % total trials 

% Plot and store handles for lines you want in the legend
hSimR = plot(concat_sim_val_r, 'r--', 'LineWidth', 1, 'DisplayName', 'Simulated Q Value Option A');
hSimB = plot(concat_sim_val_b, 'b--', 'LineWidth', 1, 'DisplayName', 'Simulated Q Value Option B');
hFitR = plot(concat_fit_val_r, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Fitted Q Value Option A');
hFitB = plot(concat_fit_val_b, 'b:', 'LineWidth', 1.5, 'DisplayName', 'Fitted Q Value Option B');

ylim([-1 1]); % Adjust if your Q values fall outside this range
xlim([1 420]);

ym      = get(gca,'ylim');
% Specify trial ranges where volatility is high and plot vertical lines
for i = 1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',"black",'linewidth',0.5);
end

% Only include specified handles in the legend
legend([hSimR, hSimB, hFitR, hFitB], 'Simulated Q Value Option A', 'Simulated Q Value Option B', 'Fitted Q Value Option A', 'Fitted Q Value Option B');

% Formatting the subplot
title('Simulated and Estimated Q Values Over Trials');
xlabel('Trial');
ylabel('Q Value');

fontsize(gcf,16,"points")
hold off;


