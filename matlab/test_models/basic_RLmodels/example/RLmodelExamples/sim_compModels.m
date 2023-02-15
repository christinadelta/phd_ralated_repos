function outData = sim_compModels(params, model, simInfo, subjData)
%% Computational models for simulations
% Simulates choices and Q-values based on apriori trial data
% Or just 'Q-values', based on actual choices and rewards observed by participants
%
% Inputs:   params      [beta1, alpha1, ... alphaN]
%           model       current model ID
%           simInfo     has simType, Qinit,
%                       & variable names for output data tables
%           subjData    data in table format
%
%%%% NS 20201117 - Optimism & confirmation bias Y3 project:
%%%%%% free choices & full feedback


%% Models
%  1 - classic RL: 1 beta, 1 alpha
%  2 - lr by outcome type X valence: 1 beta, 4 alphas
%  3 - lr by confirmation bias: 1 beta, 2 alphas
%  4 - lr by optimism bias: 1 beta, 2 alphas

%% Initialise vars

% softmax=@(x,y) 1/(1+exp( x * diff(y) )); % x:beta, y:Q

% load choice related parameters
beta1 = params(1);


% Have a continuously updated Q matrix, for 1 state with 2 options (symbols,  Best vs Worst option)
Qs = nan(1, 2);  

nTrials = height(subjData);

switch simInfo.type
    case 'data'
        % output
        actionProbs = nan(nTrials, 2);
        sim_action = nan(nTrials, 1);
        action  = nan(nTrials, 1);
        outcome_f = nan(nTrials, 1);
        outcome_cf = nan(nTrials, 1);
        correct = ones(nTrials, 1); % always correct in simulated data, just serves for the Q update later
        
    case 'Qs'
        % input
        action  = subjData.action;
        outcome_f = subjData.outcome_f;
        outcome_cf = subjData.outcome_cf;
        correct = subjData.correct;   % in case error trials (for forced choices) are included    
        %output
        sim_action = nan(nTrials, 1);
end

% common output - list for each trial (regardless of outcome condition)
out_Qs        = nan(nTrials, 2); % at start of trial
out_PE_f      = nan(nTrials, 1); % at end of trial
out_PE_cf     = nan(nTrials, 1); % at end of trial


%% Loop through the blocks
for t = 1:nTrials
    
    % if first trial in block, then Qrestart
    if subjData.trialN(t)==1
        Qs = simInfo.Qinit;   % Q-values for Best vs Worst option to be updated
    end
    % for the output record
    out_Qs(t, :) = Qs;

    %% Simulate choice & outcome
    % % % Or skip this section if just simulating Q-values

    % Get probabilities for both actions
%         actionProb(t, :) = [softmax(beta1,  Qs(thisCond, :) ),...     % Best option
%                             softmax(beta1, fliplr(Qs(thisCond, :)))]; % Worst option

    actionProbs(t, :) = exp(Qs*beta1)/sum(exp(Qs*beta1));  

    % Simulate chosen action
    sim_action(t) = randsample([1,2], 1, true, actionProbs(t, :));

    if strcmp(simInfo.type, 'data')
        
        action(t) = sim_action(t);
        
        %% Simulate outcome

        if action(t) == 1 % Best symbol chosen (in simulations, we're ignoring location of symbols, hence 1=best)
            outcome_f(t) = subjData.out_hi(t);
            outcome_cf(t) = subjData.out_lo(t);
        else
            outcome_f(t) = subjData.out_lo(t);
            outcome_cf(t) = subjData.out_hi(t);            
        end
    end


    %% Update Q-values
    % % % One Q per action that gets continually updated

    if correct(t) % But only if trial is correct (relevant if based on real data)

        out_PE_f(t) = outcome_f(t) - Qs(action(t));   % save for record, but before updating Q    
        out_PE_cf(t) = outcome_cf(t) - Qs(3-action(t));   % 3-action flips 1/2 - unchosen option index
        Qs = Q_update( Qs, params, model, action(t), outcome_f(t), outcome_cf(t));
    end
    % if error, PEout(t) = NaN;

end % for t = 1:nTrials

%% 
outData = subjData;
        
outData.Q_hi = out_Qs(:,1);
outData.Q_lo = out_Qs(:,2);
outData.p_hi = actionProbs(:,1);
outData.p_lo = actionProbs(:,2);

outData.sim_action = action;
if strcmp(simInfo.type, 'data')
    outData.sim_outcome_f = outcome_f;
    outData.sim_outcome_cf = outcome_cf;
end

outData.PE_f = out_PE_f;
outData.PE_cf = out_PE_cf;

