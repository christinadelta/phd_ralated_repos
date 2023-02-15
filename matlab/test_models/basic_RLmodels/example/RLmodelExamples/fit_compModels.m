function lik = fit_compModels(params, model, Qinit, subjData, modelInfo, prior)
%% For model optimisation, estimates all model parameters in one go
% output:  lik is negative log-likelihood or lop posterior (if MAP)
%
% NB: Error trials are assumed to have already been excluded from the input
%%%% NS 20201117 - Optimism & confirmation bias Y3 project:
%%%%%% free choices & full feedback


%% Initialise vars

% load choice related parameters
beta1 = params(1);

% Have a continuously updated Q matrix, for 1 state with 2 options (symbols,  Best vs Worst option)
Qs = nan(1, 2);  

%% Loop through the trials
nTrials = height(subjData);
p = nan(1,nTrials); % vector of probabilities for the *chosen* action

for t = 1:nTrials
    
    % if first trial in block, then Qrestart
    if subjData.trialN(t)==1
        Qs = Qinit;   % Q-values for Best vs Worst option to be updated
    end


    %% Apply softmax to get probability of chosen action - only for free choices

    % Get probability for chosen action
    p(t) = exp(Qs(subjData.action(t))*beta1)/sum(exp(Qs*beta1));  

    %% Update Q-values
    % % % One Q per action that gets continually updated
        
    Qs = Q_update( Qs, params, model, subjData.action(t), subjData.outcome_f(t), subjData.outcome_cf(t));
    
end % for t = 1:nTrials



if nargin<6
   % Negative log-likelihood
    lik  = -sum(log(p));

else % integrate lik with prior to get log posterior
    log_priors = zeros(1,length(params));
    for iP = 1:length(params)
        tmp = split(modelInfo.params{model}.labels{iP}, '_'); tmp=tmp{1};
        
        if modelInfo.priorDists.lpmf==1 % then use lpmf (Log probability mass function)
            tmpf = str2func(strcat(modelInfo.priorDists.(tmp), '_lpmf'));
            log_priors(iP) = tmpf(params(iP), prior{iP}(1), prior{iP}(2), modelInfo.tolerance);
        
        else % simple current density
            if strcmp(modelInfo.priorDists.(tmp), 'beta')        
                log_priors(iP) = log( betacdf(params(iP),prior{iP}(1),prior{iP}(2)));

            elseif strcmp(modelInfo.priorDists.(tmp), 'norm')        
                log_priors(iP) = log( normcdf(params(iP),prior{iP}(1),prior{iP}(2)));
            else  % gamma
                log_priors(iP) = log( gamcdf(params(iP),prior{iP}(1),prior{iP}(2)));
            end
        end
    end
        
    lik = -sum(log(p)) + sum(log_priors);
end

