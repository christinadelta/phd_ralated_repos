function nll = lik_modelRW_v2_sensAnalysis(params, fixedParams, modelout,m)

% objective function for version 1 of the RW model using the AL dataset
% Date created : 20/02/2024

% INPUTS: 
%       - params (1-by-2) with alpha  and beta values for the model
%       - output from model simulation (we need model choices and rewards)

% OUTPUT:
%           * nll (log likelihood) 



%% extract parameters and init variables 

% extract parameter values 
if m == 1
    alpha           = params(1); % alpha free
    beta            = fixedParams(1); % beta
    decay           = fixedParams(2); % decay
elseif m == 2
    alpha           = fixedParams(1); % alpha
    beta            = params(1); % beta free
    decay           = fixedParams(2); % decay
else
    alpha           = fixedParams(1); % alpha
    beta            = fixedParams(2); % beta
    decay           = params(1); % decay free
end


choice          = modelout.simulated_actions;
good_options    = modelout.outcome(:,1); % Good options indicating no loss
trials          = length(choice);
nll             = 0; % Initialize nll 
vA              = 0; % Initial Q-value for option A
vB              = 0; % Initial Q-value for option B

rewardvalue     = 1; % Reward for choosing the good option (no loss)
lossvalue       = -1; % Loss for choosing the bad option

% convert outcomes to 1s and 2s
good_options(find(good_options == 0)) = 2; % good option o

%% loop over trials and estimate values 

for trl = 1:trials
    PA              = exp(beta * vA) / (exp(beta * vA) + exp(beta * vB));
    PB              = 1 - PA;
    
    % Determine the probability of the choice actually made
    if choice(trl) == 1
        chosenprob  = PA;
    else
        chosenprob  = PB;
    end

    nll             = nll - log(chosenprob); % Update NLL based on the choice probability
    
    isGoodChoice    = (choice(trl) == good_options(trl));

    % is it a good option? if yes update with reward, if not update with
    % loss
    if isGoodChoice
        updateValue = rewardvalue;
    else
        updateValue = lossvalue;
    end

    
    if choice(trl) == 1
        vA = vA + alpha * (updateValue - vA);
        vB               = vB * decay; % apply decay parameter to option B
    else
        vB = vB + alpha * (updateValue - vB);
        vA               = vA * decay; % apply decay parameter to option A
    end

end % en of trials loop

end % end of function