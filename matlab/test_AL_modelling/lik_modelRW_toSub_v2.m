function nll = lik_modelRW_toSub_v2(params, modelout)

% objective function for version 1 of the RW model using the AL dataset
% Date created : 20/02/2024

% INPUTS: 
%       - params (1-by-2) with alpha  and beta values for the model
%       - output from model simulation (we need model choices and rewards)

% OUTPUT:
%           * nll (log likelihood) 



%% extract parameters and init variables 

alpha       = params(1); % alpha
beta        = params(2); % beta
decay       = params(3); % decay 

% Extract info from the data structure needed for modelling 
choice          = modelout.actions; % Actual choices made by participants
good_options    = modelout.outcome(:,1); % Good options indicating no loss
trials          = size(choice,1); % Total trials

nll             = 0; % Initialize nll 

% Initialize Q-values for both options
vA              = 0.5;
vB              = 0.5;

% Define reward/loss values 
rewardvalue = 1; % no loss (reward for choosing the good option)
lossvalue = -1; % loss for choosing the bad option

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
        vA               = vA + alpha * (updateValue - vA);
        vB               = vB * (1- decay);
        % vB               = vB - decay * (updateValue - vA); % Joe's suggestion
    else
        vB               = vB + alpha * (updateValue - vB);
        vA               = vA * (1 - decay) ; % apply decay parameter to option A
        % vA               = vA - decay * (updateValue - vB); % Joe's suggestion
    end

end % en of trials loop
end % end of function