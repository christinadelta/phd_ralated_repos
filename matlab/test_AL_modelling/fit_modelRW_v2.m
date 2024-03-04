function [modelout, nll] = fit_modelRW_v2(params, modelout)

% RW model simulation V1 

% Date created : 17/02/2024

% Rescola-Wagner model with 2 parameters to test the aversive learning task
% this v1 function fits model and computes nll for the entire task 

% INPUTS: 
%       - params (1-by-2) with alpha  and beta values for the model
%       - output from model simulation (we need model choices and rewards)

% OUTPUT:
%       - structure with model results:
%           * Qvals 
%           * probability choices
%           * ll (log likelihood) 
%           * choices (actions)
%           * reward

% feedback/outcome generated in the simualtion function and reward are also used for
% visualisation 

%% extract parameters and init variables 

% extract parameter values 
alpha           = params(1); % alpha
beta            = params(2); % beta
decay           = params(3); % decay

good_options    = modelout.outcome(:,1);
trials          = length(good_options);
nchoices        = 2;
nll             = 0;

% convert outcomes to 1s and 2s
good_options(find(good_options == 0)) = 2; % good option o

%%  Initialize parameters and estimates

vA              = 0;
vB              = 0;
values          = zeros(trials, nchoices); % column 1 for vA, column 2 for vB
p               = zeros(trials,nchoices);

% Define reward/loss values 
rewardvalue     = 1; % no loss (reward for choosing the good option)
lossvalue       = -1; % loss for choosing the bad option

%% simualte model 

for trl = 1:trials

    PA = exp(beta * vA) / (exp(beta * vA) + exp(beta * vB));
    PB = 1 - PA;
    
    simulatedChoice = rand() < PA; % Simulate choice based on PA
    
    % Use the simulated choice for NLL calculation
    if simulatedChoice
        chosenprob = PA;
    else
        chosenprob = PB;
    end
    nll = nll - log(chosenprob); % Update NLL based on simulated choice probability

    isGoodChoice        = (simulatedChoice + 1 == good_options(trl));

    % update the update value and store correct/incorrect response 
    if isGoodChoice
        updateValue     = rewardvalue;
    else
        updateValue     = lossvalue;
    end

    
    if simulatedChoice
        vA               = vA + alpha * (updateValue - vA);
        vB               = vB * (1- decay);
        % vB               = vB - decay * (updateValue - vA); % Joe's suggestion
    else
        vB               = vB + alpha * (updateValue - vB);
        vA               = vA * (1 - decay) ; % apply decay parameter to option A
        % vA               = vA - decay * (updateValue - vB); % Joe's suggestion
    end
    
    % Store simulated choices and outcomes
    modelout.simulated_actions(trl)  = simulatedChoice + 1; % Adjust for indexing
    modelout.correct(trl)            = isGoodChoice; % Store if the choice was good

    % store value estimates
    values(trl, 1)              = vA;
    values(trl, 2)              = vB;
    p(trl, 1)                   = PA;
    p(trl, 2)                   = PB;
    reward(trl)                 = updateValue;

end % end of trials loop

%% from now on the rest is for plotting 

% store output parameters in the output structure 
modelout.Qvals      = values;
modelout.allPs      = p;
modelout.reward     = reward;

end % end of function