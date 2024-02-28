function modelout = modelRW_v2(params, data)

% RW model simulation V1 

% Date created : 14/02/2024

% Rescola-Wagner model with 2 parameters to test the aversive learning task
% this v1 function simulates model for one participant and 1 condition only
% (small stochasticity/stable condition) 

% INPUTS: 
%       - params (1-by-2) with alpha  and beta values for the model
%       - data structure with simulated data

% OUTPUT:
%       - structure with model results:
%           * Qvals 
%           * probability choices
%           * ll (log likelihood) % not yet 
%           * choices (actions)
%           * reward

% feedback/outcome generated in the simualtion function and reward are also used for
% visualisation 

%% extract parameters and init variables 

% extract parameter values 
alpha           = params(1); % alpha
beta            = params(2); % beta
decay           = params(3); % decay

% extract info from the data structure that are needed for modelling 
good_options     = data.o;                      % good option outcomes (low prob of losing) 
% outcome(:,2)    = data.o;                     % bad oprion outcomes (high probability of losing)
trials          = size(good_options(:,1),1);    % total trials
nchoices        = data.cues;                    % 2 choice options 

% convert outcomes to 1s and 2s
good_options(find(good_options == 0)) = 2; % good option outcomes (low prob of losing), 1 = red, 2 = blue


%%  Initialize parameters and estimates

vA              = 0;
vB              = 0;
values          = zeros(trials, nchoices); % column 1 for vA, column 2 for vB
simulatedChoice         = zeros(1, trials);
p               = zeros(trials,nchoices);

% define the reward/loss values
rewardvalue     = 1;    % no loss (reward for choosing the good option) 
lossvalue       = -1;   % loss for choosing the bad option

%% simualte model 

for trl = 1:trials

    PA      = exp(beta * vA) / (exp(beta * vA) + exp(beta * vB)); % Since PA + PB = 1   % compute choice probabilities using the softmax function
    PB      = 1 - PA; 

    simulatedChoice = rand() < PA; % Simulate choice based on PA
    
    % Use the simulated choice for NLL calculation
    if simulatedChoice
        chosenprob = PA;
    else
        chosenprob = PB;

    end

    % is the current choice a good choice?
    isGoodChoice        = (simulatedChoice + 1 == good_options(trl));

    % update the update value and store correct/incorrect response 
    if isGoodChoice
        updateValue     = rewardvalue;
    else
        updateValue     = lossvalue;
    end

    % update the reward array
    reward(trl)         = updateValue;

    if simulatedChoice
        vA               = vA + alpha * (updateValue - vA);
        vB               = vB * decay; % apply decay parameter to option B
    else
        vB               = vB + alpha * (updateValue - vB);
        vA               = vA * decay; % apply decay parameter to option A
    end

    % Store simulated choices and outcomes
    modelout.simulated_actions(trl)  = simulatedChoice + 1; % Adjust for indexing
    modelout.correct(trl)            = isGoodChoice; % Store if the choice was good
    
    % store value estimates
    values(trl, 1)      = vA;
    values(trl, 2)      = vB;
    p(trl, 1)           = PA;
    p(trl, 2)           = PB;
  

end % end of trials loop

%% from now on the rest is for plotting 

% store output parameters in the output structure 
modelout.Qvals      = values;
modelout.allPs      = p;
modelout.reward     = reward;
modelout.outcome    = good_options;



end