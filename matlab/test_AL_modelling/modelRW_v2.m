function modelout = modelRW_v2(params, data)

% RW model simulation V1 

% Date created : 14/02/2024
% Rescola-Wagner same as version with and added parameter: DECAY 

% What does the decay parameter do?:
% The **decay factor** measures how quickly the model "forgets" or 
% devalues past information in favour of more recent observations. 
% A high decay factor (close to 1) means that forgetting is slow, and 
% past information remains relevant for longer. This slow decay is 
% beneficial in environments where changes are infrequent or non-existent, 
% as it allows for a stable buildup of knowledge.
% A low decay factor (closer to 0) accelerates the forgetting process, 
% helping the model to rapidly adjust its expectations or value estimates 
% in response to new information. This is particularly useful in 
% environments that change frequently, requiring the model to de-emphasise 
% older, potentially obsolete information to maintain flexibility in its 
% decision-making strategy.


% INPUTS: 
%       - params (1-by-2) with alpha, beta and decay values for the model
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
        vB               = vB * (1- decay);
        % vB               = vB - decay * (updateValue - vA); % Joe's suggestion
    else
        vB               = vB + alpha * (updateValue - vB);
        vA               = vA * (1 - decay) ; % apply decay parameter to option A
        % vA               = vA - decay * (updateValue - vB); % Joe's suggestion
    end

    % Store simulated choices and outcomes
    modelout.simulated_actions(trl)  = simulatedChoice + 1; % Adjust for indexing
    modelout.correct(trl)            = isGoodChoice;        % Store if the choice was good
    
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