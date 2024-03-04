function modelout = modelRW_v3(params, data)

% Adjusted from version 2 to include risk aversion parameter (rho)
% Now params is expected to be a 1-by-4 vector: [alpha, beta, decay, rho]

% The purpose of Rho is to model how different individuals value potential 
% losses versus the absence of loss, influencing their decision-making 
% process. Higher values of Rho would indicate higher sensitivity to losses, 
% potentially leading to more conservative choice behaviours to avoid loss.

%% Extract parameters and initialize variables
alpha           = params(1);    % Learning rate
beta            = params(2);    % Inverse temperature (decision noise)
decay           = params(3);    % Decay rate for non-chosen option
rho             = params(4);    % Risk aversion parameter

good_options    = data.o;               % Good option outcomes (low prob of losing)
trials          = size(good_options,1); % Total trials
nchoices        = data.cues;            % Number of choice options (2)

% Initialize value estimates and choice probabilities
vA              = 0; 
vB              = 0;                        % Initial values for options A and B
values          = zeros(trials, nchoices);  % To store value estimates
p               = zeros(trials, nchoices);  % To store choice probabilities
simulatedChoice = zeros(1, trials);         % To store simulated choices
reward          = zeros(1, trials);         % To store outcomes (rewards/losses)

% convert outcomes to 1s and 2s
good_options(find(good_options == 0)) = 2; % good option o

%% Simulation loop

for trl = 1:trials

    PA              = exp(beta * vA) / (exp(beta * vA) + exp(beta * vB)); % Compute choice probabilities
    PB              = 1 - PA;
    simulatedChoice = rand() < PA; % Simulate choice based on PA
    isGoodChoice    = (simulatedChoice + 1 == good_options(trl));

    % Adjust reward/loss value based on risk aversion for losses
    if isGoodChoice
        updateValue = 1 % No loss (reward for choosing the good option)
    else
        updateValue = -1 * rho % Adjusted loss value, scaled by rho
    end

    reward(trl)     = updateValue; % Store the outcome

    % Update values with decay for the non-chosen option
    if simulatedChoice
        vA          = vA + alpha * (updateValue - vA); % Update for choosing A
        vB          = vB * decay; % Decay value of B
    else
        vB          = vB + alpha * (updateValue - vB); % Update for choosing B
        vA          = vA * decay; % Decay value of A
    end

    % Store updates
    values(trl, :)                  = [vA, vB];
    p(trl, :)                       = [PA, PB];
    modelout.simulated_actions(trl) = simulatedChoice + 1;  % Adjust for indexing (1 or 2)
    modelout.correct(trl)           = isGoodChoice;         % Store if the choice was good
end

%% Store output parameters in the structure for plotting and analysis

modelout.Qvals                      = values;
modelout.allPs                      = p;
modelout.reward                     = reward;
modelout.outcome                    = good_options;

end
