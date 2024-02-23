function [modelout, nll] = fit_modelRW_v1(params, modelout)

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

% extract info from the data structure that are needed for modelling 
choice          = modelout.a;
outcome         = modelout.reward;
trials          = size(choice,2);    % total trials
nchoices        = 2;                    % 2 choice options 
nll             = 0;                            % init nll 

%%  Initialize parameters and estimates

vA              = 0;
vB              = 0;
values          = zeros(trials, nchoices); % column 1 for vA, column 2 for vB
p               = zeros(trials,nchoices);

%% simualte model 

for trl = 1:trials

    c       = choice(trl); % extract this trial's choice
    r       = outcome(trl); % extract this trial's reward

    PA      = exp(beta * vA) / (exp(beta * vA) + exp(beta * vB)); % Since PA + PB = 1   % compute choice probabilities using the softmax function
    PB      = 1 - PA; 
    
    % after making a choice based on probabilities...
    if c == 1
        chosenprob  = PA;
    else
        chosenprob  = PB;
    end

    nll             = nll - log(chosenprob); % accumulate the NLL

   
    % update value estimate for chosen option based on reward
    if c == 1
        vA = vA + alpha * (outcome(trl) - vA); % update vA
    else
        vB = vB + alpha * (outcome(trl) - vB); % update vB
    end

    % store value estimates
    values(trl, 1)  = vA;
    values(trl, 2)  = vB;
    p(trl, 1)       = PA;
    p(trl, 2)       = PB;

end % end of trials loop

%% from now on the rest is for plotting 

% store output parameters in the output structure 
modelout.Qvals      = values;
modelout.allPs      = p;
modelout.a          = choice;
modelout.reward     = outcome;

end % end of function