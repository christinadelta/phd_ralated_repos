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



% extract info from the data structure that are needed for modelling 
choice          = modelout.a;
outcome         = modelout.reward;
trials          = size(choice,2);    % total trials
nchoices        = 2;                 % 2 choice options 
nll             = 0;                 % init nll 

%%  Initialize parameters and estimates

vA              = 0;
vB              = 0;
values          = zeros(trials, nchoices); % column 1 for vA, column 2 for vB
p               = zeros(trials,nchoices);

%% simualte model 

for trl = 1:trials

    c       = choice(trl);  % extract this trial's choice
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
        vB = vB * decay;                    % apply decay parameter to option B
    else
        vB = vB + alpha * (outcome(trl) - vB); % update vB
        vA = vA * decay;                    % apply decay parameter to option B
    end

    % store value estimates
    values(trl, 1)  = vA;
    values(trl, 2)  = vB;
    p(trl, 1)       = PA;
    p(trl, 2)       = PB;

end % end of trials loop

end % end of function