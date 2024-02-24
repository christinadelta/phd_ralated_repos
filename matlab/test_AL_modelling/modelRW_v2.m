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
good_options     = data.oR;                     % good option outcomes (low prob of losing) 
% outcome(:,2)    = data.o;                     % bad oprion outcomes (high probability of losing)
trials          = size(good_options(:,1),1);    % total trials
nchoices        = data.cues;                    % 2 choice options 

% convert outcomes to 1s and 2s
good_options(find(good_options == 0)) = 2; % good option outcomes (low prob of losing), 1 = red, 2 = blue


%%  Initialize parameters and estimates

vA              = 0;
vB              = 0;
values          = zeros(trials, nchoices); % column 1 for vA, column 2 for vB
choices         = zeros(1, trials);
p               = zeros(trials,nchoices);

% define the reward/loss values
rewardvalue     = 1;    % no loss (reward for choosing the good option) 
lossvalue       = -1;   % loss for choosing the bad option

%% simualte model 

for trl = 1:trials

    PA      = exp(beta * vA) / (exp(beta * vA) + exp(beta * vB)); % Since PA + PB = 1   % compute choice probabilities using the softmax function
    PB      = 1 - PA; 

    % make a choice based on probabilities
    if rand < PA
        chosen = 1; % option A chosen
    else
        chosen = 2; % option B chosen
    end
    choices(trl) = chosen;

    % determine reward based on choice and good option
    if chosen == good_options(trl)
        r(trl)          = rewardvalue; % chose the good option, no loss
        correct(trl)    = 1;
    else
        r(trl)          = lossvalue; % chose the bad option, incur loss
        correct(trl)    = 0;
    end

    % update value estimate for chosen option based on reward
    if chosen == 1
        vA = vA + alpha * (r(trl) - vA);    % update vA
        vB = vB * decay;                    % apply decay parameter to option B
    else
        vB = vB + alpha * (r(trl) - vB);    % update vB
        vA = vA * decay;                    % apply decay parameter to option B
    end

    % store value estimates
    values(trl, 1) = vA;
    values(trl, 2) = vB;
    p(trl, 1) = PA;
    p(trl, 2) = PB;

end % end of trials loop

%% from now on the rest is for plotting 

% store output parameters in the output structure 
modelout.Qvals      = values;
modelout.allPs      = p;
modelout.a          = choices;
modelout.reward     = r;
modelout.correct    = correct;

%% store model stuff in table for faster reading

choices = choices'; r = r'; % transpose choices and rewards
c = correct';
modeltable          = table(choices, r, values, p,c);
% filename            = sprintf('model_table_%s.xlsx', data.volatility); % save table as xlsx file
% 
% writetable(modeltable,filename, 'Sheet', 1) 
% movefile('*.xlsx', outpath) % move file to output dir 

modelout.modeltable = modeltable;



end