function [modelout] = modelRW_v1(params, data, outpath)

% Date created : 7/1/2023
% modeified 1: 15/1/2023

% Rescola-Wagner model for the simplest version of the aversive
% learning task (simple 2-armed bandit).

% this function simulates model for one participant and 1 condition only
% (stable condition) for now

% INPUTS: 
%       - params (1-by-2) with alpha  and beta values for the model
%       - data structure obrained during parameter/data simulation

% OUTPUT:
%       - structure with model results:
%           * Qvals 
%           * probability choices
%           * ll (log likelihood)
%           * choices (actions)
%           * reward

% feedback/outcome generated in the simulate function and reward are also used for
% visualisation 

%% ------------------

% initialise Q values for the two choices 
initQ           = [.5 .5];
Qvals           = initQ;

% extract parameter values 
alpha           = params(1);
beta            = params(2); 

% extract info from the data structure that are needed for modelling 
outcome         = data.feedback;    % only to generate rewards (since feedback is 1 or 0 for each of the gabors)
trials          = data.trials;      % trials per condition
nchoices        = data.nstim;       % 2 choice options 

% init log likelihood, choice and outcome vectors
ll              = 0; % this will be needed for model fitting mainly

% run the model
for trl = 1:trials

    p           = exp(beta*Qvals) / sum(exp(beta*Qvals));   % compute choice probabilities using the softmax function
    a(trl)      = max(find([-eps cumsum(p)] < rand));       % way 1 of computing choices

    % way 2 of computing choices:
%     if rand(1) < p(1)
%         a(trl)  = 1;
%     else
%         a(trl)  = 2;
%     end

    r(trl)          = outcome(trl,a(trl));                  % select outcome based on choice made

    % store stuff
    ll              = ll + log(p(a(trl)));                  % update ll
    allQs(trl,:)    = Qvals;                                % store all Qvalues 
    allPs(trl,:)    = p;                                    % store probabilities (for visualisation)

    % update stuff
    delta           = r(trl) - Qvals(a(trl));               % compute prediction error
    Qvals(a(trl))   = Qvals(a(trl)) + alpha * delta;        % update Q values (this is what the model learns based on the feedback it receives)

end % end of trials loop

%% the rest from now on is for plotting 

% store feedback that "subject recieved" based on the choice made
a1                  = a == 1;
a2                  = a == 2;
reward(a1)          = outcome(a1,1);
reward(a2)          = outcome(a2,2);

% store output parameters in the output structure 
modelout.Qvals      = allQs;
modelout.allPs      = allPs;
modelout.a          = a;
modelout.r          = r;
modelout.reward     = reward;

%% store model stuff in table for faster reading

a = a'; r = r'; % transpose choices and rewards
modeltable          = table(a, r, allQs, allPs);
% filename            = sprintf('model_table_%s.xlsx', data.volatility); % save table as xlsx file
% 
% writetable(modeltable,filename, 'Sheet', 1) 
% movefile('*.xlsx', outpath) % move file to output dir 

modelout.modeltable = modeltable;


end % end of function