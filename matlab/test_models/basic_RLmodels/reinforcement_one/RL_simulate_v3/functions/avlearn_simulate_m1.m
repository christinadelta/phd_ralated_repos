function [data] = avlearn_simulate_m1(params, condition, probs, trials)

% Rescola-Wagner model simulation for the simplest version of the aversive
% learning task (simple 2-armed bandit).

% this function simulates data and model for one participant and 1 condition only
% (stable condition) 

% extract parameter values 
alpha           = params(1);
beta            = params(2); 

% initialise Q values for the two choices 
initQ           = [.5 .5];
Qvals           = initQ;

%% simulate task data

data            = {}; % init data structure

volatility      = 'stable';
runs            = 1; % no switch of probabilities 
prob            = probs(1); % what is the high probability?
nstim           = 2; % left and right shape
runtrials       = trials/runs;

% create feedback sequence (this will change when we'll add volatility
% component)
probtrials      = [prob 1-prob];

% simulate a trial list (a sequence with probabilities and outcomes)
counter         = 0; % init counter 

for r = 1:runs

    for x = 1:runtrials

        counter                     = counter+1; % update counter
        feedbackprob(counter,1)     = probtrials(r);
        feedback(counter,1)         = double(rand(1) <= feedbackprob(counter)); % if 1 = (high) aversive outcome, if 0 = no (for now) aversive outcome

    end % end of runtrials loop
end % end of runs loop

% if feedback(:,1) is 1 then feedback(:,2) is 0 
feedback(counter,2)                 = 1 - feedback;









end