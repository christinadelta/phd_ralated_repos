function [data] = avlearn_simulate_m1(condition, probs, trials)

% Date created : 7/1/2023
% modeified: 14/1/2023

% this function runs through the aversive_learn_part1.m script.

% It simulates data both for modelling and for visuaisation purposes 



% thsi function runs through 

%% simulate task data

data            = {};           % init data structure

volatility      = 'stable';
runs            = 1;            % no switch of probabilities 
prob            = probs(1);     % what is the high probability?
nstim           = 2;            % left and right shape
runtrials       = trials/runs;  % trials in each run

% create feedback sequence (this will change when we'll add volatility
% component)
if condition == 1
    probtrials      = [prob 1-prob];
end

% simulate a trial list (a sequence with probabilities and outcomes)
counter         = 0; % init counter 

for r = 1:runs

    for x = 1:runtrials

        counter                     = counter+1; % update counter
        feedbackprob(counter,1)     = probtrials(r);
        feedback(counter,1)         = double(rand(1) <= feedbackprob(counter)); % if 1 = (high) rewarded outcome, if 0 = no (for now) rewarded outcome

    end % end of runtrials loop
end % end of runs loop

% if feedback(:,1)--vertical column is 1, then feedback(:,2)--horizontal column is 0 
feedback(:,2)                 = 1 - feedback;

% update the data struct with the needed info
data.nstim                      = nstim;


end