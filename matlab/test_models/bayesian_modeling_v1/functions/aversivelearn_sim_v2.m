function data = aversivelearn_sim_v2(cond, probs, trials, outpath, task)

% Date created : 7/1/2023
% modified 1: 13/1/2023
% modified 2: 15/1/2023 -- version 2
% modified 3: 25/1/2023
% modified 4: 2/2/2023

% this function runs through the aversive_learn_part1.m and aversivelearn_main_v1.m (bayesian modelling) scripts.
% It simulates data both for modelling and for visualisation purposes 
% Details of the model are given in the modeling function
% details about the visualisation are given in the vis function.

% -------
% Details related to the task simulation:

% version 1 simulates only stable condition trials for a 2-armed bandit
% task:
% in Version1 of the task, participants are presented with 2 gabor stimuli
% on the left and right side of the screen. One gabor shape has vertical
% lines and the other one has horizontal lines. 

% the two gabors are associated with certain probabilities of recieving an
% reward. The probabilities are: 
% a) 75% and 25% in the stable condition. 
% b) 80% and 20% in the volatile condition
% This means that there is always a good option and a bad option. The good option is the
% one with the large probability. 
% 
% The reasoning here, is very similar to that of the reward-based version (see
% version 0 of the simulation and model). 

% For now (modification 1), I only simulate trials for the stable
% condition. As I learn the model, I will continue adding complexity to the
% simulation process and to the model. 

% NEW CHANGE (15/1/2023) included in Version 2:
% I have now added the volatile version too

% NEW CHANGE (25/1/2023):
% changed the way feedback is generated 

% This function takes as INPUTS:
%                              - condition (1-by-n array) if n = 1 then
%                               we simulate data only for the stable condition (no change of probabilities in the
%                               trials).If n = 2, then we also include high volatility trials where probailities switch frequently.

%                              - probs (a 1-by-2 array or 2-by-2 array) with probabilities. 
%                              - number of trials (i.e., 100) per condition
%                              - task version. In stable condition I have 2
%                                versions: if task = 1, then there is no
%                                switch in the 100 trials, if task = 2,
%                                there is one switch after 50 trials

% the function uses those inputs to generate trials list/sequences and
% stores them in the output structure: 
%                                       - data

% data structure cntains information that is used only for plotting
% visualisation and for modelling. 

% info used for modelling:
% 1. feedback -- a trials-by-cues vector with columns of 1 and 0 (column 1 =
% feedback for vertical gabor and column 2 = feedback for horizontal
% gabor). When column 1 row = 1, column 2 corresponding row = 0. 
% 2. number of stimuli -- this is 1-by-1 double to designate the number of
% stimuli used (i.e., 2).
% 3. number of trials (i.e., 100) per condition 

% info used for plotting:
% feedback probability (n-by-1 array), where n = number of trials (per
% condition), with the true-underlying probabilities across trials 

%% -------------------
% simulate task data

data            = {};           % init data structure

if cond == 1 % stable condition

    volatility      = 'stable'; % is it volatile or stable simulation?

    % is it stable with switch or without switch?
    if task == 1
        runs            = 1;            % no switch of probabilities 
    else
        runs            = 2;            % probabilities switch after 50 trials
    end

    prob                = probs(1);     % what is the high probability?
    nstim               = 2;            % left and right shape
    runtrials           = trials/runs;  % trials in each run
    
    % create trials sequence (this will change when we'll add volatility
    % component)
    probtrials(1,:)  = [prob 1-prob]; % feedback sequence (vertical gabor) probabilities
    probtrials(2,:)  = [1-prob prob]; % feedback sequence (horizontal gabor) probabilities

    if runs == 1
        seqtrials{1}   = trials;
    else
        seqtrials{1}   = [trials trials]/runs; % this wont probably be used but good to have it as an option for testing 
    
    end

else % if it is volatile condition 

    % in the volatile condition the probabilities are: 80 - 20 and the
    % switch evry 25 trials (e.g., 4 times if trials = 100)

    volatility          = 'volatile';   % volatile simulation!
    runtrials           = 25; 
    runs                = trials/runtrials;
    prob                = probs(2);     % for volatile condition the probs are .8 and .2

    nstim               = 2;            % left and right shape

    probtrials(1,:)     = repmat([1-prob prob], 1, runs/2); % feedback sequence (vertical gabor) probabilities
    probtrials(2,:)     = repmat([prob 1-prob], 1, runs/2); % feedback sequence (horizontal gabor) probabilities

    seqtrials{2}        = repmat(([runtrials]), 1, runs);   % trials per run for voaltile condition

end % end of condition loop

% simulate a trial list (a sequence with probabilities and outcomes
for r = 1:runs

    counter         = 0; % init counter 

    for x = 1:seqtrials{1,cond}(r)

        counter                         = counter+1; % update counter
        feedbackprob(counter,r)         = probtrials(1,r); % feedback probability for vertical gabor
        % feedback(counter,1)         = double(rand(1) <= feedbackprob(counter)); % if 1 = (high) aversive outcome, if 0 = no (for now) aversive outcome

    end % end of runtrials loop
end % end of runs loop

% There are two ways to generate feedback sequence: either randomly or it
% can be fixed. 
% If runtrials*prob(high) is  whole number, then we'll generate a fixed sequence (e.g., 75% of trials will be high rewarded-vertical and 25% low reward-vertical) 
% If runtrials*prob(high) is NOT a whole number, then we'll use rand to
% generate sequence of high and low rewarded options 

% decide whether the sequence generation will be random or fixed:
vtrials     = runtrials*prob; % is the number of high probability-vertical trials whole?

if rem(vtrials,1) == 0
    rdm     = 0; % don't generate random sequence
else
    rdm     = 1; 
end

% this will be done independently for each option/column:
for s = 1:nstim  

    for run = 1:runs

        feedback{1,s}(:,run) = computeFeedback(1:seqtrials{1,cond}(run),probtrials(s,run),rdm);

    end % end of runs loop
end % end of stimuli loop

% add feedback for vertical and horizontal gabors in one column
outcome(:,1)        = feedback{1,1}(:);
outcome(:,2)        = feedback{1,2}(:);
feedbackprob        = feedbackprob(:);
% if feedback(:,1)--vertical column is 1, then feedback(:,2)--horizontal column is 0 
%feedback(:,2)                       = 1 - feedback;

% update the data struct with the needed info
data.nstim          = nstim;
data.feedbackprob   = feedbackprob;
data.feedback       = outcome;
data.trials         = trials;
data.condition      = cond;
data.probs          = probs;
data.volatility     = volatility;

%% generate table for exporting 

% this is not required (just for faster reading (if needed) 
feedbackv           = outcome(:,1);
feedbackh           = outcome(:,2);

% add condition
if cond == 1

    condition = ones(trials,1);
else
    condition = ones(trials,1)*2;

end

% create table
simdata_onesub      = table(condition, feedbackprob, feedbackv, feedbackh);
% 
% % store table in .xlsx format
% filename = sprintf('simdata_onesub_%s.xlsx', volatility);
% writetable(simdata_onesub,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

data.datatable = simdata_onesub;

end % end of fuction
