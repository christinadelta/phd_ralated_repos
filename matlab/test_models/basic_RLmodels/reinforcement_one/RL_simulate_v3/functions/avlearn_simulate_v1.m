function [data] = avlearn_simulate_v1(condition, probs, trials, outpath, task)

% Date created : 7/1/2023
% modified 1: 13/1/2023
% modified 2: 15/1/2023
% modified 3: 25/1/2023

% this function runs through the aversive_learn_part1.m script.
% It simulates data both for modelling and for visuaisation purposes 
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
% reward; and the probabilities are 75% and 25%. This means that
% there is always a good option and a bad option. The good option is the
% one with the large probability. 
% 
% The reasoning here, is very similar to that of the reward-based version (see
% version 0 of the simulation and model). 

% For now (modification 1), I only simulate trials for the stable
% condition. As I learn the model, I will continue adding complexity to the
% simulation process and to the model. 

% NEW CHANGE (15/1/2023):
% I have now added the volatile version too

% NEW CHANGE (25/1/2023):
% changed the way feedback is generated 

% This function takes as INPUTS:
%                              - condition (1-by-n array) if n = 1 then
%                               we simulate data only for the stable condition (no change of probabilities in the
%                               trials).If n = 2, then we also include high volatility trials where probailities switch frequently.

%                              - probs (a 1-by-2 array or 2-by-2 array) with probabilities. 
%                              - number of trials (i.e., 100) per condition

% the function uses those inputs to generate trials list/sequences and
% stores them in the output structure: 
%                                       - data

% data structure cntains information that is used only for plotting
% visualisation and for modelling. 

% info used for modelling:
% 1. feedback -- a trials-by-2 vector with columns of 1 and 0 (column 1 =
% feedback for vertical gabor and column 2 = feedback for horizontal
% gabor). When column 1 row = 1, column 2 corresponding row = 0. 
% 2. number of stimuli -- this is 1-by-1 double to designitate the number of
% stimuli used (i.e., 2).
% 3. number of trials (i.e., 100) per condition 

% info used for plotting:
% feedback probability (n-by-1 array), where n = number of trials (per
% condition), with the probabilities distribution across trials 

%% -------------------
% simulate task data

data            = {};           % init data structure

if condition == 1 % stable condition only

    volatility      = 'stable'; % is it volatile or stable simulation?

    % is it stable with switch or without switch?
    if task == 1
        runs        = 1;            % no switch of probabilities 
    else
        runs        = 2;            % probabilities switch after 50 trials
    end

    prob            = probs(1);     % what is the high probability?
    nstim           = 2;            % left and right shape
    runtrials       = trials/runs;  % trials in each run
    
    % create feedback sequence (this will change when we'll add volatility
    % component)
    if runs == 1
        probtrials  = [prob 1-prob];
    else
        probtrials  = repmat([prob 1-prob], 1, runs/2);
    
    end

else % if it is stable + volatile condition 

    % in the volatile condition the probabilities are: 80 - 20 and the
    % switch evry 25 trials (e.g., 4 times if trials = 100)

    volatility      = 'volatile'; % volatile simulation!
    runtrials       = 25; 
    runs            = trials/runtrials;
    prob            = probs(2); % for volatile condition the probs are .8 and .2

    nstim           = 2;            % left and right shape
    probtrials      = repmat([prob 1-prob], 1, runs/2); % feedback sequence

end % end of condition statement 

% simulate a trial list (a sequence with probabilities and outcomes)
counter         = 0; % init counter 

for r = 1:runs

    for x = 1:runtrials

        counter                     = counter+1; % update counter
        feedbackprob(counter,1)     = probtrials(r);
        % feedback(counter,1)         = double(rand(1) <= feedbackprob(counter)); % if 1 = (high) aversive outcome, if 0 = no (for now) aversive outcome

    end % end of runtrials loop
end % end of runs loop

% if feedback(:,1)--vertical column is 1, then feedback(:,2)--horizontal column is 0 
%feedback(:,2)                       = 1 - feedback;

% one way to generate feedback is the code above (in the x for loop)
% here I will use the exact probabilities to generate feedback sequences
% (e.g., 75% and 25%) -- don't forget that feedback determines the rewarded
% action 

if condition == 1 && task == 1
    tmp_feedbackv       = cat(1,ones(runtrials * prob,1), zeros(runtrials * (1-prob),1));
    shuffled            = randperm(length(tmp_feedbackv));% shuffle the feedback 
    feedbackv           = tmp_feedbackv(shuffled);
    
    clear shuffled tmp_feedbackv
    % let's create an independent column for the second stimulus (horizontal
    % gabor) 
    tmp_feedbackh       = cat(1,zeros(runtrials * prob,1), ones(runtrials * (1-prob),1));
    shuffled            = randperm(length(tmp_feedbackh));% shuffle the feedback 
    feedbackh           = tmp_feedbackh(shuffled);
    
    % add columns two one matrix
    % feedback             = cat(2, feedbackv, feedbackh);

elseif condition == 1 && task == 2

    % generate column 1 (vertical gabor reward feedback)
    % first create the 75% trials (for vertical gabor)
    tmp_feedbackv1       = cat(1,ones(ceil(runtrials * prob),1), zeros(floor(runtrials * (1-prob)),1));
    shuffled1            = randperm(length(tmp_feedbackv1));% shuffle the feedback 
    feedbackv1           = tmp_feedbackv1(shuffled1);

    % now create the 25% trials (for vertical gabor)
    tmp_feedbackv2       = cat(1,ones(ceil(runtrials * (1-prob)),1), zeros(floor(runtrials * prob),1));
    shuffled2            = randperm(length(tmp_feedbackv2));% shuffle the feedback 
    feedbackv2           = tmp_feedbackv2(shuffled2);

    % concatinate the two arrays
    feedbackv            = cat(1, feedbackv1, feedbackv2);

    clear shuffled1 shuffled2

    % generate column 2 (horizontal gabor reward feedback)
    % first create the 75% trials (for horizontal gabor)
    tmp_feedbackh1       = cat(1,zeros(ceil(runtrials * prob),1), ones(floor(runtrials * (1-prob)),1));
    shuffled1            = randperm(length(tmp_feedbackh1));% shuffle the feedback 
    feedbackh1           = tmp_feedbackh1(shuffled1);

    % now create the 25% trials (for horizontal gabor)
    tmp_feedbackh2       = cat(1,ones(ceil(runtrials * prob),1), zeros(floor(runtrials * (1-prob)),1));
    shuffled2            = randperm(length(tmp_feedbackh2));% shuffle the feedback 
    feedbackh2           = tmp_feedbackh2(shuffled2);

    % concatinate the two arrays
    feedbackh            = cat(1, feedbackh1, feedbackh2);

    clear shuffled1 shuffled2

else % if this is volatile condition

    % first generate feedback for column 1 (vertical gabor)
    % 80% feedback probability
    tmp_feedbackv1      = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled1           = randperm(length(tmp_feedbackv1)); % shuffle the feedback
    feedbackv1          = tmp_feedbackv1(shuffled1);

    % 20% feedback probability
    tmp_feedbackv2      = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled2           = randperm(length(tmp_feedbackv2)); % shuffle the feedback
    feedbackv2          = tmp_feedbackv2(shuffled2);

    % 80% feedback probability
    tmp_feedbackv3      = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled3           = randperm(length(tmp_feedbackv3)); % shuffle the feedback
    feedbackv3          = tmp_feedbackv3(shuffled3);

    % 20% feedback probability
    tmp_feedbackv4      = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled4           = randperm(length(tmp_feedbackv4)); % shuffle the feedback
    feedbackv4          = tmp_feedbackv4(shuffled4);

    % concatinate the two arrays
    feedbackv            = cat(1, feedbackv1, feedbackv2, feedbackv3, feedbackv4);

    clear shuffled1 shuffled2 shuffled3 shuffled4

    % now generate feedback for column 2 (horizontal gabor)
    % 80% feedback probability
    tmp_feedbackh1       = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled1            = randperm(length(tmp_feedbackh1));% shuffle the feedback 
    feedbackh1           = tmp_feedbackh1(shuffled1);

    % 20% feedback probability
    tmp_feedbackh2       = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled2            = randperm(length(tmp_feedbackh2));% shuffle the feedback 
    feedbackh2           = tmp_feedbackh2(shuffled2);
    
    % 80% feedback probability
    tmp_feedbackh3       = cat(1,zeros(ceil(runtrials * prob),1), ones(ceil(runtrials * (1-prob)),1));
    shuffled3            = randperm(length(tmp_feedbackh3));% shuffle the feedback 
    feedbackh3           = tmp_feedbackh3(shuffled3);

    % 20% feedback probability
    tmp_feedbackh4       = cat(1,ones(ceil(runtrials * prob),1), zeros(ceil(runtrials * (1-prob)),1));
    shuffled4            = randperm(length(tmp_feedbackh4));% shuffle the feedback 
    feedbackh4           = tmp_feedbackh4(shuffled4);

    % concatinate the two arrays
    feedbackh            = cat(1, feedbackh1, feedbackh2, feedbackh3, feedbackh4);

    clear shuffled1 shuffled2 shuffled3 shuffled4

end

% add columns to one matrix
feedback             = cat(2, feedbackv, feedbackh);


% update the data struct with the needed info
data.nstim          = nstim;
data.feedbackprob   = feedbackprob;
data.feedback       = feedback;
data.trials         = trials;
data.condition      = condition;
data.probs          = probs;
data.volatility     = volatility;

%% generate table for exporting 

% this is not required (just for faster reading (if needed) 
feedbackv           = feedback(:,1);
feedbackh           = feedback(:,2);

% create table
simdata_onesub      = table(feedbackprob, feedbackv, feedbackh);
% 
% % store table in .xlsx format
% filename = sprintf('simdata_onesub_%s.xlsx', volatility);
% writetable(simdata_onesub,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

data.datatable = simdata_onesub;


end % end of the function

