function [data, output] = simulateRL(params, cond, sub)

% @christinadelta

% this function simulates data for every subject seperately
% subjects with odd number are included in the high volatility condition
% subjects with even number in the low volatility 

% VERSION 1 created in November 2022
% VERSION 2 modified in December 2022:
% - added task versions
% - corrected probabilities 


% ----------------------------
%% first part - prepare some initial parameters
% unpack parameters
alpha       = params(1);
beta        = params(2);

% initial value for each stimulus?
initv       = ones(1,2)*0.5;    % 0.5 for each stimulus
v           = initv;            % initialise v value
data        = struct;           % initialise output stucture

% --------------------
%% simulate some data 

volatility      = {'high', 'low'};
taskcond        = cond; % if 1=high, if 2=low

if cond == 1
    
    % create high volatility sequence data
    sequences   = [23 23 23 23 23 23 23 23];
else
    % create low volatility sequence data
    sequences   = [46 46 46 46];    
end

thisvol         = volatility{cond}; % high or low?
prb             = 85/100;           % (tone contigency) or is it outcome contigency??
nstim           = 2;                % how many stimuli?
totaltrials     = sum(sequences);   % for now 184 per condition
blocks          = length(sequences);% 8 (if high volatility, 4 if low volatility)

% generate probability and outcome sequences
trialprob       = repmat([prb 1-prb],1,blocks/2);
cnt             = 0; % counter

for block = 1:blocks
    for trl = 1:sequences(block)
        
        cnt                 = cnt+1; % update counter
        outcomeprob(cnt)    = trialprob(block);
        outcome(cnt)        = double(rand(1) <= outcomeprob(cnt)); % if random number (0 to 1) is smaller or equal to outcome probability, add 1)
        taskblock(cnt)      = block;
        
    end % end of trials loop
end % end of blocks loop

outcomeprob                 = outcomeprob(:);
outcome                     = outcome(:);
outcome(:,2)                = 1-outcome; 
taskblock                   = taskblock(:);
nchoice                     = nstim; % 2 choice options (house vs face)

% store in data struct
sim.outcomeprob         = outcomeprob;
sim.outcome             = outcome;
sim.nchoice             = nchoice;
sim.prob                = prb;
data.simdata            = sim;

% ---------------------
%% create choice and outcome vectors and store data, model values

ll                      = 0; % initiate log-likelihood
% data.choice             = nan(totaltrials,1); % may not need this now
% data.out                = nan(totaltrials,1); % may not need this now
vv                      = nan(totaltrials, nchoice); % matrix to store stimulus values for each trial
pp                      = nan(size(vv)); % matrix to store choice probabilities

for trl = 1:totaltrials
    
    % compute likelihood for each choice option
    exv                 = exp(beta * v);    % exponentiate initial value
    sexv                = sum(exv);         % calculate sum of values
    p                   = exv/sexv;         % calculate probability of each choice based on values
    
    % weighted coinflip to make a choice (e.g., use stim 1 (house) if probability for stim 1
    % is greater than tmp value).
    tmp                 = rand(1);
    if tmp < p(1)
        thischoice      = 1; % choose stimulus 1 (house)
    else % if 
        thischoice      = 2; % choose stimulus 2 (face)
    end
    
    reward              = outcome(trl,thischoice);      % select outcome -- this trial's choice
    ll                  = ll + log(p(thischoice));      % update log likelihood
    vv(trl,:)           = v;                            % start the values vector with initial v value
    pp(trl,:)           = p;                            % probability of trl choice
    
    % update values
    pe                  = reward - v(thischoice);       % compute prediction error
    v(thischoice)       = v(thischoice) + alpha * pe;   % update value for thischoice
    data.choices(trl)   = thischoice;                   % store trl choice in data structure
 
end % end of trial loop

choice1                 = data.choices == 1;            % stimulus 1 (i.e., house)
choice2                 = data.choices == 2;            % stimulus 2 (i.e., face)
data.out(choice1)       = outcome(choice1, 1);
data.out(choice2)       = outcome(choice2, 2);

% store  model values output
output.vv               = vv;       % v values
output.pp               = pp;       % probabilities of choices
output.ll               = ll;       % all lls


end % end of function