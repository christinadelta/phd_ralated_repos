function [data, output] = simulateRL_v2(params, cond, sub)

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

if cond == 1
    
    % create high volatility sequence data
    sequences   = [23 23 23 23 23 23 23 23];
else
    % create low volatility sequence data
    sequences   = [46 46 46 46];    
end

thisvol         = volatility{cond};     % high or low?
prob            = 85/100;               % (tone contigency) or is it outcome contigency??
ncues           = 2;                    % how many stimuli [1=house, 2=face]?
totaltrials     = sum(sequences);       % for now 184 per condition
blocks          = length(sequences);    % 8 (if high volatility, 4 if low volatility)
blocktrials     = totaltrials / blocks; % number of trials per block? 

% generate probability and outcome sequences
blockprob       = repmat([prob 1-prob],1,blocks/2); % is it a 0.85 0r 0.15 probability block? 
cuetrials(1)    = round(blocktrials * prob);        % within block (majority cue - 0.85%)
cuetrials(2)    = round(blocktrials * (1-prob));    % within block (minority cue - 0.15%)

for block = 1:blocks
    
    for trl = 1:sequences(block)
        
        % start with cue probabilities
        outcomeprob(trl,1)      = blockprob(block);
        
        % so, in odd number blocks:
        % 85% of blocktrials show cue 1 and the rest 15% show cue 2, 
        % and in even number blocks: 
        % 85% of blocktrials show cue 2 and the rest 15% show cue 1 
        if mod(block,2) == 1 % blocks: 1,3,5,7
            
            if trl <= cuetrials(1) % if this is trial 1-20 (85% of 23)
                cue             = 1;
                cuep = 1; % high probability (85%)
            else % if this is trial 21-23 (15% of 23)
                cue             = 2;
                cuep = 2; % low probability (15%
            end
            
        else % blocks: 2,4,6,8
            
            if trl <= cuetrials(1)
                cue             = 2; 
                cuep            = 1; %
            else
                cue             = 1;
                cuep            = 2; 
            end
            
        end
        
        % add tone based on probability
        if trl <= cuetrials(1) 
            blocktone           = 1; % high pitch tone 
        else
            blocktone           = 2; % low pitch tone
        end
        
        outcome(trl,1)          = cue;
        outcome(trl,2)          = blocktone;
        outcome(trl,3)          = cuep;
        taskblock(trl,1)        = block;
        taskcond(trl,1)         = cond;

       
    end % end of trials loop 
    
    % shuffle outcomes/cues 
    temp        = randperm(numel(outcome(:,1)));
    orderedout  = outcome(temp,:);
    
    probsperblock(:,block)  = outcomeprob; % each columns is a block
    blockcues(:,block)      = orderedout(:,1);
    blocktones(:,block)     = orderedout(:,2);
    cueprobs(:,block)       = orderedout(:,3);
    thisblock(:,block)      = taskblock;
    volatilitys(:,block)    = taskcond;
    
end % end of block loop

% flatten matricies
probabilitys                = probsperblock(:);
allcues                     = blockcues(:);
allblocks                   = thisblock(:);
alltones                    = blocktones(:);
allcueprobs                 = cueprobs(:);
allvol                      = volatilitys(:);

% store in data structure
sim.blockprob               = probabilitys;
sim.cue                     = allcues;
sim.blocknb                 = allblocks;
sim.tone                    = alltones;
sim.cueprob                 = allcueprobs;
sim.ncues                   = ncues;
data.volatile               = allvol;
data.simdata                = sim;

% ---------------------
%% simple RL model
% create choice and outcome vectors and store data, model values

ll                          = 0;                            % initiate log-likelihood
vv                          = nan(totaltrials, ncues);      % matrix to store stimulus values for each trial
pp                          = nan(size(vv));                % matrix to store choice probabilities

for trl = 1:totaltrials
    
    % compute likelihood for each choice option
    expv                = exp(beta * v);    % exponentiate initial value
    sexpv               = sum(expv);        % calculate sum of values
    p                   = expv/sexpv;       % calculate probability of each choice based on values
    
    % weighted coinflip to make a choice/prediction (e.g., use stim 1 (house) if probability for stim 1
    % is greater than tmp value).
    tmp                 = rand(1);
    if tmp < p(1)
        thischoice      = 1; % choose stimulus 1 (house)
    else % if 
        thischoice      = 2; % choose stimulus 2 (face)
    end
    
    allpredictions(trl,1)   = thischoice; % save this trial's prediction
    
    reward                  = allcues(trl,1);   % select outcome -- this trial's choice
    ll                      = ll + log(p(thischoice));      % update log likelihood
    vv(trl,:)               = v;                            % start the values vector with initial v value
    pp(trl,:)               = p;                            % probability of trl choice
    
    % update values
    pe                  = reward - v(thischoice);       % compute prediction error
    v(thischoice)       = v(thischoice) + alpha * pe;   % update value for thischoice
    
end

data.predictions        = allpredictions; % store in data structure 

allchoices              = allpredictions == allcues; % was the prediction valid? 1 =correct prediction, 0 = incorrect prediction
data.correct            = allchoices;

% store  model values output
output.vv               = vv;       % v values
output.pp               = pp;       % probabilities of choices
output.ll               = ll;       % all lls


end % end of function