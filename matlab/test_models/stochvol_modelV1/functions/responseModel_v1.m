%function [x, dp, tstable, tvolatile] = responseModel_v1(xstate, val,trialvolatile, trialstable, beta)
function [xx, mR, dp, simResp] = responseModel_v1(xstate, val,tvolatile, tstable, o, beta)

% Created May 2023
% The function generates responses for the action-learning simulated data 
% to be used for modelling 

% inputs are: Qvals , beta (used in the softmax function) the true reward
% rate and indexes for stable and volatile trials to be used to compute
% performance for each condition seperately 

% ---------------

N = size(xstate,1);
cp = find(diff(xstate(:,1))~=0);

%cp(:,2) = 1:length(cp)


% loop over trials 
for t = 1:size(xstate,1)

    mu      = xstate(t,:);
    Qval    = val(t,:);
    
    % compute choice probabilities (i.e., convert Qs into choice
    % probabiltiies)
    p           = exp(beta*Qval) / sum(exp(beta*Qval));
    a(t)        = mkchoice(p);      % make choice according to the choice probabilities
    % r(t)    = rand < mu(a(t));  
    r(t)        = o(t,a(t)); % generate reward based based on action and outcome
    correct(t)  = o(t,a(t));
 
end % end of trials loop

%% compute performance
a1          = a==1;
a2          = a==2;

reward(a1)  = o(a1,1);
reward(a2)  = o(a2,2);

reward      = reward';
correct     = correct';

% compute mean reward for each of the conditions
sR              = mean(reward(tstable,:));
vR              = mean(reward(tvolatile,:));

sPerf           = mean(correct(tstable,:));
vPerf           = mean(correct(tvolatile,:));

mR              = [sR vR];

dp              = sPerf - vPerf;
xx              = [sPerf vPerf];

simResp.r       = r;
simResp.a       = a;
simResp.reward  = reward;



end
