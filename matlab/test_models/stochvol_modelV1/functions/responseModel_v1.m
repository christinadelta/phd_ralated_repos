function [xx, mR, dp, simResp] = responseModel_v1(xstate, val,tvolatile, tstable, beta)

% Created May 2023
% The function generates responses for the action-learning simulated data 
% to be used for modelling 

% inputs are: Qvals , beta (used in the softmax function) the true reward
% rate and indexes for stable and volatile trials to be used to compute
% performance for each condition seperately 

% ---------------

N = size(xstate,1);
cp = find(diff(xstate(:,1))~=0);


% loop over trials 
for t = 1:size(xstate,1)

    mu      = xstate(t,:);
    Qval    = val(t,:);
    
    % compute choice probabilities (i.e., convert Qs into choice
    % probabiltiies)
    p       = exp(beta*Qval) / sum(exp(beta*Qval));
    a(t)    = mkchoice(p);      % make choice according to the choice probabilities
    r(t)    = rand < mu(a(t));  % generate reward based based on action

%     [~,imax]    = max(mu); % what is the contigency of this trial?
%     corr(t)     = a(t)==imax;
%     
    % is the choice correct?
    if t < cp(1)+1
        [~,imax]    = max(mu); % what is the contigency of this trial?
        corr(t)     = a(t)==imax;

    elseif t > cp(5) && t < cp(6)+1
        [~,imax]    = max(mu); % what is the contigency of this trial?
        corr(t)     = a(t)==imax;

    elseif t > cp(1) && t < cp(2)+1
        [~,imin]    = min(mu);
        corr(t)     = a(t)==imin;

    elseif t > cp(3) && t < cp(4)+1

        [~,imin]    = min(mu);
        corr(t)     = a(t)==imin;

    elseif t > cp(6) && t < cp(7)+1

        [~,imin]    = min(mu);
        corr(t)     = a(t)==imin;

    elseif t > cp(8) && t < cp(9)+1

        [~,imin]    = min(mu);
        corr(t)     = a(t)==imin;

    else 

        [~,imax]    = max(mu); % what is the contigency of this trial?
        corr(t)     = a(t)==imax;
    end

   
end % end of trials loop

% transpose
r               = r'; 
a               = a'; 
corr            = corr';

% compute p(correct) for all trials
mPerf           = mean(corr);
sPerf           = mean(corr(tstable,:));
vPerf           = mean(corr(tvolatile,:));

% difference in performance between conditions
dp              = sPerf - vPerf;
xx              = [sPerf vPerf];

% compute mean reward for each of the conditions
sR              = mean(r(tstable,:));
vR              = mean(r(tvolatile,:));

mR              = [sR vR];

simResp.r       = r;
simResp.a       = a;
simResp.corr    = corr;

end
