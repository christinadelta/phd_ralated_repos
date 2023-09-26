function [a,r,choiceProb] = respModel(vals,valsR,o,oR,beta)

% the function simulates responses for action learning (binary) outcomes

% Inputs:
% expected values (vals) - ntrials x nsims array
% outcomes (o) - ntrials x nsims array
% beta - value (>0) for inverse temperature 

% ------------------
l       = size(vals,1);
nsim    = size(vals,2);

% loop over simulations 
for s = 1:nsim

    % add vals for blue and red in one array
    allvals                 = [vals(:,s) valsR(:,s)];
    tmp_o                   = o(:,s);
    tmp_oR                  = oR(:,s);
    
    % convert zeros to twos
    outcome_zeros           = find(tmp_o(:,1)==0);
    tmp_o(outcome_zeros,1)  = 2;

    clear outcome_zeros

    outcome_zeros           = find(tmp_oR(:,1)==0);
    tmp_oR(outcome_zeros,1)  = 2;
    
    % add the feedback for blue and red shapes in one array ( for this sim)
    sim_o                   =[tmp_o tmp_oR];
    
    % loop over trials
    for t = 1:l
    
        p                   = exp(beta*allvals(t,:)) / sum(exp(beta*allvals(t,:)));   % compute choice probabilities using the softmax function
        choiceProb(t,:,s)   = p;
        a(t,s)              = max(find([-eps cumsum(p)] < rand));       % make choice based on choice probability
        tmp(t,s)            = sim_o(t,a(t,s));
    end % end of trials loop

    % store feedback that "sim subject recieved" based on the choice made
    a1                  = a(:,s) == 1; % actions blue
    a2                  = a(:,s) == 2; % actions red
    reward(a1)          = sim_o(a1,1); % 
    reward(a2)          = sim_o(a2,2);
    r(:,s)              = reward';

    clear tmp_oR tmp_o sim_o outcome_zeros allvals
end % end of simulations loop



end % end of function