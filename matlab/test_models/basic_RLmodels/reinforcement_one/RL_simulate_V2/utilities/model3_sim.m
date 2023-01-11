function [a r, updatedQvals] = model3_sim(ttrials, mu, alpha, beta)

% The rescorla wagner model has 2 free params
% Q values are updated on every trial using delta rule (line 16)

% in general this model uses the delta learning rule (to capture how fast
% the agent is learning) and the softmax decision rule to capture the way
% that the agent makes decisions (using beta as inverse temperature that
% controls the level of stochasticity or randomness in the decsion) 

initQ                   = [0.5 0.5]; % initial Q values for the two choices
Qvals                   = initQ;

for trl = 1:ttrials
    
    p                   = exp(beta*Qvals) / sum(exp(beta*Qvals));   % compute choice probabilities using softmax 
    a(trl)              = makechoice(p);                            % make choice acording to choice probabilities 
    r(trl)              = rand < mu(a(trl));                        % generate reward based on choice 
    
    % update values
    delta               = r(trl) - Qvals(a(trl));                   % prediction error
    Qvals(a(trl))       = Qvals(a(trl)) + alpha * delta;            % update Q values (only the one associated with the choice made)
    updatedQvals(trl,:) = Qvals;
    
end % end of trials loop

end % end of function