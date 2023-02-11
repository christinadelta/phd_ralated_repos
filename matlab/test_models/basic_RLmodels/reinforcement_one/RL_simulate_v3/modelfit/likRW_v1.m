function [negll output] = likRW_v1(actions, rewards, alpha, beta)

% compute negative log-likelihood

% ----------------------
% init q values
Q = [0.5 0.5];

trials = length(actions);

for trl = 1:trials

    % compute choice probabilities using the softmax function
    p                   = exp(beta*Q) / sum(exp(beta*Q));

    % compute choice probabilities for actual choices 
    choiceProb(trl)     = p(actions(trl));

    vvals(trl,:)        = Q;
    pvals(trl,:)        = p;

    % update values
    delta               = rewards(trl) - Q(actions(trl));  % pe
    Q(actions(trl))     = Q(actions(trl)) + alpha * delta; % update Q values

end % end of trials loop

% compute negative log-likelihood
negll = -sum(log(choiceProb));
% negll                   = -sum(log(pvals(actions==1,1)))+sum(log(pvals(actions==2,2)));

% store output
output.Qvals        = vvals;
output.Pvals        = pvals;
output.actions      = actions;
output.rewards      = rewards;


end % end of function