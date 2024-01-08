function [action,r] = respModel(Qvals, beta,x)

% use softmax function to simulate responses with some decisionnoise (beta
% value)
x(:,2)          = 1-x(:,1);
% ----------
for t = 1:length(Qvals)
    
    mu          = x(t,:);
    cprob       = exp(beta*Qvals(t,:)) / sum(exp(beta*Qvals(t,:)));

    % generate action based on choice probability
    action(t)   = actionChoice(cprob);
    
    % generate reward based on choice
    r(t)        = rand < mu(action(t));

end

% transpose
action  = action';
r       = r';
action  = 2 - action;


end % end of function