function [a r] = model5_sim(ttrials, mu, alpha, beta, alpha_k, beta_k)

% model 5 combines Q values from  model 3 and choice kernel from model 4

initQ   = [0.5 0.5]; 
Qvals   = initQ;
ckernel = [0 0];

for trl = 1:ttrials
    
    % compute choice probabilities using softmax
    V               = beta * Qvals + beta_k * ckernel; 
    p               = exp(V) / sum(exp(V));
    
    % make choice according to choice probabilities 
    a(trl)          = makechoice(p);
    r(trl)          = rand < mu(a(trl));
    
    % update Q values 
    delta           = r(trl) - Qvals(a(trl));                   % prediction error
    Qvals(a(trl))   = Qvals(a(trl)) + alpha * delta;            % update Q values 
    
    % update choice kernel
    ck              = (1 - alpha_k) * ckernel;                         
    ckernel(a(trl)) = ckernel(a(trl)) + alpha_k * 1; 
    
    
    
end % end of trials loop

end % end of function