function [a r] = model4_sim(ttrials, mu, alpha_k, beta_k)

ckernel = [0 0]; % this is like a Qvalue vector

for trl = 1:ttrials
    
    p               = exp(beta_k*ckernel) / sum(exp(beta_k*ckernel));   % compute choice probabilities using softmax 
    a(trl)          = makechoice(p);                                    % make choice based on the probabilities 
    r(trl)          = rand < mu(a(trl));                                % generate reward based on choice
    ck              = (1 - alpha_k) * ckernel;                          % % update choice kernel 
    ckernel(a(trl)) = ckernel(a(trl)) + alpha_k * 1; 
    
    
end % end of trials loop

end % end of function