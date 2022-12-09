function [a r] = model1_sim(trials, mu, b)

% simulates data/parameter values for the random choice model (model 1)

for trl = 1:trials
    
    % compute choice probabilities 
    p       = [b 1-b];
    
    % choose an option according to the choice probabilities computed above
    a(trl)  = makechoice(p);
    
    % generated a reward based on choice 
    r(trl)  = rand < mu(a(trl));
 
end % end of trials loop 

end 