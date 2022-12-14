function [a r] = model2_sim(ttrials, mu, epsilon)

% simulate data for model 2 
% the agent makes the same choice as long as reward (in the previous trial is 1
% if reward is 0 it swtches actions/choices 

rlast = nan; % initial previous reward
alast = nan; % initial previous action/choice

% loop over trials
for trl = 1:ttrials
    
    % compute choice probabilities
    if isnan(rlast) % if this is the first trial

        p               = [0.5 0.5];
        
    else % if this is not the first trial
        
        % choice depends on last reward 
        if rlast == 1 % if previous reward was 1 -- it's win-stay
            
            p           = epsilon / 2 * [1 1];
            p(alast)    = 1 - epsilon/2;
            
        else % if previous reward was 0 -- this is lose-switch (the probability is 1-epsilon)
            
            p           = (1 - epsilon/2) * [1 1];
            p(alast)    = epsilon/2;
            
        end
    end
    
    % make choice according to choice probabilities 
    a(trl)              = makechoice(p);
    
    % generate reward based on choice 
    r(trl)              = rand < mu(a(trl));
    
    alast               = a(trl); % update alast
    rlast               = r(trl); % update rlast
    
end % end of trials loop 


end