function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath)

% created May 2023

% 

% -------------------------
% simulate task data

data                = {};   % init data structure
NumVol              = 2;    % how many volatility conditions?
NumStoch            = 2;    % how many stochasticity conditions?

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)

% loop over stochasticity conditions 
for j = 1:NumStoch

    thisProbs = probabilities(j, :); % is it small or large stochasticity block?
    
    % loop over volatility conditions
    for i = 1:NumVol
    
        if i == 1 % if this is a stable environment
            NumSwitch = 1; % no switch of the probabilistic relationships
        else % if this is a volatile environment 
            NumSwitch = 4; % 4 switches of the prob relationships
        end

        nStim       = 2; % blue and red circles
        runTrials   = condtrials(i); 

    
    end % end of volatility loop 

end % end of stochasticity loop







end % end of function 