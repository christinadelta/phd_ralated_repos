function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath)

% created May 2023

% 

% -------------------------
% SIMULATE TASK DATA

%% generate probabilistic relationships between cues-outcomes 

data                = {};   % init data structure
NumVol              = 2;    % how many volatility conditions?
NumStoch            = 2;    % how many stochasticity conditions?
nCues               = 2;    % blue and red circles

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)

% loop over stochasticity conditions 
for j = 1:NumStoch

    thisProbs = probabilities(j, :); % is it small or large stochasticity block?
    
    % loop over volatility conditions
    for i = 1:NumVol
    
        if i == 1 % if this is a stable environment
            NumSwitch   = 1; % no switch of the probabilistic relationships
        else % if this is a volatile environment 
            NumSwitch   = 4; % 4 switches of the prob relationships
        end

        
        runTrials       = condtrials(i); 
        
        % create a sequence of the probabilistic relationships (at a trial
        % basis)
        if i == 1 % if stable

            ProbRel{j,i}(1:runTrials,1)         = thisProbs(1); % sequence of probabilities for that block
            ProbSeq(1,:)                        = thisProbs(1); % high probability of losing
            ProbSeq(2,:)                        = thisProbs(2); % low probability of losing

        else
        
            % create a sequence with the probabilities
            ProbSeq(1,:)                        = repmat(thisProbs,1,NumSwitch/2); 
            ProbSeq(2,:)                        = repmat([thisProbs(2) thisProbs(1)],1,NumSwitch/2); % feedback/loss probabilities for 

            for k = 1:NumSwitch
                ProbRel{j,i}(1:runTrials,k)     = ProbSeq(1,k);
            end

            ProbRel{j,i}                        = ProbRel{j,i}(:); % will see if this will be used 

        end % end of volatility if condition

        %% generate feedback sequences (outcome sequences) 

        % There are two ways to generate feedback sequence (outcomes): either randomly or it
        % can be fixed. 
        % If runTrials*thisProbs(1) is  whole number, then we'll generate a fixed sequence (e.g., 88% of trials will results in monetary loss and 12% will be neutral-no loss) 
        % If runTrials*thisProbs(1) is NOT a whole number, then we'll use
        % rand to generate sequence of loss/no-loss options.
        probRelTrials        = runTrials*thisProbs(1); 

        if rem(probRelTrials,1) == 0
            rdm     = 0; % generate fixed sequence
        else
            rdm     = 1; % generate random sequence
        end

        % compute feedback (independently for each cue) 
        for cue = 1:nCues

            for run = 1:NumSwitch

                % feedback{i,cue}(:,run) = computeFeedback(1:runTrials, ProbSeq(cue,run), rdm);
                feedback(:,run) = computeFeedback(1:runTrials, ProbSeq(cue,run), rdm);

            end % end of run loop

            outcome{i,cue} = feedback(:);

            clear feedback

        end % end of cues loop
        
        clear ProbSeq
        
    end % end of volatility loop 

    % add feedback for blue and red circles for low (stable-volatile) and high
    % (stable-volatile) stochasticity
    stochFeedback{1,j} = outcome;

    clear outcome

end % end of stochasticity loop

% store simulations to data struct
data.nCues      = nCues;
data.ProbRel    = ProbRel;
data.outcome    = stochFeedback;

% make tables to easily visualise the simulated data
[lowStochTable highStochTable] = makeTables(stochFeedback, ProbRel);

% store tables in data structure
data.lowTable   = lowStochTable;
data.highTable  = highStochTable;

% save tables 




end % end of function 