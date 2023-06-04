function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype)
% function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath)

% created May 2023

% To run the stochasticity/volatility model  (Piray & Daw, 2021)
% Inputs: 
% -     condition (4 conditions; double 1x1)
% -     probabilities (0.88:0.12, 0.6:0.4)
% -     trials (total number of trials
% -     condtrials (number of trials per condition)
% -     outpath (output directory)
% -     outtype (double 1x1; if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes) 

% Output:
% structure with simulated data:
% stateENV  = loss rates (1 x trials)
% outcome   = outcomes (1 x trials) 
% tables with all the simulated data (for vis)

% INFO ABOUT OUTCOMES:
% 

% -------------------------
% SIMULATE TASK DATA

%% generate probabilistic relationships between cues-outcomes 

data                = {};   % init data structure
NumVol              = 2;    % how many volatility conditions?
NumStoch            = 2;    % how many stochasticity conditions?
nCues               = 2;    % blue and red circles

if outtype == 2
    outVar          = .01;  % outcome variance (play around with the value)
end

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)

% loop over stochasticity conditions 
for j = 1:NumStoch

    thisProbs       = probabilities(j, :); % is it small or large stochasticity block?
    
    % loop over volatility conditions
    for i = 1:NumVol
    
        if i == 1 % if this is a stable environment
            NumSwitch   = 1; % no switch of the probabilistic relationships
        else % if this is a volatile environment 
            NumSwitch   = condition; % 4 switches of the prob relationships
        end

        runTrials       = condtrials(i); 
        
        % create a sequence of the probabilistic relationships (at a trial
        % basis)
        if i == 1 % if stable

            stateENV{j,i}(1:runTrials,1)         = thisProbs(1); % sequence of probabilities for that block
            ProbSeq(1,:)                        = thisProbs(1); % high probability of losing
            ProbSeq(2,:)                        = thisProbs(2); % low probability of losing

        else
        
            % create a sequence with the probabilities
            ProbSeq(1,:)                        = repmat(thisProbs,1,NumSwitch/2); 
            ProbSeq(2,:)                        = repmat([thisProbs(2) thisProbs(1)],1,NumSwitch/2); % feedback/loss probabilities for 

            for k = 1:NumSwitch
                stateENV{j,i}(1:runTrials,k)     = ProbSeq(1,k);
            end

            stateENV{j,i}                        = stateENV{j,i}(:); % will see if this will be used 

        end % end of volatility if condition

        %% generate feedback sequences (outcome sequences) 

        if outtype == 1 

            % There are two ways to generate feedback sequence (outcomes): either randomly or it
            % can be fixed. 
            % If runTrials*thisProbs(1) is  whole number, then we'll generate a fixed sequence (e.g., 88% of trials will results in monetary loss and 12% will be neutral-no loss) 
            % If runTrials*thisProbs(1) is NOT a whole number, then we'll use
            % rand to generate sequence of loss/no-loss options.
            stateTrials        = runTrials*thisProbs(1); 
    
            if rem(stateTrials,1) == 0
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

        else % if outtype == 2


        end % end of outcome type statememnt 
        
    end % end of volatility loop 

    % add feedback for blue and red circles for low (stable-volatile) and high
    % (stable-volatile) stochasticity
    stochFeedback{1,j} = outcome;

    clear outcome

end % end of stochasticity loop

% store simulations to data struct
data.nCues      = nCues;
data.state      = stateENV;
data.outcome    = stochFeedback;

% make tables to easily visualise the simulated data
[lowStochTable highStochTable] = makeTables(stochFeedback, stateENV);

% store tables in data structure
data.lowTable   = lowStochTable;
data.highTable  = highStochTable;

% save tables 

% index stable and volatile trials for each stochasticity level
for j = 1:NumStoch

    tstable                 = zeros(trials,1);
    tvolatile               = zeros(trials,1);
    tstable(1:100)          = 1;
    tvolatile(101:trials)   = 1;

    tstable                 = tstable == 1;
    tvolatile               = tvolatile == 1;

    % store in cell
    stableIndex{1,j}        = tstable;
    volIndex{1,j}           = tvolatile;

    clear tstable tvolatile

end

data.stableTrials           = stableIndex;
data.volTrials              = volIndex;

end % end of function 