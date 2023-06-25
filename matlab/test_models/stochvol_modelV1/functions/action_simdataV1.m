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

            NumSwitch       = 1;            % no switch of the probabilistic relationships
            ProbSeq(1,:)    = thisProbs(1); % high probability of losing
            ProbSeq(2,:)    = thisProbs(2); % low probability of losing

        else % if this is a volatile environment 

            NumSwitch       = condition; % 4 switches of the prob relationships
            
            % change thisProbs so that when voaltile trials start the prob
            % are [0.12 0.88 0.12 0.88]
            ProbSeq(1,:)    = repmat([thisProbs(2) thisProbs(1)],1,NumSwitch/2); % feedback/loss probabilities for 
            ProbSeq(2,:)    = repmat(thisProbs,1,NumSwitch/2); 
            

        end

        runTrials       = condtrials(i); 
        
        % create a sequence of the probabilistic relationships (at a trial
        % basis)
        for run = 1:NumSwitch

            counter = 0; % init counter 

            for trl = 1:runTrials

                counter             = counter + 1; % update counter 
                tmp(counter,run)    = ProbSeq(1,run);
            end

        end % end of runs loop
       
        % update state
        tmp                         = tmp(:);
        x(:,i)                      = tmp;
        clear tmp

        %% generate feedback sequences (outcome sequences) 

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

            outcome(:,cue) = feedback(:);

            clear feedback

        end % end of cues loop

        clear ProbSeq

        alloutcomes{1,i} = outcome;
        clear outcome

    end % end of volatility loop 

    % update state for stochasticity conditions
    x                   = x(:);
    state(:,j)          = x;
    
    % update outcomes 
    stochoutcomes{1,j} = [alloutcomes{1,1}; alloutcomes{1,2}];

    clear x alloutcomes

end % end of stochasticity loop

% update state and feedback/outcome
state           = state(:);
feedbck         = [stochoutcomes{1,1}; stochoutcomes{1,2}]; % outcomes for option 1 and option 2

% store simulations to data struct
data.nCues      = nCues;
data.state      = state;
data.feedback   = feedbck;

% make tables to easily visualise the simulated data
%[lowStochTable highStochTable] = makeTables(stochFeedback, stateENV);

% store tables in data structure
%data.lowTable   = lowStochTable;
%data.highTable  = highStochTable;

% save tables 

% index stable and volatile trials for each stochasticity level
for j = 1:NumStoch

    tstable                 = zeros(trials,1);
    tvolatile               = zeros(trials,1);
    tstable(1:100)          = 1;
    tvolatile(101:trials)   = 1;

    tstable                 = tstable == 1;
    tvolatile               = tvolatile == 1;

    tStoch{1,j}             = [tstable tvolatile];


    clear tstable tvolatile

end

%% Now let's make outcomes linear and not binary to test the model with both?

% simulate outcomes with added variance noise
totalTrials                 = length(state);
o                           = state + sqrt(outVar)*randn(totalTrials,1); % outcomes (generated based on reward rates and stochasticity??)
t                           = [tStoch{1,1}; tStoch{1,2}]; % column 1 is stable, column 2 is volatile

% index stochasticity trials
slow                        = zeros(trials*2,1);
shigh                       = zeros(trials*2,1);
slow(1:200)                 = 1;
shigh(201:end)              = 1;

slow                        = slow == 1;
shigh                       = shigh == 1;
s                           = [slow shigh]; % column 1 = low stochasticity, column 2 = high stochasticity



% update data structure
% data.stableTrials           = stableIndex;
% data.volTrials              = volIndex;
data.outcome                = o; % linear outcomes 
data.t                      = t; % index volatility condition
data.s                      = s; % index stochasticity condtion


end % end of function 