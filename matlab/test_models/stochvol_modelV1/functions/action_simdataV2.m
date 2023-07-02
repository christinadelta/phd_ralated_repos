function data = action_simdataV2(condition, probabilities, trials,condtrials, outpath, outtype, task)
% function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath)

% created July 2023

% To run the stochasticity/volatility model  (Piray & Daw, 2021)
% Inputs: 
% -     condition (4 conditions; double 1x1)
% -     probabilities (0.8:0.2, 0.68:0.32)
% -     trials (total number of trials
% -     condtrials (number of trials per condition)
% -     outpath (output directory)
% -     outtype (double 1x1; if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes) 
% -     task (if 1 = 80 stable trials, if 2 = 40 stable trials for each contigency)

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

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions?
NumStoch            = 2;        % how many stochasticity conditions?
switches            = [2 5];    % 
nCues               = 2;        % option A and option B

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
            
            NumSwitch       = switches(i);  % one switch in the contingency at 40 stable trials
            ProbSeq(1,:)    = thisProbs; 
            ProbSeq(2,:)    = 1 - thisProbs;

        else % if this is a volatile environment 

            NumSwitch       = switches(i); % 4 switches of the prob relationships
            
            % create the contigences of the volatile condition
            ProbSeq(1,:)    = [repmat(thisProbs,1,2),thisProbs(1)]; % loss probabilities for 
            ProbSeq(2,:)    = 1 - ProbSeq(1,:); % reverse

        end

        runTrials       = condtrials(i); 
        
        % create a sequence of the probabilistic relationships (at a trial
        % basis)
        for run = 1:NumSwitch

            counter = 0; % init counter 

            for trl = 1:runTrials

                counter             = counter + 1; % update counter 
                tmp(counter,run)    = ProbSeq(1,run);
                tmp2(counter,run)   = ProbSeq(2,run);
            end

        end % end of runs loop
        
        % update state(s)
        tmp                         = tmp(:);
        x(:,i)                      = tmp; % contigencies for option A
        tmp2                        = tmp2(:);
        xx(:,i)                     = tmp2; % contigencies for option B

        clear tmp tmp2

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
    xx                  = xx(:);
    state(:,j)          = x;
    state2(:,j)         = xx;
    
    % update outcomes 
    stochoutcomes{1,j} = [alloutcomes{1,1}; alloutcomes{1,2}];

    clear x xx alloutcomes

end % end of stch loop

% update state and feedback/outcome
state           = state(:);
state(:,2)      = state2(:);
feedbck         = [stochoutcomes{1,1}; stochoutcomes{1,2}]; % outcomes for option 1 and option 2

% store simulations in data struct
data.nCues      = nCues;
data.state      = state;
data.feedback   = feedbck;

% index stable and volatile trials for each stochasticity level (these will
% be used to run th PFs)
for j = 1:NumStoch

    tstable                 = zeros(trials,1);
    tvolatile               = zeros(trials,1);
    tstable(1:80)          = 1;
    tvolatile(81:trials)   = 1;

    tstable                 = tstable == 1;
    tvolatile               = tvolatile == 1;

    tStoch{1,j}             = [tstable tvolatile];


    clear tstable tvolatile

end

% index stochasticity trials
slow                        = zeros(trials*2,1);
shigh                       = zeros(trials*2,1);
slow(1:160)                 = 1;
shigh(161:end)              = 1;

slow                        = slow == 1;
shigh                       = shigh == 1;
s                           = [slow shigh]; % column 1 = low stochasticity, column 2 = high stochasticity

%% Now let's make outcomes linear and not binary to test the model with both?

% simulate outcomes with added variance noise
totalTrials                 = length(state);
o                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
oo                          = state(:,2) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option B(generated based on reward rates and stochasticity??)
t                           = [tStoch{1,1}; tStoch{1,2}]; % column 1 is stable, column 2 is volatile


end 