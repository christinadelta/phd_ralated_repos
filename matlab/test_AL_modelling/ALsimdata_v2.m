function data = ALsimdata_v2(probabilities, trials,condtrials)
% function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath)

% created August 2023
% modified January 2024

% To run the stochasticity/volatility model  (Piray & Daw, 2021)
% Inputs: 
% -     condition (6 conditions; double 1x1)
% -     probabilities (0.9:0.1, 0.8:0.2, 0.7:0.3)
% -     trials (total number of trials)
% -     condtrials (number of trials per condition)
% -     outtype (double 1x1; if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes) 

% Output:
% structure with simulated data for modeling:
% stateENV  = loss rates (1 x trials)
% outcome   = outcomes (1 x trials) (binary)
% estimated reward = saved in variable called "out" 
% tables with all the simulated data (for vis)

% state, outcomes and estimated reward are simulate for both cues/options

% INFO ABOUT OUTCOMES:
% we need hidden state, outcomes (binary and with some var), indexes for vol and stc  

% -------------------------
% SIMULATE TASK DATA

%% init variables

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions [stable volatile stable]?
NumStoch            = 3;        % how many stochasticity conditions?
switches            = [1 4];    % 
nCues               = 2;        % option A and option B
blockTrials         = 70;       %
nOut                = 2;        % number of outcomes 

% define vector of stdDev values for each stc
stdvals             = [0.15 0.20 0.25];

%% generate probabilistic relationships between cues-outcomes 

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)

% loop over stochasticity conditions 
for j = 1:NumStoch

    thisProbs       = probabilities(j, :); 

    for i = 1:NumVol

        if i == 1 % if this is a stable environment
            
            NumSwitch       = switches(i);  % one switch in the contingency at 40 stable trials
            ProbArray(1,:)  = thisProbs; 
            ProbArray(2,:)  = [thisProbs(2) thisProbs(1)];

        else % if this is a volatile environment 

            NumSwitch       = switches(i); % 4 switches of the prob relationships
            
            % create the contigences of the volatile condition
            ProbArray(1,:)  = [thisProbs(2) thisProbs(1) thisProbs(2) thisProbs(1)]; % cue-outcome probabilities
            ProbArray(2,:)  = [thisProbs(1) thisProbs(2) thisProbs(1) thisProbs(2)]; 
        end

        % number of consition trials?
        runTrials       = condtrials{1,i}; 

        for run = 1:NumSwitch

            counter         = 0; % init counter 
            thisruntrials   = runTrials(run);

            for trl = 1:thisruntrials

                counter             = counter + 1; % update counter 
                tmp(counter,run)    = ProbArray(1,run);
                tmp2(counter,run)   = ProbArray(2,run);
            end

        end % end of runs loop

        % add all runs state (probability sequence) in one matrix
        if i == 1

            x       = tmp(:); % contigencies for stable 1 option A
            x2      = tmp2(:);
           
        else
            x_vol   = nonzeros(tmp(:)); % this is  because each run has different num of trials
            x2_vol  = nonzeros(tmp2(:));
            
        end

        clear tmp tmp2

        %% generate outcome and cue sequences 

        % There are two ways to generate feedback sequence (outcomes): either randomly or it
        % can be fixed. 
        % If runTrials*thisProbs(1) is  whole number, then we'll generate a fixed sequence (e.g., 88% of trials will results in monetary loss and 12% will be neutral-no loss) 
        % If runTrials*thisProbs(1) is NOT a whole number, then we'll use
        % rand to generate sequence of loss/no-loss options.
     
        if i == 1

            cueTrials       = runTrials(run);
            stateTrials     = sum(cueTrials)*thisProbs(1); 

            if rem(stateTrials,1) == 0
                rdm         = 0; % generate fixed sequence
            else
                rdm         = 1; % generate random sequence
            end
       
            for run = 1:NumSwitch
    
                cueTrials               = runTrials(run);
    
                [f(:,run), c(:,run)]    = makefb_v2(1:cueTrials, ProbArray(1,run), rdm);
    
            end % end of run loo
    
            % 
            cues        = c(:); % all runs in one column
            outcomes    = f(:); % all runs in one column
        else

            for run = 1:NumSwitch
    
                cueTrials = runTrials(run);

                stateTrials     = sum(cueTrials)*thisProbs(1); 

                if rem(stateTrials,1) == 0
                    rdm     = 0; % generate fixed sequence
                else
                    rdm     = 1; % generate random sequence
                end
    
                [f{1,run}, c{1,run}] = makefb_v2(1:cueTrials, ProbArray(1,run),rdm);
    
            end % end of run loo
    
            % 
            cues        = [c{1,1};c{1,2};c{1,3};c{1,4}]; % all runs in one column
            outcomes    = [f{1,1};f{1,2};f{1,3};f{1,4}]; % all runs in one column
        end

        % add cues and outcomes in volatility matrix
        allcues(:,i)        = cues;
        alloutcomes(:,i)    = outcomes;
        volprobs{:,i}       = ProbArray;

        clear ProbArray cues outcomes tone_seq outcome_seq f c

    end % end of volatility loop

    % update state, cues and outcomes for stochasticity conditions
    xfinal              = [x; x_vol];
    x2final             = [x2; x2_vol];
    cues_final          = allcues(:);
    outcomes_final      = alloutcomes(:);

    state(:,j)          = xfinal;
    state2(:,j)         = x2final;
    cues_stc(:,j)       = cues_final;
    outcomes_stc(:,j)   = outcomes_final;
    stcprbs{:,j}        = volprobs;

    clear xfinal x x_vol cues_final outcomes_final alloutcomes allcues 

end % end of stochasticity loop

%% update state and feedback/outcome

state                       = state(:);
feedbck                     = outcomes_stc(:);
stateR                      = state2(:);
feedbckR                    = 1- feedbck;

%% index stable/volatile envs for each stc

% index stable and volatile trials for each stochasticity level (these will
% be used to run th PFs)
for j = 1:NumStoch

    tstable                 = zeros(trials,1);
    tvolatile               = zeros(trials,1);
    tstable(1:70)           = 1;
    tvolatile(71:140)       = 1;

    tstable                 = tstable == 1;
    tvolatile               = tvolatile == 1;

    tStoch{1,j}             = [tstable tvolatile];


    clear tstable tvolatile

end

%% index stc for each level

% index stochasticity trials
% index stochasticity trials
small                       = zeros(trials*NumStoch,1);
medium                      = zeros(trials*NumStoch,1);
large                       = zeros(trials*NumStoch,1);

small(1:140)                = 1;
medium(141:280)             = 1;
large(281:420)              = 1;

small                       = small == 1;
medium                      = medium == 1;
large                       = large == 1;

%% Now let's make outcomes linear and not binary to test the model with both?

% simulate estimated reward with added variance noise
totalTrials                 = length(state);
t                           = [tStoch{1,1}; tStoch{1,2}; tStoch{1,3}]; % column 1 is stable, column 2 is volatile
s                           = [small medium large]; % column 1 = low stochasticity, column 2 = medium stochasticity, column 3 high stochasticity

%% Add Gaussian noise to outcomes to be used by the Kalman and Particle filters 

for ss = 1:3 % three stc levels

    tmp_o           = feedbck(s(:,ss),:);
    tmp_oR          = feedbckR(s(:,ss),:);
    std_o(:,ss)     = tmp_o + stdvals(ss) * randn(size(tmp_o));
    std_oR(:,ss)    = tmp_oR + stdvals(ss) * randn(size(tmp_oR));

end 

std_o   = std_o(:);
std_oR  = std_oR(:);

%% add output in the struct

% add all output in the data structure
data.cues       = nCues;
data.x          = state;
data.xR         = stateR;
data.stcind     = s;
data.volind     = t;
data.o          = feedbckR; % the good option is the low probability option (1-outcome)
data.oR         = feedbck;  % the bad option is the high probability option (outcome)
data.std_o      = std_oR;   % the good option is the low probability option (1-outcome)
data.std_oR     = std_o;    % the bad option is the high probability option (outcome)
data.t          = tStoch{1,1}; % I think only the first one could work!
data.probs      = stcprbs;


end % end of function