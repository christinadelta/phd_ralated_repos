function data = generateData



% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
trials          = 140;                      % total trials

condtrials      = {70,[30,10,10,20]};
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes

data            = {};                       % init data structure 

NumVol          = 2;                        % how many volatility conditions [stable volatile stable]?
NumStoch        = 3;                        % how many stochasticity conditions?
switches        = [1 4];                    % 
nCues           = 2;                        % option A and option B
nOut            = 2;                        % number of outcomes 

%% generate probabilistic relationships between cues-outcomes 

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)

% loop over stochasticity conditions 
for j = 1:NumStoch

    thisProbs       = probabilities(j, :); 

    for i = 1:NumVol

        if i == 1 % if this is a stable environment
            
            NumSwitch       = switches(i);  % one switch in the contingency at 40 stable trials
            ProbSeq(1,:)    = thisProbs; 
            ProbSeq(2,:)    = [thisProbs(2) thisProbs(1)];

        else % if this is a volatile environment 

            NumSwitch       = switches(i); % 4 switches of the prob relationships
            
            % create the contigences of the volatile condition
            ProbSeq(1,:)    = [thisProbs(2) thisProbs(1) thisProbs(2) thisProbs(1)]; % cue-outcome probabilities
            ProbSeq(2,:)    = [thisProbs(1) thisProbs(2) thisProbs(1) thisProbs(2)]; 
        end

        % number of consition trials?
        runTrials       = condtrials{1,i}; 

        for run = 1:NumSwitch

            counter         = 0; % init counter 
            thisruntrials   = runTrials(run);

            for trl = 1:thisruntrials

                counter             = counter + 1; % update counter 
                tmp(counter,run)    = ProbSeq(1,run);
                tmp2(counter,run)   = ProbSeq(2,run);
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
    
                [f(:,run), c(:,run)]    = makefb_v2(1:cueTrials, ProbSeq(1,run), rdm);
    
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
    
                [f{1,run}, c{1,run}] = makefb_v2(1:cueTrials, ProbSeq(1,run),rdm);
    
            end % end of run loo
    
            % 
            cues        = [c{1,1};c{1,2};c{1,3};c{1,4}]; % all runs in one column
            outcomes    = [f{1,1};f{1,2};f{1,3};f{1,4}]; % all runs in one column
        end

        % add cues and outcomes in volatility matrix
        allcues(:,i)        = cues;
        alloutcomes(:,i)    = outcomes;

        clear ProbSeq cues outcomes tone_seq outcome_seq f c

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

    clear xfinal x x_vol cues_final outcomes_final alloutcomes allcues 

end % end of stc loop

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

%% now let's store suff in the 

t                           = [tStoch{1,1}; tStoch{1,2}; tStoch{1,3}]; % column 1 is stable, column 2 is volatile
s                           = [small medium large]; % column 1 = low stochasticity, column 2 = medium stochasticity, column 3 high stochasticity

% add all output in the data structure
data.cues       = nCues;
data.x          = state;
data.xR         = stateR;
data.stcind     = s;
data.volind     = t;
data.o          = feedbck;
data.oR         = feedbckR;
data.t          = tStoch{1,1}; % I think only the first one could work!


end 