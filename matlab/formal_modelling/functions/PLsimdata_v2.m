function data = PLsimdata_v2(probabilities, trials, condtrials, outtype)

%%% INFO WILL GO HERE


%% init variables

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions?
NumStoch            = 3;        % how many stochasticity conditions?
switches            = [1 4];    % 
nCues               = 2;        % option A and option B
blockTrials         = 70;
nOut                = 1;

if outtype == 2
    outVar          = .01;  % outcome variance (play around with the value)
end


%% generate probabilistic relationships between cues-outcomes 
% create the 6 blocks (stable - small stoch, volatile - smallstoch, stable
% - medium stoch, volatile -medium stoch, stable - large stoch, volatile - large stoch)
for j = 1:NumStoch

    thisProbs = probabilities(j, :);

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

       if i == 1
            tmp     = tmp(:);
            tmp2    = tmp2(:);
            x       = tmp; % contigencies for stable 1 option A
            x2      = tmp2;
           
        else
            tmp     = tmp(:);
            tmp2    = tmp2(:);
            x_vol   = nonzeros(tmp); % this is  because each run has different num of trials
            x2_vol  = nonzeros(tmp2);
            
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
                rdm     = 0; % generate fixed sequence
            else
                rdm     = 1; % generate random sequence
            end
       
            for run = 1:NumSwitch
    
                cueTrials               = runTrials(run);
    
                [tone(:,run), feed(:,run)]    = makeoutcome_v2(1:cueTrials, ProbSeq(1,run),j);
    
            end % end of run loo
    
            % 
            cues        = tone(:); % all runs in one column
            outcomes    = feed(:); % all runs in one column
        else

            for run = 1:NumSwitch
    
                cueTrials = runTrials(run);

                stateTrials     = sum(cueTrials)*thisProbs(1); 

                if rem(stateTrials,1) == 0
                    rdm     = 0; % generate fixed sequence
                else
                    rdm     = 1; % generate random sequence
                end
    
                [tone{1,run}, feed{1,run}] = makeoutcome_v2(1:cueTrials, ProbSeq(1,run),j);
    
            end % end of run loo
    
            % 
            cues        = [tone{1,1};tone{1,2};tone{1,3};tone{1,4}]; % all runs in one column
            outcomes    = [feed{1,1};feed{1,2};feed{1,3};feed{1,4}]; % all runs in one column
        end

        % add cues and outcomes in volatility matrix
        allcues(:,i)        = cues;
        alloutcomes(:,i)    = outcomes;

        clear ProbSeq cues outcomes tone feed

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

%% update state, cues and outcomes

state                   = state(:);
stateR                  = state2(:);
cues_array              = cues_stc(:);
outcomes_array          = outcomes_stc(:); % outcome column for high probability 
outcomesR_array         = 1 - outcomes_array; % outcome column for low probability 


%% if outcomes not binary 

if outtype == 2
    % simulate outcomes with added variance noise
    totalTrials                 = length(state);
    out                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
    outR                          = stateR(:,1) + sqrt(outVar)*randn(totalTrials,1);
end

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

%% add output in the struct

t                           = [tStoch{1,1}; tStoch{1,2}; tStoch{1,3}]; % column 1 is stable, column 2 is volatile
s                           = [small medium large]; % column 1 = low stochasticity, column 2 = medium stochasticity, column 3 high stochasticity

% add all output in the data structure
data.cues                   = nCues;
data.x                      = state;
data.xR                     = stateR;
data.stcind                 = s;
data.volind                 = t;
data.out                    = out;
data.outR                   = outR;
data.o                      = outcomes_array;
data.oR                     = outcomesR_array;
data.t                      = tStoch{1,1}; % I think only the first one could work!


end % end of function