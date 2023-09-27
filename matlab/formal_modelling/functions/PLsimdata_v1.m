function data = PLsimdata_v1(probabilities, trials, condtrials, outtype)

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

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)
for j = 1:NumStoch

    thisProbs = probabilities(j, :);

    % loop over volatility conditions
    for i = 1:NumVol

        if i == 1 % if this is a stable environment
            
            NumSwitch       = switches(i);  % one switch in the contingency at 40 stable trials
            ProbSeq(1,:)    = thisProbs; 
            ProbSeq(2,:)    = 1 - thisProbs;

        else % if this is a volatile environment 

            NumSwitch       = switches(i); % 4 switches of the prob relationships
            
            % create the contigences of the volatile condition
            ProbSeq(1,:)    = [repmat(1-thisProbs,1,2)]; % cue-outcome probabilities
            ProbSeq(2,:)    = [repmat(thisProbs,1,2)]; 
        end

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
            x       = tmp; % contigencies for stable 1 option A
           
        else
            tmp     = tmp(:);
            x_vol   = nonzeros(tmp); % contigencies for stable 1 option A
            
        end

        clear tmp tmp2

        %% generate outcome and cue sequences 

        if i == 1

            for run = 1:NumSwitch
    
                cueTrials = runTrials(run);
    
                [tone_seq(:,run), outcome_seq(:,run)] = PLmakefb(1:cueTrials, ProbSeq(1,run));
    
            end % end of run loo
            
            % 
            cues        = tone_seq(:); % all runs in one column
            outcomes    = outcome_seq(:); % all runs in one column
            
        else

            for run = 1:NumSwitch
    
                cueTrials = runTrials(run);
    
                [tone_seq{1,run}, outcome_seq{1,run}] = PLmakefb(1:cueTrials, ProbSeq(1,run));
    
            end % end of run loo
    
            % 
            cues        = [tone_seq{1,1};tone_seq{1,2};tone_seq{1,3};tone_seq{1,4}]; % all runs in one column
            outcomes    = [outcome_seq{1,1};outcome_seq{1,2};outcome_seq{1,3};outcome_seq{1,4}]; % all runs in one column

        end

        % add cues and outcomes in volatility matrix
        allcues(:,i)        = cues;
        alloutcomes(:,i)    = outcomes;

        clear ProbSeq cues outcomes tone_seq outcome_seq 

    end % end of volatility loop

    % update state, cues and outcomes for stochasticity conditions
    xfinal              = [x; x_vol];
    cues_final          = allcues(:);
    outcomes_final      = alloutcomes(:);
    
    state(:,j)          = xfinal;
    cues_stc(:,j)       = cues_final;
    outcomes_stc(:,j)   = outcomes_final;
   

    clear xfinal x x_vol cues_final outcomes_final alloutcomes allcues 

end % end of stch loop


%% update state, cues and outcomes

state                   = state(:);
cues_array              = cues_stc(:);
outcomes_array          = outcomes_stc(:); % outcome column for high probability 
outcomesR_array         = 1 - outcomes_array; % outcome column for low probability 


%% if outcomes not binary 

if outtype == 2
    % simulate outcomes with added variance noise
    totalTrials                 = length(state);
    o                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
end


%% index volatility and stochasticity 
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


%% index stochasticity trials
small                       = zeros(trials*NumStoch,1);
medium                      = zeros(trials*NumStoch,1);
large                        = zeros(trials*NumStoch,1);

small(1:140)                = 1;
medium(141:280)             = 1;
large(281:420)              = 1;

small                        = small == 1;
medium                       = medium == 1;
large                       = large == 1;
s                           = [small medium large]; % column 1 = low stochasticity, column 2 = high stochasticity
t                           = [tStoch{1,1}; tStoch{1,2}; tStoch{1,3}]; % column 1 is stable, column 2 is volatile
tt                          = tStoch{1,1};

%% store all in the data structure

data = struct('o',o, 's',s, 't',t, 'tt', tt, 'state', state, 'nCues',nCues, 'outcomes', outcomes_array, 'cues',cues_array, 'outcomesR', outcomesR_array);


end % end of function