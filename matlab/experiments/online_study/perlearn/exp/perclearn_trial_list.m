function data = perclearn_trial_list(condition, probabilities, trials,condtrials, outtype, task)

%%% INFO WILL GO HERE


%% init variables

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions?
NumStoch            = 3;        % how many stochasticity conditions?
switches            = [1 4];    % 
nCues               = 2;        % option A and option B
blockTrials         = 70;
nOut                = 2;

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
    
                [tone_seq(:,run), outcome_seq(:,run)] = makeoutcome(1:cueTrials, ProbSeq(1,run));
    
            end % end of run loo
    
            % 
            cues        = tone_seq(:); % all runs in one column
            outcomes    = outcome_seq(:); % all runs in one column
        else

            for run = 1:NumSwitch
    
                cueTrials = runTrials(run);
    
                [tone_seq{1,run}, outcome_seq{1,run}] = makeoutcome(1:cueTrials, ProbSeq(1,run));
    
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
outcomes_array          = outcomes_stc(:);

% create blocks array
blocks                  = createblocks(blockTrials);

%% if outcomes not binary 

if outtype == 2
    % simulate outcomes with added variance noise
    totalTrials                 = length(state);
    o                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
end


%% create array with stimulus names (to be used in the spreadsheet)

[stimuli_cues, stimuli_outcomes,stimuli_pitch] = writestim(outcomes_array, cues_array);

%% make table to save as spreadsheet

data_table = table(blocks,state,cues_array,outcomes_array,stimuli_cues,stimuli_outcomes,stimuli_pitch);

% store table in .xlsx format
filename = 'data_table.xlsx';
writetable(data_table,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

data.table = data_table;

end % end of function