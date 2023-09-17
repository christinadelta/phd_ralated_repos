function data = avlearn_trial_list_v2(condition, probabilities, trials,condtrials, outtype, task)

% created: September 2023

% comments will go here....


% end of preamble 
% -----------------

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

% create the 6 blocks (stable - small stoch, volatile - smallstoch, stable
% - medium stoch, volatile -medium stoch, stable - large stoch, volatile - large stoch)
for j = 1:NumStoch

    thisProbs = probabilities(j, :);

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

        % number of consition trials?
        runTrials       = condtrials{1,i}; 

        for run = 1:NumSwitch

            counter         = 0; % init counter 
            thisruntrials   = runTrials(run);

            for trl = 1:thisruntrials

                counter             = counter + 1; % update counter 
                tmp(counter,run)    = ProbSeq(1,run);
                % tmp2(counter,run)   = ProbSeq(2,run);
            end

        end % end of runs loop

        if i == 1
            tmp     = tmp(:);
            x       = tmp; % contigencies for stable 1 option A
           
        else
            tmp     = tmp(:);
            x_vol   = nonzeros(tmp); % this is  because each run has different num of trials
            
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
    cues_final          = allcues(:);
    outcomes_final      = alloutcomes(:);

    state(:,j)          = xfinal;
    cues_stc(:,j)       = cues_final;
    outcomes_stc(:,j)   = outcomes_final;

    clear xfinal x x_vol cues_final outcomes_final alloutcomes allcues 

end % end of stc loop

%% update state, cues and outcomes

state                   = state(:);
cues_array              = cues_stc(:);
outcomes_array          = outcomes_stc(:);

% create blocks array and trial randomiser array
[blocks,randomise_trials] = makeblocks(blockTrials,condtrials);

% the array cues now refferes to the high probability cue, create the low
% probability cue as 1-array
% if highProbCue = 1 p(loss|blue)
% if highProbCue = 2 p(loss|red)

highProbCue = cues_array;
lowProbCue  = 1-highProbCue;

%% convert 0s to 2s in feedback

% make the outcomes 1s & 2s
outcome_zeros                       = find(outcomes_array(:,1)==0);
outcomes_array(outcome_zeros,1)     = 2;

% make the cues 1s & 2s
cue_zeros                   = find(highProbCue(:,1)==0);
highProbCue(cue_zeros,1)    = 2;
clear cue_zeros outcome_zeros

cue_zeros                   = find(lowProbCue(:,1)==0);
lowProbCue(cue_zeros,1)    = 2;

%% add amount of loss for blue and red options 

for ii = 1:length(state)

    if outcomes_array(ii,1) == 1
        loss_blue(ii,1) = 0.1;
        loss_red(ii,1)  = 0;
    elseif outcomes_array(ii,1) == 2
        loss_blue(ii,1) = 0;
        loss_red(ii,1)  = 0.1;
    end
end

%% if outcomes not binary 

if outtype == 2
    % simulate outcomes with added variance noise
    totalTrials                 = length(state);
    o                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
end


%% create array with stimulus names (to be used in the spreadsheet)

[stimuli_left, stimuli_right, stimIdx] = writestimuli();

%% create jitter arrays

% create jitter array for fixation 
xmin    = 1000;
xmax    = 1500;
n       = 420;
xfix    = xmin+rand(1,n)*(xmax-xmin); xfix = xfix';

clear xmin xmax n

% create jitter for cues 
xmin    = 500;
xmax    = 1000;
n       = 420;
xcues   = xmin+rand(1,n)*(xmax-xmin); xcues = xcues';

clear xmin xmax n

%% make table to save as spreadsheet

leftIdx     = stimIdx(:,1);
rightIdx    = stimIdx(:,2);
data_table  = table(randomise_trials,blocks,state,highProbCue,lowProbCue,outcomes_array,stimuli_left,...
    stimuli_right,leftIdx,rightIdx,loss_blue,loss_red, xfix, xcues);

% store table in .xlsx format
filename    = 'data_table.xlsx';
writetable(data_table,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

data.table      = data_table;
data.block      = blocks;
data.s          = state;
data.cue        = cues_array;
data.outcome    = outcomes_array;
data.cueleft    = stimuli_left;
data.cueright   = stimuli_right;
data.lossblue   = loss_blue;
data.lossred    = loss_red;



end % end of function