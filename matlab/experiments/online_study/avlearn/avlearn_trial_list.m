function data = avlearn_trial_list(condition, probabilities, trials,condtrials, outtype, task)
% function data = action_simdataV1(condition, probabilities, trials,condtrials, outpath)

% created August 2023

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

%% init variables

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions [stable volatile stable]?
NumStoch            = 3;        % how many stochasticity conditions?
switches            = [1 4];    % 
nCues               = 2;        % option A and option B
blockTrials         = [40, 60]; %

if outtype == 2
    outVar          = .01;  % outcome variance (play around with the value)
end

%% generate probabilistic relationships between cues-outcomes 

% create the four blocks (stable - low stoch, volatile - low stoch, stable
% - high stoch, volatile -high stoch)

% loop over stochasticity conditions 
for j = 1:NumStoch

    thisProbs       = probabilities(j, :); 

    % loop over volatility conditions
    for i = 1:NumVol

        if i == 1 % if this is a stable environment
            
            NumSwitch       = switches(i);  % one switch in the contingency at 40 stable trials
            ProbSeq(1,:)    = thisProbs; 
            ProbSeq(2,:)    = 1 - thisProbs;

        else % if this is a volatile environment 

            NumSwitch       = switches(i); % 4 switches of the prob relationships
            
            % create the contigences of the volatile condition
            ProbSeq(1,:)    = [repmat(1-thisProbs,1,2)]; % loss probabilities for 
            ProbSeq(2,:)    = [repmat(thisProbs,1,2)]; % loss probabilities for 
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

        % update state(s)
        if i == 1
            tmp     = tmp(:);
            x       = tmp; % contigencies for stable 1 option A
            tmp2    = tmp2(:);
            xx      = tmp2; % contigencies for stable 2 option A

        else
            tmp     = tmp(:);
            x_vol   = nonzeros(tmp); % contigencies for stable 1 option A
            tmp2    = tmp2(:);
            xx_vol  = nonzeros(tmp2); 

        end

        clear tmp tmp2

        %% generate outcome and cue sequences 

        % There are two ways to generate feedback sequence (outcomes): either randomly or it
        % can be fixed. 
        % If runTrials*thisProbs(1) is  whole number, then we'll generate a fixed sequence (e.g., 88% of trials will results in monetary loss and 12% will be neutral-no loss) 
        % If runTrials*thisProbs(1) is NOT a whole number, then we'll use
        % rand to generate sequence of loss/no-loss options.
        stateTrials        = sum(runTrials)*thisProbs(1); 

        if rem(stateTrials,1) == 0
            rdm     = 0; % generate fixed sequence
        else
            rdm     = 1; % generate random sequence
        end
        
        % we need two feedbacks for stable (blue as highly predicitive and
        % then red as highly predicitive)
        if i == 1

            % compute feedback (independently for each cue) 
            for cue = 1:nCues
    
                for run = 1:NumSwitch
    
                    % feedback{i,cue}(:,run) = computeFeedback(1:runTrials, ProbSeq(cue,run), rdm);
                    feedback(:,run) = makefb(1:runTrials, ProbSeq(cue,run), rdm);
    
                end % end of run loop
    
                if cue == 1
    
                    outcome1(:,1)   = feedback(:);
                    outcome1_1(:,1) = 1-feedback;
                elseif cue == 2
                    outcome2(:,1)   = feedback(:);
                    outcome2_1(:,1) = 1-feedback;
                end
    
                clear feedback
    
            end % end of cues loop

            alloutcomes1(:,i)   = outcome1; % stable 1 blue feedback
            alloutcomes1_1(:,i) = outcome1_1; % stable 1 red feedback
            alloutcomes2(:,i)   = outcome2; % stable 2 blue feedback
            alloutcomes2_1(:,i) = outcome2_1; % stable 2 red feedback

        else

            % compute feedback (for each cue) 
            for run = 1:NumSwitch

                % feedback{i,cue}(:,run) = computeFeedback(1:runTrials, ProbSeq(cue,run), rdm);
                feedback{1,run} = makefb(1:runTrials(run), ProbSeq(1,run), rdm);

            end % end of run loop
            
            outcome1(:,1)   = [feedback{1,1}; feedback{1,2};feedback{1,3}; feedback{1,4}] ;
            outcome1_1(:,1) = 1-outcome1;

            clear feedback

            alloutcomes_vol   = outcome1; % volatile 1 blue feedback
            alloutcomes1_vol  = outcome1_1; % volatile 1 red feedback
    
         
        end % end of if statement

        clear ProbSeq clear outcome1 outcome2  outcome1_1  outcome2_1

        % create cells with stimuli names for each block and randomise
        

    end % end of vol loop

    % update state and outcomes for stochasticity conditions
    xfinal              = [x; x_vol; xx];
    outcomesfinal       = [alloutcomes1; alloutcomes_vol; alloutcomes2];
    state(:,j)          = xfinal;
    fb(:,j)             = outcomesfinal;

    clear x xx alloutcomes1

end % end of stch loop

%% update state and feedback/outcome

state                       = state(:);
feedbck                     = fb(:);

% create blocks array
blocks                      = makeblocks(blockTrials);

%% convert 0s to 2s in feedback

outcome_zeros               = find(feedbck(:,1)==0);
feedbck(outcome_zeros,1)    = 2;

%% if outcomes not binary 

if outtype == 2
    % simulate outcomes with added variance noise
    totalTrials                 = length(state);
    o                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
end

%% add amount of loss for blue and red options 

for ii = 1:length(state)

    if feedbck(ii,1) == 1
        loss_blue(ii,1) = 0.1;
        loss_red(ii,1)  = 0;
    elseif feedbck(ii,1) == 2
        loss_blue(ii,1) = 0;
        loss_red(ii,1)  = 0.1;
    end
end

%% create left and right cell arrays with stimuli names (to be used in the spreadsheet)

[stimuli_left, stimuli_right] = writestimuli(blockTrials);

%% create jitter arrays

% create jitter array for fixation 
xmin    = 700;
xmax    = 1000;
n       = 420;
xfix    = xmin+rand(1,n)*(xmax-xmin); xfix = xfix';

clear xmin xmax n

% create jitter for cues 
xmin    = 350;
xmax    = 650;
n       = 420;
xcues   = xmin+rand(1,n)*(xmax-xmin); xcues = xcues';

clear xmin xmax n

%% make table to save as spreadsheet

data_table = table(blocks,state,feedbck,stimuli_left,stimuli_right,loss_blue,loss_red, xfix, xcues);

% store table in .xlsx format
filename = 'data_table.xlsx';
writetable(data_table,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

data.table = data_table;



end 