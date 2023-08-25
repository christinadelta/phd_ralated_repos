function data = perclearn_trial_list(condition, probabilities, trials,condtrials, outtype, task)

%%% INFO WILL GO HERE


%% init variables

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions?
NumStoch            = 3;        % how many stochasticity conditions?
switches            = [1 5];    % 
nCues               = 2;        % option A and option B
blockTrials         = 50;

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
            ProbSeq(1,:)    = [thisProbs(1) thisProbs(2) thisProbs(1) thisProbs(2) thisProbs(1)]; % loss probabilities for 
            ProbSeq(2,:)    = [thisProbs(2) thisProbs(1) thisProbs(2) thisProbs(1) thisProbs(2)]; % loss probabilities for 
        end

        runTrials       = condtrials(i); 

        for run = 1:NumSwitch

            counter = 0; % init counter 

            for trl = 1:runTrials

                counter             = counter + 1; % update counter 
                tmp(counter,run)    = ProbSeq(1,run);
                tmp2(counter,run)   = ProbSeq(2,run);

            end

        end % end of runs loop

        % create blocks array
        if j == 1 && i == 1
            block(1:blockTrials,1) = 1;
        elseif j == 1 && i == 2
            block(1:blockTrials,2) = 2;
        elseif j == 2 && i == 1
            block(1:blockTrials,3) = 3;
        elseif j == 2 && i == 2
            block(1:blockTrials,4) = 4;
        elseif j == 3 && i == 1
            block(1:blockTrials,5) = 5;
        elseif j == 3 && i == 2
            block(1:blockTrials,6) = 6;
        end

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
                feedback(:,run) = makeoutcome(1:runTrials, ProbSeq(cue,run), rdm);

            end % end of run loop

            if cue == 1

                outcome1(:,1) = feedback(:);
            elseif cue == 2
                outcome2(:,1) = feedback(:);
            end

            clear feedback

        end % end of cues loop

        clear ProbSeq

        alloutcomes1(:,i) = outcome1;
        alloutcomes2(:,i) = outcome2;
        
        clear outcome1 outcome2

    end % end of volatility loop

    % update state and outcomes for stochasticity conditions
    x                   = x(:);
    xx                  = xx(:);
    alloutcomes1        = alloutcomes1(:);
    state(:,j)          = x;
    state2(:,j)         = xx;
    fb(:,j)             = alloutcomes1;

    clear x xx alloutcomes1

end % end of stch loop

%% update state and feedback/outcome

state           = state(:);
feedbck         = fb(:);
blocks          = block(:);

%% if outcomes not binary 

if outtype == 2
    % simulate outcomes with added variance noise
    totalTrials                 = length(state);
    o                           = state(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
end


%% make table to save as spreadsheet

data_table = table(blocks,state,feedbck);

% store table in .xlsx format
% filename = 'data_table.xlsx';
% writetable(data_table,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 

data.table = data_table;

%%


end % end of function