function data = avlearn_trial_list(condition, probabilities, trials,condtrials, outtype, task)
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

%% init variables

data                = {};       % init data structure
NumVol              = 2;        % how many volatility conditions?
NumStoch            = 2;        % how many stochasticity conditions?
switches            = [1 4];    % 
nCues               = 2;        % option A and option B
blockTrials         = 80;

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
                feedback(:,run) = makefb(1:runTrials, ProbSeq(cue,run), rdm);

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

    end % end of vol loop

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

%% add amount of loss for blue and red options 

for ii = 1:length(state)

    if feedbck(ii,1) == 1
        loss_blue(ii,1) = 0.025;
        loss_red(ii,1)  = 0;
    elseif feedbck(ii,1) == 2
        loss_blue(ii,1) = 0;
        loss_red(ii,1)  = 0.025;
    end
end

%% make table to save as spreadsheet

data_table = table(blocks,state,feedbck,loss_blue,loss_red);

% store table in .xlsx format
filename = 'data_table.xlsx';
writetable(data_table,filename, 'Sheet', 1)
% movefile('*.xlsx', outpath) % move file to output dir 






%% create arrays with stim names



end 