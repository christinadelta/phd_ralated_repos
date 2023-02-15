function [cfg] = setupTrials(cfg) 

%% Shuffle order of stimulus pairs

% % Group stimuli pairs for good/bad allocation - without shuffling (can
% shuffle block order in Gorilla
stimIds  = reshape(1:cfg.nSymbols, cfg.nOpts, cfg.nSymbols/cfg.nOpts)';

%% Blocks

learnBlocks = [];

for b = 1:cfg.nBlocks
    
    curStim = stimIds(b,:);
    
  % % Stimuli for this block - stimID: good, bad
    
    % But must also swap presentation of good vs. bad stim on left vs. right, 
    % hence: fullcond(w sideswap), cond(decXout), dec, out, stimID (good, bad), & stimLocRev(1 col, where 0=no, 1=yes -i.e. display bad, good)
    
    % cols: stimLocRev, stimID_hi, stimID_lo,  stimID_l, stimID_r
    bTrials = [repmat([0, curStim, curStim], cfg.nTrialsBlk/2, 1);...
                repmat([1, curStim, fliplr(curStim)], cfg.nTrialsBlk/2, 1)]; 
    
    % Shuffle order (since outcomes  already shuffled)
    bTrials = bTrials(cfg.shufflef(1:cfg.nTrialsBlk),:);
     
 
  % % Add reward info - w equal reward probability for each choice option
    rewOut = NaN(cfg.nTrialsBlk, length(cfg.rewProbs));
    for pi = 1:length(cfg.rewProbs) % For high vs. low rew prob            
        rewOut(:, pi) = getOutcomeDist(1:cfg.nTrialsBlk, cfg.rewProbs(pi), cfg.rewAssign, cfg.rewVals); % 
    end    

                
  % % Combine & add columns with outcome position reversed
    bTrialsFull = [];    
    revCol = 1;
    for iT = 1:cfg.nTrialsBlk
        if bTrials(iT,revCol)
            bTrialsFull(end+1,:) = [b, iT, bTrials(iT,:), rewOut(iT,:), fliplr(rewOut(iT,:))];
        else
            bTrialsFull(end+1,:) = [b, iT, bTrials(iT,:), rewOut(iT,:), rewOut(iT,:)];
        end
    end
    
  % % Save trial info across blocks    
    learnBlocks = [learnBlocks;...
        array2table(bTrialsFull, 'VariableNames', cfg.trialListVars)]; % convert to a table
    clear curStim bTrialsFull bTrials rewOut iT 
end

%% update
cfg.learnBlocks = learnBlocks;
cfg.stimIds = stimIds;