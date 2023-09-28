function [blocks,randomiser_col] = createblocksV2(blockTrials,condtrials)

% using the 1st input the function creates indexes for each block: 1-6 
% blocks 1,3,5 = are stable 
% blocks 2,4,6 = are volatile 

% using the 2nd input the function creates index for randomising trials
% within the blocks. So, if randomiser = 1, then all trials in block 1
% (stable - small stc condition) are randomised 
% in volatile conditions we have 4 switches and we want to randomise the
% trials within each switch, so in volatile blocks we should create 4
% randomiser indices 

% ------------------------------
%%
% create blocks array
for b = 1:6
    block(:,b)  = ones(blockTrials,1)*b;
end % end of blocks loop

blocks          = block(:);

%%
%  create within-block trial randomiser column for gorilla spreadsheet
volcond     = [1 2 1 2 1 2]; % 1=stable, 2=volatile
switches    = [1,4];
counter     = 0; % init counter 

% create stable randomisers
for i = 1:6

    % is it stable or voalitile?
    thisVol = volcond(i);
    thisTrls = condtrials{1,thisVol};


    if thisVol == 1

        counter = counter + 1;

        tmp_cols(:,i) = ones(thisTrls,1)*counter;

        clear thisTrls

    elseif thisVol == 2

        for j = 1:4

            counter         = counter + 1;
            voltrials       = thisTrls(j);
            volcols{1,j}    = ones(voltrials,1)*counter;

        end % end of j loop

        t                   = [volcols{1,1};volcols{1,2};volcols{1,3};volcols{1,4}];
        tmp_cols(:,i)       = t;

        clear t volcols voltrials thisTrls
    end % end of if
end % end of blocks loop

% flatten tmp_columns
randomiser_col = tmp_cols(:);



end % end of function