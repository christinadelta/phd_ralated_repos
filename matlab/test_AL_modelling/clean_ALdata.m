function [all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL, allsub_mrtsStc, allsub_mrtsVol,allsub_mrtsAL] = clean_ALdata(datadir,subs,nsubs)

% PILOT VERSION 2 OF THE AVERSIVE LEARNING TASK 
% This function is used to clean the AL behavioural data and prepare it for modelling 
% Author: @christinadelta
% October 2023

%% 

% loop over subejcts and import data
for sub = 1:nsubs
    
    fprintf('loading data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(datadir,subject);
    fprintf('\t reading data from subject %d\n',sub); 

    % load AL data
    ALT     = importALdata(subdir);

    trialnum                        = ALT.TrialNumber;
    subnum(1:length(trialnum),1)    = sub;
    blocknum                        = ALT.Spreadsheet_Blocks;
    state                           = ALT.Spreadsheet_State;
    feedback                        = ALT.Spreadsheet_HighProbCue; % if 1 (blue=loss), if 2 (red=loss)
    allanswers                      = ALT.Spreadsheet_Answers;
    screenum                        = ALT.Screen; % we need screen 3
    responses                       = ALT.Response;
    rts                             = ALT.ReactionTime;
    correct                         = ALT.Correct;

    gendata                         = [subnum trialnum blocknum state feedback allanswers screenum responses rts correct];

    % remove practice trials 
    gendata(isnan(gendata(:,7)),:)  = [];  

    % most of the rows are not needed. For each trial, only keep the response rows
    % (that is where screennum == 3). The results should be 420 trials 
    tmp_data                        = find(gendata(:,7) == 3);
    ALTdata                         = gendata((tmp_data),:);
    all_ALTdata{1,sub}              = ALTdata;

    % compute mean rts for all subjects
    m_rtsAL(sub,1)                  = mean(ALTdata(:,9));

    % compute accuracy for:
    % feedback all
    % feedback - per vol phase
    % feeback - per stc level for stable and volatile phases
    [m_all, m_vol, m_stc]           = getAccuracyAL(ALTdata);
    allsub_mAL(1,sub)               = m_all;
    allsub_mVol{1,sub}              = m_vol;
    allsub_mStc{1,sub}              = m_stc;
    
    % compute mean rts for each of the six conditions 
    [mrts_all, mrts_vol, mrts_stc]      = getRTsAL(ALTdata);
    allsub_mrtsAL(1,sub)                = mrts_all;
    allsub_mrtsVol{1,sub}               = mrts_vol;
    allsub_mrtsStc{1,sub}               = mrts_stc;


end % end of subjects loop



end % end of function