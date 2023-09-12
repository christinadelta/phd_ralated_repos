% load and clean data from gorilla

%%% christinadelta 11/09/2023

clear all
clc

%% define paths first 
startpath   = pwd;
datadir     = fullfile(startpath,'data');
subs        = dir(fullfile(datadir, '*data_exp_142429-v14_*'));
nsubs       = length(subs);

%%

% loop over subjects 
for sub = 1:nsubs 

    fprintf('loading data\n')  
    subject = subs(sub).name;
    subdir  = fullfile(datadir,subject);
    fprintf('\t reading data from subject %d\n',sub); 

    %% load AQ data
    
   % AQT     = importAQdata(subdir);

    %% load SPQ data

   % SPQT    = importSPQdata(subdir);

    %% load IUS data

  %  IUT     = importIUdata(subdir);

    %% load STAI data

   % STAIT   = importSTAIdata(subdir);

    %% load perceptual learning data

    PLT                             = importPLdata(subdir);

    % only extract data that we want from the table 
    eventidx                        = PLT.EventIndex; 
    subnum(1:length(eventidx),1)    = sub;
    trialnum                        = PLT.TrialNumber;
    blocknum                        = PLT.Spreadsheet_Blocks;
    allcues                         = PLT.Spreadsheet_Cues_array;
    alloutcomes                     = PLT.Spreadsheet_Outcomes_array;
    allanswers                      = PLT.Spreadsheet_Answer;
    state                           = PLT.Spreadsheet_State;
    screenum                        = PLT.Screen; % we need numbers 3 (prediction) and 5 (response)
    screencount                     = PLT.ScreenCounter; % just want to check if itis the same as screen...
    responses                       = PLT.Response;
    rts                             = PLT.ReactionTime;
    correct                         = PLT.Correct;

    % add all columns in one matrix
    gendata                = [subnum blocknum trialnum state allcues alloutcomes allanswers screenum responses rts correct];

    % I need to keep rows that correspond to screen number 3 and thrn rows
    % that correspond to screen number 5 
    predict_temp            = find(gendata(:,8) == 3); % first get all data for predictions (resp time and correct is what needed

    % extract specific trial data that correspond to  screen 3 ( that is rt
    % and response and correct
    predict_data            = gendata((predict_temp),:); % we will need columns = 7, 9, 10, 11

    % now extract rows that belong to screen number 5 (outcome response)
    outcome_temp            = find(gendata(:,8) == 5);
    outcome_data            = gendata((outcome_temp),:); 

    % merge all required columns in one matrix
    % general columns           = [1-7]
    % prediction columns        = [9-11]
    % outcome-response columns  = [9-11]
     % we mainly need cols 8, 10, 11, 13 --> predictions and
     % prediction_correct, outcome and outcome correct
    subdata = cat(2,outcome_data(:,1), outcome_data(:,2), outcome_data(:,3), outcome_data(:,4), outcome_data(:,5),... 
        outcome_data(:,6), outcome_data(:,7), predict_data(:,9), predict_data(:,10), predict_data(:,11),...
        outcome_data(:,9), outcome_data(:,10), outcome_data(:,11));

    % store sub data 
    allPL_data{1,sub} = subdata;

    % compute accuracy levels
    % get % correct for:
    % 1. predictions 
    % 2. outcomes
    % 3. preidctions (for each volatility) 
    % 4. outcomes per voaltility 
    % 5. predictions per stochasticity 
    % 6. outcomes per stochasticty 
    [m_allp(sub,:), m_allo(sub,:), m_volp(sub,:), m_volo(sub,:), m_stcp(sub,:), m_stco(sub,:),...
        m_stc1_p(sub,:), m_stc2_p(sub,:), m_stc3_p(sub,:), m_stc1_o(sub,:), m_stc2_o(sub,:), m_stc3_o(sub,:)] = getAccuracyPL(subdata);

    % for each subject compute mean rts 
    m_rtsPL(sub,1) = mean(subdata(:,9)); % mean rts for prediction
    m_rtsPL(sub,2) = mean(subdata(:,12)); % mean rts for outcome response

    
    clear subdata outcome_temp outcome_data predict_temp predict_data gendata m_all m_vol m_stc screencount rts  
    clear eventidx subnum trialnum blocknum allcues alloutcomes allanswers state screenum responses correct       
  
    %% load action learning data

    ALT                             = importALdata(subdir);
    
    trialnum                        = ALT.TrialNumber;
    subnum(1:length(trialnum),1)    = sub;
    blocknum                        = ALT.Spreadsheet_Blocks;
    state                           = ALT.Spreadsheet_State;
    feedback                        = ALT.Spreadsheet_Feedbck; % if 1 (blue=loss), if 2 (red=loss)
    allanswers                      = ALT.Spreadsheet_Answer;
    screenum                        = ALT.Screen; % we need screen 3
    responses                       = ALT.Response;
    rts                             = ALT.ReactionTime;
    correct                         = ALT.Correct;

    gendata                         = [subnum trialnum blocknum state feedback allanswers screenum responses rts correct];

    % remove rows 1-50 (practice trials)
    gendata(1:50,:)     = [];

    % change block = 0 to block = 1 if needed

    % many rows are not needed. For each trial, only keep the response rows
    % (that is where screennum == 3)
    tmp_data            = find(gendata(:,7) == 3);
    ALTdata             = gendata((tmp_data),:);
    all_ALTdata{1,sub}  = ALTdata;

    % compute mean rts for all subjects
    m_rtsAL(sub,1) = mean(ALTdata(:,9));

    % compute accuracy for:
    % feedback all
    % feedback - per vol phase
    % feeback - per stc level
    [m_all, m_vol, m_stc, m_stc1_sv, m_stc2_sv, m_stc3_sv] = getAccuracyALT(ALTdata);

    allsubAL_meanCorrect(sub,1) = m_all;
    allsubAL_meanvol(sub,:) = m_vol;
    allsubAL_meanstc(sub,:) = m_stc;
    allsubAL_meanstc1_sv(sub,:) = m_stc1_sv;
    allsubAL_meanstc2_sv(sub,:) = m_stc2_sv;
    allsubAL_meanstc3_sv(sub,:) = m_stc3_sv;

  clear subnum trialnum blocknum allanswers state screenum responses correct  rts feedback tmp_data gendata
  clear m_stc m_vol m_all

end % end of subjects loop

%% make barplots  with accuracy levels and rts


