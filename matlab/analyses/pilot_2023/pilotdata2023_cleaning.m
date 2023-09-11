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
    
    AQT     = importAQdata(subdir);

    %% load SPQ data

    SPQT    = importSPQdata(subdir);

    %% load IUS data

    IUT     = importIUdata(subdir);

    %% load STAI data

    STAIT   = importSTAIdata(subdir);

    %% load perceptual learning data

    PLT     = importPLdata(subdir);
    
    %% load action learning data

    ALT = importALdata(subdir);


end % end of subjects loop
