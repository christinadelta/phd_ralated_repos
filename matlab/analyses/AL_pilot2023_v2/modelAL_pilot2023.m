% model aversive learning data
% 

%%% comments will go here

%% prepare the behavioural data 

% get paths 


clear all
clc

%% define paths etc 

startpath   = pwd;
datadir     = fullfile(startpath,'data');
subs        = dir(fullfile(datadir, '*sub-*'));
nsubs       = length(subs);

%% get subject data
% 
[all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL] = clean_ALdata(datadir,subs,nsubs);


%%