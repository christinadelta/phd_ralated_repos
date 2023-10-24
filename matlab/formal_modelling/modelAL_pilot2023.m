% model aversive learning data
% 

%%% comments will go here

%% prepare the behavioural data 

clear all
clc

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

% model, ploting stuff
outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'stcvol');      addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);
modelfitpath    = fullfile(pwd, 'modelfit');    addpath(modelfitpath);

addpath(fullfile(pwd,'functions'))
addpath(fullfile(pwd,'codes')) % cbm


%% get subject data
% 
[all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL] = clean_ALdata(datadir,subs,nsubs);


%% get participants choice data (for testing)

for sub = 1:nsubs 
    allchoices{1,sub} = all_ALTdata{1,sub}(:,8); % get responses 
    alloutcomes{1,sub} = all_ALTdata{1,sub}(:,8); % get outcomes (high p(loss) cues) 

    alldata{sub,1}.actions = allchoices{1,sub};
    alldata{sub,1}.outcome = alloutcomes{1,sub};

end 

%% define parameters for stcvol model 

% these parameters will be used both in the particle filter (all params except the last) and response
% model (last parameter)
% 
parameters  = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1, 'beta',2);
% model_params = [.2 .2 .1 .1 2]; 

%% test model to 1 subject data/choices 

% fit to first subject for testing 
subj1.actions = allchoices{1,1};
subj1.outcome = alloutcomes{1,1};

% run model and look at ll
modelout = run_stcvol(parameters,subj1);

% check log likelihood
ll = modelout.ll

%% fit model to each subject seperately subject using laplace approximation

% fit every subject seperately in non-hierarchical way
% we a normal prior for each parameter 
% v                   = 6.25;
% prior_stc           = struct('mean', zeros(5,1),'variance',v); % mean and viarience (mean for each aprameter)
% 
% % loop over subjects 
% for i = 1:nsubs 
% 
%     % extract subject data
%     datasub     = alldata(i)
%     
%     % specify the file address 
%     fname_subj = ['lap_stcvol_sub', num2str(i), '.mat']
% 
%     % fit model
%     cbm_lap(datasub, @run_stcvol, prior_stc , fname_subj)
% 
% end % end of subjects loo
% 
% 
% 
% fname_stc           = 'lap_stcvol.mat'; % this is where output will be saved 
% 
% % fit model 1 
% cbm_lap(alldata, @run_stcvol, prior_stc, fname_stc)
% 

