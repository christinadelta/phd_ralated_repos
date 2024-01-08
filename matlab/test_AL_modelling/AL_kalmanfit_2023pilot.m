% AVERSIVE LEARNING KALMAN MODELS TO FIT 
% MODEL 1 

% christinadelta
% December 2023

% comments will go here

%%

clc 
clear all

%% dock figures first

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% get subject data
% 
[all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL, allsub_mrtsStc, allsub_mrtsVol,allsub_mrtsAL] = clean_ALdata(datadir,subs,nsubs);

%% extract participant choice data and rts

for sub = 1:nsubs 
    allchoices{1,sub}       = all_ALTdata{1,sub}(:,8); % get responses 
    alloutcomes{1,sub}      = all_ALTdata{1,sub}(:,5); % get outcomes (p(loss|blue)) 
    allrts{1,sub}           = all_ALTdata{1,sub}(:,9);

    alldata{sub,1}.actions  = allchoices{1,sub};
    alldata{sub,1}.outcome  = alloutcomes{1,sub};
    alldata{sub,1}.rts      = allrts{1,sub};

end 

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
trials          = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
stdvals         = [0.005 0.01 0.015];       % define vector of stDev values for each stc level
nCues           = 2;
beta            = 1;

data            = ALsimdata_v3(probabilities, trials,condtrials);

%% split all the required matrices for each block 

% extract all the needed arrays from data struct
stc             = data.stcind;
vol             = data.t;
x               = data.x;

% split all data in blocks 
blockdata       = splitData(stc,vol,x,alldata,stdvals);

%% test kalman model one 

% use subject 1 data to test kalman model 1
testdata = blockdata{1,1}; 

% loop over stc and vol
for i = 1:3

    for j = 1:2

        test_strcut = testdata{i,j};
        kappa       = 0.2; % free parameter of kalman model 1

        % let's test our kalman model
        NLL = kalmanModel_one(kappa,beta,test_strcut);

        block_nll(i,j) = NLL;

    end % end of vol loop

end % end of stc loop

%% now let's fit the kalman model one

% loop over subjects 
for sub = 1:nsubs

    subdat = blockdata{1,sub}; 

    for i = 1:3 % loop over stc

        for j = 1:2

            kappa                           = 1.2; % will try different values and see what works
            blockstruct                     = subdat{i,j};
            [Xfit, NLL, alpha,choice,val]   = kalmanModelOne_objective(kappa,beta,blockstruct);
            
            % store for each subject and block 
%             allXfit{sub,1}(i,j) = Xfit;
%             allNLL{sub,1}(i,j) = NLL;
%             allalpha{sub,1}{i,j} = alpha;
            allXfit1{sub,1}(i,j)    = Xfit(1);
            allXfit2{sub,1}(i,j)    = Xfit(2);
            allNLL{sub,1}(i,j)      = NLL;
            allalpha{sub,1}{i,j}    = alpha;

        end % end of volatility loop
    end % end of stochasticity loop
end % end of subjects loop

%% plot participant lr

sublr = allalpha{4,1};
[h, g, f] = plot_fitLRs(sublr)

%%




