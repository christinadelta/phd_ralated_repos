%% Data description

% simple bandit task (2-forced choices), either blue or orange machine 
% 2 levels of volatility (low and high)
% different subjects in 2 conditions (high vs low voaltility)
% subjects:
% subs 1,3,5,7: high volaitlity
% subs 2,4,6,8: low volatility 

%% play around with real data 

%load data 

% note all the data are stored in cell
% but we also have the individual subjects structures 
load('workspace.mat')


%% prepare data

clear all
clc

% this part should always be modified for different plots
subjects        = 1;
simulate        = false; 
fitData         = false; 
plotIndividual  = true; 

% specify parameters for simulating data 
% specify bounds and bins if fitting data

% set labels:
volatility      = {'high', 'low'};
label_param     = {'alpha', 'beta'};
nparams         = length(label_param);

% sort subjects per condition
highvol_s       = alldata{1,1}; % odd subejcts
lowvol_s        = alldata{1,2}; % even subjects

ntrials         = 135; % number of trials

%% prepare and plot individual data (1 subject)

% this part should always be modified for different plots
% subjects        = 1;
% simulate        = false; 
% fitData         = false; 
% plotIndividual  = true; 

% set labels:
volatility      = {'high', 'low'};
nparams         = length(volatility);

% sort subjects per condition
highvol_s       = alldata{1,1}; % odd subejcts
lowvol_s        = alldata{1,2}; % even subjects

ntrials         = 135; % number of trials

% plot individual data of one subject (1st subject)
subject         = 1;
subj            = highvol_s{1,1};
% subj            = lowvol_s{1,1};

condition       = 2-mod(subject,nparams); % if 1=high volatility, 2=low volatility

% plot data
fp              = plotInd(subject, subj, volatility{condition});


