% AL simulations with vKF 
% estmation of volatilty only

% created January 2024 @christinadelta

%% 
clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

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
nCues           = 2;
beta            = 1;

data            = ALsimdata_v3(probabilities, trials,condtrials);

%% run vKF 

% define parameter values
lambda                  = 0.1;
v0                      = 0.1;
omega                   = 0.1;
outcomes                = data.o;

% run vKF
[predictions, signals]  = vkf_bin(outcomes,lambda,v0,omega);

%% plot learning rates 

% extract the info needed to split learnig rates into conditions (stable,
% volatile)
lr          = signals.learning_rate;
ss          = data.stcind;
vv          = data.t;

% plot learning rates 
% split learning rates in stc and vol conditions
for j = 1:3

    stc_lr = lr(ss(:,j),1); % extract only lrs for p(loss|option A)

    for k = 1:2

        cond_lrs{j,k} = stc_lr(vv(:,k),:);

    end 
end % end of stc 

% plot 
g = plotLRs_vKF(cond_lrs)






