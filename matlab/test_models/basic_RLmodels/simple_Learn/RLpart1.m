% Building a reinforcement learrning model

% this is part of my PhD project PAL - PERCEPTION, ACTION LEARNING in
% Autism

% VERSION 1 - November 2022
% @christinadelta


%% RW model 

% how can we define a model that can learn about values of stimuli and
% translate these values into choices

% in the RW model the agent needs to learn the value of each option, or
% the probabiltiy that an option will pay out

% simple RW rl model:
%                       V(s,t) = V(s,t-1) + a(r(t-1) - V(s,t-1))

% here, V(s,t) is the value of the stimulus at time t, which reflects the expectation of a reward 
% r(t-1) --> the reward is received at trial t-1 (this is the outcome)
% a is the learning rate
% (r(t-1) - V(s,t-1)) is the prediction error 
% and V(s,t-1) is the prediction (based on the previous trial)

% learning rate (a) determines how much the prediction error is weighted 

%% SoftMax

% when you participate in an experiment, on every trial you are asked to decide 
% based on 2 or n number of values (options). Thus, you can either always
% pick the stimulus with the highest value or sometimes explore by choosing
% other options 

% do people use the probability maximising strategy? choose always the
% stimulus with the largest value?
%

%% Simulation description

% simple bandit task (2-forced choices), either face or house is associated with a tone 
% 2 levels of volatility (low and high)

% outcome: house or face associated with either high or low tone
% p(outcome|tone)
% 


%% Simulate one dataset and plot

% clear workspace
clear all 
clc

% add function folder to the path
addpath(genpath('functions'))

% how many subjects two simulate?
subjects        = 1;
conditions      = 2; % high volatility and low volatility

% set parameters
params(1)       = .50;  % alpha 
params(2)       = 4;    % beta 

% set labels:
volatility      = {'high', 'low'};
label_param     = {'alpha', 'beta'};
nparams         = length(volatility);

% initialise variables 
ntrials         = 184; % number of trials (per volatility condition)

% simulate data 
figh            = nan(2,subjects); % figure handle for individual trialwise plots

% simulate data for 2 subs 1 for high vol and 1 for low volatility
% condition
for sub = 1:subjects
    for cond = 1:conditions
        responses   = nan(subjects, ntrials);
        vvals       = nan(subjects, ntrials,2);
        probs       = nan(subjects, ntrials);

        % simulate data for specified number of subjects, seperately 
        [data, output]      = simulateRL(params, cond, sub);
        alldata{cond,1}     = data;
        alloutput{cond,1}   = output;
        vvals(1,:,:)        = output.vv; 
        probs(1,:)          = output.pp(:,1);
        responses(1,:)      = 2 - data.choices; % convert choices to 1 and 0 

        % let's plot the simulated results
        figh(cond,1) = figure; box off; hold on;
        set(figh(cond,1),'position', [10 60 700 400],'paperunits','centimeters','paperposition',[0 0 6 6],'Color','w');
        plot(data.simdata.outcomeprob,'k','linewidth',2) % plot stimulus probabilities
        plot(alloutput{cond,1}.pp(:,1),'b-','linewidth',2) % plot choice probabilities
        plot(alloutput{cond,1}.vv(:,1),':','color',[.4 0.4 1],'linewidth',2) % plot values for blue
        plot(alloutput{cond,1}.vv(:,2),':','color',[1 .35 0.1],'linewidth',2)% plot values for orange
        plot(responses(1,:),'k*') % plot responses

        % add legend
        legend({'p(house|tone)','p(choose house)','value(house)','value(face)','choice'},'location','northeastoutside');
        legend boxoff
        ylabel('p(house|tone)');
        ylim([-0.1 1.1]);
        xlabel('trial');
        title(sprintf('simulated data, alpha = %0.02f, beta = %02.01f',params(1),params(2)));

    end % end of simulation condition loop
end % end of subjects loop


%% simulate more participants

% how many subjects two simulate?
subjects        = 2;
conditions      = 2; % high volatility and low volatility

% set parameters
params(1)       = .50;  % alpha 
params(2)       = 4;    % beta 

% set labels:
volatility      = {'high', 'low'};
label_param     = {'alpha', 'beta'};
nparams         = length(volatility);

% initialise variables 
ntrials         = 184; % number of trials (per volatility condition)

% init variables for group plotting 
scores          = nan(conditions,subjects);
responses       = nan(subjects, ntrials);
vvals           = nan(subjects, ntrials,2); % subjects x trials x stimuli choices
probs           = nan(subjects, ntrials);   % will plot probabilities only for one stim choice

figh            = figure('Name','Data');
set(figh,'position',[10 60 900 650],'paperunits','centimeters','Color','w');


% create simulared data for subjects and task conditions
% loop over subjects 
for sub = 1:subjects
    
    for cond = 1:conditions
        
        
        
        
        
        
        
    end % end of conditions loop
    
    
    
end % end of subjects loop




