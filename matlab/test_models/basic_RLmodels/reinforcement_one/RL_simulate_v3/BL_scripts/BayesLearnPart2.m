% Bayesian Learning - Part 2 -- Bayes with & without volatility 

% Created: 12/2/2023

% 2st part of a simple bayesian learning model:
% The script goes through Bayesian modeling without volatility 

% Now let's start building a simple Bayesian model for the aversive
% learning task.

% In this simple version I will first start with the stable condition (no
% switch), then the stable condition with switch and then move on to the
% volatile condition

clear all  
clc

%% initialise variables

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);

% initialise variables 
subjects            = 1;
params              = [.25 4];                % alpha and beta values 
condition           = 1;                      % only stable condition for now (if 2 = stable and volatile)
task                = 2;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];              % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20];     % probabilities of the stable + volatile condition
end

labels              = {'alpha', 'beta'};      % for ploting
condstring          = {'stable', 'volatile'}; % for ploting 
trials              = 100;                    % per volatility condition

H                   = 1/25;                   % still not 100% sure what that does

%% simulate dataset(s)

% for now simulate one dataset only
data                = avlearn_simulate_v1(condition, probs, trials, outpath, task); % data is a structure containaining the simulation output

% what data do we need?
q                   = data.feedbackprob;    % true underlying probabilities
y                   = data.feedback(:,1);   % feedback -- 1 = vertical gabor, 0 = horizontal gabor

%% set up the state space

q_space             = [0.01:0.01:0.99];                         % possible values for q - don't allow 0 or 1, as model may explode!
prior_p_q           = NaN(length(y),length(q_space));           % initiate matrix for prior probability of each value of q: size = trials x possible values
post_p_q            = NaN(length(y),length(q_space));           % initiate matrix for posterior probability of each value of q: size = trials x possible values
prior_p_q(1,:)      = ones(1,length(q_space))./length(q_space); % on trial 1, prior probability for each value of q is uniform ie 1/number of possible values of q

%% update probabilities trial-by-trial (first version of the model)

% loop through trials
for i = 1:length(y)
    
    if y(i) == 1  % if target was vertical gabor
          p_q_given_yi(i,:)     = q_space;
    else
          p_q_given_yi(i,:)     = 1-q_space;
    end
    
    post_p_q(i,:)               = p_q_given_yi(i,:) .* prior_p_q(i,:);
    post_p_q(i,:)               = post_p_q(i,:)./sum(post_p_q(i,:),2); % normalise the posterior
    
    % set up the prior for next trial - there is a chance of (1-H) that 
    % q is the same on the next trial as on the current trial, 
    % and a chance of H that it now has a random value between 0 and 1
    prior_p_q(i+1,:)            = post_p_q(i,:)*(1-H) + H*ones(1,length(q_space))./length(q_space);
    
end

% Work out the model's best guess of pL on each trial by finding the expected value of q
est_q                           = sum(prior_p_q*q_space',2);

%% plot the model's beliefs about q (p(orange rewarded)) on the first few trials

figure; hold on;

first_trials = 5;

for c = 1:first_trials; subplot(1,first_trials,c); hold on;

    t       = 1:5; % representative trials
    plot([est_q(t(c)) est_q(t(c))],[0 1],'r-.','MarkerSize',15);
    plot(q_space,prior_p_q(t(c),:));

    if y(t(c)) == 1
        title(['prior on trial ' int2str(t(c)) ': outcome = vertical']);
    else
        title(['prior on trial ' int2str(t(c)) ': outcome = horizontal']);
    end

    xlabel('value of q'); ylabel('p(value of q)'); set(gcf,'color','w'); set(gca,'XLim',[0 1],'YLim',[0 0.025]);
    legend(['best estimate of q']);
end

%% plot the true value of q (p(orange rewarded)) and the model's estimate

figure; plot(q); hold on; plot(est_q,'r');  % plot the true values of s against the model's estimate

plot(y,'k.','MarkerSize',8);                % plot the individual data points (left and right targets)
set(gca,'YLim',[-0.2 1.2]); legend('true value of q','model estimate','Data points (1=orange, 0=blue)');
xlabel('trial number','FontSize',16); ylabel('model estimate of q','FontSize',16); set(gcf,'color','w'); set(gca,'FontSize',16);

%% plot the model's beliefs as a distribution on each trial

figure; hold on; 
imagesc(1:length(y),flipud(q_space),post_p_q'); % plot full posterior over all trials

plot(q,'w'); plot(est_q,'k');plot(est_q,'w--'); % add best estimate of p(left)
plot(0.1+y*0.8,'w.','MarkerSize',5);            % add data points (vertical/horizontal rewarded on each trial)
set(gcf,'color','w'); set(gca,'XLim',[0 length(y)],'FontSize',14);

xlabel('trial number','FontSize',14); ylabel('possible values of q','FontSize',14); 
title('colour represents probability of each possible value of q','FontSize',14);

%% plot the model's beliefs about q (p(vertical rewarded)) in the first 15 trials

figure; hold on;

for c = 1:15; subplot(3,5,c); hold on;

    t = 45:60; % 30:44; % representative trials
    plot([est_q(t(c)) est_q(t(c))],[0 1],'r-.','MarkerSize',15);
    plot(q_space,prior_p_q(t(c),:));
    
    if y(t(c)) == 1
        title(['t' int2str(t(c)) ': vertical']);
    else
        title(['t' int2str(t(c)) ': horizontal']);
    end

    xlabel('value of q'); ylabel('p(value of q)'); set(gcf,'color','w'); set(gca,'XLim',[0 1],'YLim',[0 0.05]);
end


