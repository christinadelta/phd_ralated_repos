% Building reinforcement learning models 
% Part of the PAL-ASD project 

% @christinadelta 
% Version 1 -- December 2022

% simple code to simulate data for a one-armed bandit task, with LS-WS RW and Choice kernel models and plot with different
% learning rates

% in the simulations I'll be using the "one-armed bandit" example 
% DESIGN:
% Subject makes choices T between k slot machines to try and maximise
% earnings 

% on each trial t, each machine k, pays out a reward r(t). This reward is 1
% with probability m(t)^k and 0 otherwise 
% rewards for each slot machine is different and are uknown to the subject 

% can reward ptobablities be fixed over time?

% WHAT ARE THE PARAMETERS OF THE TASK?
% 1. number of trials, T
% 2. number of slot machines, K
% 3. reward probabilities of the 2 options: machine 1 and machine 2

% win-stay-lose-switch model:
% the model repeats rewarded actions and switches away from unrewarded actions
% the model chooses randomly with ε in the noisy version. Randomness is the
% free parameter 

% RW model: 
% participants learn the expected values of each option (e.g., of each slot
% machine) based on the history of previous outcomes (Q values) and use the
% values to make decisions 

% learning rule in RW model:

%                   Q(k(t)) = Q(k(t)) + alpha(r(t) + Q(k(t)))

% Q is updated in response to reward r(t) 
% alpha = learning rate (value between 0 and 1) captures the extend to
% which the prediction error (r(t) + Q(k(t))) updates the Q (expected)
% values
% Initial values for Q are important (is it 0? is it 0.5?)

% How do participants make choices? --- softamx decision rule
% model assumes that participants most frequently choose the option with
% the highest value (but sometimes explores -- chooses low-value options)
% Softmax choice rule can model that kind of behaviour. 
% Here the probability of choosing option k (delta rule):

%               p(k(t)) = exp(beta*Q(k)) / sum(exp(beta*Q))


% while in learning rule equation alpha parameter is the learning rate, in
% the decision rule, beta is the inverse temperature parameter that
% controls the level of stochasticity/randomness in the choice (ranging
% from 0 to Inf), where 0 = completely random response and Inf=
% deterministic choice -- always choosing highest value option 

% both learning rate alpha and/or inverse temperature beta -- relates to the speed of learning but also to the
% noisiness in asymptotic behaviour 

% in all the example models, we want to find:
% how do people maximise rewards, when the most rewarding choice is unknown
% 

%% define colours for plotting

% 
clear all 
clc

% these colours will be used for all plots 
global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;

% add utilities to the path
addpath(genpath('utilities'))

%% ---------------
% EXAMPLE 1

% SIMULATING BEHAVIOUR IN THE BANDIT TASK

% first define parameters of the task:
ttrials     = 100;           % total trials 
nslots      = 2;            % number of stimuli, or options
mu          = [0.2 0.8 ];   % reward probabilities for each option
reps        = 2;            % number of repetitions for simulations 

% the two (model-independent) measures that the models should capture:
% 1. p(stay) --- the probability of repeating an action (should I change my
% behaviour in response to the feedback that I got?)
% 2. p(correct) --- the probability of choosing the correct option (have I
% learned from my action/choice?)

%% model 1: random responding

% how many times do we want to repeat the model?
% the free parameter in this model is the b (bias) that captures the bias
% of the subject in choosing one options instaed of the other
for i = 1:reps
    
    b                   = 0.5;                          % bias parameter that captures the randomness in the choices
    [a, r]              = model1_sim(ttrials, mu, b);   % simulate parameter values/estimates
    simul(1).a(:,i)     = a;                            % columns are the different repetitions of the model simulation
    simul(1).r(:,i)     = r;                            % columns are the different repetitions of the model simulation
    
end % end of random responding simulation

%% Model 2: Win-stay-lose-shift

% in this model the probability of repeating an action p(stay) should strongly depend on
% the past reward 

% the free parameter here epsilon (the level of randomness)
for i = 1:reps
    
    epsilon             = 0.1;
    [a, r]              = model2_sim(ttrials, mu, epsilon); 
    simul(2).a(:,i)     = a;                            % columns are the different repetitions of the model simulation
    simul(2).r(:,i)     = r;                            % columns are the different repetitions of the model simulation
    
end % end of Win-stay-lose-shift simulation

%% Model 3: Rescorla-Wagner 

% in this model Q values are updated using the learning (delta rule)
% and decision is made using softmax. Thus, this model has 2 free
% parameters (alpha, and beta).
for i = 1:reps
    
    % before running the model, define the free parameters:
    % alpha is the learning rate (between 0 & 1) -- captures the extent to
    % which the prediction error updates Q values 
    % Beta is the inverse temperature parameter that controls the level of
    % stochasticity in the decision 
    % beta = 0 --> response is completely random
    % beta = Inf --> agent deterministically chooses the option with the
    % highest value 
    alpha               = 0.1;
    beta                = 5;
    [a, r]              = model3_sim(ttrials, mu, alpha, beta); 
    simul(3).a(:,i)     = a;                            % columns are the different repetitions of the model simulation
    simul(3).r(:,i)     = r;                            % columns are the different repetitions of the model simulation
  
end % end of RW model

%% Model 4: Choice kernel 

% the model captures the tendency to repeat previous actions 
% very similar behaviour to the RW model however, here if choice learning 
% rate is 1 then the behaviour is very similar to win-stay-lose-switch,
% rewarded behaviour is repeated and learned faster 

% however, the way that the choice kernel is updated is different from the
% way that Qvalues vector is updated in the RW model 

for i = 1:reps
    
    % the model has 2 free parameters (alpha and beta choice kernel)
    alpha_k             = 0.1;
    beta_k              = 3;
    [a, r]              = model4_sim(ttrials, mu, alpha_k, beta_k);
    simul(4).a(:,i)     = a;
    simul(4).r(:,i)     = r;

end

%% Model 5: RW + Choice kernel 

% in this model Q values are updated according to RW learning rule but the 
% choice kernel updates according to the choice kernel learning rule 

for i = 1:reps
    
    % the model has 4 free parameters 
    alpha           = 0.1;
    beta            = 5;
    alpha_k         = 0.1;
    beta_k          = 1;
    [a, r]          = model5_sim(ttrials, mu, alpha, beta, alpha_k, beta_k);
    simul(5).a(:,i) = a;
    simul(5).r(:,i) = r;
end

%% analysis of stay-switch model

% compute stay-switch behaviour for each model
for i = 1:length(simul)
    
    for j = 1:reps
        
        simul(i).wsls(:,j) = analysis_wsls(simul(i).a(:,j)', simul(i).r(:,j)');
    end
    
    wsls(:,i) = nanmean(simul(i).wsls,2);

end % end of loop

%% plot stay-switch behaviour for all models 

figure(1); clf; hold on;

l = plot([0 1], wsls);
ylim([0 1])
set(l, 'marker', '.', 'markersize', 50, 'linewidth', 3)

legend({'M1: random' 'M2: WSLS' 'M3: RW' 'M4: CK' 'M5: RW+CK'}, ...
    'location', 'southeast')
xlabel('previous reward')
ylabel('probability of staying')

set(gca, 'xtick', [0 1], 'tickdir', 'out', 'fontsize', 18, 'xlim', [-0.1 1.1])

%% plot p(correct) over function of learning rates

% in this example we will simulate different beta values for different
% learning rates and run the RW model 

% create vectors of different alpha and beta values 
alphas  = [0.02:0.02:1];
betas   = [1 2 5 10 20];

for n = 1:1000
    n;
    for i = 1:length(alphas) % loop over different alpha rates
        
        for j = 1:length(betas) % loop over different beta values
            
            [a r]               = model3_sim(ttrials, mu, alphas(i), betas(j)); % run RW model
            [~, imax]           = max(mu);
            correctall(i,j,n)   = nanmean(a == imax);
            correctearly(i,j,n) = nanmean(a(1:10) == imax);
            correctlate(i,j,n)  = nanmean(a(end-9:end) == imax);

        end 
    end  
    
end % end 

% plot behaviour of different betas and alphas
early   = nanmean(correctearly,3);
late    = nanmean(correctlate,3);

% define dimensions of figure -- this will be a rectangle 
figure(1); clf;
set(gcf, 'Position', [284 498 750 300])
ax = easy_gridOfEqualFigures([0.2 0.1], [0.08 0.14 0.05 0.03]); % make plots (in one figure) of equal sizes

% create first plot (that is, the probability of making the same hoice
% (stay), based on previous reward
axes(ax(1)); hold on;
l = plot([0 1], wsls);
ylim([0 1])
set(l, 'marker', '.', 'markersize', 50, 'linewidth', 3)
leg1 = legend({'M1: random' 'M2: WSLS' 'M3: RW' 'M4: CK' 'M5: RW+CK'}, ...
    'location', 'southeast');
xlabel('previous reward')
% ylabel('probability of staying')
ylabel('p(stay)')
title('stay behavior', 'fontweight', 'normal')
xlim([-0.1 1.1]);
ylim([0 1.04])
set(ax(1), 'xtick', [0 1])
set(leg1, 'fontsize', 12)
set(leg1, 'position', [0.19    0.2133    0.1440    0.2617])
set(ax(1), 'ytick', [0 0.5 1])

% create the second plot, probability of correct for early trials
% (different alphas and beta values)
axes(ax(2)); hold on;
l1 = plot(alphas, early); % x and y vals 
xlabel('learning rate, \alpha')
ylabel('p(correct)')
title('early trials', 'fontweight', 'normal')

for i = 1:length(betas)
    
    leg{i} = ['\beta = ' num2str(betas(i))];
    
end % end of betas

leg2 = legend(l1(end:-1:1), {leg{end:-1:1}});

set([leg1 leg2], 'fontsize', 12)
set(leg2, 'position', [0.6267    0.6453    0.1007    0.2617]);

% create the third plot, probability of correct for late trials
% (different alphas and beta values)
axes(ax(3)); hold on;
l2 = plot(alphas, late);
xlabel('learning rate, \alpha')

% ylabel('p(correct)')
title('late trials', 'fontweight', 'normal')
for i = 1:length(l1)
    f = (i-1)/(length(l1)-1);
    set([l1(i) l2(i)], 'color', AZred*f + AZblue*(1-f));
end
set([l1 l2], 'linewidth', 3)
set(ax(3), 'yticklabel', [])

set(ax(2:3), 'ylim', [0.5 1.02])
set(ax, 'fontsize', 18, 'tickdir', 'out')
addABCs(ax(1:2), [-0.06 0.09], 32)

%% Fit parameters 
