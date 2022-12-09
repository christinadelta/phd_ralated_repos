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

%                   Q(t+1)^k = Q(t)^k + alpha(r(t) + Q(t)^k)

% Q is updated in response to reward r(t) 
% alpha = learning rate (value between 0 and 1) captures the extend to
% which the prediction error (r(t) + Q(t)^k) updates the Q (expected)
% values
% Initial values for Q are important (is it 0? is it 0.5?)

% How do participants make choices? --- softamx decision rule
% model assumes that participants most frequently choose the option with
% the highest value (but sometimes explores -- chooses low-value options)
% Softmax choice rule can model that kind of behaviour. 
% Here the probability of choosing option k (delta rule):

%               p(t)^k = exp(beta*Q(t)^k) / sum(options K) * exp(beta*Q(t)^i)


% while in learning rule equation alpha parameter is the learning rate, in
% the decision rule, beta is the inverse temperature parameter that
% controls the level of stochasticity/randomness in the choice (ranging
% from 0 to Inf), where 0 = completely random response and Inf=
% deterministic choice -- always choosing highest value option 

% both learning rate alpha and/or inverse temperature beta -- relates to the speed of learning but also to the
% noisiness in asymptotic behaviour 

%% define colours for plotting

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


%% ---------------
% EXAMPLE 1

% SIMULATING BEHAVIOUR IN THE BANDIT TASK

% first define parameters of the task:
ttrials     = 10;           % total trials 
nslots      = 2;            % number of stimuli, or options
m           = [0.2 0.8 ];   % reward probabilities 
reps        = 2;            % number of repetitions for simulations 

%% model 1: random responding

% how many times do we want to repeat the model?
for i = reps
    
    b                   = 0.5;                          % bias parameter that captures the randomness in the choices
    [a, r]              = model1_sim(ttrials, mu, b);   % simulate parameter values/estimates
    sim_rr(1).a(:,i)    = a;
    sim_rr(1).r(:,i)    = r;
    
end % end of random responding simulation












