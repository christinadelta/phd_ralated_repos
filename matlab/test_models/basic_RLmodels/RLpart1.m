% Reinforcement learning - tutorial

% tutorial on building learning model


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

%% plot the data

clear all
clc

subjects = 2;
simulate = false;









