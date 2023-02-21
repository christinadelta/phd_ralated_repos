%  VERSION 2 OF THE 2-armed bandit (aversive-learning like task) MODELlING  

% Date created: 29/12/2022 -- Version 1
% modified 1: 15/1/2023
% modified 2: 29/1/2023 (added trial-by-trial plots, and more simulated
% participants)
% modified 3: 1/2/2023 (added heatmap plots to teack performance in
% different parameter combinations)
% modified 4: 20/2/2023 -- Version 2

% In version 2 I simulate data for both conditions and plot them together (combined).
% I keep track of choice probabilities (computed using the softmax
% function) and learning rates
% The steps are:
% 1. simulate dataset(s)
% 2. model simulated dataset(s)
% 3. visualise simulation:


% Information about the task: Aversive learning task Version 2:
% A simple 2-armed bandit choice task setting is used. More details in the simulation function:
% avlearn_simulate_v2.m

clear all  
clc

%%  initialise variables

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);



