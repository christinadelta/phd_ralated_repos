%  PART 1 OF THE AVERSIVE - LEARNING VERSION. 

% In part 1 I only simulate data for the stable condition 


clear all
clc
%% simulate one dataset

% set paths
addpath(fullfile(pwd,'functions'))

% initialise variables 
subjects        = 1;
params          = [.25 4];              % alpha and beta values 
condition       = 1;                    % only stable condition for now (if 2 = stable and volatile)
probs           = [.75 .25];            % probabilities of the stable condition

labels          = {'alpha', 'beta'};    % for ploting
trials          = 100;                   % per volatility condition

% simulate the dataset
data            = avlearn_simulate_m1(params, condition, probs, trials); % data is a structure containaining the simulation output




