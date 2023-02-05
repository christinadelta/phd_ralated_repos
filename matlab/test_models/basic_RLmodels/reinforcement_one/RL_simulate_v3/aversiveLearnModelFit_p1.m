% Model fitting for the simple RW model (Part of aversive-learning
% modelling) -- PART 1

% Model fitting using fmincon for the computation of LOG-LIKELIHOOD (LL) 

% created: 04/02/2023

% STEPS:


% -----------------------------------------

clear all  
clc

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);

%% define parameters and init variables


