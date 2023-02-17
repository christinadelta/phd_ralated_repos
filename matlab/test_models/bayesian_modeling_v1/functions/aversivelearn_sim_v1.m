function data = aversivelearn_sim_v1(cond, probs, trials, outpath, task)

% Date created : 7/1/2023
% modified 1: 13/1/2023
% modified 2: 15/1/2023
% modified 3: 25/1/2023
% modified 4: 2/2/2023

% this function runs through the aversive_learn_part1.m script.
% It simulates data both for modelling and for visuaisation purposes 
% Details of the model are given in the modeling function
% details about the visualisation are given in the vis function.

% -------
% Details related to the task simulation:

% version 1 simulates only stable condition trials for a 2-armed bandit
% task:
% in Version1 of the task, participants are presented with 2 gabor stimuli
% on the left and right side of the screen. One gabor shape has vertical
% lines and the other one has horizontal lines. 

% the two gabors are associated with certain probabilities of recieving an
% reward; and the probabilities are 75% and 25%. This means that
% there is always a good option and a bad option. The good option is the
% one with the large probability. 
% 
% The reasoning here, is very similar to that of the reward-based version (see
% version 0 of the simulation and model). 

% For now (modification 1), I only simulate trials for the stable
% condition. As I learn the model, I will continue adding complexity to the
% simulation process and to the model. 

% NEW CHANGE (15/1/2023):
% I have now added the volatile version too

% NEW CHANGE (25/1/2023):
% changed the way feedback is generated 

% This function takes as INPUTS:
%                              - condition (1-by-n array) if n = 1 then
%                               we simulate data only for the stable condition (no change of probabilities in the
%                               trials).If n = 2, then we also include high volatility trials where probailities switch frequently.

%                              - probs (a 1-by-2 array or 2-by-2 array) with probabilities. 
%                              - number of trials (i.e., 100) per condition

% the function uses those inputs to generate trials list/sequences and
% stores them in the output structure: 
%                                       - data

% data structure cntains information that is used only for plotting
% visualisation and for modelling. 

% info used for modelling:
% 1. feedback -- a trials-by-2 vector with columns of 1 and 0 (column 1 =
% feedback for vertical gabor and column 2 = feedback for horizontal
% gabor). When column 1 row = 1, column 2 corresponding row = 0. 
% 2. number of stimuli -- this is 1-by-1 double to designitate the number of
% stimuli used (i.e., 2).
% 3. number of trials (i.e., 100) per condition 

% info used for plotting:
% feedback probability (n-by-1 array), where n = number of trials (per
% condition), with the probabilities distribution across trials 

%% -------------------
% simulate task data





end % end of function