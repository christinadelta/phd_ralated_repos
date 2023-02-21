% Bayesian Learning - Version 2 -- Bayes with & without volatility 

% Created: 12/2/2023 -- version 1
% Modified 21/2/2023 -- version 2

% 2nd version of a simple bayesian learning model:
% The script goes through Bayesian modeling without with and volatility 

% Now let's start building a simple Bayesian model for the aversive
% learning task.

% This version includes:
% 1. simulating datasets for volatile/stable conditions
% 2. Running simple Bayesian model where the main goal is to estimate the
% the q values (predictions of the true underlying probabilities)
% 3. adding volatility component
% 4. plotting the results 

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
condition           = 2;                      % only stable condition for now (if 2 = stable and volatile)
task                = 1;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];              % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20];     % probabilities of the stable + volatile condition
end

labels              = {'alpha', 'beta'};      % for ploting
condstring          = {'stable', 'volatile'}; % for ploting 
trials              = 100;                    % per volatility condition

H                   = 1/25;                   % used to calculate the prior (changing the value affects model predictions)

%% simulate dataset(s)

% for now simulate one dataset per condition 
for cond = 1:condition
    data{1,cond}                = aversivelearn_sim_v2(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
end

% for plotting we we need outcomes and underlying probabilities both conditions, so, concatinate the
% trials of the two conditions
feedback           = cat(1,data{1,1}.feedback, data{1,2}.feedback);
feedbackprob       = cat(1,data{1,1}.feedbackprob, data{1,2}.feedbackprob); % for ploting

% what data do we need for running the model?
q                   = feedbackprob;    % true underlying probabilities
y                   = feedback(:,1);   % feedback -- 1 = vertical gabor, 0 = horizontal gabor, is the highly rewarded

%% set up the state space

% possible values for v and h
q_candidates        = (0.01:0.01:0.99)';
H_candidates        = exp(log(0.01):(log(0.2)-log(0.01))/20:log(0.2))';

% grids
[qq,HH]             = ndgrid(q_candidates,H_candidates);

% set up the transition function
transfunc =  (reshape(repmat(eye(length(q_candidates)),1,length(H_candidates)),length(q_candidates),length(q_candidates),length(H_candidates))... % p(pL(t)| no jump occurred) * p(no jump occurred)
            .* permute(reshape(repmat(1-H_candidates,length(q_candidates),length(q_candidates)),length(H_candidates),length(q_candidates),length(q_candidates)),[2 3 1]))...
            + ((ones(length(q_candidates),length(q_candidates),length(H_candidates))./length(q_candidates))... % + p(pL(t)| jump occurred) * p(jump_occurred)
            .* permute(reshape(repmat(H_candidates,length(q_candidates),length(q_candidates)),length(H_candidates),length(q_candidates),length(q_candidates)),[2 3 1]));

prior_p_qH          = NaN(size(qq,1),size(qq,2),100); %initiate prior and post matrices
post_p_qH           = NaN(size(qq,1),size(qq,2),100);

% uniform prior to start with
prior_p_qH(:,:,1)   = ones(size(qq))./length(qq(:));


%% update probabilities trial-by-trial (first version of the model)


for i = 1:length(y)

    % 'leak' or apply transition function
    if i > 1
        unpacked_post       = permute(reshape(repmat(post_p_qH(:,:,i-1),1,length(q_candidates)),length(q_candidates),length(H_candidates),length(q_candidates)),[1 3 2]);
        unpacked_prior      = unpacked_post.*transfunc;
        prior_p_qH(:,:,i)   = squeeze(sum(unpacked_prior,1));
        prior_p_qH(:,:,i)   = prior_p_qH(:,:,i)./sum(sum(prior_p_qH(:,:,i)));
    end
        
    % p(y|q,H)
    if y(i) ==  1
        % if y(i) is a hit, p(y(i)|q) is simply q
        py_given_qH         = qq;
    else
        % if y(i) is a miss, p(y(i)|q) is 1-q
        py_given_qH         = 1-qq;
    end
    
    % Bayes' theorem: p(q,H|y_1:i)=p(y|q_i,H)p(q_1:i-1,H);
    post_p_qH(:,:,i)        = py_given_qH.*prior_p_qH(:,:,i);
    
    %normalise posterior
    post_p_qH(:,:,i)        = post_p_qH(:,:,i)./sum(sum(post_p_qH(:,:,i)));
    
    % get joint maximum of prior for q and H
    [mx,ix]                 = max(myvect(prior_p_qH(:,:,i)));
    [ixq,ixH]               = ind2sub(size(qq),ix);
    joint_mode(i,:)         = [q_candidates(ixq), H_candidates(ixH)];
    
    % get marginal expected values for q and H
    Eq(i)                   = sum(myvect(prior_p_qH(:,:,i).*qq));
    EH(i)                   = sum(myvect(prior_p_qH(:,:,i).*HH));
    
end

%% plot true q and model's estimate and H

% make the plot
h = plot_SimpleBayes_V2(q,Eq,y,H,EH);

% save the plot
% store the plot
filename = fullfile(figpath, 'SimpleBayes_singlesubPlot.fig');
saveas(h, filename)

%% plot 



