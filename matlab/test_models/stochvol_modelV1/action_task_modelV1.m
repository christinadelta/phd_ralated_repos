% first version of the aversive learning task with the
% stochasticity-volatility model 

% STEPS:

% created May 2023 

clear all
clc

%% set paths and init variables

% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'volstoch');    addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

% initialise variables 
subjects        = 1;
condition       = 4;                        % stable & volatile / small & large stochasticity
% task            = 1;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.88 .10;                 % small stochasticity probabilities
    .65 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 200;                      % total trials
condtrials      = [100 25];                 % 100 per stochasticity condition in stable env and 25 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes

% define initial learning rate inverse temperature parameters (don't think
% that this will be used)
params          = [.2 .6;
    3 6]; %alpha and beta parameteres 

% models to run:
mtorun = 3; % (first test the PF using the exact same params as Piray (2021), then run healthy model and stochasticity model

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
    alldata{sub,1}  = data;


end % end of subjects loop


%% run stoch-vol model with the Piray params 

% this will be used to define parameter values etc..
model = 1;

% define parameters and config for the particle filter 
nsim        = 10;
x           = data.state;
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);

params      = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);

% RUN PF (firts the healthy and then the lesioned) with the parameters
% defined above
[stochVolSim, vals] = runModel(data, params, config, model, condition,...
    probabilities, trials,condtrials, outpath, outtype);

% extract arrays 
nconds = 2;

xstate      = x;
xx          = nan(nsim, nconds); % I guess I wont include volatility for now 
mxx         = nan(nconds, 2); % mean correct 
exx         = nan(nconds, 2); % se correct 

% use softmax function to simulate responses for binary choices
% loop over (volatility) conditions 
for j = 1:nconds

    % loop over simulations
    for i = 1:nsim

        val = vals{j}(:,i);

        % generate responses for control and clinical simulated groups using the softmax function 
        [xx(i,:),p(:,i)] = responseModel_v1(xstate, val,tvolatile, tstable);

    end % end of simulations loop

    mxx(j,:) = median(xx);
    exx(j,:) = se_median(xx);
end % end of volatility condition


% plot the results 
fsiz        = [0 0 .45 1];

figure; 

nr          = 2;
nc          = 2;
subplots    = 1:4;

% plot estimated reward, learning rates, estimated volatility and estimated
% stochasticity 
h = plotPF(nr,nc,subplots,stochVolSim,x);

% plot performance 

col = def_actions('col');
fsy = def_actions('fsy');
subplots = 1;

labels = {'Stable','Volatile'};

h = plot_bar(1,1,subplots(1),{mxx},{exx},{'control','ASD'},{'Performance','Performance'},'',col);
set(h,'ylim',[0 1]);
legend(h,labels,'fontsize',fsy,'location','north','box','off');
title(h,'Model');


%% run healthy model only - with different parameter values

% this will be used to define parameter values etc..
model       = 2;

% define a range of parameter values:
allvols     = 0.1:0.2:1.5; % 
allstc      = 0.1:0.4:3;

% define parameters and config for the particle filter 
nsim        = 100;

for i = 1:length(allvols)

    v0 = allvols(i);

    for j = 1:length(allstc)

        s0 = allstc(j);

        % simulate dataset first
        data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);

        x               = data.state;
        tvolatile       = data.t(:,2);
        tstable         = data.t(:,1);

        params          = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',v0,'s0',s0,'s0_lesioned',0.001);
        config          = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);
        
        % RUN PF -- the healthy model only
        stochVolSim     = runModel(data, params, config, model, condition,...
            probabilities, trials,condtrials, outpath, outtype);

        % extract learning rates for each combination of vol and stoch
        ma_stable(i,j)  = stochVolSim.ma(1);
        ea_stable(i,j)  = stochVolSim.ea(1);
        ma_vol(i,j)     = stochVolSim.ma(2);
        ea_vol(i,j)     = stochVolSim.ea(2);

    end % end of stochasticities loop
end % end of volatilities loop

% plot learning rates for different parameter values 

fsiz        = [0 0 .45 1];
nr          = 1;
nc          = 2;
subplots    = 1:2;

hf = plotLR(nr, nc, subplots, ma_stable, ma_vol);


%% run the healthy and stochasticity lesion model with different parameters 

% this will be used to define parameter values etc..
model       = 1;

% define a range of parameter values:
allvols     = 0.1:0.2:1.5; % 
allstc      = 0.1:0.4:3;

% define parameters and config for the particle filter 
nsim        = 10;

for i = 1:length(allvols)

    v0 = allvols(i);

    for j = 1:length(allstc)

        s0 = allstc(j);

        % simulate dataset first
        data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);

        x               = data.state;
        tvolatile       = data.t(:,2);
        tstable         = data.t(:,1);

        params          = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',v0,'s0',s0,'s0_lesioned',0.001);
        config          = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);
        
        % RUN PF -- the healthy model only
        stochVolSim     = runModel(data, params, config, model, condition,...
            probabilities, trials,condtrials, outpath, outtype);

        % extract learning rates for each combination of vol and stoch
        % healthy model 
        hma_stable(i,j)  = stochVolSim.ma(1);
        hea_stable(i,j)  = stochVolSim.ea(1);
        hma_vol(i,j)     = stochVolSim.ma(2);
        hea_vol(i,j)     = stochVolSim.ea(2);

        % stc lesion model
        lma_stable(i,j)  = stochVolSim.ma(3);
        lea_stable(i,j)  = stochVolSim.ea(3);
        lma_vol(i,j)     = stochVolSim.ma(4);
        lea_vol(i,j)     = stochVolSim.ea(4);

    end % end of stochasticities loop
end % end of volatilities loop

%%plot learning rates for different parameter values (healthy and lesioned)

fsiz        = [0 0 .45 1];
nr          = 1;
nc          = 2;
subplots    = 1:2;

hf = plotLR_v2(nr, nc, subplots, hma_stable, hma_vol, lma_stable, lma_vol);

%% run healthy and stochasticity models with different lambda values 


