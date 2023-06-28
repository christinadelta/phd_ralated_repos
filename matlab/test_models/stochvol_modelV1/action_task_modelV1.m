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

params      = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.05);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);

% RUN PF (firts the healthy and then the lesioned) with the parameters
% defined above
[stochVolSim, vals] = runModel(data, params, config, model, condition,...
    probabilities, trials,condtrials, outpath, outtype);

%% simulate responses and compute performance
% compute responses using vals
ngroups = 2;

xstate      = x;
xx          = nan(nsim, ngroups); % I guess I wont include volatility for now 
mxx         = nan(ngroups, 2); % mean correct 
exx         = nan(ngroups, 2); % se correct 

% use softmax function to simulate responses for binary choices
% loop over (volatility) conditions 
for j = 1:ngroups

    % loop over simulations
    for i = 1:nsim

        val = vals{j}(:,i);

        % generate responses for control and clinical simulated groups using the softmax function 
        [xx(i,:)] = responseModel_v1(xstate, val,tvolatile, tstable);

    end % end of simulations loop

    mxx(j,:) = median(xx);
    exx(j,:) = se_median(xx);
end % end of volatility condition


%% plot particle filter results 
fsiz        = [0 0 .45 1];

figure; 

nr          = 2;
nc          = 2;
subplots    = 1:4;

% plot estimated reward, learning rates, estimated volatility and estimated
% stochasticity 
h = plotPF(nr,nc,subplots,stochVolSim,x);


%% plot performance 

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

        params          = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',v0,'s0',s0,'s0_lesioned',0.02);
        config          = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);
        
        % RUN PF -- the healthy model only
        [stochVolSim, vals]     = runModel(data, params, config, model, condition,...
            probabilities, trials,condtrials, outpath, outtype);

        % extract learning rates for each combination of vol and stoch
        ma_stable(i,j)      = stochVolSim.ma(1);
        ea_stable(i,j)      = stochVolSim.ea(1);
        ma_vol(i,j)         = stochVolSim.ma(2);
        ea_vol(i,j)         = stochVolSim.ea(2);
        allvals{i,j}        = vals; %(vols x stochs)

        % store estimated volatility and stochasticity for each combination
        hvols{1,i}(:,j)     = stochVolSim.m_vol;
        hevols{1,i}(:,j)    = stochVolSim.e_vol;
        hstcs{1,i}(:,j)     = stochVolSim.m_stc;
        hestcs{1,i}(:,j)    = stochVolSim.e_stc;

    end % end of stochasticities loop
end % end of volatilities loop

%% plot estimated volatility and stochasticity for the different combinations 

fsiz        = [0 0 .45 1];

figure; 

nr          = 1;
nc          = 2;
subplots    = 1:2;

for i = 1:length(allvosl)

    mvol = hvols{1,i};
    evol = hevols{1,i};
    mstc = hstcs{1,i};
    estc = hestcs{1,i};

    % plot the volatilties and stochastciities
    h = plotHealthyPF(nr,nc,subplots,mvol,evol,mstc,estc,x);




end % end of volatilities loop



%% 
ngroups         = 1;
xstate          = x;

%%% get response and store choice probabilities 
for i = 1:length(allvols)

    for j = 1:length(allstc)

        xx          = nan(nsim, 2); % I guess I wont include volatility for now 
        mxx         = nan(ngroups, 2); % mean correct 
        exx         = nan(ngroups, 2); % se correct 

        for g = 1:ngroups

            % loop over simulations
            for k = 1:nsim

                val = allvals{i,j}{g,1}(:,k);

                % generate responses for control and clinical simulated groups using the softmax function 
                [xx(k,:)] = responseModel_v1(xstate, val,tvolatile, tstable);

            end % end of simulations loop
            
            % compute mean performance and mean choice-probability (across simulations) for each group
            mxx(g,:)        = mean(xx);
            exx(g,:)        = serr(xx);
            mps(i,j,g)       = mxx(1); % mean choice probabilties (stable trials)
            mpv(i,j,g)       = mxx(2); % mean choice probabilties (volatile trials)
    
        end % end of groups loop

    end % end of stochasticities loop

end% end of volatilities loop


%%%% plot learning rates for different parameter values 

fsiz        = [0 0 .45 1];
nr          = 1;
nc          = 2;
subplots    = 1:2;



hf = plotLR(nr, nc, subplots, ma_stable, ma_vol);

% plot choice probabilities
hf = plotCP(nr, nc, subplots, mps, mpv);


%% run the healthy and stochasticity lesion model with different parameters 

% this will be used to define parameter values etc..
model       = 1;

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
        [stochVolSim, vals]     = runModel(data, params, config, model, condition,...
            probabilities, trials,condtrials, outpath, outtype);

        % extract learning rates for each combination of vol and stoch
        % healthy model 
        hma_stable(i,j)  = stochVolSim.ma(1);
        hea_stable(i,j)  = stochVolSim.ea(1);
        hma_vol(i,j)     = stochVolSim.ma(2);
        hea_vol(i,j)     = stochVolSim.ea(2);
        hvals{i,j}       = vals{1,1};

        % stc lesion model
        lma_stable(i,j)  = stochVolSim.ma(3);
        lea_stable(i,j)  = stochVolSim.ea(3);
        lma_vol(i,j)     = stochVolSim.ma(4);
        lea_vol(i,j)     = stochVolSim.ea(4);
        lvals{i,j}       = vals{1,2};

    end % end of stochasticities loop
end % end of volatilities loop

%%
ngroups         = 2;
xstate          = x;

%get response and store choice probabilities 
for i = 1:length(allvols)

    for j = 1:length(allstc)

        xx          = nan(nsim, 2); % I guess I wont include volatility for now 
        mxx         = nan(ngroups, 2); % mean correct 
        exx         = nan(ngroups, 2); % se correct 

        for g = 1:ngroups

            % loop over simulations
            for k = 1:nsim
                
                if g == 1
                    val = hvals{i,j}(:,k);
                else
                    val = lvals{i,j}(:,k);

                end
             
                % generate responses for control and clinical simulated groups using the softmax function 
                [xx(k,:)] = responseModel_v1(xstate, val, tvolatile, tstable);

            end % end of simulations loop
            
            % compute mean performance and mean choice-probability (across simulations) for each group
            mxx(g,:)            = mean(xx);
            exx(g,:)            = serr(xx);
            mps(i,j,g)          = mxx(g,1); % mean choices (stable trials)
            mpv(i,j,g)          = mxx(g,2); % mean choices (volatile trials)
    
    
        end % end of groups loop
    end % end of stochasticities loop
end % end of volatilities loop 

%% plot learning rates for different parameter values (healthy and lesioned)

fsiz        = [0 0 .45 1];
nr          = 1;
nc          = 2;
subplots    = 1:2;

hf = plotLR_v2(nr, nc, subplots, hma_stable, hma_vol, lma_stable, lma_vol);

% extract response data
h_mps = mps(:,:,1);
l_mps = mps(:,:,2);

h_mpv = mpv(:,:,1);
l_mpv = mpv(:,:,2);

hf = plotCP_v2(nr, nc, subplots, h_mps, h_mpv, l_mps, l_mpv);


%% run healthy and stochasticity models with different lambda values 

% this will be used to define parameter values etc..
model       = 1;

% define a range of parameter values:
allv_lambda      = 0.1:0.2:1; % 
alls_lambda      = 0.1:0.2:1; % 

% define parameters and config for the particle filter 
nsim        = 100;

for i = 1:length(allv_lambda)

    l_v = allv_lambda(i);

    for j = 1:length(alls_lambda)

        l_s = alls_lambda(j);

        % simulate dataset first
        data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);

        x               = data.state;
        tvolatile       = data.t(:,2);
        tstable         = data.t(:,1);

        params          = struct('nparticles',100,'x0_unc',1,'lambda_v',l_v,'lambda_s',l_s,'v0',.1,'s0',.1,'s0_lesioned',0.001);
        config          = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',nsim,'model_parameters',params);
        
        % RUN PF -- the healthy model only
        [stochVolSim, vals]     = runModel(data, params, config, model, condition,...
            probabilities, trials,condtrials, outpath, outtype);

        % extract learning rates for each combination of vol and stoch
        % healthy model 
        hma_stable(i,j)  = stochVolSim.ma(1);
        hea_stable(i,j)  = stochVolSim.ea(1);
        hma_vol(i,j)     = stochVolSim.ma(2);
        hea_vol(i,j)     = stochVolSim.ea(2);
        hvals{i,j}       = vals{1,1};

        % stc lesion model
        lma_stable(i,j)  = stochVolSim.ma(3);
        lea_stable(i,j)  = stochVolSim.ea(3);
        lma_vol(i,j)     = stochVolSim.ma(4);
        lea_vol(i,j)     = stochVolSim.ea(4);
        lvals{i,j}       = vals{1,2};


    end % end of stochasticities loop
end % end of volatilities loop

%%
ngroups         = 2;
xstate          = x;

% get response and store choice probabilities 
for i = 1:length(allv_lambda)

    for j = 1:length(alls_lambda)

        xx          = nan(nsim, 2); % I guess I wont include volatility for now 
        mxx         = nan(ngroups, 2); % mean correct 
        exx         = nan(ngroups, 2); % se correct 

        for g = 1:ngroups

            % loop over simulations
            for k = 1:nsim
                
                if g == 1
                    val = hvals{i,j}(:,k);
                else
                    val = lvals{i,j}(:,k);

                end
             
                % generate responses for control and clinical simulated groups using the softmax function 
                [xx(k,:)] = responseModel_v1(xstate, val, tvolatile, tstable);

            end % end of simulations loop
            
            % compute mean performance and mean choice-probability (across simulations) for each group
            mxx(g,:)            = mean(xx);
            exx(g,:)            = serr(xx);
            mps(i,j,g)          = mxx(g,1); % mean choices (stable trials)
            mpv(i,j,g)          = mxx(g,2); % mean choices (volatile trials)
    
    
        end % end of groups loop

    end % end of stochasticities loop

end% end of volatilities loop

%% plot learning rates for different parameter values 

fsiz        = [0 0 .45 1];
nr          = 1;
nc          = 2;
subplots    = 1:2;

% plot learning rates
hf = plotLR(nr, nc, subplots, ma_stable, ma_vol);

% extract response data
h_mps = mps(:,:,1);
l_mps = mps(:,:,2);

h_mpv = mpv(:,:,1);
l_mpv = mpv(:,:,2);

hf = plotCP_v2(nr, nc, subplots, h_mps, h_mpv, l_mps, l_mpv);


