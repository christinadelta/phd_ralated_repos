% Full (RBPF) model fitting version 1 

% Optimisation Algorithm: BADS https://github.com/acerbilab/bads
% basic optimisation algorithm is used for now. See bads_examples.m 

% free parameters: lambda_s, lambda_v, beta, v0 and s0
% 


%%

clc 
clear all

%% dock figures first

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% get subject data
% 
[all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL, allsub_mrtsStc, allsub_mrtsVol,allsub_mrtsAL] = clean_ALdata(datadir,subs,nsubs);

%% extract participant choice data and rts

for sub = 1:nsubs 
    allchoices{1,sub}       = all_ALTdata{1,sub}(:,8); % get responses 
    alloutcomes{1,sub}      = all_ALTdata{1,sub}(:,5); % get outcomes (p(loss|blue)) 
    allrts{1,sub}           = all_ALTdata{1,sub}(:,9);

    alldata{sub,1}.actions  = allchoices{1,sub};
    alldata{sub,1}.outcome  = alloutcomes{1,sub};
    alldata{sub,1}.rts      = allrts{1,sub};

end 

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
trials          = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
stdvals         = [0.1 0.15 0.2];       % define vector of stDev values for each stc level
nCues           = 2;
beta            = 1;

data            = ALsimdata_v3(probabilities, trials,condtrials);

%% test the objective function with one subject 

testsub     = alldata{1,1};

% exatact what we need from data struct and prepare data for model fitting 
state(:,1)  = data.x;
state(:,2)  = 1 - state(:,1);
ss          = data.stcind;
vv          = data.t;
bo(:,1)     = testsub.outcome; bo(find(bo==2)) = 0; % convert outcome 2 to 0
bo(:,2)     = 1 - bo(:,1);

% apply GNA to binary outcomes 
for i = 1:3

    this_stc    = bo(ss(:,i),:);
    oT(:,i)      = this_stc(:,1) + stdvals(i) * randn(size(this_stc(:,1)));
    oR(:,i)      = this_stc(:,2) + stdvals(i) * randn(size(this_stc(:,2)));

end % end of stc loop

% concatinate outcomes
o(:,1) = oT(:);
o(:,2) = oR(:);

testsub.o       = o;                                    % add linear outcomes to the data struct
parameters      = [0.4 0.5 1 0.8 1 100 1 1 0.001];  % parameters for modelling 
x               = [0.4 0.5 1 0.8 1];                           % Starting points of free parameters (lambda_s and lambda_v)

% parameters: 
% 1st: lambda_s, -- free (not used here, see x(1))
% 2nd: lambda_v, -- free (not used here, see x(2))
% 3rd: beta,
% 4rd: s0,
% 5th: v0,
% 6th: nparticles
% 7th: x0_unc
% 8th: model type (1=core, 2=stc lesioned, 3=vol lesioned)
% 9th: prior for lesion models

config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
o = testsub.o;

% I also need to reverse the actions. The output for Gorilla, saves as
% correct action the choice for the low probability cue ( because this is
% the good/no-loss option) we need to swap this and test the model 
a               = testsub.actions;
previousTwos    = a == 2;
a(a == 1)       = 2;
a(previousTwos) = 1;

% which model should we run?
[NLL,lr] = rbpf_core_full(x,o,a,parameters,config);

%% model 1: fit the RBPF using the BADS toolbox 

% fit lambda RBPF model for each subject seperately 
% exatact what we need from data struct and prepare data for model fitting 
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

% loop over subjects 
for sub = 1:nsubs 

    dat = alldata{sub,1};
    bo(:,1)     = dat.outcome; bo(find(bo==2)) = 0; % convert outcome 2 to 0
    bo(:,2)     = 1 - bo(:,1);

    % apply GNA to binary outcomes 
    for i = 1:3
    
        this_stc    = bo(ss(:,i),:);
        oT(:,i)     = this_stc(:,1) + stdvals(i) * randn(size(this_stc(:,1)));
        oR(:,i)     = this_stc(:,2) + stdvals(i) * randn(size(this_stc(:,2)));
    
    end % end of stc loop

    % concatinate outcomes 
    o(:,1)          = oT(:);
    o(:,2)          = oR(:);

    dat.o           = o;                                    % add linear outcomes to the data struct
    a               = dat.actions;
    previousTwos    = a == 2;
    a(a == 1)       = 2;
    a(previousTwos) = 1;

    parameters      = [0.4 0.2 1.5 0.3 0.3 100 1 1 0.001];  % parameters for modelling 
    config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
    % o               = testsub.o;

    %---                       ---%
    %--- fit simple BADS model ---%

    % We specify the bounds and plausible bounds that (hopefully) 
    % contain the solution. Plausible bounds represent your best guess at 
    % bounding the region where the solution might lie.
    
    x0  = [parameters(1) parameters(2) parameters(3) parameters(4) parameters(5)];  % Starting point
    lb  = [0.1 0.1 0.1 0.1 0.1];                                                    % Lower bounds
    ub  = [0.9 0.9 10 2 2];                                                             % Upper bounds
    plb = [0.2 0.3 0.5 0.2 0.2];                                                    % Plausible lower bounds
    pub = [0.9 0.8 5 2 2];                                                          % Plausible upper bounds


    % Run BADS, which returns the minimum X and the NLL.
    bFunc       = @(x) rbpf_core_full(x,o,a,parameters,config);
    [x,fval]    = bads(bFunc,x0,lb,ub,plb,pub);
   
    % now run the model with the minimised x values to get subject learning
    % rates
    [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(x,o,parameters,config);

end 

%% plot learning rates 

% split learning rates in stc and vol conditions
for j = 1:3

    stc_lr = lr(ss(:,j),1); % extract only lrs for p(loss|option A)

    for k = 1:2

        cond_lrs{j,k} = stc_lr(vv(:,k),:);

    end 
end % end of stc 

[h, g, f] = plotLRs(cond_lrs);


%% model 2: fit the RBPF using the BADS toolbox with Non-bound constraints

% fit lambda RBPF model for each subject seperately 
% extract what we need from data struct and prepare data for model fitting 
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

% loop over subjects 
for sub = 1:nsubs 

    dat = alldata{sub,1};
    bo(:,1)     = dat.outcome; bo(find(bo==2)) = 0; % convert outcome 2 to 0
    bo(:,2)     = 1 - bo(:,1);

    % apply GNA to binary outcomes 
    for i = 1:3
    
        this_stc    = bo(ss(:,i),:);
        oT(:,i)     = this_stc(:,1) + stdvals(i) * randn(size(this_stc(:,1)));
        oR(:,i)     = this_stc(:,2) + stdvals(i) * randn(size(this_stc(:,2)));
    
    end % end of stc loop

    % concatinate outcomes 
    o(:,1)          = oT(:);
    o(:,2)          = oR(:);

    dat.o           = o;                                    % add linear outcomes to the data struct
    a               = dat.actions;
    previousTwos    = a == 2;
    a(a == 1)       = 2;
    a(previousTwos) = 1;

    parameters      = [0.1 0.1 0.1 0.1 0.1 100 1 1 0.001];  % parameters for modelling 
    config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
    % o               = testsub.o;

    %---                       ---%
    %--- fit simple BADS model ---%

    % We specify the bounds and plausible bounds that (hopefully) 
    % contain the solution. Plausible bounds represent your best guess at 
    % bounding the region where the solution might lie.
    x0  = [parameters(1) parameters(2) parameters(3) parameters(4) parameters(5)];  % Starting point
    lb  = [0.1 0.1 0.1 0.1 0.1];                                                    % Lower bounds
    ub  = [1 1 10 2 2];                                                             % Upper bounds

    % Non-bound constraints are violated outside the unit circle
    nonbcon = @(x) sum(x.^2,2) > 1;


    % Run BADS, which returns the minimum X and the NLL.
    bFunc       = @(x) rbpf_core_full(x,o,a,parameters,config);
    [x,fval]    = bads(bFunc,x0,lb,ub,[],[],nonbcon);
   
    % now run the model with the minimised x values to get subject learning
    % rates
    [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(x,o,parameters,config);

end 

%% plot learning rates 

% split learning rates in stc and vol conditions
for j = 1:3

    stc_lr = lr(ss(:,j),1); % extract only lrs for p(loss|option A)

    for k = 1:2

        cond_lrs{j,k} = stc_lr(vv(:,k),:);

    end 
end % end of stc 

[h, g, f] = plotLRs(cond_lrs);

%%

