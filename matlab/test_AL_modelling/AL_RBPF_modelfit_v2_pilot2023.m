%%% fit RBPF model for each conditions seperately 

% Optimisation Algorithm: BADS https://github.com/acerbilab/bads
% basic optimisation algorithm is used for now. See bads_examples.m 

% free parameters: 
% Full model: lambda_s, lambda_v, beta, v0 and s0
% Lambdas + beta model: lambda_s, lambda_v, beta
% Lambdas model: lambda_s, lambda_v

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

% I also need to reverse the actions. The output for Gorilla, saves as
% correct action the choice for the low probability cue ( because this is
% the good/no-loss option) we need to swap this and test the model 
a               = testsub.actions;
previousTwos    = a == 2;
a(a == 1)       = 2;
a(previousTwos) = 1;

% loop over stc and vol phases 
for i = 1:3

    % extract this block data:
    stc_o = o(ss(:,i),:);
    stc_a = a(ss(:,i),:);

    for j = 1:2

        % extract this volatility outcomes and actions
        vol_o           = stc_o(vv(:,j),:);
        vol_a           = stc_a(vv(:,j),:);
        config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
        parameters      = [0.4 0.6 1 0.8 1.2 100 1 1 0.001];  % parameters for modelling 
        x               = [0.4 0.5 1 0.8 1];                % Starting points of free parameters (lambda_s and lambda_v)
        
        % which model should we run?
        [NLL{1,i}(:,j),lr{i,j}] = rbpf_core_full(x,vol_o,vol_a,parameters,config);

    end % end of volatility phases 
end % end of stc levels

% plot learning rates for testing 
[h, g, f] = plotLRs(lr);

%% fit full model to participnat data for each block seperately 

% exatact what we need from data struct and prepare data for model fitting 
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

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

    for i = 1:3

        % extract this stc outcomes and actions:
        stc_o = o(ss(:,i),:);
        stc_a = a(ss(:,i),:);

        for j = 1:2

            % extract this volatility outcomes and actions
            vol_o           = stc_o(vv(:,j),:);
            vol_a           = stc_a(vv(:,j),:);
            config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
            parameters      = [0.3 0.6 0.5 0.8 1.2 100 1 1 0.001];  % parameters for modelling 

            %---                       ---%
            %--- fit simple BADS model ---%
        
            % We specify the bounds and plausible bounds that (hopefully) 
            % contain the solution. Plausible bounds represent your best guess at 
            % bounding the region where the solution might lie.
            
            x0  = [parameters(1) parameters(2) parameters(3) parameters(4) parameters(5)];  % Starting point
            lb  = [0.1 0.1 0.1 0.1 0.1];                                                    % Lower bounds
            ub  = [0.9 0.9 5 2 2];                                                          % Upper bounds
            plb = [0.1 0.1 0.1 0.1 0.1];                                                    % Plausible lower bounds
            pub = [0.9 0.9 5 2 2];                                                          % Plausible upper bounds

            % Set BADS options
            options                         = bads('defaults');
            options.MaxObjectiveEvaluations = 50;
            options.MaxTime                 = 3600;

            % Run BADS, which returns the minimum X and the NLL.
            bFunc       = @(x) rbpf_core_full(x0,vol_o,vol_a,parameters,config);
            [x,fval]    = bads(bFunc,x0,lb,ub,plb,pub,[],options);
           
            % now run the model with the minimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(x,vol_o,parameters,config);

            % store output for each subejct, block
            allX{1,sub}{i,j}        = x;
            allNLL{1,sub}{1,i}(:,j) = fval;
            all_lr{1,sub}{i,j}      = lr;
            all_choice{1,sub}{i,j}  = choice;

        end % end of volatility loop
    end % end of stochasticity loop
end% end of subjects loop

clear a o bo dat state ss vv x fval lr choice val 

%% plot learning rates 

% loop over participants to plot learning rates
for sub = 1:nsubs
    sub_lr = all_lr{1,sub};
    subnum = sprintf('%d',sub);

    % plot learning rates for each sub seperately
    [h, g, f] = plot_fitLRs(sub_lr,subnum);
end 

%% fit lambdas_beta model to participnat data for each block seperately 

% exatact what we need from data struct and prepare data for model fitting 
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

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

    % fit model at the block level -- loop over stc levels
    for i = 1:3

        % extract this stc outcomes and actions:
        stc_o = o(ss(:,i),:);
        stc_a = a(ss(:,i),:);

        for j = 1:2

            % extract this volatility outcomes and actions
            vol_o           = stc_o(vv(:,j),:);
            vol_a           = stc_a(vv(:,j),:);
            config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
            parameters      = [0.3 0.5 0.5 0.3 1 100 1 1 0.001];  % parameters for modelling 
            % x               = [0.4 0.5 1 0.8 1];                % Starting points of free parameters (lambda_s and lambda_v)

            %---                       ---%
            %--- fit simple BADS model ---%
        
            % We specify the bounds and plausible bounds that (hopefully) 
            % contain the solution. Plausible bounds represent your best guess at 
            % bounding the region where the solution might lie.
            
            x0  = [parameters(1) parameters(2) parameters(3)];  % Starting point
            lb  = [0.1 0.1 0.1];                                % Lower bounds
            ub  = [0.9 0.9 5];                                 % Upper bounds
            plb = [0.2 0.2 0.2];                                % Plausible lower bounds
            pub = [0.9 0.9 5];                                  % Plausible upper bounds

            % Set BADS options
            options                         = bads('defaults');
            options.MaxObjectiveEvaluations = 50;
            options.MaxTime                 = 3600;

            % Run BADS, which returns the minimum X and the NLL.
            bFunc       = @(x) rbpf_lambdaBeta_full(x0,vol_o,vol_a,parameters,config);
            [x,fval]    = bads(bFunc,x0,lb,ub,plb,pub,[],options);
           
            % now run the model with the minimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_lambdaBetab_full(x,o,parameters,config);

            % store output for each subejct, block
            allX_model2{1,sub}{i,j}        = x;
            allNLL_model2{1,sub}{1,i}(:,j) = fval;
            all_lr_model2{1,sub}{i,j}      = lr;
            all_choice_model2{1,sub}{i,j}  = choice;

        end % end of volatility loop
    end % end of stochasticity loop
end% end of subjects loop

clear a o bo dat state ss vv x x0 fval lr choice val 

%% plot learning rates 

% loop over participants to plot learning rates
for sub = 1:nsubs
    sub_lr = all_lr_model2{1,sub};
    subnum = sprintf('%d',sub);

    % plot learning rates for each sub seperately
    [h, g, f] = plot_fitLRs(sub_lr,subnum);
end 

%% fit lambdas model to participnat data for each block seperately 

% exatact what we need from data struct and prepare data for model fitting 
state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

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

    % fit model at the block level -- loop over stc levels
    for i = 1:3

        % extract this stc outcomes and actions:
        stc_o = o(ss(:,i),:);
        stc_a = a(ss(:,i),:);

        for j = 1:2

            % extract this volatility outcomes and actions
            vol_o           = stc_o(vv(:,j),:);
            vol_a           = stc_a(vv(:,j),:);
            config          = struct('tvol',vv,'tstc',ss,'state',state,'nsim',100);
            parameters      = [0.3 0.5 0.5 0.3 1 100 1 1 0.001];  % parameters for modelling 
            % x               = [0.4 0.5 1 0.8 1];                % Starting points of free parameters (lambda_s and lambda_v)

            %---                       ---%
            %--- fit simple BADS model ---%
        
            % We specify the bounds and plausible bounds that (hopefully) 
            % contain the solution. Plausible bounds represent your best guess at 
            % bounding the region where the solution might lie.
            
            x0  = [parameters(1) parameters(2)];            % Starting point
            lb  = [0.1 0.1];                                % Lower bounds
            ub  = [0.9 0.9];                                % Upper bounds
            plb = [0.2 0.2];                                % Plausible lower bounds
            pub = [0.9 0.9];                                % Plausible upper bounds

            % Set BADS options
            options                         = bads('defaults');
            options.MaxObjectiveEvaluations = 50;
            options.MaxTime                 = 3600;

            % Run BADS, which returns the minimum X and the NLL.
            bFunc       = @(x) rbpf_lambdas(x0,vol_o,vol_a,parameters,config);
            [x,fval]    = bads(bFunc,x0,lb,ub,plb,pub,[],options);
           
            % now run the model with the minimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_lambdasb(x,o,parameters,config);

            % store output for each subejct, block
            allX_model3{1,sub}{i,j}        = x;
            allNLL_model3{1,sub}{1,i}(:,j) = fval;
            all_lr_model3{1,sub}{i,j}      = lr;
            all_choice_model3{1,sub}{i,j}  = choice;

        end % end of volatility loop
    end % end of stochasticity loop
end% end of subjects loop

clear a o bo dat state ss vv x x0 fval lr choice val 


%% plot learning rates 

% loop over participants to plot learning rates
for sub = 1:nsubs
    sub_lr = all_lr_model3{1,sub};
    subnum = sprintf('%d',sub);

    % plot learning rates for each sub seperately
    [h, g, f] = plot_fitLRs(sub_lr,subnum);
end 

%% for each subject and parameter, average parameter values across blocks

for sub = 1:nsubs

    this_sub_x  = allX{1,sub};
    this_sub_x2 = allX_model2{1,sub};
    this_sub_x3 = allX_model3{1,sub};

    for i = 1:3

        for j = 1:2
            
            % extract model 1 parameter values
            m1_lambda_s(i,j)    = this_sub_x{i,j}(1,1);
            m1_lambda_v(i,j)    = this_sub_x{i,j}(1,2);
            m1_beta(i,j)        = this_sub_x{i,j}(1,3);
            m1_s0(i,j)          = this_sub_x{i,j}(1,4);
            m1_v0(i,j)          = this_sub_x{i,j}(1,5);

            % extract model 2 parameter values
            m2_lambda_s(i,j)    = this_sub_x2{i,j}(1,1);
            m2_lambda_v(i,j)    = this_sub_x2{i,j}(1,2);
            m2_beta(i,j)        = this_sub_x2{i,j}(1,3);

            % extract model 3 parameter values
            m3_lambda_s(i,j)    = this_sub_x3{i,j}(1,1);
            m3_lambda_v(i,j)    = this_sub_x3{i,j}(1,2);

        end 
    end

    % average parameter values - model 1
    allsub_m1_lambda_s(sub,1)   = mean(m1_lambda_s,"all");
    allsub_m1_lambda_v(sub,1)   = mean(m1_lambda_v,"all");
    allsub_m1_beta(sub,1)       = mean(m1_beta,"all");
    allsub_m1_s0(sub,1)         = mean(m1_s0,"all");
    allsub_m1_v0(sub,1)         = mean(m1_v0,"all");

    % average parameter values - model 2
    allsub_m2_lambda_s(sub,1)   = mean(m2_lambda_s,"all");
    allsub_m2_lambda_v(sub,1)   = mean(m2_lambda_v,"all");
    allsub_m2_beta(sub,1)       = mean(m2_beta,"all");

    % average parameter values - model 3
    allsub_m3_lambda_s(sub,1)   = mean(m3_lambda_s,"all");
    allsub_m3_lambda_v(sub,1)   = mean(m3_lambda_v,"all");

end % end of subjects loop

%% 

