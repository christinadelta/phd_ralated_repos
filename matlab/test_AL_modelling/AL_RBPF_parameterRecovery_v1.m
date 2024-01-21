% AL parameter recovery of RBPF models 

% created: January 2024 @christinadelta

%% 

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

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
nCues           = 2;

data            = ALsimdata_v3(probabilities, trials,condtrials);

%% run inference model (full model)

% how many simulations?
reps = 100;

for rep = 1:reps

    % define model parameters  
    % randomly generate values for the free parameterers of the full model:
    % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
    sim_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
    sim_lambda_s    = sim_lambda(1);
    sim_lambda_v    = sim_lambda(2);
    
    % 2. softmax temperatyre beta -- sould be between 0.1 and 10
    sim_beta        = 0.1 + (5-0.1).* rand(1,1); 
    
    % 3. s0 and v0 -- should be between 0.1 and 2
    sim_init_vals   = 0.1 + (2-0.1).* rand(2,1); 
    sim_s0          = sim_init_vals(1);
    sim_v0          = sim_init_vals(2);
    
    parameters      = struct('nparticles',500,'x0_unc',1,'lambda_v',sim_lambda_v,'lambda_s',sim_lambda_s,'v0',sim_v0,'s0',sim_s0, 'beta',sim_beta);
    
    % run pf_core function which includes both the kalman and particle filters 
    output          = simModel(data, probabilities, trials,condtrials,parameters);

    % store simulated X 
    simX(rep,1) = sim_lambda_s;
    simX(rep,2) = sim_lambda_v;
    simX(rep,3) = sim_beta;
    simX(rep,4) = sim_s0;
    simX(rep,5) = sim_v0;

    % extract what we need and prepare data for model fititng 
    state(:,1)  = data.x;
    state(:,2)  = 1 - state(:,1);
    ss          = data.stcind;
    vv          = data.t;
    test_o      = output.test_o;
    test_a      = output.test_a;
    test_a(find(test_a==0)) = 2; % convert action 2 to 0
    
    config      = struct('tvol',vv,'tstc',ss,'state',state,'nsim',1);


    % recover parameters 
    % loop over stc levels
    for i = 1:3

        % extract this stc outcomes and actions:
        stc_o = test_o(ss(:,i),:);
        stc_a = test_a(ss(:,i),:);

        for j = 1:2 % loop over volatility

            % extract this volatility outcomes and actions
            vol_o           = stc_o(vv(:,j),:);
            vol_a           = stc_a(vv(:,j),:);

            % re-generate parameter values for recovery
            % 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
            lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
            lambda_s    = lambda(1);
            lambda_v    = lambda(2);
            
            % 2. softmax temperatyre beta -- sould be between 0.1 and 10
            beta        = 0.1 + (5-0.1).* rand(1,1); 
            
            % 3. s0 and v0 -- should be between 0.1 and 2
            init_vals   = 0.1 + (2-0.1).* rand(2,1); 
            s0          = init_vals(1);
            v0          = init_vals(2);
            

            % fit model to the simulated outcomes and actions 
            params          = [lambda_s lambda_v beta s0 v0 100 1 1 0.001];  % parameters for modelling 

            %---                       ---%
            %--- fit simple BADS model ---%
            % We specify the bounds and plausible bounds that (hopefully) 
            % contain the solution. Plausible bounds represent your best guess at 
            % bounding the region where the solution might lie.
            x0  = [params(1) params(2) params(3) params(4) params(5)];                      % Starting point
            lb  = [0.1 0.1 0.1 0.1 0.1];                                                    % Lower bounds
            ub  = [0.9 0.9 5 2 2];                                                          % Upper bounds
            plb = [0.1 0.1 0.1 0.1 0.1];                                                    % Plausible lower bounds
            pub = [0.9 0.9 5 2 2];                                                          % Plausible upper bounds


            % Set BADS options
            options                         = bads('defaults');
            options.MaxObjectiveEvaluations = 50;
            options.MaxTime                 = 3600;

            % Run BADS, which returns the minimum X and the NLL.
            bFunc       = @(x) rbpf_core_full(x0,vol_o,vol_a,params,config);
            [x,fval]    = bads(bFunc,x0,lb,ub,plb,pub,[],options);
           
            % now run the model with the minimised x values to get subject learning
            % rates
            [val,vol,unp,lr,unc,choice] = rbpf_coreb_full(x,vol_o,params,config);

            % store fitted parameter values
            fitX{rep}{i,j}(1)   = x(1);
            fitX{rep}{i,j}(2)   = x(2);
            fitX{rep}{i,j}(3)   = x(3);
            fitX{rep}{i,j}(4)   = x(4);
            fitX{rep}{i,j}(5)   = x(5);
            fitLrs{rep}{i,j}    = lr;
            fitaction{rep}{i,j} = choice;

        end % end of volatility loop
    end % end of stc loop    
end % end of repetition loop

%% plot parameter recovery (full model)

% first re-arrange the parameters 
[Xfit, fitParams, simParams] = rearrangeX(simX, fitX,reps);

% plot simulated vs fitted parameters for each parameter seperately 
% plot lambdas_s
y           = fitParams.lambda_s;
x           = simParams.lambda_s;
ylm         = [0 1];
colours     = 'magenta'; % let's start with magenta?
param_title = 'lambda_s';
a           = plotScatter(y, x, ylm, colours,param_title);

% plot lambda_v 
y           = fitParams.lambda_v;
x           = simParams.lambda_v;
ylm         = [0 1];
colours     = "red"; % let's start with magenta?
param_title = 'lambda_v';
b           = plotScatter(y, x, ylm, colours,param_title);

% plot beta 
y           = fitParams.betas;
x           = simParams.beta;
ylm         = [0 10];
colours     = "blue"; % let's start with magenta?
param_title = 'beta';
c           = plotScatter(y, x, ylm, colours,param_title);

% plot s0 
y           = fitParams.s0;
x           = simParams.s0;
ylm         = [0 2];
colours     = "cyan"; % let's start with magenta?
param_title = 's_0';
d           = plotScatter(y, x, ylm, colours,param_title);

% plot v0 
y           = fitParams.v0;
x           = simParams.v0;
ylm         = [0 2];
colours     = "green"; % let's start with magenta?
param_title = 'v_0';
e           = plotScatter(y, x, ylm, colours,param_title);

%%

