% MODEL FITTING OF SIMPLE RW MODEL -- VERSION 2 

% created 11/02/2023 

% 1. The script fits RW model using the fmincon optimiser for the computation of LOG-LIKELIHOOD (LL) 
% 2. Runs and plots basic parameter recovery 

% ------------------------------
%% init stuff

clear all  
clc

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);
addpath(fullfile(pwd,'modelfit'))

% define colours
global AZred AZblue

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;

%% define parameters and init variables

% initialise variables 
subjects        = 1;
params          = [.25 4];                % alpha and beta values 
condition       = 2;                      % only stable condition for now (if 2 = stable and volatile)
task            = 1;                      % stable without switch (if task = 2 then stable with one switch)

if condition == 1
    probs           = [.75 .25];          % probabilities of the stable condition
else
    probs           = [.75 .25; .80 .20]; % probabilities of the stable + volatile condition
end

labels          = {'alpha', 'beta'};      % for ploting
condstring      = {'stable', 'volatile'}; % for ploting 
trials          = 100;                    % per volatility condition

%% simulate dataset(s)

for sub = 1:subjects

    for cond = 1:condition
    % simulate dataset
        data                        = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output
        allsub_data{1,sub}{1,cond}  = data; % if many subjects add them to cell
    end

end % end of subjects loop

%% model the dataset(s)

for sub = 1:subjects
    
    for cond = 1:condition
        this_data                       = allsub_data{1,sub}{1,cond};
        [modelout]                      = modelRW_v1(params, this_data, outpath);
        allsub_modelout{1,sub}{1,cond}  = modelout;
    end

end % end of subjects loop

%% fit model (for each (simulated) subject and condition) with initial param values

% loop over subjects and conditions to fit the model
for sub = 1:subjects

    for cond = 1:condition

        this_data   = allsub_modelout{1,sub}{1,cond};

        % extract actions and rewards 
        actions     = this_data.a; actions = actions';
        rewards     = this_data.reward; rewards = rewards';
        [xfit, ll]   = modelfitRW_v1(params, actions, rewards);

        % store results
        fitted_alphas(sub,cond) = xfit(1);
        fitted_betas(sub,cond)  = xfit(2);
        lls(sub,cond)           = ll;
        
    end
end % end of subjects loop

% store them in a table
% modelfit_data = table(fitted_alphas(:,1), fitted_alphas(:,2), fitted_betas(:,1), fitted_betas(:,2), lls(:,1), lls(:,2));
modelfit_data = table(fitted_alphas, fitted_betas, lls);

%% fit model and prepare for parameter recovery

% as a first step let's run a few (hundreds) simulations with different
% alpha and beta parameter values 
rng(2);

for i = 1:1000
    
    % simulate alpha and beta values
%     tmp_alpha   = rand; 
%     tmp_beta    = exprnd(10);
%     tmp_params  = [tmp_alpha tmp_beta];

    for cond = 1:condition

        tmp_alpha   = rand; 
        tmp_beta    = exprnd(10);
        tmp_params  = [tmp_alpha tmp_beta];

        % simulate dataset 
        data{1,cond}        = avlearn_simulate_v1(cond, probs, trials, outpath, task);

        % run model
        modelout{1,cond}    = modelRW_v1(tmp_params, data{1,cond}, outpath);

        % extract actions and rewards
        actions     = modelout{1,cond}.a; actions = actions';
        rewards     = modelout{1,cond}.reward; rewards = rewards';

        % fit model 
        [xfit, ll]   = modelfitRW_v1(tmp_params, actions, rewards);
        
        % store simulated and fitted params
        Xsim{1,cond}(1,i)       = tmp_alpha;
        Xsim{1,cond}(2,i)       = tmp_beta;
        Xfitted{1,cond}(1,i)    = xfit(1);
        Xfitted{1,cond}(2,i)    = xfit(2);
        
    end % end of condition loop

end % end of iterations loop

%% plot parameter recovery

% won't probably use these
paramnames  = {'learning rate' 'softmax temperature'};
symbols     = {'\alpha' '\beta'};

% plot for each condition seperately 
for cond = 1:condition

    % make figure
    figure(1); clf;
    set(gcf, 'Position', [811   613   600   300])
    [~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);
    
    % extract condition simulated and fitted params
    cond_Xsim   = Xsim{1,cond};
    cond_Xfit   = Xfitted{1,cond};

    for j = 1:length(paramnames)

        axes(ax(j)); hold on;
        plot(cond_Xsim(j,:), cond_Xfit(j,:), 'o', 'color', AZred, 'markersize', 8, 'linewidth', 1)
        xl = get(gca, 'xlim');
        plot(xl, xl, 'k--')
    
    end % 

    % find 'bad' alpha values
    thresh = 0.25;
    ind = abs(cond_Xsim(1,:) - cond_Xfit(1,:)) > thresh;
    
    % mark the bad alphas in the plots
    for j = 1:2

        axes(ax(j));
        plot(cond_Xsim(j,ind), cond_Xfit(j,ind), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1, ...
        'markerfacecolor', [1 1 1]*0.5)
    end

    % define titles, x and y labels, etc...
    set(ax(1,2),'xscale', 'log', 'yscale' ,'log')

    axes(ax(1)); t = title('learning rate');
    axes(ax(2)); t(2) = title('softmax temperature');
    
    axes(ax(1)); xlabel('simulated \alpha'); ylabel('fit \alpha'); 
    axes(ax(2)); xlabel('simulated \beta'); ylabel('fit \beta'); 
    
    
    set(ax, 'tickdir', 'out', 'fontsize', 18)
    set(t, 'fontweight', 'normal')
    addABCs(ax(1), [-0.07 0.08], 32)
    addABCs(ax(2), [-0.1 0.08], 32, 'B')
    set(ax, 'tickdir', 'out')

        % add line to the softmax plot
    for j= 1:size(cond_Xsim,1)

        axes(ax(j));
        xl = get(gca, 'xlim');
        plot(xl, xl, 'k--')
    end

    %% scatter plot
%     figure(1); clf; hold on;
%     % plot(Xsim(1,:), Xsim(2,:),'.')
%     plot(cond_Xfit(2,:), cond_Xfit(1,:),'.')
%     set(gca, 'xscale', 'log')

end % end of cond

%% include bias to improve parameter recovery

% to be completed sometime later (if needed)


