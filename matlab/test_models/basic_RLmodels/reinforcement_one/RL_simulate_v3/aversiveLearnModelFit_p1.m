% Model fitting for the simple RW model (Part of aversive-learning
% modelling) -- PART 1

% Model fitting using fmincon for the computation of LOG-LIKELIHOOD (LL) 

% created: 04/02/2023

% STEPS:
% simulate dataset(s)
% run model 
% using actions and rewards, fit model to obtain -ll

% -----------------------------------------

clear all  
clc

% set paths
addpath(fullfile(pwd,'functions'))
outpath = fullfile(pwd, 'output'); addpath(outpath);
figpath = fullfile(pwd, 'figures'); addpath(figpath);
addpath(fullfile(pwd,'modelfit'))

%% define parameters and init variables

% initialise variables 
subjects        = 1;
params          = [.25 4];                % alpha and beta values 
condition       = 1;                      % only stable condition for now (if 2 = stable and volatile)
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
        [xfit ll]   = modelfitRW_v1(params, actions, rewards);

        % store results
        fitted_alphas(sub,cond) = xfit(1);
        fitted_betas(sub,cond)  = xfit(2);
        lls(sub,cond)           = ll;
        
    end
end % end of subjects loop

%% fit model for different combinations of parameter values 

% define range of alphas and betas
paramrange  = [20 25]; % 20 alphas and 25 betas

% initialise matrix for fitted paramters (to be used for plotting)
fitted.alpha.pdf    = nan(2,subjects,paramrange(1));
fitted.beta.pdf     = nan(2,subjects,paramrange(2));
fitted.alpha.ml     = nan(2,subjects);
fitted.beta.ml      = nan(2,subjects);
fitted.alpha.ev     = nan(2,subjects);
fitted.beta.ev      = nan(2,subjects);

% bounds for alpha and betas values
bounds              = [0 1 % alpha bounds
    0 15]; % beta bounds 

% compute binspace for alpha and beta simulated params
bin_alphas          = params(1)/bounds(1,2)*paramrange(1);
bin_betas           = params(2)/bounds(2,2)*paramrange(2);

for sub = 1:subjects 

    for cond = 1:condition

        % extract dataset to be used 
        this_data   = allsub_modelout{1,sub}{1,cond};

        % extract actions and rewards 
        actions     = this_data.a; actions = actions';
        rewards     = this_data.reward; rewards = rewards';

        % generated params ranges
        for i = 1:length(params)

            range       = linspace(bounds(i,1), bounds(i,2), paramrange(i)+1);
            simpars{i}  = range(2:end); % exlude zero

        end

        for alpha = 1:paramrange(1)

            % extract this iter alpha
            tmp_alpha = simpars{1}(alpha);

            for beta = 1:paramrange(2)

                tmp_beta = simpars{2}(beta);

                tmp_params                  = [tmp_alpha tmp_beta];
                [ll fitout]                 = likRW_v1(actions, rewards, tmp_alpha, tmp_beta);
                % [xfit ll]                   = modelfitRW_v1(tmp_params, actions, rewards);
                
                % store ll
                all_lls{sub,cond}(alpha,beta)    = ll;

            end % end of betas loop
        end % end of alphas loop
        
        % remove mnimimum in ll
        cond_ll     = all_lls{sub,cond};
        loglik      = cond_ll-min(cond_ll(:));

        % compute likelihood instead of log 
        lik         = exp(loglik);

        % compute marginal probabilities 
        for x = 1:length(params)

            tmp = sum(lik,3-x);
            marglik{x} = tmp/sum(tmp);

        end

        % plot the likelihood landschape with the maximum and true
        % value (latter for simulations only), for each individual and
        % condition
        ML(1) = myvect(simpars{1}(max(marglik{1})==marglik{1}));
        ML(2) = myvect(simpars{2}(max(marglik{2})==marglik{2}));
        EV(1) = sum(simpars{1}(:).*marglik{1}(:));
        EV(2) = sum(simpars{2}(:).*marglik{2}(:));
        
        % 
        [ll estout{sub,cond,1}] = likRW_v1(actions, rewards, ML(1), ML(2));
        [ll estout{sub,cond,2}] = likRW_v1(actions, rewards, EV(1), EV(2));
        

    end  % end of condition loop
end % end of subjects loop

%% plot model fit results 


