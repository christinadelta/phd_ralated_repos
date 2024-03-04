% Action Learning RW simulations version 2

% created February 2024

% added a decay/forgeting parameter to model the natural forgetting process
% or decrease in salience of unchosen options over time 

% 3 free parameters:
% learning rate (alpha)
% inverse temperature (beta)
% decay (lambda)


%% housekeeping commands

clear all
clc

%% set figure-docking as default 

set(0,'DefaultFigureWindowStyle','docked')

%% define paths etc 
% paths for subjects data 
startpath       = pwd;
datadir         = fullfile(startpath,'data');
subs            = dir(fullfile(datadir, '*sub-*'));
nsubs           = length(subs);

% figpath         = fullfile(pwd, 'figures');     addpath(figpath);

%% get subject data
% 
[all_ALTdata, allsub_mAL, allsub_mVol, allsub_mStc, m_rtsAL, allsub_mrtsStc, allsub_mrtsVol,allsub_mrtsAL] = clean_ALdata(datadir,subs,nsubs);


%% extract participant choice data and rts

for sub = 1:nsubs 
    allchoices{1,sub}       = all_ALTdata{1,sub}(:,9); % get responses 
    alloutcomes{1,sub}(:,1) = all_ALTdata{1,sub}(:,5); % get outcomes (low prob -- good option outcomes) 
    alloutcomes{1,sub}(:,2) = all_ALTdata{1,sub}(:,6); % get outcomes (high prob -- bad option outcomes) 
    allrts{1,sub}           = all_ALTdata{1,sub}(:,10);

    alldata{sub,1}.actions  = allchoices{1,sub};
    alldata{sub,1}.outcome  = alloutcomes{1,sub};
    alldata{sub,1}.rts      = allrts{1,sub};

end 

%% simulate some data to extract underlying loss rate stc and vol indecies

% initialise variables 
% subjects        = 1;
condition       = 6;                        % stable & volatile / small, medium & large stochasticity
task            = 2;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 70:30 or 60:40)
total_trials    = 140;                      % total trials
condtrials      = {70,[30,10,10,20]};
nCues           = 2;
params          = [0.3 1.5 0.1];            % three parameters

data            = ALsimdata_v2(probabilities, total_trials,condtrials);

%% test fitting objective function with one subject to make sure it works

state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 
dat             = alldata{1,1};

o               = dat.outcome(:,1); % extract good options only
a               = dat.actions;

for s = 1:3

    % extract this stc level outcomes and actions
    stc_o       = o(ss(:,s),:);
    stc_a       = a(ss(:,s),:);
    for v = 1:2

        test_data.outcome       = stc_o(vv(:,v),:);
        test_data.actions       = stc_a(vv(:,v),:);

        % starting parameter values 
        starting_params             = [0.2 1 0.05];
        nll(s,v)                    = lik_modelRW_toSub_v2(starting_params, test_data);

    end % end of volatility loop
end % end of stc loop

%% fit model to participnant data 

state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 

for sub = 1:nsubs

    % extract subject data
    dat             = alldata{sub,1};
    o               = dat.outcome(:,1); % extract good options only
    a               = dat.actions;

    for s = 1:3

        % extract this stc level outcomes and actions
        stc_o       = o(ss(:,s),:);
        stc_a       = a(ss(:,s),:);
        for v = 1:2
    
            sub_data.outcome            = stc_o(vv(:,v),:);
            sub_data.actions            = stc_a(vv(:,v),:);
    
            % starting parameter values 
            starting_params             = [0.2 1 0.05];

            % define the objective function
            obFunc                      = @(x) lik_modelRW_toSub_v2([x(1), x(2), x(3)],sub_data);

            X0                          = starting_params;
            LB                          = [0 0 0];
            UB                          = [1 10 1];
            [Xfit, NegLL]               = fmincon(obFunc, X0, [], [], [], [], LB, UB);
        
            allfit_alpha(s,v,sub)       = Xfit(1);
            allfit_beta(s,v,sub)        = Xfit(2);
            allfit_decay(s,v,sub)       = Xfit(3);
            allNLL(s,v,sub)             = NegLL;
        
            % fit model with participant data and optimal param values
            [sub_output, ~]             = fit_modelRW_toSub_v2(Xfit, sub_data);
            fit_PP_v2{s,v}(:,:,sub)     = sub_output.allPs;
            fit_Qval_v2{s,v}(:,:,sub)   = sub_output.Qvals;
            fit_correct_v2{s,v}(:,sub)  = sub_output.correct;

        end % end of volatility loop
    end % end of stc loop
end % end of subjects loop

%% get average of Q values and choice probs and re-arrange

count   = 1;

% for each stc/vol condition average choice probabilities over model
% fitting instances (3rd dimension)
for s = 1:3
    for v = 1:2
        all_pp{1,count}     = mean(fit_PP_v2{s,v},3);
        all_vv{1,count}     = mean(fit_Qval_v2{s,v},3);
        model_perf(:,count) = mean(fit_correct_v2{s,v},1)';
        count               = count + 1; % update counter 

    end % end of volatility loop
end % end of stochasticity loop


% loop over subject and re-arrange performance 
for sub = 1:nsubs
    % start counter 
    count = 1;
    for s = 1:3
        for v = 1:2
            sub_perf(sub,count) = allsub_mStc{1,sub}(s,v);
            count               = count + 1;

        end % end of volatility loop
    end % end of stc loop
end % end of subject loop

% concat trials of choice probabilities and q values
PPs = vertcat(all_pp{:});
VVs = vertcat(all_vv{:});

%% plot Qvals and Choice probs

% Plotting Q-values for Option A
figure; 
subplot(2,1,1); % For two subplots in one figure
plot(state(:,1), 'k-', 'LineWidth', 1); hold on; % Plot true reward rate for option A
% plot(av_Qvals(:,2), 'r--', 'LineWidth', 2);      % Plot average Q-values for option A
plot(PPs(:,1), 'r:', 'LineWidth', 2); % Plot average choice probabilities for option A
plot(PPs(:,2), 'b:', 'LineWidth', 2); % Plot average choice probabilities for option B
xlim([1 420]);
xlabel('Trial');
ylabel('Probability of Reward');
title('True Reward Rate and Choice Probabilities for Option A & B');
legend('True Reward Rate', 'Average Choice Probabilities A',  'Average Choice Probabilities B', 'Location', 'best');
fontsize(gcf,16,"points")
hold off;

% Subplot 2: Q Values
subplot(2, 1, 2); % Selects the 2nd subplot in a 2-row, 1-column grid
hold on;

% plot lines to specify where change happens 
dx          = [0; diff(state(:,1))~=0];    
ix          = (find(dx))-1;       % trials before each switch of probabilities
tt          = 1:length(state(:,1));        % total trials 

% Plot and store handles for lines you want in the legend
hSimR = plot(VVs(:,1), 'r--', 'LineWidth', 1, 'DisplayName', 'Averaged Q Value Option A');
hSimB = plot(VVs(:,2), 'b--', 'LineWidth', 1, 'DisplayName', 'Averaged Q Value Option B');

ylim([-1 1]); % Adjust if your Q values fall outside this range
xlim([1 420]);

ym      = get(gca,'ylim');
% Specify trial ranges where volatility is high and plot vertical lines
for i = 1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',"black",'linewidth',0.5);
end

% Only include specified handles in the legend
legend([hSimR, hSimB], 'Averaged Q Value Option A', 'Averaged Q Value Option B');

% Formatting the subplot
title('Averaged fitted Q Values Over Trials');
xlabel('Trial');
ylabel('Q Value');

fontsize(gcf,16,"points")
hold off;

%% look at model performance vs subject performance 

% average performance across subjcts 
% init sum matrix
sum_matrix      = zeros(3,2);

% loop through each cell to sum the matrices
for i = 1:length(allsub_mStc)
    sum_matrix  = sum_matrix + allsub_mStc{i};
end

% Calculate the average matrix
average_matrix  = sum_matrix / length(allsub_mStc);

% re-arrange 
subj_matrix     = [average_matrix(1,1), average_matrix(1,2),...
    average_matrix(2,1), average_matrix(2,2), average_matrix(3,1), average_matrix(3,2)];

% Display the result
disp(average_matrix);

% average performance across model instances 
temp            = mean(model_perf,1);
model_matrix    = [temp(1), temp(2); temp(3), temp(4); temp(5), temp(6)]; % re-arrange model matrix


disp(model_matrix);

%% plot participant vs model behaviour (per stochasticity)

% plot small vs medium vs large
stc_averages_h = mean(average_matrix, 2); % Average across columns
stc_averages_m = mean(model_matrix, 2); % Average across rows
stc_errors_h = std(average_matrix, 0, 2); % Standard deviation across rows for error bars
stc_errors_m = std(model_matrix, 0, 2); % Standard deviation across rows for error bars

% Plotting
figure;
subplot(1,2,1); % For two subplots in one figure
bar(stc_averages_h');
hold on; % Keep the bar graph and plot error bars on top

% Error bar positions and errors
errorbar(1:size(average_matrix, 1), stc_averages_h', stc_errors_h', '.');
title('Average Human Performance per Stochasticity');
xlabel('Stochasticity conditions');
ylabel('Average Performance');
xticks(1:size(average_matrix, 1));
xticklabels({'small', 'medium', 'large'})
fontsize(gcf,16,"points")
hold off;

subplot(1,2,2); % For two subplots in one figure
bar(stc_averages_m');
hold on; % Keep the bar graph and plot error bars on top

% Error bar positions and errors
errorbar(1:size(model_matrix, 1), stc_averages_m', stc_errors_m', '.');
title('Average Model Performance per Stochasticity');
xlabel('Stochasticity conditions');
ylabel('Average Performance');
xticks(1:size(model_matrix, 1));
xticklabels({'small', 'medium', 'large'})
fontsize(gcf,16,"points")
hold off;

%% plot participant vs model behaviour (per volatility)

% plot stable vs volatile 
vol_averages_h = mean(average_matrix, 1); % Average across rows
vol_averages_m = mean(model_matrix, 1); % Average across rows
vol_errors_h = std(average_matrix, 0, 1); % Standard deviation across rows for error bars
vol_errors_m = std(model_matrix, 0, 1); % Standard deviation across rows for error bars

% Plotting
figure;
subplot(1,2,1); % For two subplots in one figure
bar(vol_averages_h);
hold on; % Keep the bar graph and plot error bars on top

% Error bar positions and errors
errorbar(1:size(average_matrix, 2), vol_averages_h, vol_errors_h, '.');
title('Average Human Performance per Volatility');
xlabel('Volatility conditions');
ylabel('Average Performance');
xticks(1:size(average_matrix, 2));
xticklabels({'stable', 'volatile'})
fontsize(gcf,16,"points")
hold off;

subplot(1,2,2); % For two subplots in one figure
bar(vol_averages_m);
hold on; % Keep the bar graph and plot error bars on top

% Error bar positions and errors
errorbar(1:size(model_matrix, 2), vol_averages_m, vol_errors_m, '.');
title('Average model Performance per Volatility');
xlabel('Volatility conditions');
ylabel('Average Performance');
xticks(1:size(model_matrix, 2));
xticklabels({'stable', 'volatile'})
fontsize(gcf,16,"points")
hold off;

%% prepare parameter values for ploting 

% average estimated parameter values over model instances 
av_fit_alpha    = mean(allfit_alpha,3);
av_fit_beta     = mean(allfit_beta,3);
av_fit_decay    = mean(allfit_decay,3);

% plot parameter 1 (alpha) blockwise 
% Define the conditions for labels
stochasticity_levels = {'Small Stoch', 'Medium Stoch', 'Large Stoch'};
volatility_levels = {'Stable', 'Volatile'};

%% plot parameter values alpha

% Create a figure
figure;

% Loop through each condition
for i = 1:size(allfit_alpha,1) % Stochasticity levels
    for j = 1:size(allfit_alpha,2) % Volatility levels
        
        % Calculate subplot index
        subplotIndex = (i-1)*size(allfit_alpha,2) + j;
        
        % Select subplot
        subplot(size(allfit_alpha,1), size(allfit_alpha,2), subplotIndex);
        
        % Extract data for current condition
        thisdata = squeeze(allfit_alpha(i,j,:));
        
        % Plot histogram
        histogram(thisdata, 'BinLimits', [0, 1], 'NumBins', 20, ...
                  'Normalization', 'probability', 'FaceColor', [0.2 0.2 0.5]);
        
        % Title and labels
        title(sprintf('%s, %s', stochasticity_levels{i}, volatility_levels{j}));
        xlabel('Alpha Value');
        ylabel('Probability');
        
        % Set font size and add grid
        set(gca, 'FontSize', 14); % Adjust font size here
        grid on;
    end
end

% Super title for the figure
sgtitle('Normalized Distributions of Learning Rate \alpha');


% Adjust overall figure properties
% set(gcf, 'Position', [100, 100, 1200, 600]); % Adjust figure size

%% plot parameter values beta

% Create a figure
figure;

% Loop through each condition
for i = 1:size(allfit_beta,1) % Stochasticity levels
    for j = 1:size(allfit_beta,2) % Volatility levels
        
        % Calculate subplot index
        subplotIndex = (i-1)*size(allfit_beta,2) + j;
        
        % Select subplot
        subplot(size(allfit_beta,1), size(allfit_beta,2), subplotIndex);
        
        % Extract data for current condition
        thisdata = squeeze(allfit_beta(i,j,:));
        
        % Plot histogram
        histogram(thisdata, 'BinLimits', [0, 10], 'NumBins', 20, ...
                  'Normalization', 'probability', 'FaceColor', [0.5 0.2 0.2]);
        
        % Title and labels
        title(sprintf('%s, %s', stochasticity_levels{i}, volatility_levels{j}));
        xlabel('Beta Value');
        ylabel('Probability');
        
        % Set font size and add grid
        set(gca, 'FontSize', 14); % Adjust font size here
        grid on;
    end
end

% Super title for the figure
sgtitle('Normalized Distributions of Inverse Temperature \beta');


% Adjust overall figure properties
% set(gcf, 'Position', [100, 100, 1200, 600]); % Adjust figure size

%% plot parameter values decay

% Create a figure
figure;

% Loop through each condition
for i = 1:size(allfit_decay,1) % Stochasticity levels
    for j = 1:size(allfit_decay,2) % Volatility levels
        
        % Calculate subplot index
        subplotIndex = (i-1)*size(allfit_decay,2) + j;
        
        % Select subplot
        subplot(size(allfit_decay,1), size(allfit_decay,2), subplotIndex);
        
        % Extract data for current condition
        thisdata = squeeze(allfit_decay(i,j,:));
        
        % Plot histogram
        histogram(thisdata, 'BinLimits', [0, 1], 'NumBins', 20, ...
                  'Normalization', 'probability', 'FaceColor', [0 0.4470 0.7410]);
        
        % Title and labels
        title(sprintf('%s, %s', stochasticity_levels{i}, volatility_levels{j}));
        xlabel('Decay Value');
        ylabel('Probability');
        
        % Set font size and add grid
        set(gca, 'FontSize', 14); % Adjust font size here
        grid on;
    end
end

% Super title for the figure
sgtitle('Normalized Distributions of Decay Parameter \lambda');

% Adjust overall figure properties
% set(gcf, 'Position', [100, 100, 1200, 600]); % Adjust figure size

%%


