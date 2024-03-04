% fitting RW model to participant data 

% what does this version include:
% 1. fitting data to full task 
% 2. three free parameters 

% created February 2024

%% clc 
clear all
clc

%% dock figures first

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
subjects        = 1;
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

%% test fitting/objective function with one subject to make sure it works

state(:,1)      = data.x;
state(:,2)      = 1 - state(:,1);
ss              = data.stcind;
vv              = data.t; 
dat             = alldata{1,1};

% starting parameter values 
starting_params = [0.2 2 0.05];
nll             = lik_modelRW_toSub_v2(starting_params, dat);

%% fit model to participnant data 

for sub = 1:nsubs

    dat                 = alldata{sub,1};
    starting_params     = [0.25 2 0.05];

    % define and run the objective function
    obFunc              = @(x) lik_modelRW_toSub_v2([x(1), x(2), x(3)],dat);

    X0                  = starting_params;
    LB                  = [0 0 0];
    UB                  = [1 5 1];
    [Xfit, NegLL]       = fmincon(obFunc, X0, [], [], [], [], LB, UB);

    allfit_alpha(sub)   = Xfit(1);
    allfit_beta(sub)    = Xfit(2);
    allfit_decay(sub)   = Xfit(3);
    allNLL(sub)         = NegLL;

     % fit model with participant data and optimal param values
    [dat, ~]                = fit_modelRW_toSub_v1(Xfit, dat);
    fit_PP_v2(:,:,sub)      = dat.allPs;
    fit_Qval_v2(:,:,sub)    = dat.Qvals;
    fit_correct_v2(:,sub)   = dat.correct;

end % end of subject loop

%% get average of Q values and choice probs and plot 

av_PPs                  = mean(fit_PP_v2,3);           % average over subjects
av_Qvals                = mean(fit_Qval_v2, 3);  
fit_perf                = mean(fit_correct_v2,"all");  % model performance 
sub_perf                = mean(allsub_mAL_v2,2);       % subject performance 

% Plotting Q-values for Option A
figure; 
subplot(2,1,1); % For two subplots in one figure
plot(state(:,1), 'k-', 'LineWidth', 1); hold on; % Plot true reward rate for option A
% plot(av_Qvals(:,2), 'r--', 'LineWidth', 2);      % Plot average Q-values for option A
plot(av_PPs(:,1), 'r:', 'LineWidth', 2); % Plot average choice probabilities for option A
plot(av_PPs(:,2), 'b:', 'LineWidth', 2); % Plot average choice probabilities for option B
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
hSimR = plot(av_Qvals(:,1), 'r--', 'LineWidth', 1, 'DisplayName', 'Averaged Q Value Option A');
hSimB = plot(av_Qvals(:,2), 'b--', 'LineWidth', 1, 'DisplayName', 'Averaged Q Value Option B');

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
sum_matrix = zeros(3,2);

% Loop through each cell to sum the matrices
for i = 1:length(allsub_mStc)
    sum_matrix = sum_matrix + allsub_mStc{i};
end

% Calculate the average matrix
average_matrix = sum_matrix / length(allsub_mStc);

% Display the result
disp(average_matrix);

% now let's look at the model performance for each condition
for s = 1:3
    stc_perf                    = fit_correct_v2(ss(:,s),:);
    for v = 1:2
        model_matrix(s,v)     = mean(stc_perf(vv(:,v),:),"all");
    end % end of volatility 
end % end of 

disp(model_matrix);

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

%% plot parameter values 

% Determine bin edges for a unified scale or based on both datasets
binEdges = linspace(min([allfit_alpha'; allfit_beta']), max([allfit_alpha'; allfit_beta']), 12); % 20 bins across the combined range

% Create a figure
figure;

% First histogram (Parameter 1)
subplot(1,3,1); % 1-row, 2-column grid, 1st subplot
histogram(allfit_alpha','BinLimits',[0, 1],'NumBins', 20, 'Normalization', 'probability', 'FaceColor', [0.2 0.2 0.5]);
title('\alpha');
xlabel('Value');
ylabel('Probability');
fontsize(gcf,16,"points")
grid on; % Adds grid lines for better readability

% Second histogram (Parameter 2)
subplot(1,3,2); % 1-row, 2-column grid, 2nd subplot
histogram(allfit_beta','BinLimits',[0, 5],'NumBins', 20,'Normalization', 'probability', 'FaceColor', [0.5 0.2 0.2]);
title('\beta');
xlabel('Value');
ylabel('Probability');
fontsize(gcf,16,"points")
grid on;

% Second histogram (Parameter 2)
subplot(1,3,3); % 1-row, 2-column grid, 2nd subplot
histogram(allfit_decay','BinLimits',[0, 1],'NumBins', 20,'Normalization', 'probability', 'FaceColor', [0 0.4470 0.7410]);
title('\lambda');
xlabel('Value');
ylabel('Probability');
fontsize(gcf,16,"points")
grid on;

% Super title for the figure
sgtitle('Normalized Distributions of Learning Rate, Inverse Temperature and Decay Parameters');

%%

