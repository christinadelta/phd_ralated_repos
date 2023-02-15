%% Optimising sequence of trials
%  Creating possible sequences of outcomes & conditions
%  Testing to find most recoverable
%
%%

clearvars;
close all;

% Some of the functions are in subfolder, so adding to path
addpath('./SimFuns')
addpath('./FitFuns')
addpath('./HelperFuns')



%% run the same RL agent across 2 blocks of 2 intermixed sequences 
%  this task has 2 conditions 

dataInfo.type   = 'virtual';
dataInfo.subjN  = 100; % number of simulated agents
dataInfo.subjIDs= 1:dataInfo.subjN;

% for outputing the simulations
dataInfo.simPath = './SimResults';
dataInfo.figPath = './Figs';

dataInfo.outfname = 'simM2vAll_p70_t30_b8';

% Trial configuration information
cfg.nTrialsBlk  = 30;           % nTrials per Condition
cfg.nBlocks     = 8;
cfg.nTrialsTot  = cfg.nTrialsBlk*cfg.nBlocks;

cfg.rewProbs    = [.70, .30];   % Reward probabilities for Good vs. Bad stim
cfg.rewVals     = [1, -1];       % Reward values (good, bad)
cfg.rewAssign   = 0; % not random

cfg.nOpts       = 2; % per block, choice options
cfg.nSymbols    = cfg.nBlocks*cfg.nOpts;  % because here I am displaying some symbols to represent the choice options

cfg.shufflef = @(v)v(randperm(numel(v))); % anonymous function to randomly shuffle an array (v)

% names of variables for trial info matrix (as it's a table)
cfg.trialListVars = {'blockN', 'trialN', 'stimID_hi', 'stimID_lo', 'stimLocRev', 'stimID_l', 'stimID_r', 'out_hi', 'out_lo', 'out_l', 'out_r'};

dataInfo.cfg = cfg;

%% Models
% here I use indices to refer to models, where model = 1 is the classic RL

% fit with a priori priors
modelInfo.priors  = 0; % logical,  0 means no priors, or 1 with
if modelInfo.priors == 1
    modelInfo.priorsType = 'apriori';
else
    modelInfo.priorsType = 'noPriors';
end

%  1 - classic RL: 1 beta, 1 alpha
modelInfo.params{1}.labels      = {'b', 'a'};
modelInfo.params{1}.optStart    = [1   .5]; % plausible, but will reset in multstart
modelInfo.params{1}.optMin      = [0   0];
modelInfo.params{1}.optMax      = [50  1];


%% Simulate RL data across a range of parameters

simInfo.type = 'data';
simInfo.Qinit = repelem(mean(cfg.rewVals),1,length(cfg.rewVals));   % this will ensure the Q values are initiated to the mean of the possible reward values

modelInfo.Qinit = simInfo.Qinit;

simModel = 1; 

%% Unbiased

% draw random values for beta & alpha parameters
% assuming some priors & distributions 
inParameters = [gamrnd(2, 5, 1, dataInfo.subjN)',... % beta
    betarnd(1.1, 1.1, 1, dataInfo.subjN)']; % alpha

% store
simInfo.inParameters = inParameters;

figure;
sgtitle(['Simulated Parameters - ', dataInfo.outfname], 'Interpreter', 'none')

subplot(1,2,1);
histogram(inParameters(:,1),50);
title(modelInfo.params{simModel}.labels(1), 'Interpreter', 'none')
subplot(1,2,2)
histogram(inParameters(:,2),50);
title(modelInfo.params{simModel}.labels(2), 'Interpreter', 'none')

saveas(gcf, fullfile(dataInfo.figPath, ['SimParamHist_' dataInfo.outfname '.png']))


%% Model Simulations

% all simulated
dataInfo.allData = [];
for iSubj = 1:dataInfo.subjN
    params = inParameters(iSubj,:);
    
    cfg = setupTrials(cfg); % create the trial information
    
    subjData = cfg.learnBlocks(:, {'blockN', 'trialN', 'stimLocRev', 'out_hi', 'out_lo'});
    
    subjData = sim_compModels(params, simModel, simInfo, subjData);
    subjData.subjID(:) = dataInfo.subjIDs(iSubj); % add subjID info

    dataInfo.allData = [dataInfo.allData; subjData];   
    clear subjData params
end

fprintf('\n\n*******    Simulation Finished    *******\n');


%% Model fitting

modelInfo.whichms = 1; % indices of models to fit

modelInfo.optimOptions = optimset('Algorithm', 'interior-point', 'Display', 'notify', 'MaxIter', 10000);

outParameters = cell(1, numel(modelInfo.whichms ));
loglik = NaN(dataInfo.subjN, numel(modelInfo.whichms ));
report = loglik;
gradient = cell(dataInfo.subjN, numel(modelInfo.whichms ));
hessian = gradient;

    
for iMod = modelInfo.whichms 
    for iSubj = 1:dataInfo.subjN

        subjData = dataInfo.allData(dataInfo.allData.subjID==dataInfo.subjIDs(iSubj), :);
        
        if modelInfo.priors
            [outParameters{iMod}(iSubj,:),...
                 loglik(iSubj, iMod),...  % if no priors, it's the negative log likelihood
                 report(iSubj, iMod),...  % exit flags <= 0 imply problems
                 gradient{iSubj, iMod},...
                 hessian{iSubj, iMod}] = ...
                        fmincon(@(x) aprioriPriors(x, iMod, simInfo.Qinit, subjData, modelInfo),...
                            modelInfo.params{iMod}.optStart,... 
                            [],[],[],[],...
                            modelInfo.params{iMod}.optMin,... 
                            modelInfo.params{iMod}.optMax,...
                            [], modelInfo.optimOptions);                         

        else 
            [outParameters{iMod}(iSubj,:),...
                 loglik(iSubj, iMod),...  % if no priors, it's the negative log likelihood
                 report(iSubj, iMod),...  % exit flags <= 0 imply problems
                 gradient{iSubj, iMod},...
                 hessian{iSubj, iMod}] = ...
                        fmincon(@(x) fit_compModels(x, iMod, simInfo.Qinit, subjData, modelInfo),...
                            modelInfo.params{iMod}.optStart,... 
                            [],[],[],[],...
                            modelInfo.params{iMod}.optMin,... 
                            modelInfo.params{iMod}.optMax,...
                            [], modelInfo.optimOptions);                         
        end

        clear subjData
    end
end


% % If using priors, need to get loglikelihood & ICs, otherwise already have the ll
if modelInfo.priors && strcmp(modelInfo.priorsType, 'apriori')
    [bic, aic, ll] = get_ll_fromapriori(outParameters, modelInfo, dataInfo);

else % just Add Information Criteria
    bic = zeros(dataInfo.subjN, numel(modelInfo.whichms));
    aic = bic;

    for iSubj=1:dataInfo.subjN
        subjData = dataInfo.allData(dataInfo.allData.subj == dataInfo.subjIDs(iSubj), :);

        ll= abs(loglik(iSubj)); % this ensures values are positive, so in formula multiply by 2 to obtain positive BIC/AICs, where smaller is better(normal -2 is assuming logLik is negative)
        nTrials = height(subjData);
        df = numel(modelInfo.params{model}.labels);
        bic(iSubj, model) = 2*ll + df*log(nTrials); 
        aic(iSubj, model) = 2*ll + 2*df;    
    end        
end
fprintf('\n\n*******    Optimisation Finished    *******\n');
        

%% Save simulation and fitting output
save(fullfile(dataInfo.simPath, ['Sim&Recovery_' dataInfo.outfname '.mat']));
fprintf('\n\n*******    Results Saved!     *******\n');




%%  Heatmap of Correlations between simulated & recovered parameters 

% % % %  For M1

simModel=1;

% save
recovery = NaN(numel(modelInfo.params{simModel}.labels), numel(modelInfo.params{simModel}.labels));

for iSimP = 1: numel(modelInfo.params{simModel}.labels)
    for iRecP = 1: numel(modelInfo.params{simModel}.labels)    
        recovery(iSimP, iRecP) = corr(inParameters(:,iSimP), outParameters{simModel}(:,iRecP));
    end
end


% % Plot correlation matrices
colRange = [-.3, 1]; % revising given results
if min(recovery,[],'all')< colRange(1)
    warning('Colour scale needs adjusted - minimum smaller than range')
end


figure;
imagesc(recovery, colRange);
xlabel('Simulated')
xticks(1:numel(modelInfo.params{simModel}.labels))
xticklabels(modelInfo.params{simModel}.labels)                
ylabel('Recovered')
yticks(1:numel(modelInfo.params{simModel}.labels))
yticklabels(modelInfo.params{simModel}.labels)
colormap('jet')
colorbar('EastOutside')
title(['Correlations between Simulated & Recovered Parameters - ', dataInfo.outfname], 'Interpreter', 'none')

saveas(gcf, fullfile(dataInfo.figPath, ['SimVRecCorrMat_' dataInfo.outfname '.png']))



%% Dist of sim v rec params
figure;
set(gcf, 'Position', [200   800   1000   400])

sgtitle(['Simulated vs Recovered Parameters - ', dataInfo.outfname], 'Interpreter', 'none')

for iP = 1: numel(modelInfo.params{simModel}.labels)

    subplot(1, numel(modelInfo.params{simModel}.labels), iP)
    plot(inParameters(:,iP), outParameters{simModel}(:,iP), 'o', 'markersize', 8, 'linewidth', 1);
    hold on;
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
    title(modelInfo.params{simModel}.labels(iP))
    xlabel('Simulated')
    ylabel('Recovered')
    hold off;
end

saveas(gcf, fullfile(dataInfo.figPath, ['SimVRecCorrPlots_' dataInfo.outfname '.png']))




%% Performance levels


pHiXtrialN=NaN(dataInfo.subjN, cfg.nTrialsBlk);
for iSubj = 1:dataInfo.subjN
    for iT = 1:cfg.nTrialsBlk
        % choice=1 is high value choice
        pHiXtrialN(iSubj, iT) = sum((dataInfo.allData.subjID==iSubj & dataInfo.allData.trialN==iT) & dataInfo.allData.action==1) / sum(dataInfo.allData.subjID==iSubj & dataInfo.allData.trialN==iT);
    end
end

mean_line = mean(pHiXtrialN);
sds = std(pHiXtrialN);
min_line = mean_line - sds;
max_line = mean_line + sds; max_line(max_line>1)=1;

x=1:cfg.nTrialsBlk;
px=[x,fliplr(x)]; % make closed patch
py=[min_line, fliplr(max_line)];


figure;
hold on

patch(px,py,1, 'FaceColor', 'b','EdgeColor','none','FaceAlpha', .3);
plot(mean_line,'--ks','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[230 230 250]/255,'MarkerSize',7)
hold off
ylim([.2,1]);

sgtitle(['Learning Curve (SD) - ', dataInfo.outfname], 'Interpreter', 'none')

saveas(gcf, fullfile(dataInfo.figPath, ['LearnCurve_' dataInfo.outfname '.png']))





