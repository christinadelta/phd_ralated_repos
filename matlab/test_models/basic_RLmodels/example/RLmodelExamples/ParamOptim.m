% This function finds the best fitting model/parameters
% considering priors (in paramPriors)
% thus obtaining the log of the posterior probability of the model

%%%% NS 20201117/20210204 - Optimism & confirmation bias Y3 project:
%%%%%% free choices & full feedback


%% Models
%  1 - classic RL: 1 beta, 1 alpha
%  2 - lr by outcome type X valence: 1 beta, 4 alphas
%  3 - lr by confirmation bias: 1 beta, 2 alphas
%  4 - lr by optimism bias: 1 beta, 2 alphas


%% Population definition

clearvars;
close all

addpath('./FitFuns')

dataInfo.type       = 'real';
dataInfo.dataPath   = '../data/';
dataInfo.outPath    = './OptimResults';

dataInfo.infname    = 'taskData_NoErr.csv' ;
dataInfo.subj2excl  = [2664319, 2700344, 2700384]; % Ss with poor learning
modelInfo.Qinit     = zeros(1,2);    % middle of outcome range [-1,1]



dataInfo.exclProbSs = 0; % logical, exclude problematic subjects before MAP?
if dataInfo.exclProbSs == 1
    dataInfo.outfnameSs = 'goodSs';
else
    dataInfo.outfnameSs = 'allSs';
end

dataInfo.exclErrs   = 1; % logical, ignore errors? - if not, should record - normally already removed
dataInfo.inDataVars  = {'subj', 'blockN', 'trialNblk_fin', 'stimLocRev', 'hiResp', 'out_f', 'out_cf'};
dataInfo.modelDataVars = {'subj', 'blockN', 'trialN', 'stimLocRev', 'action', 'outcome_f', 'outcome_cf'}; % rename vars & convert hiResp[1,0] to action[1,2]


%% Set up modelInfo vars

% which models to run?
modelInfo.whichms = 1:4;
modelInfo.outfnameMds = 'm1-4'; % sprintf('m%02d', modelInfo.whichms);

% Include priors? 
modelInfo.priors  = 1; % logical
modelInfo.priorsType = 'apriori';

if modelInfo.priors
    if strcmp(modelInfo.priorsType, 'apriori')
        dataInfo.outfname = 'ParamOptim_aprioriPriors';
    else
        error('Prior Type not specified');
    end
else
    dataInfo.outfname = 'ParamOptim_noPriors';
end


% filename for output data
switch dataInfo.type
    case 'real'
        dataInfo.outfnameXtra = [dataInfo.outfnameSs '_' modelInfo.outfnameMds];
    otherwise        
        dataInfo.outfnameXtra = [dataInfo.fakeDataName '_' modelInfo.outfnameMds];
end


modelInfo.nStartPoints = 1; % if 1, use what's in optStart, else randsample
% modelInfo.nStartPoints = 3; % if 1, use what's in optStart, else randsample
modelInfo.maxsearch = 5;

% standard, untransformed
modelInfo.transform = 0;
modelInfo.priorDists.lpmf = 0;


%  1 - classic RL: 1 beta, 1 alpha
modelInfo.params{1}.labels      = {'Beta', 'Alpha'};
modelInfo.params{1}.optStart    = [1   .5]; % plausible, but will reset in multstart
modelInfo.params{1}.optMin      = [0   0];
modelInfo.params{1}.optMax      = [Inf  1];


%  2 - lr by outcome type X valence: 1 beta, 4 alphas
modelInfo.params{2}.labels      = {'Beta', 'Alpha_Factual_Positive', 'Alpha_Factual_Negative',...
                                    'Alpha_Counterfactual_Positive', 'Alpha_Counterfactual_Negative'};
modelInfo.params{2}.optStart    = [1   repmat(.5, 1, 4)]; % plausible, but will reset in multstart
modelInfo.params{2}.optMin      = [0   zeros(1, 4)];
modelInfo.params{2}.optMax      = [Inf ones(1, 4)];


%  3 - lr by confirmation bias: 1 beta, 2 alphas
modelInfo.params{3}.labels      = {'Beta', 'Alpha_Confirmatory', 'Alpha_Disconfirmatory'};
modelInfo.params{3}.optStart    = [1   repmat(.5, 1, 2)]; % plausible, but will reset in multstart
modelInfo.params{3}.optMin      = [0   zeros(1, 2)];
modelInfo.params{3}.optMax      = [Inf ones(1, 2)];


%  4 - lr by confirmation bias: 1 beta, 2 alphas
modelInfo.params{4} = modelInfo.params{3};
modelInfo.params{4}.labels = {'Beta', 'Alpha_Positive', 'Alpha_Negative'};



%% Load Data

allData = readtable(fullfile(dataInfo.dataPath, dataInfo.infname));
allData = allData(:,dataInfo.inDataVars);
% convert vars for modelling
allData.Properties.VariableNames{3} = 'trialN';
allData.Properties.VariableNames(end-1:end) = {'outcome_f', 'outcome_cf'};
allData.action = nan(height(allData),1);
allData.action(allData.hiResp==1) = 1; % 1 is best option, 2 is worse
allData.action(allData.hiResp==0) = 2;
allData=allData(:,dataInfo.modelDataVars);

if dataInfo.exclProbSs 
    allData = allData(~ismember(allData.subj, dataInfo.subj2excl), :);
end

dataInfo.allData = allData;
dataInfo.subjIDs = unique(allData.subj);
dataInfo.subjN   = numel(dataInfo.subjIDs);

%% Fit models

% % Save command window output in a "diary" file:
diary(fullfile(dataInfo.outPath, sprintf('%s_%s_Diary_%s.txt', dataInfo.outfname, dataInfo.outfnameSs, datestr(now,'yyyymmdd-HHMM'))));
fprintf('\n   %s, run on %s\n\n', dataInfo.outfname, datestr(now,'yyyymmdd-HHMM'));

%%
modelInfo.optimOptions = optimset('Algorithm', 'interior-point', 'Display', 'notify', 'MaxIter', 10000);

fullRes = struct;

for model = modelInfo.whichms

    fprintf('\nFitting model %d\n\n', model);

    for iSubj = 1:dataInfo.subjN    

        subjData = allData(allData.subj == dataInfo.subjIDs(iSubj), :);
    

        tmp_loglik = []; tmp_report = tmp_loglik; 
        tmp_parameters = cell(1, modelInfo.nStartPoints); tmp_gradient = tmp_parameters; tmp_hessian = tmp_parameters;
                
        for run = 1:modelInfo.nStartPoints

            if modelInfo.nStartPoints > 1
                start_points = [];
                for iP = 1:length(modelInfo.params{model}.optStart)
                    if modelInfo.params{model}.optMax(iP) > modelInfo.maxsearch
                        start_points(iP) = modelInfo.params{model}.optMin(iP) + (modelInfo.maxsearch-modelInfo.params{model}.optMin(iP))*rand;
                    else
                        start_points(iP) = rand;
                    end                    
                end                
            else
                start_points = modelInfo.params{model}.optStart;
            end
            fullRes.start_points{iSubj,model}(run,:) = start_points;
                 
            if modelInfo.priors
            
                [tmp_parameters{run},...
                 tmp_loglik(run),...  % if priors, it's the log of the posterior probability
                 tmp_report(run),...  % exit flags <= 0 imply problems
                 tmp_gradient{run},...
                 tmp_hessian{run}] = ...
                        fmincon(@(x) aprioriPriors(x, model, modelInfo.Qinit, subjData, modelInfo),...
                            start_points,... 
                            [],[],[],[],...
                            modelInfo.params{model}.optMin,... 
                            modelInfo.params{model}.optMax,...
                            [], modelInfo.optimOptions);      
                                       
            else % no priors, obtain negative log likelihood, then calculate bic

                [tmp_parameters{run},...
                 tmp_loglik(run),...  % if no priors, it's the negative log likelihood
                 tmp_report(run),...  % exit flags <= 0 imply problems
                 tmp_gradient{run},...
                 tmp_hessian{run}] = ...
                        fmincon(@(x) fit_compModels(x, model, modelInfo.Qinit, subjData, modelInfo),...
                            start_points,... 
                            [],[],[],[],...
                            modelInfo.params{model}.optMin,... 
                            modelInfo.params{model}.optMax,...
                            [], modelInfo.optimOptions);                         
                        
            end % if modelInfo.priors
            
        fullRes.parameters{iSubj,model}(run,:) = tmp_parameters{run};
        fullRes.loglik{iSubj,model}(run) = tmp_loglik(run);
        fullRes.report{iSubj,model}(run) = tmp_report(run);
        
        end % for run = 1:modelInfo.nStartPoints
                        
        [~, winRun] = min(tmp_loglik);
        
        parameters{model}(iSubj,:)  = tmp_parameters{winRun};
        loglik(iSubj,model)         = tmp_loglik(winRun); 
        report(iSubj,model)         = tmp_report(winRun); 
        gradient{iSubj,model}       = tmp_gradient{winRun}; 
        hessian{iSubj,model}        = tmp_hessian{winRun}; 
        
        fullRes.winRun(iSubj,model) = winRun;
    end
end
clear iSubj model subjData tmp_parameters tmp_loglik tmp_report tmp_gradient tmp_hessian winRun run start_points iP

fprintf('\n\n*******    Model Optim Finished!     *******\n');


%% get loglikelihood & ICs
if modelInfo.priors && strcmp(modelInfo.priorsType, 'apriori')
    [bic, aic, ll] = get_ll_fromapriori(parameters, modelInfo, dataInfo);
end

%%
save(fullfile(dataInfo.outPath, [dataInfo.outfname, '_', dataInfo.outfnameXtra]));
fprintf('\n\n*******    Results Saved!     *******\n');



%% Close the diary
diary off;  


%% save model summaries csvs 
if modelInfo.priors && strcmp(modelInfo.priorsType, 'apriori')
    sumPath    = fullfile(dataInfo.outPath, 'summaries');
    sum_outname = 'aprioriPriors_%s_sum_m%d.csv';

    for model = modelInfo.whichms

       temp = array2table([parameters{model}, loglik(:, model), report(:,model),...
                    dataInfo.subjIDs, ll(:, model), bic(:, model), aic(:, model)],...
                        'VariableNames', horzcat(modelInfo.params{model}.labels, 'logposterior', 'exitflag', 'subj', 'likelihood','bic','aic'));
       writetable(temp, fullfile(sumPath,  sprintf(sum_outname, dataInfo.outfnameSs, model)));

    end
    clear temp sum_outname
    fprintf('\n\n*******    Model summary tables saved!     *******\n');    
end