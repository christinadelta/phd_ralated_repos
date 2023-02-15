function [bic, aic, ll] = get_ll_fromapriori(parameters, modelInfo, dataInfo)
% Fits model with optimised parameters to obtain the corresponding
% negative log likelihood, rather than the log posterior output 

bic = zeros(dataInfo.subjN, numel(modelInfo.whichms));
aic = bic;
ll  = bic;

for model = modelInfo.whichms

    m_results=parameters{model};
    
    likelihood=zeros(dataInfo.subjN,1);
    
    for iSubj=1:dataInfo.subjN
        subjData = dataInfo.allData(dataInfo.allData.subj == dataInfo.subjIDs(iSubj), :);
        
        subjParams = m_results(iSubj, 1:numel(modelInfo.params{model}.labels));
        likelihood(iSubj) = fit_compModels(subjParams, model, modelInfo.Qinit, subjData, modelInfo);
        
        %% Add Information Criteria
        logLik = abs(likelihood(iSubj)); % this ensures values are positive, so in formula multiply by 2 to obtain positive BIC/AICs, where smaller is better(normal -2 is assuming logLik is negative)
        nTrials = height(subjData);
        df = numel(modelInfo.params{model}.labels);
        bic(iSubj, model) = 2*logLik + df*log(nTrials); 
        aic(iSubj, model) = 2*logLik + 2*df;    
        ll(iSubj, model)  = logLik;
                
    end        
end

