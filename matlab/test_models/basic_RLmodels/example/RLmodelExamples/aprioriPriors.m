function [post] = aprioriPriors(params, model, Qinit, subjData, modelInfo)
%% Combines priors with the posterior likelihood of the estimated parameters 
% Input: parameters, model info, and
% % data for 1 subject (1 trial per entry in important vars (choice, action, outcome, Qrestart))
% Output: posterior likelihood

%%%% NS 20201117 - Optimism & confirmation bias Y3 project:
%%%%%% free choices & full feedback


%% Models

%  1 - classic RL: 1 beta, 1 alpha
%  2 - lr by outcome type X valence: 1 beta, 4 alphas
%  3 - lr by confirmation bias: 1 beta, 2 alphas
%  4 - lr by optimism bias: 1 beta, 2 alphas

%% log prior of parameters

if model == 1
    beta1   = params(1);
    alpha1  = params(2);

    pBeta1  = log(gampdf(beta1,   1.2 , 5.0));
    pAlpha1 = log(betapdf(alpha1, 1.1 , 1.1));

    prior = [pBeta1  pAlpha1];       
            
elseif model == 2 % 1 beta, 2 alphas
    beta1   = params(1);
    alpha1  = params(2);
    alpha2  = params(3);
    alpha3  = params(4);
    alpha4  = params(5);

    pBeta1  = log(gampdf(beta1,   1.2 , 5.0));
    pAlpha1 = log(betapdf(alpha1, 1.1 , 1.1));
    pAlpha2 = log(betapdf(alpha2, 1.1 , 1.1));
    pAlpha3 = log(betapdf(alpha3, 1.1 , 1.1));
    pAlpha4 = log(betapdf(alpha4, 1.1 , 1.1));    

    prior = [pBeta1  pAlpha1  pAlpha2  pAlpha3  pAlpha4];       
        
    
elseif ismember(model, 3:4) % 1 beta, 2 alphas
        
    beta1   = params(1);
    alpha1  = params(2);
    alpha2  = params(3);

    pBeta1  = log(gampdf(beta1,   1.2 , 5.0));
    pAlpha1 = log(betapdf(alpha1, 1.1 , 1.1));
    pAlpha2 = log(betapdf(alpha2, 1.1 , 1.1));

    prior = [pBeta1  pAlpha1  pAlpha2];       
    

else
    error( 'Error!!  in paramPriors  --  No Prior settings for model %d!!\n', model)
end

prior = -sum(prior);

%% Estimate model

lik = fit_compModels(params, model, Qinit, subjData, modelInfo);

%% Combine prior & likelihood
post = prior + lik;

