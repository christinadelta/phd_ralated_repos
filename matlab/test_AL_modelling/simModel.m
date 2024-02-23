function output = simModel(data, probabilities, trials,condtrials,parameters)

% created: Januray 2024
% @christinadelta

% MAIN CHARACTERISTICS OF THE MODEL:

% Part of Parameter Recovery Process. 

% For inferences about volatility and stochasticity and about the true
% reward/loss rate we are using the stc-vol model that jointly estimates
% stochasticity and volatility using the kalman (to estimate x) and
% particle (to estimate stc and vol) filters.

% Mathematical Formulation: The model uses specific applications of 
% Bayes' rule, representing beliefs about the reward rate at each step as 
% a Gaussian distribution with a mean and variance. This formulation leads 
% to a rationale for the error-driven update (a prominent concept in 
% neuroscience and psychology) and accounts for the learning rate, 
% which depends on the agent's uncertainty
% 
% Kalman Filter Foundation: The model begins with a Kalman filter, which 
% is optimal for learning from noisy data. It assumes that inferences are 
% drawn from observations corrupted by two types of noise: process noise 
% (volatility) and outcome noise (stochasticity). Volatility refers to 
% how quickly the true value being estimated changes, while stochasticity 
% describes additional noise in the observation of each outcome
% 
% Probabilistic Model Development: To handle the complexities of 
% simultaneous estimation of volatility and stochasticity, the model uses 
% particle filtering, a standard approximation approach. This method 
% allows for tractable inference for the reward rate, replacing the true 
% values of stochasticity and volatility with their corresponding samples
% 
% Learning Rate Dependence: The model underscores that volatility 
% increases the learning rate while stochasticity reduces it. Learning 
% both parameters simultaneously is critical for efficient learning due 
% to their opposite effects on the learning rate​
% 
% --------------

%% prepare data 

% unpack data structure
stc_small   = data.stcind(:,1);
stc_medium  = data.stcind(:,2);
stc_large   = data.stcind(:,3);
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);
x           = data.x;
xR          = data.xR;
nCues       = data.cues;
o           = data.std_o;
oR          = data.std_oR;
beta        = parameters.beta;

% define configuration 
config          = struct('tvolatile',tvolatile,'tstable',tstable,'stc_small', stc_small,'stc_medium',stc_medium,...
    'stc_large',stc_large,'state',x,'rng_id',0,'nsim',1,'model_parameters',parameters);

rng(config.rng_id); 
nsim        = config.nsim;

%% prepare variable for the inference model

mdl         = 1; % only core model to run 
% initialise variables
% for pf
N           = length(o);
outcome     = nan(N,nsim);
outcomeR    = nan(N,nsim);
vol         = nan(N,nsim);
stc         = nan(N,nsim);
lr          = nan(N,nsim);
val         = nan(N,nsim);  
lnames      = {'Healthy'};  

%% run the kalman and particle filters

for sim = 1:nsim

    simdata             = ALsimdata_v3(probabilities, trials, condtrials);
    outcome(:,1,sim)    = simdata.std_o; % 
    outcome(:,2,sim)    = simdata.std_oR;
    o_binary(:,1,sim)   = simdata.o;
    o_binary(:,2,sim)   = simdata.oR;

    % loop over cues 
    for cue = 1:nCues
        % run rbpf model 
        [vol(:,cue,sim),stc(:,cue,sim),lr(:,cue,sim),val(:,cue,sim)] = model_pf(outcome(:,cue,sim) ,config.model_parameters,lnames{1});

    end % end of cues loop  

    % run response model to simulate actions 
    Qvals(:,1)                  = val(:,1,sim);
    Qvals(:,2)                  = val(:,2,sim);
    [action(:,sim),r(:,sim)]    = respModel(Qvals,beta,x);
    
end % end of simulations loop

% for testing model fitting 
test_a          = action(:,1);
test_o(:,1)     = outcome(:,1);
test_o(:,2)     = outcome(:,2);


%% store model results to output structure

output.outcomes         = outcome;
output.outcomesR        = outcomeR;
output.vvals            = Qvals;
output.volatility       = vol;
output.stochasticity    = stc;
output.binary_o         = o_binary;
output.rewrad           = r;
output.lrs              = lr;
output.test_a           = test_a;
output.test_o           = test_o;

end % end of fucntion