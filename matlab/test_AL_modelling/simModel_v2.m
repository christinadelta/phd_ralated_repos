function output = simModel_v2(data, probabilities, trials,condtrials,parameters,i,j)

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
lnames      = {'Healthy'};  

t           = [tstable tvolatile];
s           = [stc_small stc_medium stc_large];

%% run the kalman and particle filters

% loop over simulations 
for sim = 1:nsim

    simdata             = ALsimdata_v2(probabilities, trials, condtrials);
    outcome(:,1,sim)    = simdata.std_o; % 
    outcome(:,2,sim)    = simdata.std_oR;
    o_binary(:,1,sim)   = simdata.o;
    o_binary(:,2,sim)   = simdata.oR;

    % extract only ith and jth outcomes and state
    stc_o = outcome(s(:,i),:,:); vol_o = stc_o(t(:,j),:,:);
    stc_x = x(s(:,i),:,:); vol_x = stc_x(t(:,j),:,:);
    stc_binary_o = o_binary(s(:,i),:,:); vol_binary_o = stc_binary_o(t(:,j),:,:);

    % loop over cues 
    for cue = 1:nCues
        % run rbpf model 
        [vol(:,sim,cue),stc(:,sim,cue),lr(:,sim,cue),val(:,sim,cue)] = model_pf(vol_o(:,cue,sim),config.model_parameters,lnames{1});

        all_outcomes(:,cue,sim) = vol_o(:,cue,sim);
        all_binary_outcomes(:,cue,sim) =vol_binary_o(:,cue,sim);

    end % end of cues loop  

    % run response model to simulate actions 
    Qvals(:,1)                  = val(:,sim,1);
    Qvals(:,2)                  = val(:,sim,2);
    [action(:,sim),r(:,sim)]    = respModel(Qvals,beta,vol_x);
    
end % end of simulations loop

%% store output 
output.o        = all_outcomes;
output.binary_o = all_binary_outcomes;
output.vals     = Qvals;
output.vol      = vol;
output.stc      = stc;
output.lr       = lr;
output.a        = action;
output.reward   = r;
output.x        = vol_x;

end % end of fucntion