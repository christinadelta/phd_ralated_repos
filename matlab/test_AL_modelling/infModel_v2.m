function output = infModel_v2(data, probabilities, trials,condtrials,beta)

% created: January 2024
% @christinadelta

% MAIN CHARACTERISTICS OF THE MODEL:

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

%% define model parameters and config 

% define model parameters  
% randomly generate values for the free parameterers of the full model:
% 1. lambda_s and lambda_v -- should be at unit ramge (between 0 and 1) 
sim_lambda      = 0.1 + (0.9-0.1).* rand(2,1); 
sim_lambda_s    = sim_lambda(1);
sim_lambda_v    = sim_lambda(2);

% 2. softmax temperatyre beta -- sould be between 0.1 and 10
sim_beta        = 0.1 + (5-0.1).* rand(1,1); 

% 3. s0 and v0 -- should be between 0.1 and 2
sim_init_vals   = 0.1 + (2-0.1).* rand(2,1); 
sim_s0          = sim_init_vals(1);
sim_v0          = sim_init_vals(2);

% define parameters for th particle filter
%parameters  = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
parameters  = struct('nparticles',500,'x0_unc',1,'lambda_v',sim_lambda_v,'lambda_s',sim_lambda_s,'v0',sim_v0,'s0',sim_s0);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'stc_small', stc_small,'stc_medium',stc_medium,...
    'stc_large',stc_large,'state',x,'rng_id',0,'nsim',100,'model_parameters',parameters);

rng(config.rng_id); 
nsim        = config.nsim;

%% prepare variable for the inference model

mdl         = 1; % only core model to run 
% initialise variables
% for pf
N           = length(o);
outcome     = nan(N,nCues,nsim);
outcomeR    = nan(N,nCues,nsim);
vol         = nan(N,nCues,nsim);
stc         = nan(N,nCues,nsim);
lr          = nan(N,nCues,nsim);
val         = nan(N,nCues,nsim);    

% init varaibles that will be averaged across pf simulations
vols        = cell(1,nCues); % 
stcs        = cell(1,nCues); % 
lrs         = cell(1,nCues); % 
vals        = cell(1,nCues); % 
valsR       = cell(1,nCues); % 
v_example   = cell(1,nCues); % 
% glabels     = {'neurotypical','neurodiverse'};
lnames      = {'Healthy'};  

%% run the kalman and particle filters


% loop over simulations 
for i = 1:nsim

    simdata         = ALsimdata_v2(probabilities, trials, condtrials);
    outcome(:,1,i)  = simdata.std_o; % 
    outcome(:,2,i)  = simdata.std_oR;
    o_binary(:,1,i) = simdata.o;
    o_binary(:,2,i) = simdata.oR;

    % loop over cues 
    for cue = 1:nCues
        % run rbpf model 
        [vol(:,i,cue),stc(:,i,cue),lr(:,i,cue),val(:,i,cue)] = model_pf(outcome(:,cue,i),config.model_parameters,lnames{1}); 
    end % end of cues loop  

    % run response model to simulate actions 
    Qvals(:,1)              = val(:,i,1);
    Qvals(:,2)              = val(:,i,2);
    [action(:,i),r(:,i)]    = respModel(Qvals, sim_beta,x);

end % end of simulations loop

% extract values for option A and option B
for j = 1:nCues
    vols{j}     = vol(:,:,j);
    stcs{j}     = stc(:,:,j);
    lrs{j}      = lr(:,:,j);
    vals{j}     = val(:,:,j);
end

% for testing model fitting 
test_a          = action(:,1);
test_o(:,1)     = outcome(:,1,1);
test_o(:,2)     = outcome(:,2,1);

%% compute mean of estimated values over simulations 

t           = [tstable tvolatile];
s           = [stc_small stc_medium stc_large];
volcols     = {'Stable','Volatile'};
stccols     = {'Small', 'Medium', 'Large'};

% estimate volatility, stochasticity and learning rates (means and
% error of means accross simulations)
% loop over cues/options
for j = 1:nCues

    % compute means and SEs over simulations for vol, stc, lr
    m_vol(:,j)  = mean(vols{j},2);
    m_stc(:,j)  = mean(stcs{j},2);
    m_lr(:,j)   = mean(lrs{j},2);
    m_val(:,j)  = mean(vals{j},2);
    e_vol(:,j)  = serr(vols{j},2);
    e_stc(:,j)  = serr(stcs{j},2);
    e_lr(:,j)   = serr(lrs{j},2);  
    e_val(:,j)  = serr(vals{j},2);

    % for each cue split estimated values per stc and vol levels (
    % this could be used for visualization???)
    for ss = 1:3 % 3 stc levels 
        tmp_vol     = vols{j}(s(:,ss),:); % 
        tmp_stc     = stcs{j}(s(:,ss),:); % 
        tmp_lrs     = lrs{j}(s(:,ss),:); % 
        tmp_vals    = vals{j}(s(:,ss),:); % 

        % and the averaged and serr values
        tmp_mvol     = m_vol(s(:,ss),j);
        tmp_mstc     = m_stc(s(:,ss),j);
        tmp_mlr      = m_lr(s(:,ss),j);
        tmp_mval     = m_val(s(:,ss),j);

        tmp_evol = e_vol(s(:,ss),j);
        tmp_estc = e_stc(s(:,ss),j);
        tmp_elr = e_lr(s(:,ss),j);
        tmp_eval = e_val(s(:,ss),j);

        % now that we have the stc estimates split per voaltility
        for k = 1:2 % stable - volatile

            allvols{1,j}{ss,k}      = tmp_vol(t(:,k),:); % split trials based on volatility
            allstcs{1,j}{ss,k}      = tmp_stc(t(:,k),:); % 
            allvals{1,j}{ss,k}      = tmp_vals(t(:,k),:); % split trials based on volatility
            alllrs{1,j}{ss,k}       = tmp_lrs(t(:,k),:); % split trials based on volatility

            all_mvols{1,j}{ss,k}    = tmp_mvol(t(:,k),:);
            all_mstc{1,j}{ss,k}     = tmp_mstc(t(:,k),:);
            all_mlr{1,j}{ss,k}      = tmp_mlr(t(:,k),:);
            all_mval{1,j}{ss,k}     = tmp_mval(t(:,k),:);

            all_evols{1,j}{ss,k}    = tmp_evol(t(:,k),:);
            all_estc{1,j}{ss,k}     = tmp_estc(t(:,k),:);
            all_elr{1,j}{ss,k}      = tmp_elr(t(:,k),:);
            all_evals{1,j}{ss,k}    = tmp_eval(t(:,k),:);

        end % end of volatility (k) loop
    end % end of stc loop
end % end of options loop

%% split actions per stc and volatility

for stc = 1:3

    tmp_actions     = action(s(:,stc),:);
    tmp_outcomes    = o_binary(s(:,stc),:,:);

    for vol = 1:2
        allactions{stc,vol}     = tmp_actions(t(:,vol),:);
        allbinary_o{stc,vol}    = tmp_outcomes(t(:,vol),:,:);

    end % end of volatility condition
end % end of stc level

%% store model results to output structure

sim_data    = struct('allvols',allvols,'allstcs',allstcs, 'allvals',allvals,'alllrs',alllrs);
mean_data   = struct('mvols',all_mvols,'mstcs',all_mstc,'mlr',all_mlr,'mvals',all_mval);
se_data     = struct('evols',all_evols,'estcs',all_estc,'elr',all_elr,'evals',all_evals);

output.simdata      = sim_data;
output.mdata        = mean_data;
output.sedata       = se_data;
output.outcomes     = outcome;
output.actions      = allactions;
output.binary_o     = allbinary_o;
output.rewrad       = r;
output.test_a       = test_a;
output.test_o       = test_o;

end % end of function