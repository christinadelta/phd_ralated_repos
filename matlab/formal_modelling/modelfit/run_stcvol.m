function modelout = run_stcvol(parameters,subj1)


% this function runs particle filters for healthy model only 
% comments will go here 
% and here 

%% define initial variables to simulate some data

probabilities   = [.90 .10;                 % small stochasticity probabilities
    .80 .20;                                % medium stochasticity
    .70 .30];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 140;                      % total trials per stc condition
condtrials      = {70,[30,10,10,20]};
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes
nCues           = 2;

% simulate data ( we just want stochasticity and volatility indexes for
% now)
data            = ALsimdata_v2(probabilities, trials, condtrials, outtype);

%% unpack data structure 

stc_small   = data.stcind(:,1);
stc_medium  = data.stcind(:,2);
stc_large   = data.stcind(:,3);
tvolatile   = data.t(:,2);      % volatilitile trials indexed with 1
tstable     = data.t(:,1);      % stable trials indexed with 1
x           = data.x;           % true loss rate
xR          = data.xR;

% unpack subject structure 
actions     = subj1.actions;
outcome     = subj1.outcome;
totalTrials = length(outcome);
mdl         = 1;                % how many models to run? if 1 = only healthy model 

% outcome variance
outVar      = .01;


%% configure params for particle filter

config      = struct('tvolatile',tvolatile,'tstable',tstable,'stc_small', stc_small,'stc_medium',stc_medium,...
    'stc_large',stc_large,'state',x,'rng_id',0,'nsim',1, 'model_parameters', parameters);

nsim        = config.nsim;
rng(config.rng_id); 

%% initialise variable to store estimated parameters 

% for pf
N           = length(outcome);
vol         = nan(N,nsim);
stc         = nan(N,nsim);
lr          = nan(N,nsim);
val         = nan(N,nsim);    

lnames      = {'Healthy'};  

%% 
out                         = x(:,1) + sqrt(outVar)*randn(totalTrials,1); % outcomes for option A(generated based on reward rates and stochasticity??)
outR                        = xR(:,1) + sqrt(outVar)*randn(totalTrials,1);

%% normalise free parameters

% % transform lambdas, V0 and S0 to bebetween 0 and 1 
% tmp_lambdav = parameters(1);
% lambda_v    = 1/(1+exp(-tmp_lambdav));
% tmp_lambdas = parameters(2);
% lambda_s    = 1/(1+exp(-tmp_lambdas));
% tmp_v0      = parameters(3);
% v0          = 1/(1+exp(-tmp_v0));
% tmp_s0      = parameters(4);
% s0          = 1/(1+exp(-tmp_s0));
% 
% % transform beta 
% tmp_beta    = parameters(5);
% beta        = exp(tmp_beta);

%% add optimised paramteres to the configuration file 

% mparams                 = [lambda_v lambda_s v0 s0];
% config.model_parameters = mparams;

%% run particle filter with subject choices? 
for i = 1:nsim
    for cue = 1:nCues
        if cue == 1
            
            [vol(:,i),stc(:,i),lr(:,i),val(:,i,cue)]    = model_parfilter(out(:,i),config.model_parameters,lnames{1});  

        else
            [vol(:,i),stc(:,i),lr(:,i),val(:,i,cue)]    = model_parfilter(outR(:,i),config.model_parameters,lnames{1}); 

        end

        
    end
end % end of simulations loop 

vals        = val(:,:,1); % expected values for option 1
valsR       = val(:,:,2); % expected values for option 2

%% compute means for all estimated values above 

t           = [tstable tvolatile];
s           = [stc_small stc_medium stc_large];
volcols     = {'Stable','Volatile'};
stccols     = {'Small', 'Medium', 'Large'};

% these are not needed (lets leave it as it is for now)
allvols     = cell(1,mdl);
allstcs     = cell(1,mdl);
allvals     = cell(1,mdl);
alllrs      = cell(1,mdl);

% loop over stc levels
for ss = 1:3

    tmp_vol             = vol(s(:,ss),:);   % 
    tmp_stc             = stc(s(:,ss),:);   % 
    tmp_lrs             = lr(s(:,ss),:);    % 
    tmp_vals            = vals(s(:,ss),:);  % 

    for k = 1:2 % stable - volatile

        allvols{ss,k}  = tmp_vol(t(:,k),:);     % split trials based on volatility
        allstcs{ss,k}  = tmp_stc(t(:,k),:);     % 
        allvals{ss,k}  = tmp_vals(t(:,k),:);    % split trials based on volatility
        alllrs{ss,k}   = tmp_lrs(t(:,k),:);     % split trials based on volatility

    end % end of volatility loop
end % end of stc loop

%% compute ll using response model 

beta                    = parameters.beta; % last parameter
xpvals                  = [vals valsR];
ll                      = 0; 

% convert action values to choice probablities and compute ll
for t = 1:totalTrials

    p                   = exp(beta*xpvals(t,:)) / sum(exp(beta*xpvals(t,:)));   % compute choice probabilities using the softmax function
    choiceProbs(t,:)    = p;
    a                   = actions(t); % get choice of that trial
    o                   = outcome(t); % get outcome of that trial

    % store probability of the chosen action
    if a == 1
        choiceProb(t)   = p(1);
    elseif a == 2
        choiceProb(t)   = p(2);

    end

   % ll = ll - log(choiceProb(t)+eps);

end % end of trials loop

% update ll
ll              = sum(log(choiceProb+eps));

%% update model output

modelout.ll     = ll;
modelout.cb     = choiceProbs;
modelout.vols   = allvols;
modelout.vals   = allvals;
modelout.lrs    = alllrs;
modelout.stc    = allstcs;
modelout.t      = t; % volatile/stable trials 
modelout.s      = s; % all stc trials indexes 


end % end of function