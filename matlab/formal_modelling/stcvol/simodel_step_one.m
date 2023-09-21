function mout = simodel_step_one(data, probabilities, trials,condtrials, outtype)


% this function simulates data (outcomes) to run particle filtering 
% it computes vol stc, lr and action values on a trial-bytrial basis

% -------------------------
%%

% which model?
mdl         = 1; % if 1 = only healthy model, if 2 = healthy and stc, if 3 healthy, stc and vol models

% unpack data structure
stc_small   = data.stcind(:,1);
stc_medium  = data.stcind(:,2);
stc_large   = data.stcind(:,3);
tvolatile   = data.volind(:,2);
tstable     = data.volind(:,1);
x           = data.x;
o           = data.o;
out         = data.out;

% define parameters for th particle filter
parameters  = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'stc_small', stc_small,'stc_medium',stc_medium,...
    'stc_large',stc_large,'state',x,'rng_id',0,'nsim',2,'model_parameters',parameters);

rng(config.rng_id); 
nsim = config.nsim;

%%

% initialise variables
% for pf
N           = length(o);
outcome     = nan(N,nsim);
vol         = nan(N,nsim);
stc         = nan(N,nsim);
lr          = nan(N,nsim);
val         = nan(N,nsim);    

% init varaibles that will be averaged across pf simulations
outcomes    = cell(1,mdl);
vols        = cell(1,mdl); % for model group
stcs        = cell(1,mdl); % for model group
lrs         = cell(1,mdl); % 
vals        = cell(1,mdl); % vals --> action values? 
v_example   = nan(N,mdl);
glabels     = {'neurotypical','neurodiverse'};

if mdl == 1
    lnames  = {'Healthy'};  
elseif mdl == 2
    lnames  = {'Healthy', sprintf('%s lesion',def('stc'))};      
else
    lnames  = {'Healthy', sprintf('%s lesion',def_actions('stc')), sprintf('%s lesion',def_actions('vol'))};      
end


%% run the particle filter 

% loop over models (healthy, stc, vol) 
for j = 1:mdl

    % loop over simulations 
    for i = 1:nsim
        simdata                                 = ALsimdata_v1(probabilities, trials, condtrials, outtype);
        outcome(:,i)                            = simdata.o; % run pf model fith binary outcomes 
        [vol(:,i),stc(:,i),lr(:,i),val(:,i)]    = model_parfilter(outcome(:,i),config.model_parameters,lnames{j});   % is that the inference model?       

    end % end of simulations loop 

    if j==1
        o_example   = o;
    end
    v_example(:,j)  = val(:,1);
    
    outcomes{j}     = outcome; % I don't think this is needed 
    vols{j}         = vol;
    stcs{j}         = stc;
    lrs{j}          = lr;
    vals{j}         = val;

end % end of model loop

specs(1,:)          = glabels;

%% compute means for all values estimated above

t           = [tstable tvolatile];
s           = [stc_small stc_medium stc_large];
strcols     = {'Stable','Volatile'};

 % estimate volatility, stochasticity and learning rates (for mean and
 % error accross simulations)
m_vol       = nan(N,mdl);
m_stc       = nan(N,mdl);
m_lr        = nan(N,mdl);
e_vol       = nan(N,mdl);
e_stc       = nan(N,mdl);
e_lr        = nan(N,mdl);

% these are not needed (lets leave it as it is for now)
allvols        = cell(N,2);
allstcs        = cell(N,2);
alllrs         = cell(N,4);

for j = 1:mdl

    % compute means and SEs over simulations for vol, stc, lr
    m_vol(:,j)  = mean(sim.vols{j},2);
    m_stc(:,j)  = mean(sim.stcs{j},2);
    m_lr(:,j)   = mean(sim.lrs{j},2);
    e_vol(:,j)  = serr(sim.vols{j},2);
    e_stc(:,j)  = serr(sim.stcs{j},2);
    e_lr(:,j)   = serr(sim.lrs{j},2);    

end % end of model groups 


end % end of model function 