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
nsim        = config.nsim;

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
volcols     = {'Stable','Volatile'};
stccols     = {'Small', 'Medium', 'Large'};

 % estimate volatility, stochasticity and learning rates (for mean and
 % error accross simulations)
m_vol       = nan(N,mdl);
m_stc       = nan(N,mdl);
m_lr        = nan(N,mdl);
e_vol       = nan(N,mdl);
e_stc       = nan(N,mdl);
e_lr        = nan(N,mdl);

% these are not needed (lets leave it as it is for now)
allvols     = cell(1,mdl);
allstcs     = cell(1,mdl);
allvals     = cell(1,mdl);
alllrs      = cell(1,mdl);

for j = 1:mdl

    % compute means and SEs over simulations for vol, stc, lr
    m_vol(:,j)  = mean(vols{j},2);
    m_stc(:,j)  = mean(stcs{j},2);
    m_lr(:,j)   = mean(lrs{j},2);
    m_val(:,j)  = mean(vals{j},2);
    e_vol(:,j)  = serr(vols{j},2);
    e_stc(:,j)  = serr(stcs{j},2);
    e_lr(:,j)   = serr(lrs{j},2);  
    e_val(:,j)  = serr(vals{j},2);

    % for each model group split estimated values per stc and vol levels
    for ss = 1:3 % 3 stc levels 
        tmp_vol     = vols{j}(s(:,ss),:); % 
        tmp_stc     = stcs{j}(s(:,ss),:); % 
        tmp_lrs     = lrs{j}(s(:,ss),:); % 
        tmp_vals    = vals{j}(s(:,ss),:); % 

        % now that we have the stc estimates split per voaltility
        for k = 1:2 % stable - volatile
        allvols{1,j}{ss,k} = tmp_vol(t(:,k),:); % split trials based on volatility



        end % end of volatility (k) loop
    end % end of stc loop


end % end of model groups 

end % end of model function 