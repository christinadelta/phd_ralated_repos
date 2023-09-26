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
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);
x           = data.x;
xR          = data.xR;
o           = data.o;
oR          = data.oR;
out         = data.out;
nCues       = data.cues;

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
outcomes    = cell(1,mdl); % blue feedback 
outcomesR   = cell(1,mdl); % red feedback
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
        for cue = 1:nCues
            simdata                                                 = ALsimdata_v1(probabilities, trials, condtrials, outtype);
            outcome(:,i)                                            = simdata.o; % run pf model fith binary outcomes 
            outcomeR(:,i)                                           = simdata.oR; % we also want to estimate vals for the red feedback outcomes
            if cue == 1
                [vol(:,i,cue),stc(:,i,cue),lr(:,i,cue),val(:,i,cue)]    = model_parfilter(outcome(:,i),config.model_parameters,lnames{j});   
            else
                [vol(:,i,cue),stc(:,i,cue),lr(:,i,cue),val(:,i,cue)]    = model_parfilter(outcomeR(:,i),config.model_parameters,lnames{j}); 
            end
        end
    end % end of simulations loop 

    if j==1
        o_example   = o;
    end
    v_example(:,j)  = val(:,1);
    
    outcomes{j}     = outcome; % I don't think this is needed 
    outcomesR{j}    = outcomeR;
    vols{j}         = vol(:,:,1);
    stcs{j}         = stc(:,:,1);
    lrs{j}          = lr(:,:,1);
    vals{j}         = val(:,:,1);
    valsR{j}        = val(:,:,2);

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

    % for each model group split estimated values per stc and vol levels (
    % this could be used for visualization???
    for ss = 1:3 % 3 stc levels 
        tmp_vol             = vols{j}(s(:,ss),:); % 
        tmp_stc             = stcs{j}(s(:,ss),:); % 
        tmp_lrs             = lrs{j}(s(:,ss),:); % 
        tmp_vals            = vals{j}(s(:,ss),:); % 

        % now that we have the stc estimates split per voaltility
        for k = 1:2 % stable - volatile

            allvols{1,j}{ss,k}  = tmp_vol(t(:,k),:); % split trials based on volatility
            allstcs{1,j}{ss,k}  = tmp_stc(t(:,k),:); % 
            allvals{1,j}{ss,k}  = tmp_vals(t(:,k),:); % split trials based on volatility
            alllrs{1,j}{ss,k}   = tmp_lrs(t(:,k),:); % split trials based on volatility

        end % end of volatility (k) loop
    end % end of stc loop
end % end of model groups 

%% compute mean lrs for each stc level
% compute averaged learning rates for each of the 6 conditions for each
% model group (heathy, stc, etc..)
% loop over group
for g = 1:mdl

    group_lrs = alllrs{1,g};

    for i = 1:3 % loop over stochasticities
        for ii = 1:2 % loop over volatilies 
            stc_lrs                 = group_lrs{i,ii};
            tmp_m                   = mean(stc_lrs,2);  % average over simulations
            mean_alphas{g}(i,ii)    = mean(tmp_m,1);    % average over trials
            mstc_lrs{ii}            = tmp_m;

            clear stc_lrs tmp_m
        end % end of voalitlity loop

        mean_stc_lrs(:,i) = [mstc_lrs{1}; mstc_lrs{2}];

    end % end of stcs loop
end % end of model groups loop

%% compute mean vols for each stc level

% loop over group
for g = 1:mdl

    group_vols = allvols{1,g};

    for i = 1:3 % loop over stochasticities
        for ii = 1:2 % loop over volatilies 
            stc_vols            = group_vols{i,ii};
            tmp_m               = mean(stc_vols,2);  % average over simulations
            mean_vols{g}(i,ii)  = mean(tmp_m,1);    % average over trials
            mstc_vol{ii}        = tmp_m;

            clear stc_vols tmp_m
        end % end of voalitlity loop
        mean_stc_vols(:,i) = [mstc_vol{1}; mstc_vol{2}];
    end % end of stcs loop
end % end of model groups loop

%% compute mean stc for each stc level
% loop over group
for g = 1:mdl

    group_stcs                  = allstcs{1,g};

    for i = 1:3 % loop over stochasticities
        for ii = 1:2 % loop over volatilies 
            stc_stcs            = group_stcs{i,ii};
            tmp_m               = mean(stc_stcs,2);  % average over simulations
            mean_stcs{g}(i,ii)  = mean(tmp_m,1);    % average over trials
            mstc_stc{ii}        = tmp_m;
            
            clear stc_stcs tmp_m
        end % end of voalitlity loop
        mean_stc_stcs(:,i) = [mstc_stc{1}; mstc_stc{2}];
    end % end of stcs loop
end % end of model groups loop

%% compute mean expected values for each stc level

% loop over group
for g = 1:mdl

    group_vals                  = allvals{1,g};

    for i = 1:3 % loop over stochasticities
        for ii = 1:2 % loop over volatilies 
            stc_vals            = group_vals{i,ii};
            tmp_m               = mean(stc_vals,2);  % average over simulations
            mean_vals{g}(i,ii)  = mean(tmp_m,1);    % average over trials
            mstc_val{ii}        = tmp_m;

            clear stc_vals tmp_m
        end % end of voalitlity loop
        mean_stc_vals(:,i) = [mstc_val{1}; mstc_val{2}];
    end % end of stcs loop
end % end of model groups loop

%% store all needed data to output struct

sim_data                = struct('vols',vols, 'vals', vals,'valsR',valsR, 'stcs', stcs,'lrs', lrs, 'o', outcomes, 'oR', outcomesR);


mean_data.m_vol         = m_vol;
mean_data.m_val         = m_val;
mean_data.m_stc         = m_stc;
mean_data.m_lr          = m_lr;
mean_data.mean_vol      = mean_vols;
mean_data.mean_val      = mean_vals;
mean_data.mean_stc      = mean_stcs;
mean_data.mean_alphas   = mean_alphas;
mean_data.mstc_vols     = mean_stc_vols;
mean_data.mstc_vals     = mean_stc_vals;
mean_data.mstc_stcs     = mean_stc_stcs;
mean_data.mstc_alphas   = mean_stc_lrs;


mout.sim_data           = sim_data;
mout.mean_data          = mean_data;

end % end of model function 