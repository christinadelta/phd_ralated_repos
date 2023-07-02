function [stochVolSim, vals] = runModel(cue, config, model, condition, probabilities, trials,condtrials, outpath, outtype)

% runs stochasticity-volatility model (main, lesioned)
% Created in June 2023

%% unpack data , params and config structures

rng(config.rng_id); 
nsim        = config.nsim;

% init params to be estimated 
N           = length(config.state);
outcome     = nan(N,nsim);
vol         = nan(N,nsim);
stc         = nan(N,nsim);
lr          = nan(N,nsim);
val         = nan(N,nsim);  

% which model is it?
if model == 1 % run the simulations using the Piray 2021 parameter values

    g           = 2; % 2 groups --> healthy, ASD
    glabels     = {'Control','ASD'};
    lnames      = {'Healthy', sprintf('%s lesion',def_actions('stc'))};      
elseif model == 2
    g = 1;
    glabels     = {'Control'};
    lnames      = {'Healthy'};  

elseif model == 3

    g = 3;
    glabels     = {'ASD'};
    lnames      = {sprintf('%s lesion',def_actions('stc'))};

end

% init cells to store the above parameters
outcomes    = cell(1,g);
vols        = cell(1,g); % stable, volatile
stcs        = cell(1,g); % small, large
lrs         = cell(1,g); % 
vals        = cell(1,g); % vals --> action values? 
v_example   = nan(N,g);

%% run the particle filter 

for j=1:g % this loop is for the group (control, anxiety)
    for i=1:nsim            
        % [o] = timeseries;
        % simulate dataset(s)
        d                                       = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
        % o                                       = d.outcome(:,cue);
        % o                                       = d.feedback(:,cue);
        o                                       = d.outcome(:,cue);
        [vol(:,i),stc(:,i),lr(:,i),val(:,i)]    = model_parfilter(o,config.model_parameters,lnames{j});   % is that the inference model?       
    end 

    if j==1
        o_example       = o;
    end
    v_example(:,j)      = val(:,1);
    
    outcomes{j}         = outcome; % I don't think this is needed 
    vols{j}             = vol;
    stcs{j}             = stc;
    lrs{j}              = lr;
    vals{j}             = val;
end

% specs(1,:)              = glabels;
t                       = [config.tstable config.tvolatile]; % all trials
strcols                 = {'Stable','Volatile'};

if model == 1
    ncond                   = 4;
else
    ncond                   = 2;
end


% define empty arrays for volatility, stochasticity and learning rates
m_vol       = nan(N,g);
m_stc       = nan(N,g);
m_lr        = nan(N,g);
e_vol       = nan(N,g);
e_stc       = nan(N,g);
e_lr        = nan(N,g);

% this is not needed for now
% if model == 1
%     all_vols    = cell(ncond,2);
%     all_stcs    = cell(ncond,2);
%     all_lrs     = cell(ncond,4);
% else
%     all_vols    = cell(ncond,1);
%     all_stcs    = cell(ncond,1);
%     all_lrs     = cell(ncond,2);
% 
% end
% 
%% loop over groups and compute means and SEs for stochasticity, volatility and learning rates 

for j=1:g

    % compute means and SEs over simulations 
    m_vol(:,j)      = mean(vols{j},2);
    m_stc(:,j)      = mean(stcs{j},2);
    m_lr(:,j)       = mean(lrs{j},2);
    e_vol(:,j)      = serr(vols{j},2);
    e_stc(:,j)      = serr(stcs{j},2);
    e_lr(:,j)       = serr(lrs{j},2);   
    
    % split the trials based on low-high volatility and low-high
    % stochasticity
    for k = 1:2
        l           = (j-1)*2 + k; % have no idea what that does

        all_vols{l,1} = vols{j}(t(:,k),:);
        all_stcs{l,1} = stcs{j}(t(:,k),:);
        all_lrs{l,1}  = lrs{j}(t(:,k),:);

        all_vals{l,1} =  vals{j}(t(:,k),:);

        specs{1,l}  = glabels{j};
        specs{2,l}  = strcols{k};            
        specs{3,l}  = sprintf('%s-%s',glabels{j},strcols{k});            
    end
end 

% compute mean learning rate
ma          = nan(1,ncond);
ea          = nan(1,ncond);


for j = 1:ncond
    a       = mean(all_lrs{j},1)';    

    ma(j)   = mean(a);
    ea(j)   = serr(a);
end

vval        = v_example; % estimated reward

% store output from PF simulations 
stochVolSim = struct('v_example', vval, 'specs',{specs}, 'm_vol', m_vol, 'm_stc',...
    m_stc, 'm_lr', m_lr,'e_vol', e_vol, 'e_stc', e_stc, 'e_lr', e_lr, 'ma', ma, 'ea', ea);


end % end of function 