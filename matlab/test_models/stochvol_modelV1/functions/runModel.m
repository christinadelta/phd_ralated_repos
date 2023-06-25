function stochVolSim = runModel(data, params, config, model, condition, probabilities, trials,condtrials, outpath, outtype)

% runs stochasticity-volatility model (main, lesioned)
% Created in June 2023

%% unpack data , params and config structures

rng(config.rng_id); 
nsim        = config.nsim;
o           = data.outcome; 

% init params to be estimated 
N           = length(o);
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
else
    g = 1;
    glabels     = {'Control'};
    lnames      = {'Healthy'};   

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
        o                                       = d.outcome;
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














end % end of function 