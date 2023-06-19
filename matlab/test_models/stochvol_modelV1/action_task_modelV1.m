% first version of the aversive learning task with the
% stochasticity-volatility model 

% STEPS:

% created May 2023 

clear all
clc

%% set paths and init variables

% set paths
addpath(fullfile(pwd,'functions'))

outpath         = fullfile(pwd, 'output');      addpath(outpath);
figpath         = fullfile(pwd, 'figures');     addpath(figpath);
modelpath       = fullfile(pwd, 'volstoch');    addpath(modelpath);
plotpath        = fullfile(pwd, 'plotting');    addpath(plotpath);

% initialise variables 
subjects        = 1;
condition       = 4;                        % stable & volatile / small & large stochasticity
% task            = 1;                        % stable without switch (if task = 2 then stable with one switch)
probabilities   = [.88 .12;                 % small stochasticity probabilities
    .60 .40];                               % large stochasticity probabilities (either 60:40 or 64:36)
trials          = 200;                      % total trials
condtrials      = [100 25];                 % 100 per stochasticity condition in stable env and 25 trials in volatile condition;
outtype         = 2;                        % if 1 = outcomes are binary [0,1], if 2 = outcome variance [0.01] is added to outcomes

% define initial learning rate inverse temperature parameters (don't think
% that this will be used)
params          = [.2 .6;
    3 6]; %alpha and beta parameteres 

%% simulate dataset

for sub = 1:subjects

    % simulate dataset(s)
    data            = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
    alldata{sub,1}  = data;


end % end of subjects loop

%% run particle filter 

% unpack data struct 
o           = data.outcome; 
x           = data.state;
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);

% define parameters and config for the particle filter 
params      = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',2,'model_parameters',params);

rng(config.rng_id); 
nsim = config.nsim;

% init params to be estimated 
N           = length(o);
outcome     = nan(N,nsim);
vol         = nan(N,nsim);
stc         = nan(N,nsim);
lr          = nan(N,nsim);
val         = nan(N,nsim);  

% init cells to store the above parameters
outcomes    = cell(1,2);
vols        = cell(1,2); % stable, volatile
stcs        = cell(1,2); % small, large
lrs         = cell(1,2); % 
vals        = cell(1,2); % vals --> action values? 
v_example   = nan(N,2);
glabels     = {'Control','Anxious'};
lnames      = {'Healthy', sprintf('%s lesion',def_actions('stc'))};      

% run the particle filter 
for j=1:2 % I think this loop is for the group (control, anxiety)
    for i=1:nsim            
        % [o] = timeseries;
        % simulate dataset(s)
        data                                    = action_simdataV1(condition, probabilities, trials,condtrials, outpath, outtype);
        o                                       = data.outcome;
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

% estimate volatility, stochasticity and learning rates
m_vol       = nan(N,2);
m_stc       = nan(N,2);
m_lr        = nan(N,2);
e_vol       = nan(N,2);
e_stc       = nan(N,2);
e_lr        = nan(N,2);

all_vols    = cell(N,2);
all_stcs    = cell(N,2);
all_lrs     = cell(N,4);

% loop over groups?
for j=1:2

    % compute means and SEs
    m_vol(:,j)  = mean(vols{j},2);
    m_stc(:,j)  = mean(stcs{j},2);
    m_lr(:,j)   = mean(lrs{j},2);
    e_vol(:,j)  = serr(vols{j},2);
    e_stc(:,j)  = serr(stcs{j},2);
    e_lr(:,j)   = serr(lrs{j},2);   

    for k = 1:2
        l = (j-1)*2 + k; % have no idea what that does

        all_vols{l} = vols{j}(t(:,k),:);
        all_stcs{l} = stcs{j}(t(:,k),:);
        all_lrs{l} = lrs{j}(t(:,k),:);

        all_vals{l} =  vals{j}(t(:,k),:);

        specs{1,l} = glabels{j};
        specs{2,l} = strcols{k};            
        specs{3,l} = sprintf('%s-%s',glabels{j},strcols{k});            
    end
end 

% compute mean learning rate
ncond       = 4;
ma          = nan(1,ncond);
ea          = nan(1,ncond);

for j = 1:ncond
    a       = mean(all_lrs{j},1)';    

    ma(j)   = mean(a);
    ea(j)   = serr(a);
end

vval        = v_example; % estimated reward

%% plot first results

fsiz = [0 0 .45 1];

figure; 

nr          = 2;
nc          = 2;
subplots    = 1:4;
labels      = specs(1:2,:);
clabels     = specs(2,[1 2]);
xstr        = {def_actions('lr'), def_actions('vol'), def_actions('stc')};

alf = .3; % I think this is oppacity? 
col = def_actions('col_br');
fsy = def_actions('fsy');

% plot the estimated reward (expected values computed with the particle
% filter) and the true state of the environment
N       = size(vval,1); % trials
Nline   = nan;
ii      = [2 1];
[hx]    = plot_signal(nr,nc,subplots(1),{vval(:,ii)},{zeros(N,2)},{'Estimated reward'},'',Nline,[],'',col(ii,:));
hold on;
plot(hx, x,'color',.6*ones(1,3),'linewidth',2);

h(1) = hx;

% plot learning rates 
h(2) = plot_bar(nr,nc,subplots(2),{ma},{ea},labels,xstr(1));
lg = legend(h(2),clabels,'fontsize',fsy,'location','north','box','off');

% plot volatility and stochasticity 
ii = [1 2];
[hx, hp] = plot_signal(nr,nc,subplots(3:4),{ m_vol(:,ii), m_stc(:,ii)},{e_vol(:,ii), e_stc(:,ii)},xstr(2:3),'',nan,[],'',col);
lg = legend(hp(1,:),glabels,'fontsize',fsy,'location','northwest','box','off','autoupdate','off');

h(3:4) = hx;

%% generate responses

xstate = x;
% generate responses for control and clinical simulated groups using the softmax function 
[xx, dp, trialstable, trialvolatile] = responseModel_v1(xstate, vval);

%% plot second results

%% run stoch-vol

%% plot third results 


