%% run particle filter 

% unpack data struct 
o           = data.outcome; 
x           = data.state;
tvolatile   = data.t(:,2);
tstable     = data.t(:,1);

% define parameters and config for the particle filter 
params      = struct('nparticles',100,'x0_unc',1,'lambda_v',.2,'lambda_s',.2,'v0',.1,'s0',.1,'s0_lesioned',0.001);
config      = struct('tvolatile',tvolatile,'tstable',tstable,'state',x,'rng_id',0,'nsim',10,'model_parameters',params);

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
glabels     = {'Control','ASD'};
lnames      = {'Healthy', sprintf('%s lesion',def_actions('stc'))};      

% run the particle filter 
for j=1:2 % this loop is for the group (control, anxiety)
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
ncond                   = 4;

% estimate volatility, stochasticity and learning rates
m_vol       = nan(N,2);
m_stc       = nan(N,2);
m_lr        = nan(N,2);
e_vol       = nan(N,2);
e_stc       = nan(N,2);
e_lr        = nan(N,2);

all_vols    = cell(ncond,2);
all_stcs    = cell(ncond,2);
all_lrs     = cell(ncond,4);

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

alf         = .3; % I think this is oppacity? 
col         = def_actions('col_br');
fsy         = def_actions('fsy');

% plot the estimated reward (expected values computed with the particle
% filter) and the true state of the environment
N           = size(vval,1); % trials
Nline       = nan;
ii          = [2 1];
[hx]        = plot_signal(nr,nc,subplots(1),{vval(:,ii)},{zeros(N,2)},{'Estimated reward'},'',Nline,[],'',col(ii,:));
hold on;
plot(hx, x,'color',.6*ones(1,3),'linewidth',2);

h(1)        = hx;

% plot learning rates 
h(2)        = plot_bar(nr,nc,subplots(2),{ma},{ea},labels,xstr(1));
lg          = legend(h(2),clabels,'fontsize',fsy,'location','north','box','off');

% plot volatility and stochasticity 
ii          = [1 2];
[hx, hp]    = plot_signal(nr,nc,subplots(3:4),{ m_vol(:,ii), m_stc(:,ii)},{e_vol(:,ii), e_stc(:,ii)},xstr(2:3),'',nan,[],'',col);
lg          = legend(hp(1,:),glabels,'fontsize',fsy,'location','northwest','box','off','autoupdate','off');

h(3:4)      = hx;

%% generate responses

xstate      = x;
xx          = nan(nsim, ncond/2); % I guess I wont include volatility for now 
mxx         = nan(ncond/2, 2); % mean correct 
exx         = nan(ncond/2, 2); % se correct 

% loop over (volatility) conditions 
for j = 1:ncond/2

    % loop over simulations
    for i = 1:nsim

        thisVal = vals{j}(:,i);

        % generate responses for control and clinical simulated groups using the softmax function 
        [xx(i,:)] = responseModel_v1(xstate, thisVal);

    end % end of simulations loop

    mxx(j,:) = median(xx);
    exx(j,:) = se_median(xx);
end % end of volatility condition


%% plot second results

col = def_actions('col');
fsy = def_actions('fsy');
subplots = 1;

labels = {'Stable','Volatile'};

h = plot_bar(1,1,subplots(1),{mxx},{exx},{'control','ASD'},{'Performance','Performance'},'',col);
set(h,'ylim',[0 1]);
legend(h,labels,'fontsize',fsy,'location','north','box','off');
title(h,'Model');

%% run stoch-vol particle filter

% simulate data for the lesioned models (compute learning rates, true
% voaltility and true stochasticity)
factors.stc     = {'Small','Small','Large','large'};
factors.vol     = {'Small','Large','Large','Small'};

parameters      = struct('nparticles',100,'x0_unc',100,'lambda_v',.1,'lambda_s',.1,'v0',1,'s0',2);
config          = struct('true_vol',[.5 1.5 .5 1.5],'true_stc',[1 1 3 3],'factors',factors,'N',400,...
    'rng_id',0,'nsim',2,'model_parameters',parameters);
            
rng(config.rng_id);
true_vol        = config.true_vol;
true_stc        = config.true_stc;
N               = config.N;
nsim            = config.nsim;    

% init cells to store estimated volatiltiies, stochasticities and learning
% rates
outcomes        = cell(1,4);
vols            = cell(1,4);
stcs            = cell(1,4);
lrs             = cell(1,4);
vals            = cell(1,4);
specs           = cell(3,4);

for j=1:4 % 4 conditions

    outcome = nan(N,nsim);        
    vol     = nan(N,nsim);
    stc     = nan(N,nsim);
    lr      = nan(N,nsim);
    val     = nan(N,nsim);  

    for i=1:nsim
       
        data                                    = action_simdataV1(condition, probabilities, trials, condtrials, outpath, outtype);
        outcome(:,i)                            = data.outcome;  
        [vol(:,i),stc(:,i),lr(:,i),val(:,i)]    = model_parfilter(outcome(:,i),config.model_parameters);
    end

    % store for each condition
    outcomes{j} = outcome;        
    specs{1,j}  = sprintf('%s %s',factors.stc{j},lower(def_actions('stc')));
    specs{2,j}  = sprintf('%s %s',factors.vol{j},lower(def_actions('vol')));
    specs{3,j}  = sprintf('true_vol=%0.2f, true_stc=%0.2f',true_vol(j),true_stc(j));        
    
    vols{j}     = vol;
    stcs{j}     = stc;
    lrs{j}      = lr;
    vals{j}     = val;
end

config.specs    = specs;    
sim             = struct('config',config,'specs',{config.specs},'outcomes',{outcomes},...
    'vols',{vols},'stcs',{stcs},'lrs',{lrs},'vals',{vals}); %#ok<NASGU>

clear outocme vol stc lr val

% compute means and SEs
N       = sim.config.N;
nsim    = sim.config.nsim;

vol     = nan(N,4);
stc     = nan(N,4);
lr      =  nan(N,4);
e_vol   = nan(N,4);
e_stc   = nan(N,4);
e_lr    =  nan(N,4);

tend    = (.9*N+1):N;
a       = nan(nsim,4); % lr 
v       = nan(nsim,4); % vol
u       = nan(nsim,4); % stoch

for j=1:4            
    vol(:,j)    = mean(sim.vols{j},2);
    stc(:,j)    = mean(sim.stcs{j},2);
    lr(:,j)     = mean(sim.lrs{j},2);

    e_vol(:,j)  = serr(sim.vols{j},2);
    e_stc(:,j)  = serr(sim.stcs{j},2);
    e_lr(:,j)   = serr(sim.lrs{j},2);

    a(:,j)      = mean(sim.lrs{j}(tend,:),1);
    v(:,j)      = mean(sim.vols{j}(tend,:),1);
    u(:,j)      = mean(sim.stcs{j}(tend,:),1);

end

mv = mean(v); % mean volatility 
ms = mean(u); % mean stochasticity
ma = mean(a); % mean learning rate
ev = serr(v);
es = serr(u);
ea = serr(a);    

% split the conditions
ii1     = [1 3];
ii2     = [2 4];
ma      = [ma(ii1)' ma(ii2)'];
mv      = [mv(ii1)' mv(ii2)'];
ms      = [ms(ii1)' ms(ii2)'];

ea      = [ea(ii1)' ea(ii2)'];
ev      = [ev(ii1)' ev(ii2)'];
es      = [es(ii1)' es(ii2)'];

%% plot third results 

% init ploting 
colstrs     = {'Small','Large'};
glbl        = {sprintf('Small true %s',lower(def_actions('vol'))),sprintf('Large true %s',lower(def_actions('vol')))};
levels      = {'Small','Large'};
lgtitle     = sprintf('True %s',lower(def_actions('stc')));
xltitle     = sprintf('True %s',lower(def_actions('vol')));

xstr        = {def_actions('lr'), def_actions('vol'), def_actions('stc')};
ylb         = {[0 .75],[0 0.05],[0 0.05]};
yls         = {[0.2 .81],[0 1.8],[0 3.5]};
fsy         = def_actions('fsy');
abc         = def_actions('abc');

nr          = 1;
nc          = 3;     
subplots    = 1:3;    
fsiz        = [0 0 .6 .25];
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

h(1:3)      = plot_bar(nr,nc,subplots(1:3),{ma',mv',ms'},{ea',ev',es'},colstrs,xstr,abc,[],ylb);
loc         = {'northwest','northwest','north'};

for i=1:length(h)
    xlabel(h(i),xltitle,'fontsize',fsy);
    if i==2
    lg = legend(h(i),levels,'fontsize',fsy,'location',loc{i},'box','off');
    title(lg,lgtitle);    
    end
end