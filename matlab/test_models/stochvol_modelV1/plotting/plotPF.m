function h = plotPF(nr,nc,subplots,stochVolSim,x)

% created in June 2023

% plots particle filter results:
% - estimated reward 
% - learning rates
% - volatility 
% - stochasticity 

%% define some params 

% extract data from struct 
vval        = stochVolSim.v_example;
ma          = stochVolSim.ma;
ea          = stochVolSim.ea;
m_vol       = stochVolSim.m_vol;
m_stc       = stochVolSim.m_stc;
e_vol       = stochVolSim.e_vol;
e_stc       = stochVolSim.e_stc;
specs       = stochVolSim.specs;
glabels     = {'Control','ASD'};

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


end