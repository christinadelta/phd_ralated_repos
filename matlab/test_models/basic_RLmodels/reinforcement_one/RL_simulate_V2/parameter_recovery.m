% parameter recovery in RL 
% I'll be using this in the RW model

% in Parameter recovery first simulate data (actions and rewards), fit the data and then recover the
% simulated parameters 



%% set colours 

clear all
clc

global AZred AZblue AZcactus AZsky AZriver AZsand AZmesa AZbrick

AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;

% add utilities to the path
addpath(genpath('utilities'))

%% define experimental parameters 

ttrials     = 100; % number of trials
mu          = [0.2 0.8]; % mean reward of bandits

rng(2);

% run the RW model 
for count = 1:1000

    alpha                   = rand; % choose a random number betwen 0 and 1
    beta                    = exprnd(10); % generate a random number from the exponential distributuin with mean value of 10
    [a,r,updatedQvals]      = model3_sim(ttrials, mu, alpha, beta); % simulate RW data
    [Xf, LL, BIC]           = model3_fit(a, r);        % fit RW model

    Xsim(1,count)           = alpha;
    Xsim(2,count)           = beta;
    Xfit(1,count)           = Xf(1); % fitted alpha
    Xfit(2,count)           = Xf(2); % fitted beta

end


%% basic parameter recovery plots

names       = {'learning rate' 'softmax temperature'};
symbols     = {'\alpha' '\beta'};

figure(1); clf;
set(gcf, 'Position', [811   613   600   300])

[~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);

for i= 1:size(Xsim,1)
    
    axes(ax(i)); hold on;
    plot(Xsim(i,:), Xfit(i,:), 'o', 'color', AZred, 'markersize', 8, 'linewidth', 1)
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
    
end

% find 'bad' alpha values
thresh  = 0.25;
ind     = abs(Xsim(1,:) - Xfit(1,:)) > thresh;

for i = 1:2
    axes(ax(i));
    plot(Xsim(i,ind), Xfit(i,ind), 'o', 'color', AZblue, 'markersize', 8, 'linewidth', 1, ...
    'markerfacecolor', [1 1 1]*0.5)
end

set(ax(1,2),'xscale', 'log', 'yscale' ,'log')

axes(ax(1)); t      = title('learning rate');
axes(ax(2)); t(2)   = title('softmax temperature');

axes(ax(1)); xlabel('simulated \alpha'); ylabel('fit \alpha'); 
axes(ax(2)); xlabel('simulated \beta'); ylabel('fit \beta'); 


set(ax, 'tickdir', 'out', 'fontsize', 18)
set(t, 'fontweight', 'normal')
addABCs(ax(1), [-0.07 0.08], 32)
addABCs(ax(2), [-0.1 0.08], 32, 'B')
set(ax, 'tickdir', 'out')

for i = 1:size(Xsim,1)
    axes(ax(i));
    xl = get(gca, 'xlim');
    plot(xl, xl, 'k--')
end


%%

figure(1); clf; hold on;
plot(Xsim(1,:), Xsim(2,:),'.')
% plot(Xfit(2,:), Xfit(1,:),'.')
set(gca, 'xscale', 'log')
