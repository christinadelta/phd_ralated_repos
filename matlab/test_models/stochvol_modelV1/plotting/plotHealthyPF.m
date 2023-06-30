function [h,f] = plotHealthyPF(mvol,evol,mstc,estc,x)


% plots the estimated volatility and stochasticity on a trial by trial
% basis 

% used to plot the different patamter value combinations

%% define some params 


xstr        = {def_actions('lr'), def_actions('vol'), def_actions('stc')};
glabels     = {'stc=0.1','stc=0.7','stc=1.3','stc=1.9', 'stc=2.5'};

alf         = .3; % I think this is oppacity? 
col         = def_actions('col_br');
fsy         = def_actions('fsy');

% plot volatility and stochasticity 
ii          = size(mvol,2);
counter     = 1;

yll = {[0 0.1] [0 0.2]};

% plot first 4 combs
fsiz = [0 0 .45 1];
figure(1); 
nr          = 4;
nc          = 2;
subplots    = 1:8;

% row 1
[hx, hp]    = plot_signal(nr,nc,subplots(1:2),{mvol(:,1), mstc(:,1)},{evol(:,1), estc(:,1)},xstr(2:3),'',nan,yll,'',col);
lg          = legend(hp(1,:),glabels{1},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
h(1:2)      = hx;

% row 2
[hx2, hp]    = plot_signal(nr,nc,subplots(3:4),{mvol(:,2), mstc(:,2)},{evol(:,2), estc(:,2)},xstr(2:3),'',nan,yll,'',col);
lg          = legend(hp(1,:),glabels{2},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
h(3:4)      = hx2;

% row 3
[hx3, hp]    = plot_signal(nr,nc,subplots(5:6),{mvol(:,3), mstc(:,3)},{evol(:,3), estc(:,3)},xstr(2:3),'',nan,yll,'',col);
lg          = legend(hp(1,:),glabels{3},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
h(5:6)      = hx3;

% raw 4
[hx4, hp]    = plot_signal(nr,nc,subplots(7:8),{mvol(:,4), mstc(:,4)},{evol(:,4), estc(:,4)},xstr(2:3),'',nan,yll,'',col);
lg          = legend(hp(1,:),glabels{4},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
h(7:8)      = hx4;

%%% plot last 4 combs
ffsiz = [0 0 .45 1];
figure(2); 
nr          = 4;
nc          = 2;
subplots    = 1:8;

% row 1
[hx5, hp]    = plot_signal(nr,nc,subplots(1:2),{mvol(:,5), mstc(:,5)},{evol(:,5), estc(:,5)},xstr(2:3),'',nan,yll,'',col);
lg          = legend(hp(1,:),glabels{5},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
f(1:2)      = hx5;

% % row 2
% [hx6, hp]    = plot_signal(nr,nc,subplots(3:4),{mvol(:,6), mstc(:,6)},{evol(:,6), estc(:,6)},xstr(2:3),'',nan,yll,'',col);
% lg          = legend(hp(1,:),glabels{6},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
% f(3:4)      = hx6;
% 
% % row 3
% [hx7, hp]    = plot_signal(nr,nc,subplots(5:6),{mvol(:,7), mstc(:,7)},{evol(:,7), estc(:,7)},xstr(2:3),'',nan,yll,'',col);
% lg          = legend(hp(1,:),glabels{7},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
% f(5:6)      = hx7;
% 
% % raw 4
% [hx8, hp]    = plot_signal(nr,nc,subplots(7:8),{mvol(:,8), mstc(:,8)},{evol(:,8), estc(:,8)},xstr(2:3),'',nan,yll,'',col);
% lg          = legend(hp(1,:),glabels{8},'fontsize',fsy,'location','northwest','box','off','autoupdate','off');
% f(7:8)      = hx8;


end % end of function
