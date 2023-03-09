function h = plotHGF_bin(plt_vars)

% Date created: 9/03/2023

% plot stuff:
% volatility(:,1)
% learning rate(:,1)
% predected state(:,1) -- m1

%% ----------
% define plotting parameters, plotting font-size, font types, etc..
fnt         = 'Helvetica';
fn          = 'Arial';
abc         = 'ABC';
fsl         = 14; 
fsA         = 18;
xsA         = -.15;
ysA         = 1.1;
yst         = 1.15; %0.15;
fst         = 18;
fsy         = 14;
% colours     = [1 .4 0;0 .4 1]; % orange and blue 
colours     = [0.6350 0.0780 0.1840; 0 0.4470 0.7410]; % red blue 
lcol        = [.8 .8 .8];

%% unpack structure

fp                  = plt_vars.fp;
u                   = plt_vars.u;
actions             = plt_vars.y;
rho                 = plt_vars.rho;
kappa               = plt_vars.kappa;
omega               = plt_vars.omega;
lr                  = plt_vars.lr;
state_predictions   = plt_vars.sp;
vol                 = plt_vars.volatility;

%% 










end