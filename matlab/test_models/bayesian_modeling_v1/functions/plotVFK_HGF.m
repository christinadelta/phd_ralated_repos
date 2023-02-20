function plotVFK_HGF(probability,outc,allsignals)

% Date created: 19/02/2023

% unpack params
m1      = allsignals.m1;
v1      = allsignals.v1;
lr1     = allsignals.lr1;
m2      = allsignals.m2;
v2      = allsignals.v2;
lr2     = allsignals.lr2;

%% ----- 
% define plotting parameters, plotting font-size, font types, etc..
fnt         = 'Helvetica';
fn          = 'Arial';
abc         = 'ABCDEFG';
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

nr          = 3; % number of rows
nc          = 2; % number of columns
fpos0       = [0.2    0.0800    .55*1.0000    .65*0.8133];

%---------
figure;
set(gcf,'units','normalized');
set(gcf,'position',fpos0);

yl          = [0 1.5; .2 1.2;-2.5 2.5];
sub_plts    = [1 3 5];
hl          = sim_plot(nr,nc,sub_plts,probability,outc,v1,lr1,m1,yl);

text(.5,yst,'VKF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hl(1),'HorizontalAlignment','Center','fontweight','bold');

for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hl(i));
end

%---------
sub_plts    = [2 4 6];
hr          = sim_plot(nr,nc,sub_plts,probability,outc,v2,lr2,m2,yl);

text(.5,yst,'HGF','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',hr(1),'HorizontalAlignment','Center','fontweight','bold');

for i=1:2
    text(xsA,ysA,abc(i+2),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',hr(i));
end



end % end of function