function h = plot_SimpleBayes_V2(q,Eq,y,H,EH)
% function h = plot_SimpleBayes_V2(q,Eq,y,H,EH)

% plot simple Bayesian model results:
% q (true probabilities
% estimated q
% outcomes for highly rewarded p(vertical)
% H % estimated H

%% ----- 
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
fsy2        = 10;
% colours     = [1 .4 0;0 .4 1]; % orange and blue 
% coloursEV   = [0.6350 0.0780 0.1840; 0 0.4470 0.7410];
colours     = [0.6350 0.0780 0.1840;    % red
    0 0.4470 0.7410;                    % blue 
    0 0 0];                             % black
lcol        = [.8 .8 .8];

figure;
set(gcf,'units','normalized');

% define a few more ploting parameters
nr          = 2; % number rows 
nc          = 1; % number columns 
sub_plts    = [1 2];
ylm         = [0 0.1];
ylm2        = [.08 .16; .35 .5; -.2 1.2];
j           = 1;

dx          = [0; diff(q)~=0];    % deal with feedback probabilities (when 1 = switch of probabilties occured)
ix          = (find(dx))-1;                 % trials before each switch of probabilities
tt          = 1:length(q);        % total trials 

%% make the subplots 

% plot expected values as dotted lines 
h(1)    = subplot(nr, nc, sub_plts(1)); 
plot([0 length(y)],[1/25 1/25],'color',colours(1,:),'linewidth',2); hold on;   % plot Q value for vertical gabor in red
plot(EH,':','color',colours(2,:),'linewidth',2); hold on;   % plot Q value for horizontal gabor in blue
if ~any(isnan(ylm(1,:))), ylim(ylm(1,:)); end

ym      = get(gca,'ylim');

% specify trial ranges where volatility is high and plot vertical lines
for i = 1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end

% plot H and Expected H
plot([0 length(y)],[1/25 1/25],'color',colours(1,:),'linewidth',2); hold on;
plot(EH,':','color',colours(2,:),'linewidth',2); hold on;
set(gca,'fontname',fn); 

ylabel('Hazard rate (H)','fontsize',fsy);
set(gca,'ticklength', [0 0]);

% plot q and predicted q
h(2) = subplot(nr,nc,sub_plts(2)); 
plot(Eq,'color',colours(1,:),'linewidth',2); hold on;
plot(q,'color',colours(3,:),'linewidth',1); hold on;

if all(~isnan(y(:,1)))
    plot(y(:,1),'.','color','k');
end

set(gca,'fontname',fn);
ylabel('Estimated q','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsl);

if ~any(isnan(ylm2(3,:))), ylim(ylm2(3,:)); end

%% deal with legends

% add plot number 
for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
%     set(h(i),'ylim',yl(i,:));
end

% add legend to plot A
hlg = legend(h(1),{'true H','estimated H'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.2*pos(2);
set(hlg,'Position',pos);

% add legend to plot B 
hlg2 = legend(h(2),{'true q', 'estimated q'},'fontsize',fsl,'location','east');
pos2 = get(hlg2,'Position');
pos2(2) = 1.2*pos2(2);
set(hlg2,'Position',pos2);

% i =1;
% text(0.5,yst,'Binary','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');


end % end of function