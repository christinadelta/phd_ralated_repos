function h = plotRW_V2(probability,choiceProb,Vvals,actions)
% function h = plotRW_V2(probability,choiceProb,Vvals,actions)

% plot RW results:
% choices/actions
% expected values(:,:)
% choice probabilities (:,1)

%% ----- 
% define plotting parameters, plotting font-size, font types, etc..
fnt         = 'Helvetica';
fn          = 'Arial';
abc         = 'ABC';
fsl         = 14; 
fsA         = 18;
xsA         = -.15;
ysA         = 1.1;
yst         = 1.05; %0.15;
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
ylm         = [-0.1 1.1];
ylm2        = [.08 .16; .35 .5; -.2 1.2];
j           = 1;

dx          = [0; diff(probability)~=0];    % deal with feedback probabilities (when 1 = switch of probabilties occured)
ix          = (find(dx))-1;                 % trials before each switch of probabilities
tt          = 1:length(probability);        % total trials 

%% make the subplots 

% plot expected values as dotted lines 
h(1)    = subplot(nr, nc, sub_plts(1)); 
plot(Vvals(:,1),':','color',colours(1,:),'linewidth',2); hold on;   % plot Q value for vertical gabor in red
plot(Vvals(:,2),':','color',colours(2,:),'linewidth',2); hold on;   % plot Q value for horizontal gabor in blue
if ~any(isnan(ylm(1,:))), ylim(ylm(1,:)); end

ym      = get(gca,'ylim');

% specify trial ranges where volatility is high and plot vertical lines
for i = 1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end

% plot expected values lines 
plot(Vvals(:,1),':','color',colours(1,:),'linewidth',2); hold on;
plot(Vvals(:,2),':','color',colours(2,:),'linewidth',2); hold on;
set(gca,'fontname',fn); 

ylabel('Expected Values','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsl);

% plot choice probabilities
h(2) = subplot(nr,nc,sub_plts(2)); 
plot(choiceProb(:,1),'color',colours(1,:),'linewidth',2); hold on;
plot(probability,'color',colours(3,:),'linewidth',1); hold on;

if all(~isnan(actions(:,1)))
    plot(actions(:,1),'.','color','k');
end

set(gca,'fontname',fn);
ylabel('Choice Probabilities','fontsize',fsy); hold on;
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
hlg = legend(h(1),{'EV vertical','EV horizontal'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.2*pos(2);
set(hlg,'Position',pos);

% add legend to plot B 
hlg2 = legend(h(2),{'p(choose vertical)', 'p(reward|verical)'},'fontsize',fsl,'location','east');
pos2 = get(hlg2,'Position');
pos2(2) = 1.2*pos2(2);
set(hlg2,'Position',pos2);

i =1;
text(0.5,yst,'alpha = 0.25, beta = 4','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');


end % end of function 