function h = plotRW_v1(x,Ps,Qvals,a)
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
ylm         = [-1.1 1.1];
ylm2        = [.08 .16; .35 .5; -.2 1.2];
j           = 1;

dx          = [0; diff(x)~=0];    % deal with feedback probabilities (when 1 = switch of probabilties occured)
ix          = (find(dx))-1;                 % trials before each switch of probabilities
tt          = 1:length(x);        % total trials 

%% make the subplots 

% plot expected values as dotted lines 
h(1)    = subplot(nr, nc, sub_plts(1)); 
plot(Qvals(:,1),':','color',colours(1,:),'linewidth',2); hold on;   % plot Q value for high probability of loss
plot(Qvals(:,2),':','color',colours(2,:),'linewidth',2); hold on;   % plot Q value for low probability of loss
if ~any(isnan(ylm(1,:))), ylim(ylm(1,:)); end

ym      = get(gca,'ylim');
xlim([1 420]);

% specify trial ranges where volatility is high and plot vertical lines
for i = 1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end

% plot expected values lines 
plot(Qvals(:,1),':','color',colours(1,:),'linewidth',2); hold on;
plot(Qvals(:,2),':','color',colours(2,:),'linewidth',2); hold on;
set(gca,'fontname',fn); 

ylabel('Expected Values','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsl);

% plot choice probabilities
h(2) = subplot(nr,nc,sub_plts(2)); 
plot(Ps(:,1),'color',colours(1,:),'linewidth',2); hold on; % probability values for choosing the low prob/good option
plot(Ps(:,2),':','color',colours(2,:),'linewidth',2); hold on;
plot(x,'color',colours(3,:),'linewidth',1); hold on;

if all(~isnan(a(1,:)))
    plot(a(1,:),'.','color','k');
end

set(gca,'fontname',fn);
ylabel('Choice Probabilities','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsl);

if ~any(isnan(ylm2(3,:))), ylim(ylm2(3,:)); end
xlim([1 420]);

%% deal with legends

% add plot number 
for i=1:2
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
%     set(h(i),'ylim',yl(i,:));
end

% add legend to plot A
hlg = legend(h(1),{'EV for option A','EV for option B'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.2*pos(2);
set(hlg,'Position',pos);

% add legend to plot B 
hlg2 = legend(h(2),{'p(option A)','p(option B)', 'reward rate'},'fontsize',fsl,'location','east');
pos2 = get(hlg2,'Position');
pos2(2) = 1.2*pos2(2);
set(hlg2,'Position',pos2);

i =1;
text(0.5,yst,'alpha = 0.3, beta = 1.5','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');


end % end of function 