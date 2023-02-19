function h = plotVFK_bin(volatility,learning_rate,m1,probability,outc)
% function h = plotVFK_bin(volatility,learning_rate,m1,probability,outc)

% Date created: 17/02/2023

% plot stuff:
% volatility(:,1)
% learning rate(:,1)
% predected state(:,1) -- m1

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
% colours     = [1 .4 0;0 .4 1]; % orange and blue 
colours     = [0.6350 0.0780 0.1840; 0 0.4470 0.7410]; % red blue 
lcol        = [.8 .8 .8];

figure;
set(gcf,'units','normalized');

nr          = 3; % number rows ?
nc          = 1; % number columns ?
sub_plts    = [1 2 3];
ylm         = [.08 .16; .35 .5; -.2 1.2];
j           = 1;

dx          = [0; diff(probability)~=0]; % deal with feedback probabilities
ix          = (find(dx))-1;
tt          = 1:length(probability);

%% make the subplots 

% plot volatility 
h(1)    = subplot(nr, nc, sub_plts(1)); 
plot(volatility(:,1),'color',colours(1,:),'linewidth',2); hold on;
if ~any(isnan(ylm(1,:))), ylim(ylm(1,:)); end

ym      = get(gca,'ylim');

% specify trial ranges where volatility is high
for i = 1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end

% plot volatility line
plot(volatility(:,1),'color',colours(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn); 

ylabel('Volatility','fontsize',fsy);
set(gca,'ticklength', [0 0]);

% plot learning rate
h(2)    = subplot(nr,nc,sub_plts(2)); 
plot(learning_rate(:,1),'color',colours(1,:),'linewidth',2); hold on;
if ~any(isnan(ylm(2,:))), ylim(ylm(2,:)); end

ym      = get(gca,'ylim');

% specify trial ranges where volatility is high
for i=1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end

plot(learning_rate(:,1),'color',colours(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn);
ylabel('Learning rate','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);

% plot predictions
h(3) = subplot(nr,nc,sub_plts(3)); 
plot(m1(:,1),'color',colours(1,:),'linewidth',2); hold on;
plot(probability,'color',colours(2,:),'linewidth',1); hold on;

if all(~isnan(outc(:,1)))
    plot(outc(:,1),'.','color','k');
end

% deal with labels
set(gca,'fontname',fn);
ylabel('Predictions','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsy);

if ~any(isnan(ylm(3,:))), ylim(ylm(3,:)); end

% deal with legend 
for i=1:3
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
%     set(h(i),'ylim',yl(i,:));
end

hlg = legend(h(3),{'Predicted','True'},'fontsize',fsl,'location','east');
pos = get(hlg,'Position');
pos(2) = 1.2*pos(2);
set(hlg,'Position',pos);

i =1;
text(0.5,yst,'Binary','fontsize',fst,'Unit','normalized','fontname',fnt,'parent',h(i),'HorizontalAlignment','Center','fontweight','bold');

end % end of function

