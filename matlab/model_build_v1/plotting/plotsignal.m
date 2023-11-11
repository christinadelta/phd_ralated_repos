function f = plotsignal(col,row,m1,vol,x,y)

% plot volatility and predicted signal trial-by-trial

%% ----
figure;
set(gcf,'units','normalized');
% set(gcf,'position',fpos0);

colours = [1.0000 0.4000 0; 0 0.4000 1.0000];
lcol    = [.8 .8 .8];
yl      = [0 .4;.5 .85;-1.9 1.9];

% mark where the switches start
if isrow(x)
    deltax  = [0 diff(x)~=0];
else
    deltax  = [0; diff(x)~=0];
end

% at which trials do the switches start?
switchx     = (find(deltax))-1;
trls        = 1:length(x);

%% plot signals 

% plot volatility 
f(1) = subplot(row,col,1);
plot(vol,'color',colours(1,:),'linewidth',2); hold on;

if ~any(isnan(yl(1,:))), ylim(yl(1,:)); end
ym = get(gca,'ylim');

% define the switches 
for i=1:length(switchx)
    plot([trls(switchx(i)); trls(switchx(i))],ym','color',lcol,'linewidth',2);
end

plot(vol,'color',colours(1,:),'linewidth',2); hold on;
set(gca,'fontname', 'Arial');
ylabel('Volatility','fontsize', 14);
set(gca,'ticklength', [0 0]);


% plot predicted signal
h(2) = subplot(row,col,2); 
plot(m1,'color',colours(1,:),'linewidth',2); hold on;
plot(x,'color',colours(2,:),'linewidth',1); hold on;

% plot outcomes 
if all(~isnan(y))
    plot(y,'.','color','k');
end

set(gca,'fontname','Arial');
ylabel('Predictions','fontsize',14); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',14);
if ~any(isnan(yl(3,:))), ylim(yl(3,:)); end







end % end of function loop