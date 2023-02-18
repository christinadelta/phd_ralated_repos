function h = sim_plot(nr,nc,sub_plts,x,v,m,yl)
fn = getdefaults('fn');
fsy = getdefaults('fsy');

cols = getdefaults('colmap');
lcol = [.8 .8 .8];
%---------

if isrow(x)
    dx = [0 diff(x)~=0];
else
    dx = [0; diff(x)~=0];
end
ix = (find(dx));
tt = 1:length(x);

%---------

h(1) = subplot(nr,nc,sub_plts(1)); 
plot(v,'color',cols(1,:),'linewidth',2); hold on;
if ~any(isnan(yl(1,:)))
    ylim(yl(1,:)); 
else
    ym = get(gca,'ylim');
    ym(2) = 1.01*ym(2);
    set(gca,'ylim',ym);
end
ym = get(gca,'ylim');
for i=1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end
plot(v,'color',cols(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn);
ylabel('Volatility','fontsize',fsy);
set(gca,'ticklength', [0 0]);

h(2) = subplot(nr,nc,sub_plts(2)); 
plot(m,'color',cols(1,:),'linewidth',2); hold on;
if ~any(isnan(yl(2,:))), ylim(yl(2,:)); end
ym = get(gca,'ylim');
for i=1:length(ix)
    plot([tt(ix(i)); tt(ix(i))],ym','color',lcol,'linewidth',2);
end
plot(m,'color',cols(1,:),'linewidth',2); hold on;
set(gca,'fontname',fn);
ylabel('State predictions','fontsize',fsy); hold on;
set(gca,'ticklength', [0 0]);
xlabel('Trial','fontsize',fsy);
end