function hf = plotLR_v2(nr, nc, subplots, hma_stable, hma_vol, lma_stable, lma_vol)

% plot learning rates for different parameter values (main model and lesioned)
% created in June 2023

%%

% plot stable learning rates - healthy
hf = figure(1); clf;
set(gcf, 'Position', [400   405   900   800]);
ax = easy_gridOfEqualFigures([0.02 0.23 0.17], [0.1 0.14 0.03]);
axes(ax(1)); 
t = imageTextMatrix(hma_stable);
set(t(hma_stable'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(hma_stable);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% plot volatile learning rates - healthy
axes(ax(2)); 
t = imageTextMatrix(hma_vol);
set(t(hma_vol'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(hma_vol);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% title 
a = annotation('textbox', [0 0.94 1 0.06]);
set(a, 'string', 'Learning rates for different parameter values of volatility and stochasticity (healthy model)', 'fontsize', 28, ...
    'horizontalalignment', 'center', ...
    'linestyle', 'none', 'fontweight', 'normal')

% plot stable learning rates - lesioned
axes(ax(3)); 
t = imageTextMatrix(lma_stable);
set(t(lma_stable'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(lma_stable);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% plot volatile learning rates - lesioned
axes(ax(4)); 
t = imageTextMatrix(lma_vol);
set(t(lma_vol'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(lma_vol);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

a = annotation('textbox', [0 0.94-0.52 1 0.06]);
set(a, 'string', 'Learning rates for different parameter values of volatility and stochasticity (lesioned model)', 'fontsize', 28, ...
    'horizontalalignment', 'center',...
    'linestyle', 'none', 'fontweight', 'normal')




end % end of function