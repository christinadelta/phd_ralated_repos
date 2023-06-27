function hf = plotCP_v2(nr, nc, subplots, h_mps, h_mpv, l_mps, l_mpv)

% plot learning rates for different parameter values (main model and lesioned)
% created in June 2023

%%

% plot stable learning rates - healthy
hf = figure(1); clf;
set(gcf, 'Position', [400   405   900   800]);
ax = easy_gridOfEqualFigures([0.02 0.23 0.17], [0.1 0.14 0.03]);
axes(ax(1)); 
t = imageTextMatrix(h_mps);
set(t(h_mps'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(h_mps);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% plot volatile learning rates - healthy
axes(ax(2)); 
t = imageTextMatrix(h_mpv);
set(t(h_mpv'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(h_mpv);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% title 
a = annotation('textbox', [0 0.94 1 0.06]);
set(a, 'string', 'P(choose option A) for different parameter values of volatility and stochasticity (healthy model)', 'fontsize', 28, ...
    'horizontalalignment', 'center', ...
    'linestyle', 'none', 'fontweight', 'normal')

% plot stable learning rates - lesioned
axes(ax(3)); 
t = imageTextMatrix(l_mps);
set(t(l_mps'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(l_mps);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% plot volatile learning rates - lesioned
axes(ax(4)); 
t = imageTextMatrix(l_mpv);
set(t(l_mpv'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(l_mpv);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

a = annotation('textbox', [0 0.94-0.52 1 0.06]);
set(a, 'string', 'P(choose option A) for different parameter values of volatility and stochasticity (lesioned model)', 'fontsize', 28, ...
    'horizontalalignment', 'center',...
    'linestyle', 'none', 'fontweight', 'normal')




end % end of function