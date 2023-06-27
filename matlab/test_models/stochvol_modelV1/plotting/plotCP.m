function hf = plotCP(nr, nc, subplots, mps, mpv)

% plot learning rates for different parameter values (main model)
% created in June 2023

%% plot the matrix

% plot stable learning rates
hf = figure(1); clf;
set(gcf, 'Position', [811   613   600   300])
[~,~,~,ax] = easy_gridOfEqualFigures([0.2  0.1], [0.1 0.18 0.04]);
axes(ax(1)); 
t = imageTextMatrix(mps);
set(t(mps'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(mps);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% plot volatile learning rates
axes(ax(2)); 
t = imageTextMatrix(mpv);
set(t(mpv'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(mpv);
set(t, 'fontsize', 14)
xlabel('stochasticity')
ylabel('volatility')

% title 
a = annotation('textbox', [0 0.94 1 0.06]);
set(a, 'string', 'P(choose option A) for different parameter values of volatility and stochasticity', 'fontsize', 28, ...
    'horizontalalignment', 'center', ...
    'linestyle', 'none', 'fontweight', 'normal')


% change the scales
% set(ax, 'xtick', [0.1:0.4:3], 'ytick', [0.1:0.2:1.5], 'fontsize', 14, ...
%     'xaxislocation', 'top', 'tickdir', 'out')
% addABCs(ax, [-0.05 0.08], 40)



end % end of function