function t = plotheatmap(avhm)
% function will go here


% this function creates heatmap for performance using different
% combinations of alpha and beta parameters

%----------------------

% make figure
figure(1); clf;
set(gcf, 'Position', [400   405   900   800]);
ax = easy_gridOfEqualFigures([0.02 0.23], [0.1 0.14]); % one figure for each condition
axes(ax(1));
t = imageTextMatrix(avhm);
set(t(avhm'<0.3), 'color', 'w')
hold on;
[l1, l2] = addFacetLines(avhm);
set(t, 'fontsize', 18)
xlabel('alpha values')
ylabel('beta values')


a = annotation('textbox', [0 0.94 1 0.06]);
set(a, 'string', 'performance (%) using different combinations of parameter values', 'fontsize', 38, ...
    'horizontalalignment', 'center', ...
    'linestyle', 'none', 'fontweight', 'normal')


set(ax, 'xtick', [1:5], 'ytick', [1:5], 'fontsize', 24, ...
    'xaxislocation', 'top', 'tickdir', 'out')
addABCs(ax, [-0.05 0.08], 40)

end % end of function