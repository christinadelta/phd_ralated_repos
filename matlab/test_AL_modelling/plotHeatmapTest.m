function h = plotHeatmapTest(Hopt)

% normalise matrix
normH = Hopt / 10^8;
xvalues = {'lambda_s','lambda_v','beta','s0','v0'};
yvalues = {'lambda_s','lambda_v','beta','0','v0'};

h = heatmap(xvalues,yvalues,normH);

h.Title = 'Hessian Matrix Full Model';
h.XLabel = 'Parameters';
fontsize(gcf,16,"points")

end % end of function