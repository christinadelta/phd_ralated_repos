function h = plotHeatmapFull(Hmatrix)

% plot hessian matrix from parameter recovery 

%% define a few params
cols = 3;
rows = 2;

%% extract the matrices and average over repetitions 

% average condition 1
tmp_cond1 = Hmatrix{1,1};

cond1_hmat = mean(tmp_cond1,3);
cond1_mat = cond1_hmat/10^9;

% average condition 2
tmp_cond2 = Hmatrix{1,2};
cond2_hmat = mean(tmp_cond2,3);
cond2_mat = cond2_hmat/10^9;

% average condition 3
tmp_cond3 = Hmatrix{2,1};
cond3_hmat = mean(tmp_cond3,3);
cond3_mat = cond3_hmat/10^9;

% average condition 4
tmp_cond4 = Hmatrix{2,2};
cond4_hmat = mean(tmp_cond4,3);
cond4_mat = cond4_hmat/10^9;

% average condition 5
tmp_cond5 = Hmatrix{3,1};
cond5_hmat = mean(tmp_cond5,3);
cond5_mat = cond5_hmat/10^9;


% average condition 5
tmp_cond6 = Hmatrix{3,2};
cond6_hmat = mean(tmp_cond6,3);
cond6_mat = cond6_hmat/10^9;

%% plot matrices as heatmaps

xvalues = {'lambda_s','lambda_v','beta','s_0','v_0', 'lambda reg'};
yvalues = {'lambda_s','lambda_v','beta','s_0','v_0','lambda reg'};

figure 

% plot column 1 (blocks 1-3)
% block 1
subplot(rows,cols,1)
h(1)            = heatmap(xvalues,yvalues,cond1_mat);
heatmapPosition = h(1).Position; % Get the position of the heatmap
axes('Position', heatmapPosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Create invisible axes on top of the heatmap
text(0.5, 1.1, 'Block 1', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % Add text at a desired location

% block 2
subplot(rows,cols,2)
h(2)            = heatmap(xvalues,yvalues,cond2_mat);
heatmapPosition = h(2).Position; % Get the position of the heatmap
axes('Position', heatmapPosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Create invisible axes on top of the heatmap
text(0.5, 1.1, 'Block 2', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % Add text at a desired location

% block 3
subplot(rows,cols,3)
h(3)            = heatmap(xvalues,yvalues,cond3_mat);
heatmapPosition = h(3).Position; % Get the position of the heatmap
axes('Position', heatmapPosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Create invisible axes on top of the heatmap
text(0.5, 1.1, 'Block 3', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % Add text at a desired location

% plot column 2 (blocks 4-6)
% block 4
subplot(rows,cols,4)
h(4)            = heatmap(xvalues,yvalues,cond4_mat);
heatmapPosition = h(4).Position; % Get the position of the heatmap
axes('Position', heatmapPosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Create invisible axes on top of the heatmap
text(0.5, 1.1, 'Block 4', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % Add text at a desired location

% block 5
subplot(rows,cols,5)
h(5)            = heatmap(xvalues,yvalues,cond5_mat);
heatmapPosition = h(5).Position; % Get the position of the heatmap
axes('Position', heatmapPosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Create invisible axes on top of the heatmap
text(0.5, 1.1, 'Block 5', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % Add text at a desired location

% block 6
subplot(rows,cols,6)
h(6)            = heatmap(xvalues,yvalues,cond6_mat);
heatmapPosition = h(6).Position; % Get the position of the heatmap
axes('Position', heatmapPosition, 'Color', 'none', 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Create invisible axes on top of the heatmap
text(0.5, 1.1, 'Block 6', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top'); % Add text at a desired location

% annotation('textbox',[0.05 0.9 0.1 0.1],'String','Hessian Matrix Full Model','EdgeColor',[105/255 105/255 105/255])
%t = sgtitle('Hessian Matrix Full Model');
fontsize(gcf,14,"points")

%%

end % end of function