function [h, g, f] = plotPerf(stable_corr,vol_corr)

% visulaise performance for:
% A) different stochasticity and volatility conditions in box plots 
% B) Stable vs Volatile phases with all stc levels averaged in box-plots
% C) Small vs medium vs large stc levels with all phases averaged

% -------
%% define a few params
cols = 2;
rows = 2;

colours = [0 0.4470 0.7410; % blue
    0.8500 0.3250 0.0980;   % orange
    0.9290 0.6940 0.1250];  % yellow

%% extract actions for each condition and level 

% stable vs volatile performance
all_corr(:,1) = stable_corr(:);
all_corr(:,2) = vol_corr(:);

% small vs medium vs large stc performance
stc_corr(:,1) = [stable_corr(:,1); vol_corr(:,1)];
stc_corr(:,2) = [stable_corr(:,2); vol_corr(:,2)];
stc_corr(:,3) = [stable_corr(:,3); vol_corr(:,3)];

%% plot the performance

% plot stable phase
x = {'small','medium','large'}; % how many box plots per subplot?
figure 
subplot(rows,cols,1)

boxplot(stable_corr,x)
h(:,1) = findobj(gca,'Tag','Box');
for j=1:length(x)
    patch(get(h(j,1),'XData'),get(h(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Stable Phase')
ylabel('% correct')
text(-0.03,1.12,'A','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?

% plot volatile phase 
subplot(rows,cols,2)
boxplot(vol_corr,x)
h(:,2) = findobj(gca,'Tag','Box');
for j=1:length(x)
    patch(get(h(j,2),'XData'),get(h(j,2),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Volatile Phase')
ylabel('% correct')
% text(-0.03,1.12,'B','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?


% --------------
% plot stable vs volatile phase for all stc levels averaged 
xx = {'stable','volatile'}; % how many box plots per subplot?
subplot(rows,cols,3)

boxplot(all_corr,xx)
g(:,1) = findobj(gca,'Tag','Box');
for j=1:length(xx)
    patch(get(g(j,1),'XData'),get(g(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Phases')
ylabel('% correct')
text(-0.03,1.12,'B','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?

% ------------
% plot small vs medium vs large stc levels

x = {'small','medium','large'}; % how many box plots per subplot?
subplot(rows,cols,4)

boxplot(stc_corr,x)
f(:,1) = findobj(gca,'Tag','Box');
for j=1:length(x)
    patch(get(f(j,1),'XData'),get(f(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('STC Levels')
ylabel('% correct')
text(-0.03,1.12,'C','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?
fontsize(gcf,16,"points")

% title handle
t = sgtitle('Performance of model over 100 simulations');
t.FontSize = 20;


end % end of function