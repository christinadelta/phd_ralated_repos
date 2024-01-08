function [h, g, f] = plot_fitLRs(sublr)

% plot fit lrs 

% visulaise learning rates for:
% A) different stochasticity and volatility conditions in box plots 
% B) Stable vs Volatile phases with all stc levels averaged in box-plots
% C) Small vs medium vs large stc levels with all phases averaged

%% define a few params
cols = 2;
rows = 2;

colours = [0 0.4470 0.7410; % blue
    0.8500 0.3250 0.0980;   % orange
    0.9290 0.6940 0.1250];  % yellow

%% make the arrays

% extract stable and volatile learning rates for small, medium and large
% stc levels
stable_lr(:,1) = sublr{1,1}(:,1);
stable_lr(:,2) = sublr{2,1}(:,1);
stable_lr(:,3) = sublr{3,1}(:,1);

vol_lr(:,1) = sublr{1,2}(:,1);
vol_lr(:,2) = sublr{2,2}(:,1);
vol_lr(:,3) = sublr{3,2}(:,1);

% stable vs volatile lrs
alllrs(:,1)   = stable_lr(:);
alllrs(:,2)   = vol_lr(:);

% small vs medium vs large stc lrs 
stc_lrs(:,1) = [stable_lr(:,1); vol_lr(:,1)];
stc_lrs(:,2) = [stable_lr(:,2); vol_lr(:,2)];
stc_lrs(:,3) = [stable_lr(:,3); vol_lr(:,3)];


%% plot the lrs

% plot stable phase
x = {'small','medium','large'}; % how many box plots per subplot?
figure 
subplot(rows,cols,1)

boxplot(stable_lr,x)
h(:,1) = findobj(gca,'Tag','Box');
for j=1:length(x)
    patch(get(h(j,1),'XData'),get(h(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Stable Phase')
ylabel('Learning Rates')
text(-0.03,1.12,'A','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?

% plot volatile phase 
subplot(rows,cols,2)
boxplot(vol_lr,x)
h(:,2) = findobj(gca,'Tag','Box');
for j=1:length(x)
    patch(get(h(j,2),'XData'),get(h(j,2),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Volatile Phase')
ylabel('Learning Rates')
% text(-0.03,1.12,'B','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?

% --------------
% plot stable vs volatile phase for all stc levels averaged 
xx = {'stable','volatile'}; % how many box plots per subplot?
subplot(rows,cols,3)

boxplot(alllrs,xx)
g(:,1) = findobj(gca,'Tag','Box');
for j=1:length(xx)
    patch(get(g(j,1),'XData'),get(g(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Phases')
ylabel('Learning Rates')
text(-0.03,1.12,'B','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?

% ------------
% plot small vs medium vs large stc levels

x = {'small','medium','large'}; % how many box plots per subplot?
subplot(rows,cols,4)

boxplot(stc_lrs,x)
f(:,1) = findobj(gca,'Tag','Box');
for j=1:length(x)
    patch(get(f(j,1),'XData'),get(f(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('STC Levels')
ylabel('Learning Rates')
text(-0.03,1.12,'C','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?
fontsize(gcf,16,"points")

% title handle
t = sgtitle('Participant 4 Learning Rates');
t.FontSize = 20;


%%








end 