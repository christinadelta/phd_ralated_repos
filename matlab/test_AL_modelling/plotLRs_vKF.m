function g = plotLRs_vKF(cond_lrs)

% visulaise learning rates for:
% Stable vs Volatile phases with all stc levels averaged in box-plots

%% define a few params
cols = 1;
rows = 1;

colours = [0 0.4470 0.7410; % blue
    0.8500 0.3250 0.0980;   % orange
    0.9290 0.6940 0.1250];  % yellow

%% make the arrays

% extract stable and volatile learning rates for small, medium and large
% stc levels
stable_lr(:,1) = cond_lrs{1,1};
stable_lr(:,2) = cond_lrs{2,1};
stable_lr(:,3) = cond_lrs{3,1};

vol_lr(:,1) = cond_lrs{1,2};
vol_lr(:,2) = cond_lrs{2,2};
vol_lr(:,3) = cond_lrs{3,2};

% stable vs volatile lrs
alllrs(:,1)   = stable_lr(:);
alllrs(:,2)   = vol_lr(:);

%% plot the lrs

% plot stable vs volatile phase for all stc levels averaged 
xx = {'stable','volatile'}; % how many box plots per subplot?
subplot(rows,cols,1)

boxplot(alllrs,xx)
g(:,1) = findobj(gca,'Tag','Box');
for j=1:length(xx)
    patch(get(g(j,1),'XData'),get(g(j,1),'YData'),colours(j,:),'FaceAlpha',.5);
end

xlabel('Phases')
ylabel('Learning Rates')
text(-0.03,1.12,'B','Units', 'Normalized', 'VerticalAlignment', 'Top') % add label?

% title handle
t = sgtitle('Simulated sub learning rates with the volatile Kalman Filter');
t.FontSize = 20;

%%

end