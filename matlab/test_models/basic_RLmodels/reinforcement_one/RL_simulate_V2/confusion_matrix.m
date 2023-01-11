% this script creates a confusuon matrix with the model comparisons 


%%

clear all
clc

% add utilities to the path
addpath(genpath('utilities'))

%% initialise variables

CM  = zeros(5);

T   = 1000;
mu  = [0.2 0.8];

for count = 1:100
    count;
    
    figure(1); clf;
    FM      = round(100*CM/sum(CM(1,:)))/100;
    t       = imageTextMatrix(FM);

    set(t(FM'<0.3), 'color', 'w')
    hold on;
    [l1, l2] = addFacetLines(CM);
    set(t, 'fontsize', 22)
    title(['count = ' num2str(count)]);
    set(gca, 'xtick', [1:5], 'ytick', [1:5], 'fontsize', 28, ...
        'xaxislocation', 'top', 'tickdir', 'out')
    xlabel('fit model')
    ylabel('simulated model')
    
    
    drawnow
    % Model 1
    b                   = rand;
    [a, r]              = model1_sim(T, mu, b);
    [BIC, iBEST, BEST]  = fit_all_v1(a, r);
    CM(1,:)             = CM(1,:) + BEST;
    
    % Model 2
    epsilon             = rand;
    [a, r]              = model2_sim(T, mu, epsilon);
    [BIC, iBEST, BEST]  = fit_all_v1(a, r);
    CM(2,:)             = CM(2,:) + BEST;
    
    % Model 3
    
    alpha               = rand;
    beta                = 1+exprnd(1);
    [a, r]              = model3_sim(T, mu, alpha, beta);
    [BIC, iBEST, BEST]  = fit_all_v1(a, r);
    CM(3,:)             = CM(3,:) + BEST;
    
    % Model 4
    alpha_c             = rand;
    beta_c              = 1+exprnd(1);
    [a, r]              = model4_sim(T, mu, alpha_c, beta_c);
    [BIC, iBEST, BEST]  = fit_all_v1(a, r);
    CM(4,:)             = CM(4,:) + BEST;
    
    % Model 5
    alpha               = rand;
    beta                = 1+exprnd(1);
    alpha_c             = rand;
    beta_c              = 1+exprnd(1);
    [a, r]              = model5_sim(T, mu, alpha, beta, alpha_c, beta_c);
    [BIC, iBEST, BEST]  = fit_all_v1(a, r);
    CM(5,:)             = CM(5,:) + BEST;
    
end


%%
figure(1); 
title('')
set(gcf, 'Position', [811   417   500   400])
set(gca, 'fontsize', 28);
% saveFigurePdf(gcf, '~/Figures/Figure5b')
%
%
% [Xf, LL, BIC] = fit_M1random_v1(a, r);
% [Xf, LL, BIC] = fit_M2WSLS_v1(a, r);
% [Xf, LL, BIC] = fit_M3RescorlaWagner_v1(a, r);
% [Xf, LL, BIC] = fit_M4CK_v1(a, r);
% [Xf, LL, BIC] = fit_M5RWCK_v1(a, r);