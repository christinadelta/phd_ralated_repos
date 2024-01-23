function loss_corr = plotCorrLoss(sum_simloss, sum_fitloss)

%% first correlate time series

[r,p]            = corrcoef(sum_fitloss, sum_simloss);
pvals            = p(1,2);
loss_corr  = r(1,2);


% Display the correlation coefficient and p-value
fprintf('Correlation coefficient (r): %f\n', loss_corr);
fprintf('P-value: %f\n', pvals);
 
% Interpret the result
if pvals < 0.05
    fprintf('The correlation is statistically significant.\n');
else
    fprintf('The correlation is not statistically significant.\n');
end


end 
