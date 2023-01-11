function [BIC, iBEST, BEST] = fit_all_v1(a, r)

[~, ~, BIC(1)]  = model1_fit(a, r);
[~, ~, BIC(2)]  = model2_fit(a, r);
[~, ~, BIC(3)]  = model3_fit(a, r);
[~, ~, BIC(4)]  = model4_fit(a, r);
[~, ~, BIC(5)]  = model5_fit(a, r);

[M, iBEST]      = min(BIC);
BEST            = BIC == M;
BEST            = BEST / sum(BEST);

return