function [avCorrs, blockCorrs] = runCorr(simX, fitParams)


% 1. for each fit parameter average values across blocks and correlate with
% corresponding sim parameters
% lambda_s
sim_lambda_s    = simX(:,1);
fit_lambda_s    = mean(fitParams.lambda_s,2);
[R,P] = corrcoef(sim_lambda_s,fit_lambda_s);
R_lambda_s(:,1) = R(2,1);
P_lambda_s(:,1) = P(2,1);

clear R P 

% lambda_v
sim_lambda_v = simX(:,2);
fit_lambda_v = mean(fitParams.lambda_v,2);
[R,P] = corrcoef(sim_lambda_v,fit_lambda_v);
R_lambda_v(:,1) = R(2,1);
P_lambda_v(:,1) = P(2,1);
clear R P 

% beta
sim_beta = simX(:,3);
fit_beta = mean(fitParams.betas,2);
[R,P] = corrcoef(sim_beta,fit_beta);
R_beta(:,1) = R(2,1);
P_beta(:,1) = P(2,1);
clear R P 

% s0
sim_s0    = simX(:,4);
fit_s0    = mean(fitParams.s0,2);
[R,P] = corrcoef(sim_s0,fit_s0);
R_s0(:,1) = R(2,1);
P_s0(:,1) = P(2,1);
clear R P 

% v0
sim_v0    = simX(:,5);
fit_v0    = mean(fitParams.v0,2);
[R,P] = corrcoef(sim_v0,fit_v0);
R_v0(:,1) = R(2,1);
P_v0(:,1) = P(2,1);
clear R P 

avCorrs.R_lambda_s = R_lambda_s;
avCorrs.P_lambda_s = P_lambda_s;
avCorrs.R_lambda_v = R_lambda_v;
avCorrs.P_lambda_v = P_lambda_v;
avCorrs.R_beta = R_beta;
avCorrs.P_beta = P_beta;
avCorrs.R_s0 = R_s0;
avCorrs.P_s0 = P_s0;
avCorrs.R_v0 = R_v0;
avCorrs.P_v0 = P_v0;

%% correlate for each block seperately

for i = 1:6

    % correlate lambda_s
    sim_lambda_s        = simX(:,1);
    fit_lambda_s        = fitParams.lambda_s(:,i);
    [R,P]               = corrcoef(sim_lambda_s,fit_lambda_s);
    R_alllambda_s(:,i)  = R(2,1);
    P_alllambda_s(:,i)  = P(2,1);

    % correlate lambda_v
    sim_lambda_v        = simX(:,2);
    fit_lambda_v        = fitParams.lambda_v(:,i);
    [R,P]               = corrcoef(sim_lambda_v,fit_lambda_v);
    R_alllambda_v(:,i)  = R(2,1);
    P_alllambda_v(:,i)  = P(2,1);

    % correlate betas
    sim_beta            = simX(:,3);
    fit_beta            = fitParams.betas(:,i);
    [R,P]               = corrcoef(sim_beta,fit_beta);
    R_allbeta(:,i)     = R(2,1);
    P_allbeta(:,i)     = P(2,1);

    % correlate betas
    sim_s0            = simX(:,4);
    fit_s0            = fitParams.s0(:,i);
    [R,P]               = corrcoef(sim_s0,fit_s0);
    R_alls0(:,i)     = R(2,1);
    P_alls0(:,i)     = P(2,1);

    sim_v0            = simX(:,5);
    fit_v0            = fitParams.v0(:,i);
    [R,P]               = corrcoef(sim_v0,fit_v0);
    R_allv0(:,i)     = R(2,1);
    P_allv0(:,i)     = P(2,1);

    blockCorrs.R_alllambda_s = R_alllambda_s;
    blockCorrs.P_alllambda_s = P_alllambda_s;
    blockCorrs.R_alllambda_v = R_alllambda_v;
    blockCorrs.P_alllambda_v = P_alllambda_v;
    blockCorrs.R_beta = R_allbeta;
    blockCorrs.P_beta = P_allbeta;
    blockCorrs.R_s0 = R_alls0;
    blockCorrs.P_s0 = P_alls0;
    blockCorrs.R_v0 = R_allv0;
    blockCorrs.P_v0 = P_allv0;

end 



end % end of function