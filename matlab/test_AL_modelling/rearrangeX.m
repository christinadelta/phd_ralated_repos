function [Xfit, fitParams, simParams] = rearrangeX(simX, fitX,reps)

% split X fit for each block 
for i = 1:size(fitX,2)
    Xfit.fit1(i,:) = fitX{i}{1,1}; % block 1 - small/stable 
    Xfit.fit2(i,:) = fitX{i}{1,2}; % block 2 - small/volatile
    Xfit.fit3(i,:) = fitX{i}{2,1}; % block 3 - medium/stable
    Xfit.fit4(i,:) = fitX{i}{2,2}; % block 4 - medium/volatile
    Xfit.fit5(i,:) = fitX{i}{3,1}; % block 5 - large/stable
    Xfit.fit6(i,:) = fitX{i}{3,2}; % block 6 - large/volatile

    

    all_lambda_s(i,1) =  fitX{i}{1,1}(1,1); 
    all_lambda_v(i,1) =  fitX{i}{1,1}(1,2); 
end

% for each block store each estimated parameter in different arrays
for i = 1:size(fitX,2)

    % extract lambda_s
    lambdas_s(i,1) =  Xfit.fit1(i,1);
    lambdas_s(i,2) =  Xfit.fit2(i,1);
    lambdas_s(i,3) =  Xfit.fit3(i,1);
    lambdas_s(i,4) =  Xfit.fit4(i,1);
    lambdas_s(i,5) =  Xfit.fit5(i,1);
    lambdas_s(i,6) =  Xfit.fit6(i,1);

    % extract lambda_v
    lambdas_v(i,1) =  Xfit.fit1(i,2);
    lambdas_v(i,2) =  Xfit.fit2(i,2);
    lambdas_v(i,3) =  Xfit.fit3(i,2);
    lambdas_v(i,4) =  Xfit.fit4(i,2);
    lambdas_v(i,5) =  Xfit.fit5(i,2);
    lambdas_v(i,6) =  Xfit.fit6(i,2);

    % extract betas
    betas(i,1) =  Xfit.fit1(i,3);
    betas(i,2) =  Xfit.fit2(i,3);
    betas(i,3) =  Xfit.fit3(i,3);
    betas(i,4) =  Xfit.fit4(i,3);
    betas(i,5) =  Xfit.fit5(i,3);
    betas(i,6) =  Xfit.fit6(i,3);

    % extract s0
    s0(i,1) =  Xfit.fit1(i,4);
    s0(i,2) =  Xfit.fit2(i,4);
    s0(i,3) =  Xfit.fit3(i,4);
    s0(i,4) =  Xfit.fit4(i,4);
    s0(i,5) =  Xfit.fit5(i,4);
    s0(i,6) =  Xfit.fit6(i,4);

     % extract v0
    v0(i,1) =  Xfit.fit1(i,5);
    v0(i,2) =  Xfit.fit2(i,5);
    v0(i,3) =  Xfit.fit3(i,5);
    v0(i,4) =  Xfit.fit4(i,5);
    v0(i,5) =  Xfit.fit5(i,5);
    v0(i,6) =  Xfit.fit6(i,5);

    % extract simulated parameters for each block
    sim_lambda_s(i,1) = simX(i,1);
    sim_lambda_v(i,1) = simX(i,2);
    sim_beta(i,1) = simX(i,3);
    sim_s0(i,1) = simX(i,4);
    sim_v0(i,1) = simX(i,5);

end

fitParams.lambda_s  = lambdas_s;
fitParams.lambda_v  = lambdas_v;
fitParams.betas     = betas;
fitParams.s0        = s0;
fitParams.v0        = v0;

simParams.lambda_s  = sim_lambda_s;
simParams.lambda_v  = sim_lambda_v;
simParams.beta      = sim_beta;
simParams.s0        = sim_s0;
simParams.v0        = sim_v0;

end % end of function