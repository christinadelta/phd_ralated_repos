function [Xfit, fitParams, simParams] = rearrangeX_v2(simX, fitX,reps)

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

% split X sim for each block 
for i = 1:size(fitX,2)
    Xsim.sim1(i,:) = simX{i}{1,1}; % block 1 - small/stable 
    Xsim.sim2(i,:) = simX{i}{1,2}; % block 2 - small/volatile
    Xsim.sim3(i,:) = simX{i}{2,1}; % block 3 - medium/stable
    Xsim.sim4(i,:) = simX{i}{2,2}; % block 4 - medium/volatile
    Xsim.sim5(i,:) = simX{i}{3,1}; % block 5 - large/stable
    Xsim.sim6(i,:) = simX{i}{3,2}; % block 6 - large/volatile

    

    all_lambda_s(i,1) =  simX{i}{1,1}(1,1); 
    all_lambda_v(i,1) =  simX{i}{1,1}(1,2); 
end

%% for each block store each estimated parameter in different arrays

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

end

fitParams.lambda_s  = lambdas_s;
fitParams.lambda_v  = lambdas_v;
fitParams.betas     = betas;
fitParams.s0        = s0;
fitParams.v0        = v0;

clear lambdas_v lambdas_s betas s0 v0

%% for each block store each estimated parameter in different arrays (simulated)

for i = 1:size(fitX,2)

    % extract lambda_s
    lambdas_s(i,1) =  Xsim.sim1(i,1);
    lambdas_s(i,2) =  Xsim.sim2(i,1);
    lambdas_s(i,3) =  Xsim.sim3(i,1);
    lambdas_s(i,4) =  Xsim.sim4(i,1);
    lambdas_s(i,5) =  Xsim.sim5(i,1);
    lambdas_s(i,6) =  Xsim.sim6(i,1);

    % extract lambda_v
    lambdas_v(i,1) =  Xsim.sim1(i,2);
    lambdas_v(i,2) =  Xsim.sim2(i,2);
    lambdas_v(i,3) =  Xsim.sim3(i,2);
    lambdas_v(i,4) =  Xsim.sim4(i,2);
    lambdas_v(i,5) =  Xsim.sim5(i,2);
    lambdas_v(i,6) =  Xsim.sim6(i,2);

    % extract betas
    betas(i,1) =  Xsim.sim1(i,3);
    betas(i,2) =  Xsim.sim2(i,3);
    betas(i,3) =  Xsim.sim3(i,3);
    betas(i,4) =  Xsim.sim4(i,3);
    betas(i,5) =  Xsim.sim5(i,3);
    betas(i,6) =  Xsim.sim6(i,3);

    % extract s0
    s0(i,1) =  Xsim.sim1(i,4);
    s0(i,2) =  Xsim.sim2(i,4);
    s0(i,3) =  Xsim.sim3(i,4);
    s0(i,4) =  Xsim.sim4(i,4);
    s0(i,5) =  Xsim.sim5(i,4);
    s0(i,6) =  Xsim.sim6(i,4);

     % extract v0
    v0(i,1) =  Xsim.sim1(i,5);
    v0(i,2) =  Xsim.sim2(i,5);
    v0(i,3) =  Xsim.sim3(i,5);
    v0(i,4) =  Xsim.sim4(i,5);
    v0(i,5) =  Xsim.sim5(i,5);
    v0(i,6) =  Xsim.sim6(i,5);


end

simParams.lambda_s  = lambdas_s;
simParams.lambda_v  = lambdas_v;
simParams.betas     = betas;
simParams.s0        = s0;
simParams.v0        = v0;


%%





end % end of function