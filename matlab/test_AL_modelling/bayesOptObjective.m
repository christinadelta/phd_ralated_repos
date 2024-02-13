function nll = bayesOptObjective(params, o, a, otherParams, config)
    % Unpack the parameters from the BayesOpt structure
    lambda_s = params.lambda_s;
    lambda_v = params.lambda_v;
    beta = params.beta;
%     s0 = params.s0;
%     v0 = params.v0;
    lambda = params.lambda;
    
    % Construct the parameter array for the RBPF model
    % parameterArray = [lambda_s, lambda_v, beta, s0, v0, lambda];
    parameterArray = [lambda_s, lambda_v, beta, lambda];
    
    % Call your RBPF core function with all necessary inputs
    % nll = rbpf_core_full_reg(parameterArray, o, a, otherParams, config);
    nll = rbpf_lambdaBeta_reg(parameterArray, o, a, otherParams, config);

    
end
