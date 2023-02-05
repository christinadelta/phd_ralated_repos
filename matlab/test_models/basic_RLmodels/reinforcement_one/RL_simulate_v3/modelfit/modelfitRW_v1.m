function [xfit ll]  = modelfitRW_v1(params, actions, rewards)

% this function runs fmincon.m to get optimal X parameter values
% is used to compute -ll

%--------------------
llFunc  = @(x) likRW_v1(actions, rewards, x(1), x(2));

% what are the initial param values?
x0              = params;
lb              = [0 0];    %lower bound
ub              = [1 inf];  % upper bound

% search for optimal param values and compute negative ll
[xfit negll]    = fmincon(llFunc, x0, [], [], [], [], lb, ub);

ll = -negll;

end  % end of function 