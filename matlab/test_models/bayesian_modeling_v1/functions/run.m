function [ckv, tx] = run(modelname, simy,simx)

% how many simulations?
nsim = 500;

% choose model to run
switch modelname
    case 'hgf'
        model = @model_hgf;
        d = 3;
    case 'vkf'
        model = @model_vkf;
        d = 3;
    otherwise
        error('!');
end

rng(0);
Ys = simy;
Xs = simx;

% variance of the fit?
vfit = 15.23; % variance 

% get prior and variance 
config = struct('range',[-5*ones(1,d);5*ones(1,d)],'numinit',7*d,'prior_for_bads',1);
prior = struct('mean',zeros(d,1),'variance',vfit);
cbm  = cbm_lap({Ys}, model, prior, [], config); %#ok<NASGU>

% evaluate model
[ckv, tx] = eval_model(simy, model, cbm);