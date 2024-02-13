function NLL_reg = rbpf_lambdaBeta_reg(x,o,a,parameters,config)
% function [val,vol,unp,lr,unc] = rbpf_core(o,x0_unc,lambda_v,lambda_u,v0,u0,nparticles)

%------

% Fitting the RB particle filter to participant data 
% Inputs: 
          % 1) parameters structure: 
                % 1st: lambda_s, -- free
                % 2nd: lambda_v, -- free
                % 3rd: beta, -- free
                % 4rd: s0, -- fixed
                % 5th: v0, -- fixed
                % 6th: nparticles
                % 7th: x0_unc
                % 8th: model type (1=core, 2=stc lesioned, 3=vol lesioned)
                % 9th: prior for lesion models
          % 2) o: transformed outcomes ntrials x ncues array (column 1 = p(loss|blue cue), column 2 = p(loss|red cue)) 
          % 3) a: participant choices (binary) ntrials x 1 array
          % 4) configuration structure

% Output: negative log-likelihood    

% ------

% unpack parameters array
lambda_s    = x(1);           % lambda stc
lambda_v    = x(2);           % lambda vol
beta        = x(3); 
lambda      = x(4);
s0          = parameters(4);
v0          = parameters(5);
nparticles  = parameters(6);
x0_unc      = parameters(7);

z_rng = [v0 v0].^-1;
y_rng = [s0 s0].^-1;

state_model = @(particles)pf_state_transition(particles, lambda_s,lambda_v);
measurement_model = @pf_measurement_likelihood;

pf          = pf_initialize(state_model, measurement_model, nparticles, [z_rng; y_rng]);
N           = length(o);
estimated   = nan(N,2);
val         = nan(N,2); % two cues/options 
unc         = nan(N,1);
lr          = nan(N,1);
m           = zeros(1,nparticles); 
w           = x0_unc*ones(1,nparticles);
loglik      = 0; % init ll
state       = config.state;

%------                                ------%
%---- RUN THE ACTUAL RB PARTICLE FILTER -----%
%------                                ------%

% run the actual model. This starts with the particle filter:
% step 1: prediction stage 
% step 2: correction of prediction based on actual observation 
% step 3: run the kalman filter 
% particles: row 1: volatility 
% particles: row 2: stochasticity 

% loop over outcomes ( for both cues/options)
for t=1:size(o(:,1),1)   

    % loop over cues
    for c = 1:size(o,2)   
        
        pf              = pf_predict(pf);
        estimated(t,:)  = pf.weights*pf.particles';
        
        [pf, idx]       = pf_correct(pf,o(t,c),m,w);
        [m,w,k]         = kalman(pf.particles,o(t,c),m(idx),w(idx));
        val(t,c)        = pf.weights*m';
        unc(t)          = pf.weights*w';
        lr(t)           = pf.weights*k';

    end % end of cues/options loop

    % now run the softmax function and compute ll
    mu          = state(t,:);                       % probabilities for this trial
    p           = exp(beta*val(t,:)) / sum(exp(beta*val(t,:)));
    cprobs(t,:) = p;
    choice      = a(t); % what action did the participant take on that trial?

    % choose action prob based on action
    if choice == 1
        cprob(t)    = p(1);
    elseif choice == 2
        cprob(t)    = p(2);
    end

    % update log-likelihood
    loglik          = loglik + (log(cprob(t)));

end % end of trials loop

NLL                 = -loglik;

% Implement Ridge Regularization:
% calculating Ridge Penalty - only include parameters you're optimizing 
ridge_penalty       = lambda * (lambda_s^2 + lambda_v^2 + beta^2);

% Update NLL to include Ridge penalty
NLL_reg             = NLL + ridge_penalty;

vol                 = estimated(:,1).^-1;
unp                 = estimated(:,2).^-1;

end

%------------------------------
function particles = pf_state_transition(particles, lambda_s,lambda_v)

% what are the free parameters??
% lambda_s = x(1);
% lambda_v = x(2);

% estimate z (the inverse of volatility) 
z = particles(1,:);
eta = 1-lambda_v;
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(z)) + eps; % rescaled beta distr to account for indipendent random noise in the 
e = (eta.^-1)*epsil;
z = z.*e;

% estimate s (the inverse of stochasticity)
y = particles(2,:);
eta = 1-lambda_s; % eta is between 0 and 1 
nu = .5/(1-eta);
epsil = betarnd(eta*nu,(1-eta)*nu, size(y)) + eps; % random variable 
e = (eta.^-1)*epsil;
y = y.*e;

particles = [z; y];

end

function likelihood = pf_measurement_likelihood(particles, measurement, m, w)

z = particles(1,:);
y = particles(2,:);
v = z.^-1;
u = y.^-1;

likelihood = normpdf(measurement,m,sqrt(w+v+u));

end

function [m, w, k]=kalman(particles,outcome,m,w)
z = particles(1,:);
y = particles(2,:);
v = z.^-1;
u = y.^-1;

k = (w+v)./(w+v+u);
m = m + k.*(outcome-m);
% w = (1-k).*(w+v);

w = u./(w+v+u).*(w+v);

end