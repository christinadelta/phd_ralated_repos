function [signals]  = kf_v1_1(params)

% comments will go here

%% init params 

% upack params
y                   = params.y;             % outcomes 
omega               = params.sigma;         % noise used to estimate the kalman gain (constant variance - fixed param) 
lambda              = params.lambda_v; 
v0                  = params.v0;

%% 
if lambda <= 0 || lambda >= 1
    error('lambda should be in the unit range');
end
if omega <= 0
    error('omega should be positive');
end
if v0 <= 0
    error('v0 should be positive');
end  

%%
% init posterior variance 
w0                  = omega; 
[trls,cues]         = size(y);

m                   = zeros(1,cues);       % estimated posterior mean
w                   = w0*ones(1,cues);     % estimated posterior variance
v                   = v0*ones(1,cues);     % estimated volatility

p                   = nan(trls,cues); % predicted signal (posterior mean)
lr                  = nan(trls,cues); % learning rates
vol                 = nan(trls,cues); % volatility 
pe                  = nan(trls,cues); % prediction error 
ve                  = nan(trls,cues); % volatility error 

sigmoid             = @(x)1./(1+exp(-x));

%% update parameters on every trial

for t = 1:trls
    
    % update outcome, volatility and posterior
    o               = y(t,:);
    p(t,:)          = m;    
    vol(t,:)        = v;   

    mpre            = m;
    wpre            = w;

    delta_m         = o - sigmoid(m);                       % PE  
    k               = (w+v)./(w+v+ omega);                  % kalman gain
    alpha           = sqrt(w+v);                            % lr
    m               = m + alpha.*delta_m;                   % posterior 
    w               = (1-k).*(w+v);                         % posterior variance
    
    wcov            = (1-k).*wpre;                          % posterior covariance :)
    delta_v         = (m-mpre).^2 + w + wpre - 2*wcov - v;  % volatility error   
    v               = v +lambda.*delta_v;                   % volatility update

    % update paramters of interest
    lr(t,:)         = alpha;
    pe(t,:)         = delta_m;
    ve(t,:)         = delta_v;    

end 

signals = struct('predictions',p,'volatility',vol,'learning_rate',lr,...
                 'prediction_error',pe,'volatility_prediction_error',ve);

end 