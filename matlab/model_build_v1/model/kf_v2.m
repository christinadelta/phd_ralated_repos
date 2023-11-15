function [signals] = kf_v2(params)

% comments will go here

%% init params 

% upack params
y                   = params.y;                 % outcomes 
lambda_s            = params.lambda_s;          % noise used to estimate the kalman gain (constant variance - fixed param) 
lambda_v            = params.lambda_v; 
v0                  = params.v0;
s0                  = params.s0;

%% 
if lambda_v <= 0 || lambda_v >= 1
    error('lambda vol should be in the unit range');
end
if lambda_s <= 0 || lambda_s >= 1
    error('lambda stc should be in the unit range');
end
if s0 <= 0
    error('initial stc should be positive');
end
if v0 <= 0
    error('initial vol should be positive');
end  

%%
% init posterior variance 
w0                  = s0; 
[trls,cues]         = size(y);

m                   = zeros(1,cues);       % estimated posterior mean
w                   = w0*ones(1,cues);     % estimated posterior variance
v                   = v0*ones(1,cues);     % estimated volatility (dependent to process noise)
s                   = s0*ones(1,cues);     % estimated stochasticity (dependent to observation noise)

p                   = nan(trls,cues); % predicted signal (posterior mean)
lr                  = nan(trls,cues); % learning rates
vol                 = nan(trls,cues); % volatility 
stc                 = nan(trls,cues); % stochasticity 
pe                  = nan(trls,cues); % prediction error 
ve                  = nan(trls,cues); % volatility error
se                  = nan(trls,cues); % stochasticity error

sigmoid             = @(x)1./(1+exp(-x));

%% update parameters on every trial

for t = 1:trls
    
    % update outcome, volatility and posterior
    o               = y(t,:);
    p(t,:)          = m;    
    vol(t,:)        = v; 
    stc(t,:)        = s;

    mpre            = m;
    wpre            = w;

    delta_m         = o - sigmoid(m);                       % PE  
    k               = (w+v)./(w+v+s);                       % kalman gain
    %alpha           = sqrt(w+v+s);                          % lr
    alpha           = k;
    m               = m + alpha.*delta_m;                   % posterior 
%    w               = (1-k).*(w+v);                         % posterior variance
    w               = s./(w+v+s).*(w+v);
    
    wcov            = (1-k).*wpre;                          % posterior covariance :)
    delta_v         = (m-mpre).^2 + w + wpre - 2*wcov - v;  % volatility error  
    delta_s         = (m-mpre).^2 + w + wpre - 2*wcov - s;  % stochasticity error
    v               = v + lambda_v.*delta_v;                % volatility update
    s               = s + lambda_s.*delta_s;                % stochasticity update    

    % update paramters of interest
    lr(t,:)         = alpha;
    pe(t,:)         = delta_m;
    ve(t,:)         = delta_v;  
    se(t,:)         = delta_s;

end 

signals = struct('predictions',p,'volatility',vol,'stochasticity',stc,'learning_rate',lr,...
                 'prediction_error',pe,'volatility_prediction_error',ve, 'stochasticity_prediction_error',se);

end % end of function 