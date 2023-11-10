function [predict, vol, lr, sigm]  = kf_v1(params)

% comments will go here

%% init params 

% upack params
y                   = params.y;             % outcomes 
convar              = params.sigma;         % noise used to estimate the kalman gain (constant variance - fixed param) 
lambda              = params.lambda_v; 
v0                  = params.v0;

[trls,cols]         = size(y);

m                   = 0*ones(1,cols);       % posterior mean
w                   = convar*ones(1,cols);   % posterior variance 
v                   = v0*ones(1,cols);      % volatility

predict             = nan(trls,cols);
lr                  = nan(trls,cols);
vol                 = nan(trls,cols);
sigm                = nan(trls,cols);
sigmoid             = @(x)1./(1+exp(-x));  % sigmoid function to deal with binary outcomes

%% loop over trials and estimate parameters 

for t = 1:trls
    
    % update predicted signal, lr, and volatility of trial t-1
    predict(t,:)    = m(1,:);
    lr(t,:)         = sqrt(w + v); % this will change when we add stochasticity 
    vol(t,:)        = v;
    wpre            = w; 

    % now estimate parameters at trial t
    delta           = sqrt(w+v).*(y(t,:) - sigmoid(m)); % we use the sigmoid function to transform the normally distributed psoterior mean into unit range
    m               = m + delta;                        % update mean
    k               = (w+v)./(w+v + convar);            % update the kalman gain
    w               = (1-k).*(w+v);                     % update posterior variance

    % update volatility acording to expected value of [m - m(t-1)] which is the delta value and lambda 
    v               = v +lambda.*(delta.^2 + k.*wpre - k.*v); 
    sigm(t,:)       = sigmoid(m);


end % end of trials loop 


end % end of fucntion