function [predictions, signals] = vkf_v1(data,model_params)

% first version of the volatile kalman filter model
% adapted from Piray & Daw (2020); Piray et al. (2019)

% From the Piray & Daw (2020) paper, equations 14 to 19 are used

% created: 15/02/2023

%--------

% INPUTS: 
%           - data -- structure that contains simulation data(feedback/outcome is needed
%           - model_params - structure with parameters needed to run the model: lambda (volatility learning rate), initial
%             volatility, variance, omega (noise parameter)

% OUTPUTS:
%           - predictions: predicted state
%           - signals: predictions, voalitity, learning rate, prediction
%             error, volatility prediction error

%%% INFO ABOUT PARAMETER VALUES:
% lAMBDA:               SHOULD BE BETWEEN 0 AND 1 
% OMEGA:                SHOULD BE POSITIVE (>0) ---- NOTE: OMEGA IS ASSUMED TO BE INITIAL VARIANCE HERE
% INITIAL VOLATILITY:   SHOULD BE POSITIVE (>0)

% output depends on initial variance (omega)

% ---------------------------------

%% unpack params and init variables

% unpack parameter values 
lambda      = model_params.lambda;
v0          = model_params.init_vol; % initial volatility
omega       = model_params.omega;

outcomes    = data; % 2 columns 

if lambda <= 0 || lambda >= 1 % should be between 0 and 1 (high lambda = fast diffusion)
    error('lambda should be in the unit range');
end
if omega <= 0 % should be greater than 0
    error('omega should be positive');
end
if v0 <= 0 % should be greater than 0
    error('v0 should be positive');
end  

[ntrials, ncues]        = size(outcomes);       % T: number of trials, C: number of cues
w0                      = omega;                % initial noise param
m                       = zeros(1,ncues);       % initial mean of the posterior
w                       = w0*ones(1,ncues);
v                       = v0*ones(1,ncues);

predictions             = nan(ntrials,ncues);   % predicted state of the underlying probabilities 
learning_rate           = nan(ntrials,ncues);
volatility              = nan(ntrials,ncues);
prediction_error        = nan(ntrials,ncues);   % difference between expected and true value 
volatility_error        = nan(ntrials,ncues);

sigmoid                 = @(x)1./(1+exp(-x));   % sigmoid function to map normally distributed X to the unit range

% loop over trials and update parameter values
for t = 1:ntrials

    o                       = outcomes(t,:);        % this-trial outcome
    predictions(t,:)        = m;                    % posterion mean 
    volatility(t,:)         = v;    
    mpre                    = m;
    wpre                    = w;

    delta_m                 = o - sigmoid(m);       % delta rule or prediction error update 
    k                       = (w+v)./(w+v+ omega);  % kalman gain or learning update rate 
    alpha                   = sqrt(w+v);            % speed of learning
    m                       = m + alpha.*delta_m;
    w                       = (1-k).*(w+v);         % posterior variance 
    wcov                    = (1-k).*wpre;  
    delta_v                 = (m-mpre).^2 + w + wpre - 2*wcov - v;  % volatility error 
    v                       = v +lambda.*delta_v;                   % volatility

    learning_rate(t,:)      = alpha;
    prediction_error(t,:)   = delta_m;
    volatility_error(t,:)   = delta_v;    

end % end of trials loop

signals = struct('predictions',predictions,'volatility',volatility,'learning_rate',learning_rate,...
                 'prediction_error',prediction_error,'volatility_prediction_error',volatility_error);


end % end of function 