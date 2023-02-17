function [predictions, signals] = vkf_v1(data,model_params)

% first version of the volatile kalman filter model
% adapted from Piray & Daw (2020); Piray et al. (2019)

% created: 16/02/2023 

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
v0          = model_params.init_vol;
omega       = model_params.omega;

outcomes    = data.feedback; % 2 columns 

if lambda <= 0 || lambda >= 1 % should be between 0 and 1
    error('lambda should be in the unit range');
end
if omega <= 0 % should be greater than 0
    error('omega should be positive');
end
if v0 <= 0 % should be greater than 0
    error('v0 should be positive');
end  

w0          = omega; 

[T,C]       = size(outcomes); % T: number of trials, C: number of cues









end % end of function 