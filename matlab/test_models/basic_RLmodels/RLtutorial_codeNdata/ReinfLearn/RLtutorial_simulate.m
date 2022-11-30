function [data pout]= RLtutorial_simulate(params,sID)
% FUNCTION [data pout] = RLtutorial_simulate(params,sID)
% 
% params = vector of parameter values for alpha and beta
% sID = subject number in [1 - 999] range
% 
% The learning model used to generate choices is a simple Rescorla-Wagner
% (Rescorla & Wagner 1972) function, using a softmax
% choice/link/observation function:  
%   learning model (RW):     vA <- vA + alpha*(r-vA)
%   choice function (softmax):    p(A) = exp(beta*vA)/(exp(beta*vA)_exp(beta*vB))
% 
% feedback sequence is generated from the reversal learning task calling
% the revlParams function
% 
% ------------------------------------------------------------------------
% Written by Hanneke den Ouden 2015 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------


% get parameter values
alpha   = params(1);% sensitivity to reward and punishment
beta    = params(2); %learning rate

% Define initial value of each stimulus.
v0      = .5*ones(1,2); % initial value 
% v0      = [0.8 0.5];
v       = v0; % initialise v

% initialise output structure
data = {};

% get the task information
data.prep   = revlParams(sID,false);
outcome     = data.prep.feedback;
nt          = data.prep.nt;
nc          = data.prep.nStim; % 2 choice options

% initialise choice and outcome vectors to store the data.
loglik      = 0;
data.choice = nan(nt,1);
data.outcome= nan(nt,1);
VV          = nan(nt,nc); % matrix to store value over trials
PP          = nan(size(VV)); % matrix to store choice probability over trials

for t = 1:nt
    
    % Compute likelihood of the each choice option 
    ev  = exp(beta*v);  % exponentiate the value
    sev = sum(ev);		% compute the sum of the values
    p   = ev/sev; 		% probability each choice i.e. ratio of each of the values and their sum. 
                        % This code effectively does: 
                        % p(1) = ev(1)/sev;
                        % p(2) = ev(2)/sev;

    % Do a weighted coinflip to make a choice: choose stim 1 if random
    % number is in the [0 p(1)] interval, and 2 otherwise
    if rand(1)<p(1)
        c = 1;
    else
        c = 2;
    end
    r = outcome(t,c); % select the outcome for this choice
    
    % update the log likelihood with the p(choice|model):
    % note that this line can be done more efficiently after the trial
    % loop, see below
    loglik = loglik+log(p(c));
    
    % store value of V and choice probability
    VV(t,:) = v; % store value of V
    PP(t,:) = p;
    
    % update values
    dv  = r-v(c); % compute prediction error
    v(c)= v(c) + alpha*dv;    % update value
    
    data.choice(t) = c;
end

% store the feedback the subject received
c1 = data.choice==1;
c2 = data.choice==2;
data.outcome(c1) = outcome(c1,1);
data.outcome(c2) = outcome(c2,2);

% alternative, more efficient way of computing the likelihood
% loglik = sum(log(PP(c1,1)))+sum(log(PP(c2,2)));


pout = struct('VV',VV,'PP',PP, 'loglik',loglik,'params',params,'data',data);
end


