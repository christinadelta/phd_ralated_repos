% softmax function (observation equation) simple demonstration 

% softmax function is used to convert values associated with stimuli into
% probability choices 

% this short script is used to show the effect of the beta value 

% the code computes choice probabilities based on the value difference
% between A & B using the softmax function:

%   p(A) = exp(beta*vA)/(exp(beta*vA)+exp(beta*vB))

% where, vA is 1-vB
% the code then generates a choice based on theses probabilities 

% the beta value in the softmax function controls the level of
% stochasticity (or randomness) in the choice. 
% when beta = 0, the choice is completely random, 
% when choice beta = inf, choice is deterministic (always choosing the
% highest value option)

%% -------------------

clear all 
clc

% define a range of betas
betaArray = [0 9 10];

if any(abs(betaArray)>700)
    error('Please pick beta values between -700 and +700, as matlab cannot deal with such small/big numbers!');
end

% init figure
h = figure; set(h,'position',[10 60 400 400 ],'Color','w'); 
hold on; box off; 

va      = 0:0.02:1;
vb      = 1-va; % value of B is [1 - value of A]
pa      = nan(length(betaArray),length(va));
data    = nan(length(betaArray),length(va));


ct      = 0;
x       = .02;

for beta = betaArray

    ct          = ct+1;
    pa(ct,:)    = exp(beta*va)./(exp(beta*va)+exp(beta*vb)); % probability of choosing A
    leg{ct}     = sprintf('beta = %d',beta);
    tmp         = double(rand(1,length(va))< pa(ct,:))*(1+x*ct); % generate choices based on this data
    tmp(tmp==0) = tmp(tmp==0)-x*ct;
    data(ct,:)  = tmp;
end

%% visualise softmax function

plot(va-vb,pa);
legend(leg,'location','best'); legend boxoff;

plot(va-vb,data,'*');
xlabel('Value(A) - Value (B)');
ylabel('p(choice = A)');
title(sprintf('softmax function'));
ylim([-x*ct 1+x*ct]);


