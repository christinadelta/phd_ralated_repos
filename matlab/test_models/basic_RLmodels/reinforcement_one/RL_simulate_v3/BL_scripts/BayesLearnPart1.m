 % Bayesian Learning - Part 1 -- Likelihoods and uncertainty 

 % Created: 11/2/2023

 % 1st part of a simple bayesian learning model:
 % The script goes through likelihood and uncertainty using heads/tails as
 % example 

 %%%%%%%%% SOME GENERAL INFORMATION: %%%%%%%%%%
 % Let's say that probability of vertical gabor -- p(vertical) = q
 % and probability of horizontal gabor -- p(horizontal) = 1-q
 % since each observation is independent, we can compute the overall probability
 % as:
 % p(vvhhv) = q^2*(1-q)^2*q

 % If q = 0.6 then: 
 % p(vvhhv) = 0.6^2*0.4^2*0.6 = 0.0346

 % Something to remember:
 % there is a mathematical relationship between:
 % 1. the probability of the data given the parameters of the environment: p(vvhhv|q)
 % 2. the likelihood of the parameters, given the data observed: p(q|vvhhv)

 % so...: p(q|data) = p(data|q)

 % 

 %% section 1: vvhhv

 clear all
 clc

 qVals          = 0:0.01:1; % values of q in steps of 0.01

 % compute the probability of data given q (for each value of q)
 pdata_givrn_q  = qVals.^2 .* (1-qVals).^2 .*qVals;

 % plot the likelihood function
figure(1);
plot(qVals,pdata_givrn_q, 'b');
xlabel('value of q','FontSize',14);
ylabel('p(VVHHV) given q','FontSize',14);
set(gcf,'Color','w');

%% section 2: let's look at 50 observations 

% p(data|q) = qxqxqx.....xq x (1-q)x(1-q)...(1-q) 
%           = q^30 (30 vertical) x (1-q)^20 (20 horizontal)

clear;
qVals           = 0:0.01:1; %possible values of q
pdata_given_q   = qVals.^30 .* (1-qVals).^20;

% plot new likelihood function on top of old graph
figure(2); hold on;
plot(qVals,pdata_given_q,'k');
xlabel('value of q','FontSize',14);
ylabel('p(30xV, 20xH) given q','FontSize',14);
set(gcf,'Color','w');

%% section 3: sequential learning

% instead of revealing the entire sequence all at once, the observations
% occur one at a time. After each obsrevation is presented the observer UPDATES their
% belief and we can model this using Bayes theorem.
% Sequence of observations: y1,y2,y3....yi

%                               p(q|Y1:i) ∝ p(Yi|q) p(q|Y1:i-1)

% p(q|Y1:i) is the posterior probability of some value of q given all the observed data on trials 1-i
% p(Yi|q) is the likelihood function for some value of q given most recent observation
% p(q|Y1:i-1) is the prior probability of some value of q given all the
% observations up until the most recent ones 

clear;
% data = [1 1 0 1 0 ]; %VVHVH
% data = [ 1 0 1 0 1 0 1 1 0 1 1 1 1 1 0 1 1 1 0 1 0 0 0 0 1 1 0 0 1 1 1 0 1 1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 1 1 ];

data        = [1 0 0 1 1 ]; %VHHVV
qVals       = [0:0.01:1]'; % candidate values of q
prior(:,1)  = ones(size(qVals)); prior(:,1) = prior./sum(prior(:,1));

for i = 1:length(data) % loop over trials

    if data(i) == 1
        L(:,i)      = (0:0.01:1)'; % Likelihood function p(q|vertical) for values of q 0-1 in steps of 0.01
    elseif data(i) == 0
        L(:,i)      = (1:-0.01:0)'; % Likelihood function p(q|horizontal) for values of q 0-1 in steps of 0.01
    end

    posterior(:,i)  = L(:,i).*prior(:,i);
    posterior(:,i)  = posterior(:,i)./sum(posterior(:,i)); % normalise so probabilities add up to 1
    prior(:,i+1)    = posterior(:,i);

    % plot the posterior
    subplot(ceil(length(data)/5),5,i);
    plot(qVals,posterior(:,i),'k'); set(gca,'YLim',[0 0.03]); set(gca,'yTick',[]);
    
    xlabel('q','FontSize',14); ylabel(['p(q|y_{1:' int2str(i) '})'],'FontSize',14);
    if data(i) == 1
        title(['Trial ' int2str(i) ': Vertical'],'FontSize',14);
    else
        title(['Trial ' int2str(i) ': Horizontal'],'FontSize',14);
    end
    set(gcf,'Color','w');

end % end of trials loop


