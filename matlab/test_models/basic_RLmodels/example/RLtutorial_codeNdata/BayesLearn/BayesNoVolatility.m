% Code by Jill O'Reilly <jill.oreilly@ndcn.ox.ac.uk >
% Oxford Centre for Functional MRI of the Brain (FMRIB)
% Oxford University 
% version: 11 August 2015
% ------------------------------------------------------------------------

clear; close all; 

%% user options
H=1/25; %probability of side reversal


%% First load the data set that we will fit the model to
q=load('data/schedule_volatile.txt');
y=load('data/trials_volatile.txt');


%%
% set up the state space

q_space=[0.01:0.01:0.99]; % possible values for q - don't allow 0 or 1, as model may explode!
prior_p_q=NaN(length(y),length(q_space)); % initiate matrix for prior probability of each value of q: size = trials x possible values
post_p_q=NaN(length(y),length(q_space)); % initiate matrix for posterior probability of each value of q: size = trials x possible values

prior_p_q(1,:)=ones(1,length(q_space))./length(q_space); % on trial 1, prior probability for each value of q is uniform ie 1/number of possible values of q


%% go through data trial by trial updating probabilities
for i=1:length(y)
    
    if y(i)==1  % if target was on Left on this trial
          p_q_given_yi(i,:) = q_space;
    else
          p_q_given_yi(i,:) = 1-q_space;
    end;
    
    post_p_q(i,:) = p_q_given_yi(i,:) .* prior_p_q(i,:);
    post_p_q(i,:) = post_p_q(i,:)./sum(post_p_q(i,:),2); % normalise the posterior
    
    % set up the prior for next trial - there is a chance of (1-H) that 
    % q is the same on the next trial as on this trial, 
    % and a chance of H that it now has a random value between 0 and 1
    prior_p_q(i+1,:) = post_p_q(i,:)*(1-H) + H*ones(1,length(q_space))./length(q_space);
    
end

% Work out the model's best guess of pL on each trial by finding the
% expected value of q

est_q = sum(prior_p_q*q_space',2);

%% Plot the model's beliefs about q (p(orange rewarded)) on the first few trials
figure; hold on;

for c=1:5; subplot(1,5,c); hold on;
    t=1:5; % representative trials
    plot([est_q(t(c)) est_q(t(c))],[0 1],'r-.','MarkerSize',15);
    plot(q_space,prior_p_q(t(c),:));
    if y(t(c))==1
        title(['prior on trial ' int2str(t(c)) ': outcome = Orange']);
    else
        title(['prior on trial ' int2str(t(c)) ': outcome = Blue']);
    end
    xlabel('value of q'); ylabel('p(value of q)'); set(gcf,'color','w'); set(gca,'XLim',[0 1],'YLim',[0 0.025]);
    legend(['best estimate of q']);
end;

%% Plot the true value of q (p(orange rewarded)) and the model's estimate

figure; plot(q); hold on; plot(est_q,'r'); % plot the true values of s against the model's estimate
plot(y,'k.','MarkerSize',8); % plot the individual data points (left and right targets)
set(gca,'YLim',[-0.2 1.2]); legend('true value of q','model estimate','Data points (1=orange, 0=blue)');
xlabel('trial number','FontSize',16); ylabel('model estimate of q','FontSize',16); set(gcf,'color','w'); set(gca,'FontSize',16);



%% Plot the model's beliefs as a distribution on each trial
figure; hold on; 
imagesc(1:length(y),flipud(q_space),post_p_q'); % plot full posterior over all trials
plot(q,'w'); plot(est_q,'k');plot(est_q,'w--'); % add best estimate of p(left)
plot(0.1+y*0.8,'w.','MarkerSize',5); % add data points (Orange/Blue rewarded on each trial)
set(gcf,'color','w'); set(gca,'XLim',[0 length(y)],'FontSize',14);
xlabel('trial number','FontSize',14); ylabel('possible values of q','FontSize',14); 
title('colour represents probability of each possible value of q','FontSize',14);


%% Plot the model's beliefs about q (p(orange rewarded)) around a reversal
figure; hold on;

for c=1:15; subplot(3,5,c); hold on;
    t=30:44; % representative trials
    plot([est_q(t(c)) est_q(t(c))],[0 1],'r-.','MarkerSize',15);
    plot(q_space,prior_p_q(t(c),:));
    if y(t(c))==1
        title(['t' int2str(t(c)) ': Orange']);
    else
        title(['t' int2str(t(c)) ': Blue']);
    end
    xlabel('value of q'); ylabel('p(value of q)'); set(gcf,'color','w'); set(gca,'XLim',[0 1],'YLim',[0 0.05]);
end;


