% Code by Jill O'Reilly <jill.oreilly@ndcn.ox.ac.uk >
% Oxford Centre for Functional MRI of the Brain (FMRIB)
% Oxford University 
% version: 11 August 2015
% ------------------------------------------------------------------------

%% Section 1: HHTHT
clear;
qVals = 0:0.01:1; %possible values of q

% work out probability of data given q for each value of q
% observations are HHTHT
pData_given_q = qVals.*qVals.*(1-qVals).*qVals.*(1-qVals);

% plot probability
figure(1); 
plot(qVals,pData_given_q,'b');
xlabel('value of q','FontSize',14);
ylabel('p(HHTHT) given q','FontSize',14);
set(gcf,'Color','w');

%% Section 2: 50 coin tosses
% p(data|q) = qxqxqx.....xq x (1-q)x(1-q)...(1-q) 
%           = q^30 (30 heads) x (1-q)^20 (20 tails)
clear;
qVals = 0:0.01:1; %possible values of q
pData_given_q = qVals.^30 .* (1-qVals).^20;

% plot new likelihood function on top of old graph
figure(2); hold on;
plot(qVals,pData_given_q,'k');
xlabel('value of q','FontSize',14);
ylabel('p(30xH, 20xT) given q','FontSize',14);
set(gcf,'Color','w');


%% Section 3: Sequential learning
% Heads are 1's, tails are 0'1
clear;
% data = [1 1 0 1 0 ]; %HHTHT
data = [1 0 0 1 1 ]; %HHTHT

qVals=[0:0.01:1]'; % candidate values of q
prior(:,1) = ones(size(qVals)); prior(:,1)=prior./sum(prior(:,1));

for i=1:length(data)
    if data(i)==1
        L(:,i)=(0:0.01:1)'; % Likelihood function p(q|head) for values of q 0-1 in steps of 0.01
    elseif data(i)==0
        L(:,i)=(1:-0.01:0)'; % Likelihood function p(q|tail) for values of q 0-1 in steps of 0.01
    end;
    posterior(:,i)=L(:,i).*prior(:,i);
    posterior(:,i)=posterior(:,i)./sum(posterior(:,i)); % normalise so probabilities add up to 1
    prior(:,i+1)=posterior(:,i);
    
    % plot the posterior
    subplot(ceil(length(data)/5),5,i);
    plot(qVals,posterior(:,i),'k'); set(gca,'YLim',[0 0.03]); set(gca,'yTick',[]);
    xlabel('q','FontSize',14); ylabel(['p(q|y_{1:' int2str(i) '})'],'FontSize',14);
    if data(i)==1
        title(['Trial ' int2str(i) ': HEADS'],'FontSize',14);
    else
        title(['Trial ' int2str(i) ': TAILS'],'FontSize',14);
    end
    set(gcf,'Color','w');
end;



% Heads are 1s, tails are 0s - 50 trials = 30 heads, 20 tails



data = [ 1 0 1 0 1 0 1 1 0 1 1 1 1 1 0 1 1 1 0 1 0 0 0 0 1 1 0 0 1 1 1 0 1 1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 1 1 ];
