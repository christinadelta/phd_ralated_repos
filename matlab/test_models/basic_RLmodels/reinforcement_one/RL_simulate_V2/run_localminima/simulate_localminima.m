% running local minima part two

% still I have no idea what that does

clear all
clc

%% define parameters 

% define a banch of initial parameter values for alpha, rho and beta:
alphas      = [.06:.01:.5]; % learning rate
betas       = [1 4:2:20];   % inverse temperature
rhos        = [.5:.01:.98]; % WM memory weight
Ks          = 2:6;          % capacity

% also define real simulation parametrs 
realalpha   = .1;
realbeta    = 8;
realrho     = .9;
realK       = 4;

% simulate one dataset
[stim, update, choice, rew, setsize] = simulate_data(realalpha,realbeta,realrho,realK);

%% fit simulated data using fmincon

% set fmincon options
options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%
        
% run optimization over 10 starting points
for init = 1:10
    
    % random starting point
    x0 = rand(1,3); % (rho, alpha, beta)
    
    % optimize
    [pval,fval,bla,bla2] = fmincon(@(x) computell(x,realK,stim,update,choice,rew,setsize),x0,[],[],[],[],...
        [0 0 0], [1 1 1], [],options); %
    
    % fval that is returned is the function value (the ll value) 
    
    % store optimization result
    pars(init,:) = [pval,fval];
    
end

%% find the global best

[mf, i] = min(pars(:,end)); % get minimum of all fvals (last column)

pars = pars(i,:); % optimised rho, alpha, beta for minimum ll

%% brute force fitting

i1 = 0;
for alpha = alphas
    
    i1      = i1+1;
    i2      = 0;
    
    for beta = betas
        
        i2  = i2+1;
        j1  = 0;
        
        for rho = rhos
            
            j1 = j1+1;
            j2 = 0;
            
            for K = realK
                
                j2  = j2+1;
                p   = [rho,alpha,beta/50];
                
                % store likelihood over parameters in a mesh
                llh(i1,i2,j1,j2)=-computell(p,K,stim,update,choice,rew,setsize);
            end
        end
    end
end

%% plot the results - in 2d

figure;
subplot(2,2,1)

llh2    = squeeze(max(squeeze(llh),[],2));
mi      = min(llh2(:));
ma      = max(llh2(:));
x       = repmat(1:length(alphas),length(rhos),1)';
y       = repmat(1:length(rhos),length(alphas),1);
[mb,i]  = max(llh2(:));
imagesc(alphas(1:end),rhos(1:end),llh2',[mi,ma])

colorbar
hold on
plot(alphas(x(i)),rhos(y(i)),'ok')
plot(realalpha,realrho,'xr')
plot(pars(2),pars(1),'*k')
xlabel('alpha')
ylabel('rho')
    set(gca,'fontsize',16)
    
%% iterate simulation and fitting
    
options = optimset('MaxFunEval',100000,'Display','off','algorithm','active-set');%

% number of random starting points for optimizer
ninitialpoints = 10;

% for 100 simulations
for iter = 1:100
    
    disp(['simulation #',num2str(iter)])
    % generate data
    [stim,update,choice,rew,setsize] = simulate_data(realalpha,realbeta,realrho,realK);
    pars = [];
    
    % fit simulated data with ninitialpoints random starting points
    for init=1:ninitialpoints
        
        x0                      = rand(1,3);
        [pval,fval,bla,bla2]    = fmincon(@(x) computell(x,realK,stim,update,choice,rew,setsize),x0,[],[],[],[],...
            [0 0 0],[1 1 1],[],options);
        
        pars(init,:)            = [pval,fval];
        [m,i]                   = min(pars(:,end));
        bestllh(iter,init)      = m;
        bestpars(iter,init,:)   = pars(i,1:end-1);
    end
    
    % find global best fit
    [mf,i]          = min(pars(:,end));
    % find at which random starting point it was found
    when(iter,1)    = i;
    
    % find at which random starting point a likelihood within .01 of the
    % global best was found
    i               = find(bestllh(iter,:)<bestllh(iter,end)+.01);
    when(iter,2)    = i(1);
    
    % find at which random starting point a likelihood within .1 of the
    % global best was found
    i               = find(bestllh(iter,:)<bestllh(iter,end)+.1);
    when(iter,3)    = i(1);
    
end

% compute what the best log-likelihood found was up to random starting
% point i, substracting the final best log-likelihood (putative global
% best)
bestllh = bestllh(:,1:end-1)-repmat(bestllh(:,end),[1,ninitialpoints-1]);

subplot(2,2,2)
errorbar(mean(bestllh),std(bestllh)/sqrt(iter),'linewidth',1)
set(gca,'fontsize',14)
xlabel('starting point iteration')
ylabel('local-global best nlh')

% compute distance to "gloabl" best parameters as a function of optimizer
% iteration over random starting points.
bestpars    = bestpars(:,1:end-1,:)-repmat(bestpars(:,end,:),[1,ninitialpoints-1,1]);
bestpars    = sum(bestpars.^2,3);
subplot(2,2,3)
errorbar(mean(bestpars),std(bestpars)/sqrt(iter),'linewidth',1)
set(gca,'fontsize',14)
xlabel('starting point iteration')
ylabel('d(local-global best param)')

% plot when
subplot(2,2,4)
hold on

for j=1:2
    plot(sort(when(:,j)),'o-','linewidth',1)
end

set(gca,'fontsize',14)
legend('global = best','global = |llh-best|<.01')
ylabel('iteration where global llh first reached')
xlabel('sorted simulation number')
