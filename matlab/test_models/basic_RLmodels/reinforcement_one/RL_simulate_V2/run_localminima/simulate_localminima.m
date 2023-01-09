% this script is generated based on the paper: 10 simple rules for the
% computational modelling of behavioural data

% in this script I will be running local/global minima for the simulated WM
% task

% 1. I will first simulate some data for the task (one dataset)
% 2. then I'll fit the free parameters (alpha, rho, beta) using the simulated
% dataset with 10 different starting points (for the 3 params)
% 3. and find the (global) best fit of the 10 (the minimum ll)
% 4. the next step is to brute force compute maximum ll (using the range of different parameter values) and plot 
% 5. simulate 100 datasets, for each dataset run step 2 and 3
% 6. and plot the results from step 5

% the 4 plots:
% 1. plot the brute force max lls using one simulated dataset
% 2. plot best ll of the 100 simulated datasets (across the 10 startings
% points for each dataset)
% 3. plot the distance to global best params 
% 4. 



%% define initial parameters 

% define a banch of different parameters values:
alphas      = [0.06:0.01:0.5];  % range of alphas: 0.06 to 0.5
betas       = [1 4:2:20];       % range of betas: from 1 to 20
rhos        = [0.5:0.01:0.98];  % wm weight ranges: from 0.5 to 0.98
Ks          = 1:6;              % capacity: ranges from 1 to 6

% define real simulation parameters (like initial parameters, that will be used to simulate data)
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
    pars(init,:) = [pval,fval]; % pval stands for parameter value (for the three parameters), fvall is the function value (ll)
    
end

%% find the global best

[mf, i] = min(pars(:,end)); % get minimum of all fvals (last column)

pars = pars(i,:); % optimised rho, alpha, beta for minimum ll

%% brute force fitting 

% for this we'll use the range of alphas, rhos and betas specified in the
% beginning (the parameter space) 

% this is the simplest approach fot finding the maximum likelihood 
i1 = 0;

for alpha = alphas
    
    i1      = i1+1; % increment alpha
    i2      = 0;
    
    for beta = betas
        
        i2  = i2+1; % increment beta
        j1  = 0;
        
        for rho = rhos
            
            j1 = j1+1; % increment rho
            j2 = 0;
            
            for K = realK
                
                j2  = j2+1; % increment k 
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


