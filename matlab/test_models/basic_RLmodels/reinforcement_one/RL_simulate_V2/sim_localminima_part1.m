% PART 1 OF FINDING THE LOCAL MINIMA 

% MORE COMMENTS SHOULD GO HERE

clear all 
clc

%% initialise parameters 

% first define a range of learning rates
alphas      = [0.05:.05:1]; % parameter 1

% define a range of softmax parameter values to test
betas       = [1 5:5:50]; % parameter 2

% define a range of wm reliance to test
rhos        = [0:.05:1]; % parameter 3

% define a range of capacities to test
ks          = 2:6;

% define simulation parameters
realalpha   = .1;
realbeta    = 10;
realrho     = .9;
realk       = 4;

%% simulate the wm task

b = 0;
t = 0;

for rep = 1:3 % three repetitions for each set size
    
    for ns = 2:6 % block of set size
        
        b           = b + 1; % update b
        update(t+1) = 1;
        w           = realrho * (min(1, realk/ns)); % initialise wm mixture weights
        
        % init RL and wm agents
        Q           = .5 + zeros(ns,3);
        wm          = .5 + zeros(ns,3);
        trials      = repmat(1:ns,1,15); %define a sequence of trials with 15 iterations of each stimulus
        
        % loop over trials
        for trl = trials
            
            t = t + 1; % update t
            
            % rl policy
            softmax1 = exp(realbeta*Q(trl,:))/sum(exp(realbeta*Q(trl,:)));
            
            % wm policy 
            softmax2 = exp(50*wm(trl,:))/sum(exp(50*wm(trl,:)));
            
            % mixture policy
            pr = (1-w)*softmax1 + w*softmax2;
            
            % action choice
            r = rand;
            
            if r<pr(1)
                choice(t)=1;
            elseif r<pr(1)+pr(2)
                choice(t)=2;
            else
                choice(t)=3;
            end
            
            % reward correct action
            rew(t) = choice(t) == (rem(trl,3)+1);
            
            % update rl wm agents
            Q(trl,choice(t))=Q(trl,choice(t))+realalpha*(rew(t)-Q(trl,choice(t)));
            wm(trl,choice(t))=rew(t);
            
            % store information
            stim(t) = trl;
            setsize(t) = ns;

        end % end of trials loop
    end % end of ns loop
end % end of reps loop

% check that it worked by making sure that performance in higher set sizes
% is lower than in high set sizes
update(t) = 0;
for ns = 2:6
    [ns mean(rew(setsize==ns))]
end

%% compute ll for multiple parameters 

i1 = 0;

% loop over alphas
for alpha = alphas 
    
    i1                              = i1+1; % update i1
    i2                              = 0;
    
    % loop over betas now
    for beta = betas
        
        i2                          = i2+1; % update i2
        j1                          = 0;
        
        % loop over wm reliance
        for rho = rhos
            
            j1                      = j1+1; % update j1
            j2                      = 0;
            
            for k = ks
                
                j2                  = j2+1; % update j2
                l                   = 0; 
                
                for t = 1:length(stim)
                    
                    s = stim(t);
                    if update(t)
                        
                        Q           = 0.5 + zeros(setsize(t),3);
                        wm2         = 0.5 + zeros(setsize(t),3);
                        
                    end
                    
                    w               = rho*(min(1,k/setsize(t))); % wm mixture weight
                    softmax1        = exp(beta*Q(s,:))/sum(exp(beta*Q(s,:)));
                    softmax2        = exp(50*wm(s,:))/sum(exp(50*wm(s,:)));
                    pr              = (1-w)*softmax1 + w*softmax2; % mixture policy
                    l               = l+log(pr(choice(t)));
                    Q(s,choice(t))  = Q(s,choice(t))+alpha*(rew(t)-Q(s,choice(t)));
                    wm(s,choice(t)) = rew(t);
                    
                end
                
                llh(i1,i2,j1,j2)=l;
                
                
            end
  
        end

    end
end

%% plot ll projected on dimensions alpha and rho

figure;

% project
llh2    = squeeze(max(max(llh,[],4),[],2));

% take out extreme values, because they make plotting less visually
% readable
llh2    = llh2(1:end-1,1:end-1);

% find min and max extrema
mi      = min(llh2(:));
ma      = max(llh2(:));
x       = repmat(1:length(alphas)-1,length(rhos)-1,1)';
y       = repmat(1:length(rhos)-1,length(alphas)-1,1);
[~,i]   = max(llh2(:));

% plot the surface
imagesc(alphas(1:end-1),rhos(1:end-1),llh2',[mi,ma])
colorbar
hold on
plot(alphas(x(i)),rhos(y(i)),'ok')
plot(realalpha,realrho,'xr')
xlabel('alpha')
ylabel('rho')
    set(gca,'fontsize',16)


    %% plot 1d versions of the likelihood

ps{1} = alphas; na{1} = 'alpha';
ps{2} = betas;  na{2} = 'beta';
ps{3} = rhos;   na{3} = 'rho';
ps{4} = ks;     na{4} = 'K';

figure

for i = 1%:4
    out     = setdiff(1:4,i);
    llh1    = squeeze(max(max(max(llh,[],out(3)),[],out(2)),[],out(1)));
    [v,w]   = max(llh1);
    
    %subplot(1,4,i)
    plot(ps{i},llh1,'o-','linewidth',2)
    hold on
    plot(ps{i}(w),v,'rx','linewidth',2)
    set(gca,'fontsize',16)
    title(na{i})
    ylabel(llh)
end



