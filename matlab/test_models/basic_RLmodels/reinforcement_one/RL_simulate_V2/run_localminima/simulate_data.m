function [stim, update, choice, rew, setsize] = simulate_data(realalpha,realbeta,realrho,realK)

b = 0;
t = 0;

% 3 iterations
for rep = 1:3

    % of blocks of set sizes 2 through 6
    for ns = 2:6

        b           = b+1;
        update(t+1) = 1;

        % WM weight
        w           = realrho*(min(1,realK/ns));

        % initialize RL and WM
        Q           = (1/3)+zeros(ns,3);
        WM          = (1/3)+zeros(ns,3);
        trials      = repmat(1:ns,1,15); % starts with 30 trials, increments by 15 for each ns(30, 45, 60, 75, 90)

        for s = trials
            t           = t+1;
            stim(t)     = s; % 3 options (can be either 1,2 or 3)
            setsize(t)  = ns; % set-size can be 2 to 6

            % RL policy
            softmax1    = exp(realbeta*Q(s,:))/sum(exp(realbeta*Q(s,:)));
            
            % WM policy (high beta=50 captures perfect 1-trial memory)
            softmax2    = exp(50*WM(s,:))/sum(exp(50*WM(s,:)));
            
            % mixture policy
            pr          = (1-w)*softmax1 + w*softmax2;
            
            % make choice stochastically
            r           = rand;
            if r<pr(1)
                choice(t) = 1;
            elseif r<pr(1)+pr(2)
                choice(t) = 2;
            else
                choice(t) = 3;
            end

            % feedback
            rew(t)          = choice(t)==(rem(s,3)+1);

            % RL learning
            Q(s,choice(t))  = Q(s,choice(t))+realalpha*(rew(t)-Q(s,choice(t)));

            % WM update
            WM(s,choice(t)) = rew(t);
        end
    end
end

update(t)                   =0;

end