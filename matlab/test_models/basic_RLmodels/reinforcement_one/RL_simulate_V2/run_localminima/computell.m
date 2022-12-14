function llh = computell(p,K,stim,update,choice,rew,setsize)

global ppath;
ppath   = [ppath;p];
rho     = p(1);
alpha   = p(2);
beta    = 50*p(3);
l       = 0;

for t = 1:length(stim)
    s               = stim(t);
    
    if update(t)
        Q           = (1/3)+zeros(setsize(t),3);
        WM          = (1/3)+zeros(setsize(t),3);
    end
    
    w               = rho*(min(1,K/setsize(t)));
    softmax1        = exp(beta*Q(s,:))/sum(exp(beta*Q(s,:)));
    softmax2        = exp(50*WM(s,:))/sum(exp(50*WM(s,:)));
    pr              = (1-w)*softmax1 + w*softmax2;
    l               = l+log(pr(choice(t)));
    Q(s,choice(t))  = Q(s,choice(t))+alpha*(rew(t)-Q(s,choice(t)));
    WM(s,choice(t)) = rew(t);
    
end

llh                 = -l;

end