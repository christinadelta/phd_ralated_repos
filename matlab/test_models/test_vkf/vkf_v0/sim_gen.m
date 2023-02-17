function [Ys, Xs] = sim_gen(fsim,nsim)
if exist(fsim,'file')
    data = load(fsim); 
    Ys = data.y;
    Xs = data.x;    
    return;
end

n  = 20;
p0 = .8;

p  = [p0*ones(1,2) repmat([1-p0 p0],1,2) (1-p0)*ones(1,4)];
p  = [p p];

nb = length(p);  
N  = nb*n;

Xs = nan(N/2,nsim);
Ys = nan(N/2,nsim);
for j=1:nsim
    x  = nan(N,1);    
    y  = zeros(N,1);

    t0 = 0;
    for i=1:nb
        ii = t0 + (1:n);
        x(ii) = p(i);
        ni = randperm(n);
        ni = ni(1: round(p(i)*n));
        y(ii(ni))  = 1;
        t0 = t0 + n;
    end


    y = y(1:N/2,:);
    x = x(1:N/2,:);
    
    Xs(:,j) = x;
    Ys(:,j) = y;
end

end