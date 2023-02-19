function [ckv, tx] = eval_model(fsim, model, cbm)
data = load(fsim);
Ys = data.y;

nsim = size(Ys,2);

fx = cbm.output.parameters;

ckv = nan(nsim,1);
for n=1:nsim
    y = Ys(:,n);
    [~, tx, cc]= model(fx,y);   % model_hgf.m and model_vkf.m                   
    ckv(n) = cc;        
end

end