function [mc, nc, tx] = mainf(simy,simx)

models = {'hgf','vkf'};

for i=1:length(models)
    [ckv(:,i), tx(i,:)] = run(models{i},simy,simx); %#ok<AGROW>
end

ckv = fisher(ckv);
mc = invfisher(mean(ckv));
nc = mean(ckv<0);

end