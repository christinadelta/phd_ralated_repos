function [mc, nc, tx] = main

models = {'hgf','vkf'};

for i=1:length(models)
    [ckv(:,i), tx(i,:)] = run(models{i}); %#ok<AGROW>
end

ckv = fisher(ckv);
mc = invfisher(mean(ckv));
nc = mean(ckv<0);

end
