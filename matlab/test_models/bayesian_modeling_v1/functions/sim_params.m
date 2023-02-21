function [tx_vkf, tx_hgf, mc, nc, nc_supp] = sim_params(simy,simx,outc)


[mc, nc, tx] = mainf(simy,simx);
nc_supp = supp(outc);

tx_vkf = tx(2,:);
tx_hgf = tx(1,:);


end