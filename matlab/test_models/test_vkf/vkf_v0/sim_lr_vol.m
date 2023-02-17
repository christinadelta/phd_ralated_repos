function [tx_vkf, tx_hgf, mc, nc, nc_supp] = sim_lr_vol


[mc, nc, tx] = main;
nc_supp = supp;

tx_vkf = tx(2,:);
tx_hgf = tx(1,:);


end