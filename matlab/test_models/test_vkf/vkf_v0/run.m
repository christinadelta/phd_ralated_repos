function [ckv, tx] = run(modelname)
simcat = 'comp_switch';
nsim = 500;


switch modelname
    case 'hgf'
        model = @model_hgf;
        d = 3;
    case 'vkf'
        model = @model_vkf;
        d = 3;
    otherwise
        error('!');
end

pipedir = getdefaults('pipedir');
fsim1    = fullfile(pipedir,simcat, sprintf('data_train.mat') );
fsim2    = fullfile(pipedir,simcat, sprintf('data_test.mat') );
ffit  = fullfile(pipedir,simcat, sprintf('fit_%s.mat',modelname) );
makedir(fullfile(pipedir,simcat));

rng(0);
Ys = sim_gen(fsim1,nsim);
sim_gen(fsim2,nsim);

do_fit = ~exist(ffit,'file');


vfit = 15.23;
if do_fit
    config = struct('range',[-5*ones(1,d);5*ones(1,d)],'numinit',7*d,'prior_for_bads',1);
    prior = struct('mean',zeros(d,1),'variance',vfit);
    cbm  = cbm_lap({Ys}, model, prior, [], config); %#ok<NASGU>
    save(ffit, 'cbm');    
end

cbm = load(ffit);
cbm = cbm.cbm;
[ckv, tx] = eval_model(fsim1, model, cbm);

end