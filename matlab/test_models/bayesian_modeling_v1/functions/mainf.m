function [mc nc tx] = mainf(xs,ys)

% models to be used:
models              = {'hgf', 'vkf'};

% run models 
for i = 1:length(models)
    
    % which model to run?
    modelname       = models{i};

    switch modelname
        case 'hgf'
            model   = @model_hgf;
            d = 3;
        case 'vkf'
            model   = @model_vkf;
            d = 3;
        otherwise
            error('!');
    end
    
    % 
    vfit            = 15.23; % variance 

    % run cbm to get the params 
    config          = struct('range',[-5*ones(1,d);5*ones(1,d)],'numinit',7*d,'prior_for_bads',1);
    prior           = struct('mean',zeros(d,1),'variance',vfit);
    cbm             = cbm_lap({ys}, model, prior, [], config); %#ok<NASGU>
    
    % run model evaluation to get parameter values
    [ckv(:,i), tx(i,:)]       = eval_model(ys, model, cbm);


end % end of models loop

% hmmm don't understand these yet
ckv             = fisher(ckv);
mc              = invfisher(mean(ckv));
nc              = mean(ckv<0);



end 