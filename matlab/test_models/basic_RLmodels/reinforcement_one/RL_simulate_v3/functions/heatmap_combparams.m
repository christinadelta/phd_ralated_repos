function [performance] = heatmap_combparams(tmp_params,cond,probs,trials,outpath,task)

% the function is used to simulate, data run the RW model and average
% choice probabilities and performance across all combinations of alpha and
% beta parameter values

% --------------------------------------
% simulate data 
data        = avlearn_simulate_v1(cond, probs, trials, outpath, task); % data is a structure containaining the simulation output

% run model
modelout    = modelRW_v1(tmp_params, data, outpath);

if cond == 1
    if task == 1
        [~,imax]    = max(probs(1,:));
        tmptrial    = modelout.a == imax;
    end
else % if cond == 2

    runs            = 4; 
    trlruns         = trials/runs;
    counter         = 1;

    for run = 1:runs

        % work with p(correct) for the volatile condition
        tmp                     = modelout.a(counter:trlruns*run);
    
        if mod(run,2)
            [~,imax]            = max(probs(2,:));
            trls{run}           = tmp == imax;
        else
            [~,imin]            = min(probs(2,:));
            trls{run}           = tmp == imin;
        end
    
        counter                 = counter + trlruns;

    end % end of for loop

    tmptrial                    = [trls{1} trls{2} trls{3} trls{4}]; % concatenate
    
end % end of conditions loop

performance                     = nanmean(tmptrial); % get mean performance

end % end of function