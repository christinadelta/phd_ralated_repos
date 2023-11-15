function h = plotlrsV1(lr, volindex)

% plot stable and volatile learning rates

%% first re-arrange the data 

colours = ['#EDB120'; "#D95319"];

% extract learning rates for stable and voltile conditions
stbl_lr = lr(volindex(:,1),:);
vol_lr  = lr(volindex(:,2),:);

% add vol/stbl string
for i = 1:2
    for ii = 1:length(vol_lr)

        if i == 1 % if stable 
            env(ii,i) = "stable";
        else % if volatile
            env(ii,i) = "volatile";
        end
    end
end

env = env(:); 
env = categorical(env);
lrs = [stbl_lr; vol_lr];

clear stbl_lr vol_lr i ii

%% now create box plots

h = boxchart(lrs,'GroupByColor',env);




end % end of function