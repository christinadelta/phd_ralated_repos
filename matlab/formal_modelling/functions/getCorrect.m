function [all_stc, corr_stc] = getCorrect(a,o,ss,tt)

% compute performance 
nsim = size(a,2);

for i = 1:nsim
    choice = a(:,i);
    choice = 2 - choice; % convert to 0 and 1
    sim_o = o(:,i);

    for t = 1:length(choice)

        % for each trial, is the simulated response correct or not?
        if choice(t,1) == sim_o(t,1)
            correct(t,i) = 1;
        else 
            correct(t,i) = 0;
        end
    end % end of trials loop
end % end of simulation loop

%%
% split responses and correctness per stc and vol
smallstc    = correct(ss(:,1),:);
mediumstc   = correct(ss(:,2),:);
largestc    = correct(ss(:,3),:);

% split per voalitlity 
small_stbl  = smallstc(tt(:,1),:);
small_vol   = smallstc(tt(:,2),:);

medium_stbl = mediumstc(tt(:,1),:);
medium_vol  = mediumstc(tt(:,2),:);

large_stbl  = largestc(tt(:,1),:);
large_vol   = largestc(tt(:,2),:);

%% compute averaged performance 

all_stc(1,:) = mean(smallstc,"all");
all_stc(2,:) = mean(mediumstc,"all");
all_stc(3,:) = mean(largestc,"all");

% compute performance for small (stable and volatile)
corr_stc(1,1) = mean(small_stbl,"all");
corr_stc(1,2) = mean(small_vol,"all");

% compute performance for medium (stable and volatile)
corr_stc(2,1) = mean(medium_stbl,"all");
corr_stc(2,2) = mean(medium_vol,"all");

% compute performance for large (stable and volatile)
corr_stc(3,1) = mean(large_stbl,"all");
corr_stc(3,2) = mean(large_vol,"all");



end % end of function