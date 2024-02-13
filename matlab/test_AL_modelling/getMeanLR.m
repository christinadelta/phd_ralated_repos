function [meanLR,allrepsLR] = getMeanLR(all_lr)

%% HOW MANY REPETITIONS?
nreps = size(all_lr,2);

for rep = 1:nreps

    templr = all_lr{1,rep};
    % loop over stc and vol conditions
    for i = 1:3
        for j = 1:2
            allrepsLR{i,j}(:,rep) = templr{i,j}(:,1);
        end 
    end
end % end of reps loop

%% average lr

for i = 1:3

    for j = 1:2
        tmp         = allrepsLR{i,j};
        meanLR{i,j} = mean(tmp,2);

    end 
end 


end % end of function