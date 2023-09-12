function [m_all, m_vol, m_stc] = getAccuracyPL(subdata)

% get % correct for:
% 1. predictions 
% 2. outcomes
% 3. preidctions (for each volatility) 
% 4. outcomes per voaltility 
% 5. predictions per stochasticity 
% 6. outcomes per stochasticty 

% -------------

%% compute mean for prediction correct and outcome correct (all data)

p               = subdata(:,10); % extract preidctions correct
o               = subdata(:,13); % extract outcome response correct

% look at proportion correct for outcomes and preidctions all blocks
% and conditions
m_all(1,1) = mean(p);
m_all(1,2) = mean(o);

%% compute mean for prediction correct and outcome correct for stable vs volatile 

% look at differences between stable and volatile phases (all
% stochasticities) -- split stable (1,3,5) and volatile (2,4,6) blocks
tmp_stbl1       = find(subdata(:,2) == 1);
tmp_stbl2       = find(subdata(:,2) == 3);
tmp_stbl3       = find(subdata(:,2) == 5);

tmp_vol1        = find(subdata(:,2) == 2);
tmp_vol2        = find(subdata(:,2) == 4);
tmp_vol3        = find(subdata(:,2) == 6);

stbl1           = subdata((tmp_stbl1),:);
stbl2           = subdata((tmp_stbl2),:);
stbl3           = subdata((tmp_stbl3),:);

allstable       = cat(1,stbl1, stbl2, stbl3);

vol1            = subdata((tmp_vol1),:);
vol2            = subdata((tmp_vol2),:);
vol3            = subdata((tmp_vol3),:);

allvol          = cat(1, vol1, vol2, vol3);

% extarct predictions and outcomes for stable and volatile 
p_stable        = allstable(:,10);
o_stable        = allstable(:,13);
p_vol           = allvol(:,10);
o_vol           = allvol(:,13);

% compute means 

m_vol(1,1)     = mean(p_stable);
m_vol(1,2)     = mean(o_stable);

m_vol(2,1)     = mean(p_vol);
m_vol(2,2)     = mean(o_vol);

%% compute mean for prediction correct and outcome correct per stochasticity level

% difine stc levels matrices
stc1            = cat(1,stbl1,vol1);
stc2            = cat(1,stbl2,vol2);
stc3            = cat(1,stbl3,vol3);

% extract predictions and outcomes
p_stc1          = stc1(:,10);
p_stc2          = stc2(:,10);
p_stc3          = stc3(:,10);

o_stc1          = stc1(:,13);
o_stc2          = stc2(:,13);
o_stc3          = stc3(:,13);

% get means
m_stc(1,1)      = mean(p_stc1);
m_stc(1,2)      = mean(o_stc1);

m_stc(2,1)      = mean(p_stc2);
m_stc(2,2)      = mean(o_stc2);

m_stc(3,1)      = mean(p_stc3);
m_stc(3,2)      = mean(o_stc3);

end % end of function