function [m_allp, m_allo, m_volp, m_volo, m_stcp, m_stco, m_stc1_p, m_stc2_p, m_stc3_p, m_stc1_o, m_stc2_o, m_stc3_o] = getAccuracyPL(subdata)

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
m_allp(1,1) = mean(p);
m_allo(1,1) = mean(o);

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

% compute means per volatility
m_volp(1,1)     = mean(p_stable);
m_volp(1,2)     = mean(p_vol); 

m_volo(1,1)     = mean(o_stable);
m_volo(1,2)     = mean(o_vol);

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
m_stcp(1,1)      = mean(p_stc1);
m_stcp(1,2)      = mean(p_stc2);
m_stcp(1,3)      = mean(p_stc3);

m_stco(1,1)      = mean(o_stc1);
m_stco(1,2)      = mean(o_stc2);
m_stco(1,3)      = mean(o_stc3);

%% compute mean of stable-volatile accuracy % for each stch

% extract responses for stc 1 
stc1_stblp       = stbl1(:,10);
stc1_volp        = vol1(:,10);
stc1_stblo       = stbl1(:,13);
stc1_volo        = vol1(:,13);

% extract responses for stc 2 
stc2_stblp       = stbl2(:,10);
stc2_volp        = vol2(:,10);
stc2_stblo       = stbl2(:,13);
stc2_volo        = vol2(:,13);

% extract responses for stc 3 
stc3_stblp       = stbl3(:,10);
stc3_volp        = vol3(:,10);
stc3_stblo       = stbl3(:,13);
stc3_volo        = vol3(:,13);

% compute means for predictions
m_stc1_p(1,1) = mean(stc1_stblp); % stable
m_stc1_p(1,2) = mean(stc1_volp); % volatile

m_stc2_p(1,1) = mean(stc2_stblp); % stable
m_stc2_p(1,2) = mean(stc2_volp); % volatile

m_stc3_p(1,1) = mean(stc3_stblp); % stable
m_stc3_p(1,2) = mean(stc3_volp); % stable

% compute means for outcomes
m_stc1_o(1,1) = mean(stc1_stblo); % stable
m_stc1_o(1,2) = mean(stc1_volo); % volatile

m_stc2_o(1,1) = mean(stc2_stblo); % stable
m_stc2_o(1,2) = mean(stc2_volo); % volatile

m_stc3_o(1,1) = mean(stc3_stblo); % stable
m_stc3_o(1,2) = mean(stc3_volo); % stable

end % end of function