function [mrts_all, mrts_vol, mrts_stc] = getRTsAL(ALTdata)

% compute mean RTs for:
% 1. all data
% 2. Stable & volatile data
% 3. Stable & volatile data per stochasticity 

% extract accuracy
rt_data    = ALTdata(:,10);

%% 1. compute mean rts for all data
mrts_all       = mean(rt_data);

%% 2. compute mean rts for stable & volatile 

% get stable trials
stable1 = find(ALTdata(:,3)==1);
stable2 = find(ALTdata(:,3)==3);
stable3 = find(ALTdata(:,3)==5);

stb1    = ALTdata((stable1),:);
stb2    = ALTdata((stable2),:);
stb3    = ALTdata((stable3),:);

% concatenate them
all_stable = cat(1,stb1,stb2,stb3);

% get volatile trials
vol1        = find(ALTdata(:,3)==2);
vol2        = find(ALTdata(:,3)==4);
vol3        = find(ALTdata(:,3)==6);

vol1        = ALTdata((vol1),:);
vol2        = ALTdata((vol2),:);
vol3        = ALTdata((vol3),:);
all_vol     = cat(1,vol1,vol2,vol3);

% compute means: all stable and volatile
mrts_vol(1,1)  = mean(all_stable(:,10));
mrts_vol(1,2)  = mean(all_vol(:,10));

%% 3. compute mean rts for stable & volatile per stochasticity 

mrts_stc(1,1)  = mean(stb1(:,10));
mrts_stc(1,2)  = mean(vol1(:,10));
mrts_stc(2,1)  = mean(stb2(:,10));
mrts_stc(2,2)  = mean(vol2(:,10));
mrts_stc(3,1)  = mean(stb3(:,10));
mrts_stc(3,2)  = mean(vol3(:,10));






end % end of function