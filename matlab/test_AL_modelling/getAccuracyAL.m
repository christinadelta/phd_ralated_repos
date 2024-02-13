function [m_all, m_vol, m_stc] = getAccuracyAL(ALTdata)

% compute accuracy for:
% 1. all data
% 2. Stable & volatile data
% 3. Stable & volatile data per stochasticity 

% extract accuracy
acc_data    = ALTdata(:,11);

% 

%% 1. compute mean accuracy for all data
m_all       = mean(acc_data);

%% 2. compute mean accuracy for stable & volatile 

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
m_vol(1,1)  = mean(all_stable(:,11));
m_vol(1,2)  = mean(all_vol(:,11));

%% 3. compute means for table & volatile per stochasticity 

m_stc(1,1)  = mean(stb1(:,11));
m_stc(1,2)  = mean(vol1(:,11));
m_stc(2,1)  = mean(stb2(:,11));
m_stc(2,2)  = mean(vol2(:,11));
m_stc(3,1)  = mean(stb3(:,11));
m_stc(3,2)  = mean(vol3(:,11));


end % end of function