function [m_all, m_vol, m_stc, m_stc1_sv, m_stc2_sv, m_stc3_sv] = getAccuracyALT(ALTdata)

% extract accuracy
acc_data = ALTdata(:,10);
m_all = mean(acc_data);

% split data per volatility
s1 = find(ALTdata(:,3)==0);
s2 = find(ALTdata(:,3)==1);
s3 = find(ALTdata(:,3)==3);
s4 = find(ALTdata(:,3)==4);
s5 = find(ALTdata(:,3)==6);
s6 = find(ALTdata(:,3)==7);
s7 = find(ALTdata(:,3)==9);

stb1 = ALTdata((s1),:);
stb2 = ALTdata((s2),:);
stb3 = ALTdata((s3),:);
stb4 = ALTdata((s4),:);
stb5 = ALTdata((s5),:);
stb6 = ALTdata((s6),:);
stb7 = ALTdata((s7),:);

stbl    = cat(1,stb1,stb2,stb3,stb4,stb5,stb6,stb7);

v1      = find(ALTdata(:,3)==2);
v2      = find(ALTdata(:,3)==5);
v3      = find(ALTdata(:,3)==8);

vol1 = ALTdata((v1),:);
vol2 = ALTdata((v2),:);
vol3 = ALTdata((v3),:);
vol     = cat(1,vol1,vol2,vol3);

% compute means all stable and volatile
acc_stable = stbl(:,10);
m_vol(1,1) = mean(acc_stable);

acc_vol = vol(:,10);
m_vol(1,2) = mean(acc_vol);


% split data per stochasticity
stc1 = cat(1,stb1,stb2,vol1,stb3);
stc2 = cat(1,stb4,vol2,stb5);
stc3 = cat(1,stb6,vol3,stb7);

acc_stc1 = stc1(:,10);
acc_stc2 = stc2(:,10);
acc_stc3 = stc3(:,10);

% compute mean stochasticity
m_stc(1,1) = mean(acc_stc1);
m_stc(1,2) = mean(acc_stc2);
m_stc(1,3) = mean(acc_stc3);

%% for rach stochasticity get means for stable and voaltile envs

% extract stable trials for stc 1,2,3
stc1_s = cat(1,stb1,stb2,stb3);
stc2_s = cat(1,stb4,stb5);
stc3_s = cat(1,stb6,stb7);

acc_stc1_s = stc1_s(:,10);
acc_stc2_s = stc2_s(:,10);
acc_stc3_s = stc3_s(:,10);

% extract volatile trials for stc 1,2,3
stc1_v = vol1;
stc2_v = vol2;
stc3_v = vol3;

acc_stc1_v = stc1_v(:,10);
acc_stc2_v = stc2_v(:,10);
acc_stc3_v = stc3_v(:,10);

% compute means for each stc level
m_stc1_sv(1,1) = mean(acc_stc1_s);
m_stc1_sv(1,2) = mean(acc_stc1_v);

m_stc2_sv(1,1) = mean(acc_stc2_s);
m_stc2_sv(1,2) = mean(acc_stc2_v);

m_stc3_sv(1,1) = mean(acc_stc3_s);
m_stc3_sv(1,2) = mean(acc_stc3_v);



end % end of function