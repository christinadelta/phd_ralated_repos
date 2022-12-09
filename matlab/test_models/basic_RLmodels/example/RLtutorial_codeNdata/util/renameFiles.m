% rename odd subject numbers
ct = 1; 
for s = 1:2:60
    try
        load(sprintf('tutorialRevLearn_high_s%03.0f_data.mat',s));
        data.prep.sID = ct;
        data.subTag = sprintf('tutorialRevLearn_high_s%03.0f',ct);
        save(sprintf('tutorialRevLearn_high_s%03.0f_data.mat',ct));
        ct = ct+2;
        
    catch
        disp(['sub ' num2str(s) ' does not exist']);
    end
end

% rename even subject numbers
ct = 2; 
for s = 2:2:60
    try
        load(sprintf('tutorialRevLearn_low_s%03.0f_data.mat',s));
        data.prep.sID = ct;
        data.subTag = sprintf('tutorialRevLearn_low_s%03.0f',ct);
        save(sprintf('tutorialRevLearn_low_s%03.0f_data.mat',ct));
        ct = ct+2;
        
    catch
        disp(['sub ' num2str(s) ' does not exist']);
    end
end
