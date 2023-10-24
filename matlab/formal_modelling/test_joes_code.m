parfor i = 1:k 
    data_subj = sb(i); 
    fname_subj = fullfile('MatlabFittng/lap_subjects_USPilot',['lap_USPilot_', num2str(i), '.mat']); 
    cbm_lap(data_subj, @ABA_shift_Gen, prior_USPilot , fname_subj); 
end