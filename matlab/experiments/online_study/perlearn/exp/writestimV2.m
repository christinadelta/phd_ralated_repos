function [stimuli_cues, stimuli_outcomes, pitch_pngs] = writestimV2(outcomes_array, cues_array)

% create cells with stimuli names for cues (low and high tone) and outcomes
% (face and house)

%% create cell array for cues 

for i = 1:length(cues_array)

    if cues_array(i,1) == 1
        stimuli_cues{i} = 'low.mp3';
        stimuli_pitch{i} = 'low.png';

    else 
        stimuli_cues{i} = 'high.mp3';
        stimuli_pitch{i} = 'high.png';
    end


end % end of cues loop

stimuli_cues    = stimuli_cues';
pitch_pngs      = stimuli_pitch'; % in case we will have pngs

%% create cell array for outcomes

for i = 1:length(outcomes_array)

    if outcomes_array(i,1) == 1
        stimuli_outcomes{i} = 'house.png';
   
    else 
        stimuli_outcomes{i} = 'face.png';
    end


end % end of cues loop


stimuli_outcomes = stimuli_outcomes';

%%

end