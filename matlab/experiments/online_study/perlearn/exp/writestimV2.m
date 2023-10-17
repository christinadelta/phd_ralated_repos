function [stimuli_cues, stimuli_outcomes, pitch_pngs] = writestimV2(outcomes_array, cues_array)

% create cells with stimuli names for cues (low and high tone) and outcomes
% (face and house)

%% create cell array for cues 

for i = 1:length(cues_array)

    if cues_array(i,1) == 1
        stimuli_cues{i} = 'low.mp3';
        stimuli_pitch{i} = 'low.png';

    elseif cues_array(i,1) == 2
        stimuli_cues{i} = 'high.mp3';
        stimuli_pitch{i} = 'high.png';
    end


end % end of cues loop

stimuli_cues    = stimuli_cues';
pitch_pngs      = stimuli_pitch'; % in case we will have pngs

%% create cell array for outcomes

for j = 1:length(outcomes_array)
    
    % we have 2 face and 2 house stimuli, so make sure that half stimuli are
    % stimuli 1 and the other half are stimuli 2 (both for faces and
    % houses)
    if rem(j,2) == 1

        if outcomes_array(j,1) == 1
            stimuli_outcomes{j} = 'house1.png';
    
        elseif outcomes_array(j,1) == 2 
            stimuli_outcomes{j} = 'face1.png';
        end

    elseif rem(j,2) == 0

        if outcomes_array(j,1) == 1
            stimuli_outcomes{j} = 'house2.png';
    
        elseif outcomes_array(j,1) == 2
            stimuli_outcomes{j} = 'face2.png';
        end
    end
end % end of cues loop


stimuli_outcomes = stimuli_outcomes';

%%

end