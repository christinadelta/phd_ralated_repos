function wd = revlDrawStim(wd,prep,cond,data,t,img)

%==========================================================================
% Code for drawing the different stimulus events. 
% - cond = 1 draw stimulus
% - cond = 2 draw stimulus + feedback
%
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 11-08-2015
%==========================================================================

% draw fixation cross.
drawfix

% draw stimuli at specified locations.
for iCue = 1:prep.nStim
    location = prep.locs(t,iCue);
    Screen('DrawTexture',wd,img.stim{iCue},[],prep.draw.rect.stim{location});
end

% draw choice.
if cond > 0 
    choiceLoc = prep.locs(t,data.choice(t));
    Screen('DrawTexture',wd,img.select{data.choice(t)},[],prep.draw.rect.stim{choiceLoc});
end

% draw feedback: for 0, get punish feedback, for 1, get reward feedback
if cond > 1 
    fbimg = img.feedback{data.outcome(t)+1};
    Screen('DrawTexture',wd,fbimg,[],prep.draw.rect.stim{5});
end
