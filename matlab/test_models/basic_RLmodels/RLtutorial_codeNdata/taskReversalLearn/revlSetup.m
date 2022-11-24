% SCRIPT revlSetup
% Script to unclutter the main Run script. Sets up all variables needed to 
% run the task, allocates storage, load image files. Here screensize can be
% specified for debugging. 
%
% Hanneke den Ouden 
% Donders Institute for Brain, Cognition and Behaviour
% h.denouden@gmail.com
% 
% version 11-08-2015
%==========================================================================

% A. Set stuff up start PTB - Screen.
%--------------------------------------------------------------------------
% file and directory names
subTag      = sprintf('tutorialRevLearn_%s_s%03.0f',prep.volatility,sID);
dataFile = fullfile(dataPath,sprintf('%s_data.mat', subTag));

% create output directory if it doesn't exist yet
if ~exist(dataPath,'dir')
    mkdir(dataPath);
end

today       = date;
nt          = prep.nt;
data     = {};
data.subTag = subTag;

% allocate storage
data.choice  = nan(nt,1);
data.outcome = NaN(nt,1);
data.RT      = NaN(nt,1);
tm.stim         = NaN(nt,1);
tm.choice       = NaN(nt,1);

% Load all the image files.
tmp.feedback{1}    = importdata(fullfile(prep.dir.imgs,'noeuro.jpg'));  % feedback
tmp.feedback{2}    = importdata(fullfile(prep.dir.imgs,'euro.jpg')); % feedback
for t = 1:length(prep.cue)
    tmp.cue{t}    = importdata(prep.cue{t});
end

% B.   Start Psychtoolbox screen, show instructions
%--------------------------------------------------------------------------
if debug == 1 % Use when in debugging mode
    Screen('Preference', 'Verbosity', 4);
    % Uncomment to test for smaller screen
    wd          = Screen('OpenWindow', 0, [0 0 0], [0 0 800 600]);
    % Uncomment to test for MR screen res.
%     wd          = Screen('OpenWindow', 0, [0 0 0], [0 0 1024 768]);
else % % Use when in run mode
    Screen('Preference', 'Verbosity', 0);
     wd = Screen('OpenWindow', 0, [0 0 0], []);
end

Screen(wd,'TextStyle',1); % Bold
Screen(wd,'TextFont',prep.draw.font);
Screen(wd,'TextSize',prep.draw.txtsize);
black       = BlackIndex(wd);
white       = WhiteIndex(wd);
HideCursor; 

% get the size of the screen and figure out the locations of the stimuli
[wdw, wdh]      = Screen('WindowSize', wd);
centre          = round([wdw/2, wdh/2]);
dist            = prep.stimDistance;
loc1            = [centre(1)-dist, centre(2)]; % centre of pie 1: centre left.
loc2            = [centre(1)+dist, centre(2)]; % centre of pie 2: centre right.
loc5            = centre; % centre of screen: feedback.

rect.stim{1}          = [(loc1(1)-prep.draw.rStim) (loc1(2)-prep.draw.rStim) ...
    (loc1(1)+prep.draw.rStim) (loc1(2)+prep.draw.rStim)]; %rectangle defining image left
rect.stim{2}          = [(loc2(1)-prep.draw.rStim) (loc2(2)-prep.draw.rStim) ...
    (loc2(1)+prep.draw.rStim) (loc2(2)+prep.draw.rStim)]; %rectangle defining image right
rect.stim{5}          = [(loc5(1)-prep.draw.rSmile) (loc5(2)-prep.draw.rSmile) ...
    (loc5(1)+prep.draw.rSmile) (loc5(2)+prep.draw.rSmile)]; %rectangle defining feedback
rect.frame{1}         = [(loc1(1)-prep.draw.rFrame) (loc1(2)-prep.draw.rFrame) ...
    (loc1(1)+prep.draw.rFrame) (loc1(2)+prep.draw.rFrame)];
rect.frame{2}         = [(loc2(1)-prep.draw.rFrame) (loc2(2)-prep.draw.rFrame) ...
    (loc2(1)+prep.draw.rFrame) (loc2(2)+prep.draw.rFrame)]; 

prep.draw.loc{1}    = loc1;
prep.draw.loc{2}    = loc2;
prep.draw.loc{5}    = loc5;
prep.draw.rect.stim{1}   = rect.stim{1};
prep.draw.rect.stim{2}   = rect.stim{2};
prep.draw.rect.stim{5}   = rect.stim{5};
prep.draw.rect.frame{1}   = rect.frame{1};
prep.draw.rect.frame{2}   = rect.frame{2};
prep.draw.black     = black;
prep.draw.white     = white;
prep.draw.wdw       = wdw;
prep.draw.wdh       = wdh;
prep.draw.x0        = round(wdw/2);
prep.draw.y0        = round(wdh/2);
top                 = round(wdh/12); % location of the top of the text
prep.top            = top;

% make the textures for the feedback
img.feedback{1}   = Screen('MakeTexture',wd,tmp.feedback{1});
img.feedback{2}   = Screen('MakeTexture',wd,tmp.feedback{2});

% make the textures for the stimuli
for j = 1:2
    img.stim{j}   = Screen('MakeTexture',wd,tmp.cue{j});
    img.select{j}      = Screen('MakeTexture',wd,tmp.cue{j+2});
end

% present black screen for 1 sec before starting
Screen('FillRect',wd,black);
Screen('Flip', wd);
WaitSecs(1);
