
% the script creates high and low frequency sounds to be used as cues in
% the perceptual learning task

% a) create high pitch and low pich audios
% b) store as wav files 

% files should be converted to mp3 in different apps as I can't create mp3
% files with the "audiowrite" function (removed since r2015b).

% created in September 2023
% @christinadelta

% -------------

clear all 
clc

%% general variables
amp         = 10;
fs          = 20500; % sampling frequency
duration    = 1; % duration in sec
values      = 0:1/fs:duration; % 1 second
fhigh       = 500; % frequency for high pitch
flow        = 50;        
high_pitch  = amp*sin(2*pi*fhigh*values);
low_pitch   = amp*sin(2*pi*flow*values);

%% display
sound(high_pitch,fs); % high pitch
sound(low_pitch,fs); % low pitch

%% save as .wav (mp3 is not supported with audiowrite) 

filename_high   = 'high_pitch.wav';
filename_low    = 'low_pitch.wav';
audiowrite(filename_high,high_pitch,fs);
audiowrite(filename_low,low_pitch,fs);


