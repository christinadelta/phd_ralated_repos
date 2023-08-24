% ADD GAUSSIA NOISE TO IMAGES FOR THE PERCEPTUAL LEARNING TASK (PAL
% PROJECT ONLINE)
% GAUSSIAN NOISE AT 3 LEVELS: 


% christinadelta - August 2023

%% clear all

clc
clear all

%% load images

face1   = imread("face1.png");
% imshow(face1) % display image

face2   = imread("face2.png");
house1  = imread("house1.png");
house2  = imread("house2.png");

%% add gaussian noise (3 levels)

% add noise level 1
face1_l1 = imnoise(face1,'salt & pepper', 0.30);
face1_l2 = imnoise(face1,'salt & pepper', 0.60);

face2_l1 = imnoise(face2,'salt & pepper', 0.30);
face2_l2 = imnoise(face2,'salt & pepper', 0.60);

house1_l1 = imnoise(house1,'salt & pepper', 0.30);
house1_l2 = imnoise(house1,'salt & pepper', 0.60);

house2_l1 = imnoise(house2,'salt & pepper', 0.30);
house2_l2 = imnoise(house2,'salt & pepper', 0.60);


% display image
imshow(house2_l2)
