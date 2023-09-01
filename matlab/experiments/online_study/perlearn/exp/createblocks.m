function blocks = createblocks(blockTrials)

% create blocks array
block1(1:blockTrials,1)  = 1;
block2(1:blockTrials,1)  = 2;
block3(1:blockTrials,1)  = 3;
block4(1:blockTrials,1)  = 4;
block5(1:blockTrials,1)  = 5;
block6(1:blockTrials,1)  = 6;

blocks                      = [block1;block2;block3;block4;block5;block6];

end % end of function