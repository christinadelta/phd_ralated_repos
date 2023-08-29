function blocks = makeblocks(blockTrials)

% create blocks array
block1(1:blockTrials(1),1)  = 1;
block2(1:blockTrials(2),1)  = 2;
block3(1:blockTrials(1),1)  = 3;
block4(1:blockTrials(1),1)  = 4;
block5(1:blockTrials(2),1)  = 5;
block6(1:blockTrials(1),1)  = 6;
block7(1:blockTrials(1),1)  = 7;
block8(1:blockTrials(2),1)  = 8;
block9(1:blockTrials(1),1)  = 9;

blocks                      = [block1;block2;block3;block4;block5;block6;block7;block8;block9];

end % end of function