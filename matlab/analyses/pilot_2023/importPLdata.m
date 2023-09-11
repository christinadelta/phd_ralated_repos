function PLT = importPLdata(subdir)

% function imports perceptual learning trials store in xls files from gorilla

fprintf('\t\t loading perceptual data file \n\n');
PLfile = dir(fullfile(subdir,'*data_exp_142429-v14_task-7lag*'));

% take only the name 
PLname      = PLfile.name;
PLtempdir   = append(subdir,'/');
PLfilename  = append(PLtempdir,PLname); % combine the name with the direcory 
PLT         = readtable(PLfilename);

end % end of fucntion 