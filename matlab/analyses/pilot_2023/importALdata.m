function ALT = importALdata(subdir)

% function imports action learning learning trials store in xls files from gorilla

fprintf('\t\t loading action learning data file \n\n');
ALfile = dir(fullfile(subdir,'*data_exp_142429-v14_task-hu2y*'));

% take only the name 
ALname      = ALfile.name;
ALtempdir   = append(subdir,'/');
ALfilename  = append(ALtempdir,ALname); % combine the name with the direcory 
ALT         = readtable(ALfilename);

end % end of fucntion 