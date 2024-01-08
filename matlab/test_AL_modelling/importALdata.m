function ALT = importALdata(subdir)

% function imports action learning learning trials stored in csv files from gorilla

fprintf('\t\t loading action learning data file \n\n');
ALfile = dir(fullfile(subdir,'*sub_*'));

% take only the name 
ALname      = ALfile.name;
ALtempdir   = append(subdir,'/');
ALfilename  = append(ALtempdir,ALname); % combine the name with the direcory 
ALT         = readtable(ALfilename);

end % end of fucntion 