function SPQT = importSPQdata(subdir)

% function imports AQ xls files from gorilla

fprintf('\t\t loading SPQ file \n\n');
SPQfile = dir(fullfile(subdir,'*data_exp_142429-v14_questionnaire-jeno*'));

% take only the name 
SPQname      = SPQfile.name;
SPQtempdir   = append(subdir,'/');
SPQfilename  = append(SPQtempdir,SPQname); % combine the name with the direcory 
SPQT         = readtable(SPQfilename);

end % end of fucntion 