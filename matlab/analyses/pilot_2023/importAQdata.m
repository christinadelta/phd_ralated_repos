function AQT = importAQdata(subdir)

% function imports AQ xls files from gorilla

fprintf('\t\t loading AQ file \n\n');
AQfile = dir(fullfile(subdir,'*data_exp_142429-v14_questionnaire-wlgi*'));

% take only the name 
AQname      = AQfile.name;
AQtempdir   = append(subdir,'/');
AQfilename  = append(AQtempdir,AQname); % combine the name with the direcory 
AQT         = readtable(AQfilename);










end % end of fucntion 