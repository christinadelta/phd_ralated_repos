function IUT = importIUdata(subdir)

% function imports perceptual learning trials store in xls files from gorilla

fprintf('\t\t loading IUS data file \n\n');
IUfile = dir(fullfile(subdir,'*data_exp_142429-v14_questionnaire-qpdd*'));

% take only the name 
IUname      = IUfile.name;
IUtempdir   = append(subdir,'/');
IUfilename  = append(IUtempdir,IUname); % combine the name with the direcory 
IUT         = readtable(IUfilename);

end % end of fucntion 