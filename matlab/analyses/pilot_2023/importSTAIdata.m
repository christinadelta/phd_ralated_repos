function STAIT = importSTAIdata(subdir)

% function imports perceptual learning trials store in xls files from gorilla

fprintf('\t\t loading STAI data file \n\n');
STAIfile = dir(fullfile(subdir,'*data_exp_142429-v14_questionnaire-q9my*'));

% take only the name 
STAIname      = STAIfile.name;
STAItempdir   = append(subdir,'/');
STAIfilename  = append(STAItempdir,STAIname); % combine the name with the direcory 
STAIT         = readtable(STAIfilename);

end % end of fucntion 