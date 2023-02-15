%rating files from .csv to .mat 
subjects = 16:34;

for subject = subjects
    
    dinfo = dir('*.csv');
    filenames = {dinfo.name};
    
    for K = 1 : length(filenames)
        thisfile = filenames{K};
        
        T = readtable(thisfile,'ReadVariableNames',false);
        output = T;
      
        [~, basename, ~] = fileparts(thisfile);
        newfile = [basename '.mat'];
        
        save(newfile,'output')
        
        %S = regexp( fileread(thisfile), '\r?\n', 'split' );
        %S( cellfun(@isempty, S) ) = [];   %delete empties, especially likely right at end
        %clear TS
        %TS.(basename) = S;
        %output = S;
        %save(newfile);
    end
    
end

%T=readtable('mycsv.csv');
     %       ^^^^^^^^^------ your csv filename
%p=T{:,1};
%q=T{:,2};
%save('mymat.mat','p','q')

