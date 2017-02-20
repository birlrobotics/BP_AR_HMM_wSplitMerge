% Extract a given data file from multiple folders
% Input:
%   - ParentDir: main directory path to data results
%   - dataTyep:  name of data file that we want to extract
%   - nTrail:    number of time series we want to extract to build model
function [data, folderNames] = extract_file(base_dir, dataSource, numTrials)

global DB_FILTER;

% Initialization of local structures
data(numTrials).obs;
folderNames = cell(numTrials,1); % Preallocate some folder size

% Randomly select folder data
cd(base_dir);
folders = dir(base_dir);
if(numTrials>length(folders))
    numTrials=length(folders);
end
id = randperm(length(folders(1:numTrials))); 

for i=1:length(id) 
    % Ensure it's a desired folder
    if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git'))
        
        % Iteratively append folder names to this cell
        folderNames = [folderNames, { folders(id(i)).name } ];
        
        % Switch to a folder
        fn = strcat(base_dir,folders(id(i)).name,'/');
        cd(fn);
        
        % Load the data not including time.
        for j=1:length(dataSource)
            file = dir([char(dataSource(j)) '.dat*']); % Assumes a list of data sources
            filename = strcat(fn,file(1).name);
            
            % Data
            d = load(filename);
            d = d(:,2:end);  %delete the time column
            
            if(DB_FILTER)
                if strcmp(dataSource(j), 'R_Torques')
                    d = smoothWrenchFilter(d); 
                end
            end
            
            % Insert data into an observation element in the data_structure.
            data(i).obs = d;
            
        end % Go through all data sources        
    end     % Validate folders
end         % Go through all data trials
end         % End function