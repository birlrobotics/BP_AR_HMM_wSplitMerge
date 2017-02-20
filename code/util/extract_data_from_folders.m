%=========================================================================================
% data = extract_data_from_folders(base_dir, dataSource, data_params)
%=========================================================================================
% Extract a given data file from multiple folders
%=========================================================================================
% Input:
%   - ParentDir:    main directory path to data results
%   - dataTyep:     name of data file that we want to extract
%   - data_params:  data structure containing VAR order, #States,#Dims,#Trials,timeSteps
% Output:
%   - data:         contains observations (N,dim*numTrials).
%   - data_params:  updated data_params
%=========================================================================================
function [data,data_params] = extract_data_from_folders(base_dir, dataSource, data_params)

global DB_FILTER;

%% Initialization of local structures
numTrials=data_params.N;
timeSteps = 0;

data = struct('obs',{});
data(100).obs=[];
folderNames = cell(numTrials,1); % Preallocate some folder size

%% Randomly select folder data
cd(base_dir);
folders = dir(base_dir);
folder_len=length(folders);
if(numTrials>folder_len)
    numTrials=folder_len;
end
id = randperm(length(folders)); % may need to change for cross-validation
id = id(1:numTrials);
%% Extract data
for i=1:numTrials 
    % Ensure it's a desired folder
    if(folders(id(i)).isdir == 1 && ~strcmp(folders(id(i)).name,'..') && ~strcmp(folders(id(i)).name,'.') && ~strcmp(folders(id(i)).name,'.git'))
        
        % Iteratively append folder names to this cell
        folderNames(i) = { folders(id(i)).name };
        
        % Switch to a folder
        fn = strcat(base_dir,folders(id(i)).name,'/');
        cd(fn);
        
        % Load the data not including time.
        for j=1:length(dataSource)
            file = dir([char(dataSource(j)) '.dat*']); % Assumes a list of data sources
            filename = strcat(fn,file(1).name);
            
            % Data
            d = load(filename);
            d = d(:,2:end);             % Do not include the time column
            len_d = length(d);
            % Record largest number of time-steps in all series
            if len_d>timeSteps
                data_params.T=len_d;        % # of time-steps
                timeSteps=len_d;
            end
            
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