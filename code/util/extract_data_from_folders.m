%=========================================================================================
% data = extract_data_from_folders(base_dir, dataSource, data_params)
%=========================================================================================
% Extract a given data file from multiple folders
%=========================================================================================
% Input:
%   - base_dir:     main directory path to data results
%   - dataSource:   the type of data we want to include (wrenches, joing angles, etc). From the possible terms: [R_Torques.dat, R_Angles.dat, R_CartPos.dat]
%   - data_params:  data structure containing VAR order, #States,#Dims,#Trials,timeSteps
% Output:
%   - data:         contains observations (N,dim*numTrials).
%   - data_params:  updated data_params
%=========================================================================================
function [data,data_params] = extract_data_from_folders(base_dir, dataSource, data_params)

global DB_FILTER;

%% Initialization of local structures
numTrials=data_params.N;
timeStepsVec=zeros(1,numTrials);

data = struct('obs',{});
data(numTrials).obs=[];
folderNames = cell(numTrials,1); % Preallocate some folder size

%% Randomly select folder data
folders = dir(base_dir);
folder_len=length(folders);

% Remove non-folders
for k = folder_len:-1:1
    if ~folders(k).isdir || folders(k).name(1) == '.'    
        folders(k) = [ ];
        continue
    end
end

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
        cd(fn); % TODO: remove cd, just refer to its path.
        
        % Load the data not including time.
        for j=1:length(dataSource)
            file = dir([char(dataSource(j)) '.dat*']); % Assumes a list of data sources
            filename = strcat(fn,file(1).name);
            
            % Data
            d = load(filename);
            if( isempty(d) )
                fprintf('Data from the %ith trial is empty',j)
                exit;
            end
            d = d(:,2:end);             % Do not include the time column
            len_d = length(d);
            
            if(DB_FILTER)
                if strcmp(dataSource(j), 'R_Torques')
                    d = smoothWrenchFilter(d); 
                end
            end
            
            % Insert data into an observation element in the data_structure in a row-wise fashion.
            data(i).obs = d';        
            timeStepsVec(i)=len_d;
            
        end % Go through all data sources
    end     % Validate folders
end         % Go through all data trials

% Save the number of time-steps per trial
data_params.T=timeStepsVec;

% Get the total number of time-steps for all trials.
%for i=1:data_params.N
    % data_params.T=data_params.T+size(data(i).obs,2); 
%end
end         % End function