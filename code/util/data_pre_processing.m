% ============================================================
% [data]=data_pre_processing(data,dp)
% Convert's the data's mean to 0 can computes the covariance.
% Inputs:
%   data:           structure that holds the observations for all included trials.
%                   this is a T*(D*N) struc, where T:time_steps, D:dimensions, N:numTrials
%   data_params_p:  hold parameters that reflect properties of the data, pointer.
%------------------------------------------------------------------------------------------
% Outputs:
% data:             The observation mean will have been shifted to be equal to zero.
% meanSigma:        The first diff of the covariance set to 1.0
% ============================================================
function [data,meanSigma]=data_pre_processing(data,dp)

% Local variables
DB_FILTER=1;

% Large column vectors
Ybig1 = zeros(dp.T,(dp.nDim*dp.N));         % Holds smoothed observations
Ybig2 = Ybig1;                              % Holds variance
Ybig3 = Ybig1;                              % Holds the first diff. of the covariance 
%% Smooth the Data: requires you to set a radius parameter for smoothing
% 
if(DB_FILTER)
    radius = 20; % Default value 10

    for i=1:length(data)
        [dims,len] = size(data(i).obs); % # rows not even for all trials
        smoothed = zeros(dims,len);

        for j=1:dims
            for k=1:len
                low = max(1,k-radius);
                high = min(len,k+radius);
                smoothed(j,k) = mean(data(i).obs(j,low:high));
            end
        end
        data(i).obs = smoothed;
        Ybig1 = [Ybig1 data(i).obs];
    end
end
% Append data for each of the elements in the data strucuture in column-wise format
% for ii=1:length(data)
%     Ybig1 = [Ybig1 data(ii).obs];
% end

% Make copies of orig data for YBig2 and 3
Ybig2=Ybig1;
Ybig3=Ybig1;
%% Adjust each dimesnion to Mean=0

mY = mean((Ybig1')); % computes mean of each matrix-column
for i=1:length(data)
    for j=1:length(mY)
        data(i).obs(j,:) = data(i).obs(j,:) - mY(j);
    end
end

%% Adjust 1st diff variance of the observations to set them to 1.0
vY = var(diff(Ybig2'));
for i=1:length(data)
    for j=1:length(vY)
        data(i).obs(j,:) = data(i).obs(j,:) ./ sqrt(vY(j));
    end
end

%% Compute the Covariance matrix

kernel = 5.0;                                                       % Scale coefficient. The larger the value the more vague it makes the covariance. Typical range 0.5-5.
meanSigma = kernel * cov(diff(Ybig3'));  
for i=1:size(meanSigma,1)
    for j=1:size(meanSigma,2)
        if(i~=j)
            meanSigma(i,j) = 0;
        end
    end
end

end