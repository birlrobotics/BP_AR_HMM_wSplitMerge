% ============================================================
% [data]=data_pre_processing(data,dp)
% Convert the data mean to 0 and compute the covariance.
% Inputs:
%   data:           structure that holds the observations for all included trials.
%                   this is a D*(\Sum_i^N T_i) struc, where T_i:time_steps for ith trial, D:dimensions, N:numTrials
%   data_params_p:  hold parameters that reflect properties of the data, pointer.
%------------------------------------------------------------------------------------------
% Outputs:
% data:             The observation mean will have been shifted to be equal to zero.
% meanSigma:        The first diff of the covariance set to 1.0
% ============================================================
function [mY,meanSigma]=data_pre_processing(data,dp)
%function [full_data_set,mY,meanSigma]=data_pre_processing(data,dp)

% Local variables
DB_FILTER=1;
cidx1=1; % column-wise indeces indicating start/finish of observations in each trial.
cidx2=0;

% Large row-vectors (6,dp.T)
full_data_set = zeros(dp.nDim,sum(dp.T));         
                                       % Ybig1 - Holds smoothed observations
                                       % Ybig2 - Holds variance
                                       % Ybig3 - Holds the first diff. of the covariance 
%% Smooth the Data: requires you to set a radius parameter for smoothing
if(DB_FILTER)
    radius = 10; % Default value 10

    for i=1:dp.N
        [r,c] = size(data(i).obs); % # rows not even for all trials
        smoothed = zeros(r,c);

        for j=1:r
            for k=1:c
                low  = max(1,k-radius);
                high = min(c,k+radius);
                smoothed(j,k) = mean(data(i).obs(j,low:high));
            end
        end
        data(i).obs = smoothed;
        
        cidx2                        =cidx2+length(data(i).obs);
        full_data_set(:,cidx1:cidx2) = data(i).obs;
        cidx1                        =cidx2+1;
    end
end
% Append data for each of the elements in the data strucuture in column-wise format
% for ii=1:length(data)
%     Ybig1 = [Ybig1 data(ii).obs];
% end

% Make copies of orig data for YBig1,2,3
Ybig1=full_data_set;
Ybig2=Ybig1;
Ybig3=Ybig1;
%% Adjust each dimesnion to Mean=0

mY = mean((Ybig1')); % computes mean of each matrix-column
for i=1:dp.N
    for j=1:dp.nDim
        data(i).obs(j,:) = data(i).obs(j,:) - mY(j);
    end
end

%% Adjust 1st diff variance of the observations to set them to 1.0
% In general, it seems that this will help 
vY = var(diff(Ybig2'));
for i=1:length(data)
    for j=1:length(vY)
        data(i).obs(j,:) = data(i).obs(j,:) ./ sqrt(vY(j));
    end
end

%% Compute the Covariance matrix

kernel = 0.5;                                                       % Scale coefficient is the gain by which we multiple the empirical covariance of data set. If covariance computed from pooling of all data overestimates the mode-specific covariance, it is good to downscale and vice-versa. See p. 156 of Fox PhD Thesis.
meanSigma = kernel * cov(diff(Ybig3'));                             % Covariance currently 6x6. Covariance of Fx-Mz across trials? Covariance of diff almost all zero. Covariance of orig data has stronger correlations.
for i=1:size(meanSigma,1)
    for j=1:size(meanSigma,2)
        if(i~=j)
            meanSigma(i,j) = 0; %% SHOULD WE TAKE THIS DEPENDENCY OFF??
        end
    end
end

end