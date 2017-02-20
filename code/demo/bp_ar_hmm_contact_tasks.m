% Modified EasyDemo to run BP_AR_HMM for contact tasks. 
% This demo first loads manipulation tasks data for the BP-HMM dataset,
%  and then runs fast BP-HMM inference and visualizes the results!
%
% Make sure you've done these simple things to run this script:
%   -- install Eigen C++ library
%   -- install Lightspeed toolbox
%   -- Compile the MEX routines for fast sampling (./CompileMEX.sh)
%   -- Create local directories for saving results (./ConfigToolbox.sh)
% See QuickStartGuide.pdf in doc/ for details on configuring the toolbox
function bp_ar_hmm_contact_tasks(StrategyType)

% ============================================================   GLOBALS
global DB_FILTER;
global DB_VISUALIZE;

DB_FILTER   =0;
DB_VISUALIZE=0;

% ============================================================   LOAD DATA!
fprintf( 'Load task data\n' );

% Set Path, Strategy, and get folders and data. 
path='/media/vmrguser/hdd/school/research/AIST/Results/'; % The path at which you want to save the main body of results. Folders will be created within this folder for different strategyTypes.
strat=AssignDir(StrategyType);
if(exist(path,'dir')~=7)
    fprintf('Data path does not exist. Please check your path.')
    return;
end
base_dir=strcat(path,strat);
folders = dir(base_dir);

% Extract data according to data types you want to use. Possible data:
dataSource = {'R_Torques'}; % {'R/L_Torques, R/L_CartPos_Corrected','R/L_Torques'};

% Ground Truth data
data_params=struct('flag',      1,...
                   'R',         1,...
                   'nStates',   4,...
                   'nDim',      6,...
                   'N',         10,...
                   'T',         0);
data_params.flag=1;
if(data_params.flag)
    % data_params_p = libpointer('struct',data_params);

    % Compile the data from all experiments for which we will compute parameters \theta.    
    [data, data_params] = extract_data_from_folders(base_dir, dataSource, data_params);
    % TODO: still need to recover data_params.T either through output or ptr
    
% ============================================================   DATA PRE-PROCESSING  
    % Change the mean of all trial data to zero and compute the covariance.
    [data,meanSigma]=data_pre_processing(data,data_params_p);

% ============================================================   GEN AR-HMM PARAMS    
    [data,TruePsi] = gen_ARGaussian(data,meanSigma,data_params);
    
end

% ============================================================   VISUALIZE  
% Visualize the raw data time series
%   with background colored by "true" hidden state
figure( 'Units', 'normalized', 'Position', [0.1 0.25 0.75 0.5] );
subplot(2, 1, 1 );
plotData( data, 1 );
subplot(2, 1, 2 );
plotData( data, 3 );

% Visualize the "true" generating parameters
% Feat matrix F (binary 5 x 4 matrix )
figure('Units', 'normalized', 'Position', [0 0.5 0.3 0.5] );
plotFeatMat( TruePsi.F );
title( 'True Feature Matrix', 'FontSize', 20 );

% Emission parameters theta (Gaussian 2D contours)
if(~flag)
    figure('Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5] );
    plotEmissionParams( TruePsi.theta, data );
    title( 'True Emission Params (with all data points)', 'FontSize', 20 );
end
pause;

% -------------------------------------------------   RUN MCMC INFERENCE!
modelP = {'bpM.gamma', 2}; % Set the model as a BetaProcess model. 2nd arg is dimension
algP   = {'Niter'               , 100, ...
          'HMM.doSampleHypers'  ,   0, ... 
          'BP.doSampleMass'     ,   0, ...
          'BP.doSampleConc'     ,   0}; 
% Start out with just one feature for all objects
initP  = {'F.nTotal', 1}; 
CH = runBPHMM( data, modelP, {1, 1}, algP, initP );
% CH is a structure that captures the "Chain History" of the MCMC
%  it stores both model config at each each iteration (in Psi field)
%             and diagnostic information (log prob, sampler stats, etc.)


% -------------------------------------------------   VISUALIZE RESULTS!
% Remember: the actual labels of each behavior are irrelevent
%   so there won't in general be direct match with "ground truth"
% For example, the true behavior #1 may be inferred behavior #4

% So we'll need to align recovered parameters (separately at each iter)
% Let's just look at iter 90 and iter 100

Psi90 = CH.Psi( CH.iters.Psi == 90 );
alignedPsi90 = alignPsiToTruth_OneToOne( Psi90, data );

Psi100 = CH.Psi( CH.iters.Psi == 100 );
alignedPsi100 = alignPsiToTruth_OneToOne( Psi100, data );


% Estimated feature matrix F
figure( 'Units', 'normalized', 'Position', [0 0.5 0.5 0.5] );
subplot(1,2,1);
plotFeatMat( alignedPsi90 );
title( 'F (@ iter 90)', 'FontSize', 20 );
subplot(1,2,2);
plotFeatMat( alignedPsi100 );
title( 'F (@ iter 100)', 'FontSize', 20 );

% Estimated emission parameters
figure( 'Units', 'normalized', 'Position', [0.5 0.5 0.5 0.5] );
subplot(1,2,1);
plotEmissionParams( Psi90 );
title( 'Theta (@ iter 90)', 'FontSize', 20 );
subplot(1,2,2);
plotEmissionParams( Psi100 );
title( 'Theta (@ iter 100)', 'FontSize', 20 );

% Estimated state sequence
plotStateSeq( alignedPsi100, [1 3] );
set( gcf, 'Units', 'normalized', 'Position', [0.1 0.25 0.75 0.5] );
subplotHandles = findobj(gcf,'type','axes');
title(min(subplotHandles), 'Est. Z : Seq 1', 'FontSize', 20 );
title(max(subplotHandles), 'Est. Z : Seq 3', 'FontSize', 20 );

fprintf( 'Remember: actual labels for behaviors are *irrelevant* from model perspective\n');
fprintf( '  what matters: *aligned* behaviors consistently assigned to same datapoints as ground truth\n' );
end