% ============================================================
% [data,PsiTrue] = gen_ARGaussian(data,meanSigma,dp)
% INPUTS ----------------------------------------------------------
% data:      structure that holds the observations for all included trials.
%            this is a T*(D*N) struc, where T:time_steps, D:dimensions, N:numTrials
% meanSigma: first diff of covariance set to 1.0
% dp:        hold parameters that reflect properties of the data, pointer.
%            dp structure, that includes the following information:
%            dp.nStates  = # of available Markov states
%            dp.nDim     = number of observations at each time instant
%            dp.N        = number of time series objects
%            dp.T        = length of each time series
%            dp.R        = order of the autoregressive process
% OUTPUT ----------------------------------------------------------
%    data  :  ARSeqData object
%               note that each sequence will *actually* have ar.T-ar.R
%               observations, since we need "ar.R" to properly
%                define the likelihood of the first "kept" observation
%   PsiTrue : consists of existing features data
%             could have mean/sigma.
% ------------------------------- Remember old state to use again afterward
% ============================================================
function [data,PsiTrue] = gen_ARGaussian(data,meanSigma,dp)

curStream = RandStream.getGlobalStream();
entryState = curStream.State;

% Reset PRNG state to default value with SEED 0
%       so that we always get same synth data regardless of when called
reset( RandStream.getGlobalStream(), 0);

if ar.T < 0
    doVaryLength = 1;
    ar.T = abs(ar.T);
else
    doVaryLength = 0;
end

% Set the probability that any given feature is included in the initial
% feature matrix.
pIncludeFeature = 0.75;

% Usually high prob. Look at ratio of StateNum with Time-steps.
pSelfTrans = 1-(2*ar.nStates)/ar.T;

% Create initial state distribution (uniform)
Pz_init = ones(1, ar.nStates);

% Create state-to-state transition matrix Pz
Pz = zeros(  ar.nStates, ar.nStates );
for k = 1:ar.nStates
    Pz(k,k) = pSelfTrans;
    Pz(k, [1:k-1 k+1:end] ) = (1-pSelfTrans)/(ar.nStates-1);
end

% Check valididit of Pz
doRowsSumToOne = ( sum(Pz,2) - ones(size(Pz,1),1) ) <= 1e-10;
assert( all(doRowsSumToOne), 'ERROR: Not a valid transition distr.' );

% Create state-specific emission params Px
%   Means are evenly spaced around the unit circle
%   Covariances are aligned so major diagonal of ellipsoid
%        points toward the origin

as = linspace( -0.9, 0.9, ar.nStates );
A  = zeros( ar.nDim, ar.nDim*ar.R, ar.nStates);

nu = ar.nDim+2;                % inverse Wishart degrees of freedom, nu
% meanSigma = 0.5*eye(ar.nDim);  % inverse Wishart mean covariance matrix
nu_delta = (nu-ar.nDim-1)*meanSigma; % inverse wishart Sigma

% Covariance amtrix
My2DSigs(:,:,1) = [15 0; 0 0.5];
My2DSigs(:,:,2) = [10 5; 5 3];
My2DSigs(:,:,3) = [0.5 0; 0 15];
My2DSigs(:,:,4) = [3 -5; -5 10];

% Variables that are later used to create the Sigma matrix in the VAR process
MyVarCoefs = repmat( [0.001 0.1 1 2 3 4 5 10 15], 1, ceil(ar.nDim/5) );

if ar.R > 1
    muteVec = logspace( 0, -ar.R/2, ar.R);
else
    muteVec = 1;
end

% Create the A matrix
for kk = 1:ar.nStates
    
    for dd = 1:ar.nDim
         A(dd, [dd:ar.nDim:ar.nDim*ar.R], kk ) = as(kk) .* muteVec;
    end
    
    if ar.nDim == 1
        % IW(0.5, 3)
        Sigma(:,:,kk) = iwishrnd( nu_delta,  nu  );
    else
        Sigma(1:2, 1:2, kk) = My2DSigs(:,:, mod(kk-1,4)+1 );
        Sigma(3:ar.nDim, 3:ar.nDim, kk) = diag( randsample( MyVarCoefs, ar.nDim-2 ) );
    end
end

Px.A = A;
Px.Sigma = Sigma;

% Instantiate an ARSeqData object 
data = ARSeqData( ar.R );
F = zeros( ar.N, ar.nStates );
for i = 1:ar.N
    if doVaryLength
        Ti = poissrnd(ar.T);
    else
        Ti = ar.T;
    end
    
    % Draw subset of states that this time-series exhibits
    mask = rand( 1, ar.nStates ) < pIncludeFeature;
    % Ensure mask isn't all zeros
    if sum( mask ) < 1
        kk = randsample( ar.nStates, 1);
        mask(  kk  ) = 1;
    end
    F(i,:) = mask;
    
    % Ensure mask isn't all zeros
    if sum( mask ) < 1
        kk = randsample( ar.nStates, 1);
        %kk = multinomial_single_draw( ones( Ktrue,1 ) );
        mask(  kk  ) = 1;
    end
    F(i,:) = mask;
    
    
    zTrue = zeros(1,Ti);
    X = zeros( ar.nDim, Ti );
    xprev = zeros( ar.nDim*ar.R, 1 );
    for t = 1:Ti
        % ---------------------------------------------- Assign true label
        if t == 1
            zcur = multinomial_single_draw( mask.*Pz_init );
        else
            zprev = zTrue(t-1);
            zcur = multinomial_single_draw( mask.*Pz(zprev,:) );
        end
        zTrue(t) = zcur;
        
        % ---------------------------------------- Assign emitted data
 
        X(:,t) = mvnrnd( Px.A(:,:,zcur)*xprev,  Px.Sigma(:,:,zcur) );
        xprev = [X(:,t); xprev( 1:(end-ar.nDim) ) ];
        
    end
    data = data.addSeq( X, num2str(i), zTrue );
end

% ---------------------------------------------------------  Reset stream
curStream = RandStream.getGlobalStream();
curStream.State = entryState;

% ---------------------------------------------------------  Compute Psi??
% This computation was nor originally in AR. Should we include it?
PsiTrue.F = zeros(ar.N, ar.nStates);
for ii = 1:ar.N
    PsiTrue.F(ii, unique( data.zTrue(ii) ) ) = 1;
end
% Below section further commented out b/c we do not have a Px struct.
% for kk = 1:ar.nStates
%     PsiTrue.theta(kk).mu = Px.Mu(kk,:); 
%     PsiTrue.theta(kk).invSigma = inv( Px.Sigma(:,:,kk) );
% end
% PsiTrue.Pz = Pz;
% PsiTrue.z = zTrue;

end % main function



