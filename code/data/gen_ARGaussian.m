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
%            dp.T        = length of each/multiple time series
%            dp.R        = order of the autoregressive process
% OUTPUT ----------------------------------------------------------
%    data  :  ARSeqData object
%               note that each sequence will *actually* have dp.T-dp.R
%               observations, since we need "dp.R" to properly
%                define the likelihood of the first "kept" observation
%   PsiTrue : consists of existing features data
%             could have mean/sigma.
%             sigma is computed from nu_delta which is created from
%             meanSigma
% ------------------------------- Remember old state to use again afterward
% ============================================================
function [data,PsiTrue] = gen_ARGaussian(raw_data,mu,meanSigma,dp)

curStream = RandStream.getGlobalStream();
entryState = curStream.State;

% Reset PRNG state to default value with SEED 0
%       so that we always get same synth data regardless of when called
reset( RandStream.getGlobalStream(), 0);

%% Check Length of Time Series
if min(dp.T) < 0
    doVaryLength = 1;
    dp.T = abs(min(dp.T));
else
    doVaryLength = 0;
end

total_time_steps=sum(dp.T);

%% Self-Transition Probability
% Set the probability that any given feature is included in the initial
% feature matrix.
pIncludeFeature = 1.0; % We assume all of our features will happen in this case. Orig value 0.75
pSelfTrans      = 1-(2*dp.nStates)/min(dp.T); % Usually high prob. Look at ratio of StateNum with Time-steps.

% Create initial state distribution 
Pz_init = [1 zeros(1,dp.nStates-1)]; % Original state distribution was Pz_init ones(1, dp.nStates);
Pz      = zeros(  dp.nStates, dp.nStates ); % Create state-to-state transition matrix Pz

%% Original Transition Matrix
% for k = 1:dp.nStates
%     Pz(k,k) = pSelfTrans;
%     Pz(k, [1:k-1 k+1:end] ) = (1-pSelfTrans)/(dp.nStates-1);
% end

% Creating a left-to-right transition matrix
for k = 1:dp.nStates
    Pz(k,k) = pSelfTrans;
end
for k = 1:dp.nStates
    Pz(k,k) = pSelfTrans;
    if k~=dp.nStates
        Pz(k, k+1) = (1-pSelfTrans);
    else
        Pz(k,k)=1;
    end
end

% Check validity of Pz
doRowsSumToOne = ( sum(Pz,2) - ones(size(Pz,1),1) ) <= 1e-10;
assert( all(doRowsSumToOne), 'ERROR: Not a valid transition distr.' );

%--------------------------------------------------------------------------
%% Placeholder for state-specific emission params Px
%--------------------------------------------------------------------------
%   - Means:        obtained from raw data. 
%   - Covariances   obtained from raw data. 
A    = zeros( dp.nDim, dp.nDim*dp.R, dp.nStates);
Sigma=zeros( dp.nDim, dp.nDim*dp.R, dp.nStates);

nu       = dp.nDim+2;                   % inverse Wishart degrees of freedom, nu: typically 2.
nu_delta = (nu-dp.nDim-1)*meanSigma;    % used to compute the inverse wishart Sigma

% Covariance amtrix
% My2DSigs(:,:,1) = [15 0; 0 0.5];
% My2DSigs(:,:,2) = [10 5; 5 3];
% My2DSigs(:,:,3) = [0.5 0; 0 15];
% My2DSigs(:,:,4) = [3 -5; -5 10];

% Variables that are later used to create the Sigma matrix in the VAR process
MyVarCoefs = repmat( [0.001 0.1 1 2 3 4 5 10 15], 1, ceil(dp.nDim/5) );

if dp.R > 1
    muteVec = logspace( 0, -dp.R/2, dp.R);
else
    muteVec = 1;
end

% Instantiate the A matrix
m=dp.nDim*dp.R;
for kk = 1:dp.nStates
    A(:,:,kk)     = inv(diag( 5*ones(1,m) ));   % How might one know how to set up these coefficients to start of?
    Sigma(:,:,kk) = iwishrnd( nu_delta,  nu  ); % Sigma depends on nu_delta which in turn depends on meanSigma, which was derived from the raw data.
%     for dd = 1:dp.nDim
%          A(dd, (dd:dp.nDim:dp.nDim*dp.R), kk) = as(kk)*muteVec;         
%     end
%     Sigma(:,:,kk) = iwishrnd( nu_delta,  nu  );
%     if dp.nDim == 1
%         % IW(0.5, 3)
%         Sigma(:,:,kk) = iwishrnd( nu_delta,  nu  );
%     else
%         Sigma(1:2, 1:2, kk) = My2DSigs(:,:, mod(kk-1,4)+1 );
%         Sigma(3:dp.nDim, 3:dp.nDim, kk) = diag( randsample( MyVarCoefs, dp.nDim-2 ) );
%     end
end

%--------------------------------------------------------------------------
%% MNIW Parameters 
%--------------------------------------------------------------------------
% Mean and covariance for VAR dynamic parameters A,\Sigma. Both are
% uncertain. MNIW is the conjugate prior to the AR likelihood. 
% The prior will allow for tractable posterior inference. 
% A is a linear function on Yt_1, that updates to the mean of the signal
% observed at Yt. Sigma diagonal??
Px.A     = A;
Px.Sigma = Sigma;

% Instantiate an ARSeqData object 
data = ARSeqData( dp.R );
F = zeros( dp.N, dp.nStates );

for i = 1:dp.N
    if doVaryLength
        dp.T(i) = poissrnd(dp.T(i));
    else
        T = dp.T(i);
    end

    % Create structures for X, Xprev, and zTrue
    zTrue = zeros(1, T);
    X     = zeros( dp.nDim, T );
    xprev = zeros( dp.nDim*dp.R, 1 );
    
    % Draw subset of states that this time-series exhibits
    mask = rand( 1, dp.nStates ) < pIncludeFeature;
    % Ensure mask isn't all zeros
    if sum( mask ) < 1
        kk = randsample( dp.nStates, 1);
        mask(  kk  ) = 1;
    end
    F(i,:) = mask;
    
    % Ensure mask isn't all zeros
    if sum( mask ) < 1
        kk = randsample( dp.nStates, 1);
        %kk = multinomial_single_draw( ones( Ktrue,1 ) );
        mask(  kk  ) = 1;
    end
    F(i,:) = mask;
       
    for t = 1:T
        % ---------------------------------------------- Assign true label
        if t == 1
            zcur = multinomial_single_draw( mask.*Pz_init );
        else
            zprev = zTrue(t-1);
            zcur = multinomial_single_draw( mask.*Pz(zprev,:) );
        end
        zTrue(t) = zcur;
%--------------------------------------------------------------------------------   
%         Another way to assign zTrue??        
%         Ninit = 4; %2
%         init_blocksize = floor(T/Ninit);
%         z_init = [];
%         for i=1:Ninit
%             z_init = [z_init i*ones(1,init_blocksize)];      % Z is as long as observations. But, it is indexed by a state number on a per-bloc basis.
%         end
%         z_init(Ninit*init_blocksize+1:T) = Ninit;        
%--------------------------------------------------------------------------------        
        % ---------------------------------------- Assign emitted data
        X(:,t) = mvnrnd( Px.A(:,:,zcur)*xprev,  Px.Sigma(:,:,zcur) );
        xprev = [X(:,t); xprev( 1:(end-dp.nDim) ) ];
        
    end

    % The "addSeq" function appends current X and zTrue to the
    % corresponding elements of the data object
    data = data.addSeq( X, num2str(i), zTrue );
end

% ---------------------------------------------------------  Reset stream
curStream = RandStream.getGlobalStream();
curStream.State = entryState;

% ---------------------------------------------------------  Compute Psi??
% This computation was nor originally in AR. Should we include it?
PsiTrue.F = zeros(dp.N, dp.nStates);
for ii = 1:dp.N
    PsiTrue.F(ii, unique( data.zTrue(ii) ) ) = 1;
end
% Below section further commented out b/c we do not have a Px struct.
% for kk = 1:dp.nStates
%     PsiTrue.theta(kk).mu = Px.Mu(kk,:); 
%     PsiTrue.theta(kk).invSigma = inv( Px.Sigma(:,:,kk) );
% end
% PsiTrue.Pz = Pz;
% PsiTrue.z = zTrue;

end % main function



