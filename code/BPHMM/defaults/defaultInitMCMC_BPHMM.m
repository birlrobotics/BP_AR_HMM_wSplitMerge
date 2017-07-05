function initParams = defaultInitMCMC_BPHMM()

% ---------------------------------------------  Initializations
initParams.InitFunc = @initBPHMMFresh; % See more choices in /HOME/code/BPHMM/init
% Initialize BP-HMM hidden variables for MCMC sampling as follows
%   1) fix F to user-specified configuration each object has U unique features
%   2) sample TS from prior given F
%   3) sample state asgns Z given F, TS
%             object i has U_i features in F
%             divide z_i into U_i equal sized contiguous blocks
%                 assign each to a distinct element of the U_i features
%   4) sample emit params Theta from posterior given stateSeq, data

initParams.F.nUniquePerObj  = 4; % How many features are shared across all data. 1 is the lowest number. 

initParams.z.doPartition    = 1; % 0 := init from prior

