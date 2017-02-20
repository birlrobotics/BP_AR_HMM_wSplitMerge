function initParams = defaultInitMCMC_BPHMM()

% ---------------------------------------------  Initializations
initParams.InitFunc = @initBPHMMFresh; % See more choices in /HOME/code/BPHMM/init

initParams.F.nUniquePerObj  = 4; % How many features are shared across all data. 1 is the lowest number. 

initParams.z.doPartition    = 1; % 0 := init from prior

