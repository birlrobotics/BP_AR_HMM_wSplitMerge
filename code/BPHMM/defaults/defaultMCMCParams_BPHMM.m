function algP = defaultMCMCParams_BPHMM()
% Creates a struct encoding the default settings for MCMC inference

algP.Niter = 500;

algP.doSampleFShared = 1;
algP.doSampleFUnique = 1;
algP.doSampleUniqueZ = 0;
algP.doSplitMerge    = 1; % (Split-Merge Move - high performance)
algP.doSplitMergeRGS = 0; 
algP.doSMNoQRev = 0; % ignore proposal prob in accept ratio... not valid!


algP.SM.featSelectDistr = 'splitBias+margLik';
algP.SM.doSeqUpdateThetaHatOnMerge = 0;

algP.nSMTrials = 5;
algP.nSweepsRGS = 5;

algP.doSampleZ      = 1;
algP.doSampleTheta  = 1;
algP.doSampleEta    = 1;

algP.BP.doSampleMass = 1;
algP.BP.doSampleConc = 1;
algP.BP.Niter        = 10;
algP.BP.var_c        = 2; % Part of BetaProcess. Affects generations of weights \omega_i

algP.HMM.doSampleHypers = 1;
algP.HMM.Niter          = 20;
algP.HMM.var_alpha      = 2;
algP.HMM.var_kappa      = 10;

% Reversible Jump Proposal settings
algP.RJ.doHastingsFactor = 1;
algP.RJ.birthPropDistr   = 'DataDriven'; % Great results in combination with SM in 2014paper.
algP.RJ.minW = 15; 
algP.RJ.maxW = 100;

algP.doAnneal = 0;
algP.Anneal.T0 = 100;
algP.Anneal.Tf = 10000;
