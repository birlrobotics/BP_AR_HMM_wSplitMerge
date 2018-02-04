function [] = RunMocap6Experiment( jobID, taskID, infName, initName, TimeLimit )
%INPUT
%  jobID : integer/name of job
%  taskID : integer/name of task
%  infName  : {'Prior','DD','SM'} indicates type of sampler     
%  initName : {'one', 'unique5'}  indicates initialization
%  TimeLimit : # of seconds to allow for MCMC


dataP = {'Mocap6'};
modelP = {};

TimeLimit = force2double( TimeLimit );
T0 = 50;
Tf = 10000;
switch infName
    case 'Prior'
        algP = {'doSampleFUnique', 1, 'doSplitMerge', 0, 'RJ.birthPropDistr', 'Prior'};         
    case {'DD', 'DataDriven'}
        algP = {'doSampleFUnique', 1, 'doSplitMerge', 0, 'RJ.birthPropDistr', 'DataDriven'};  
    case {'zDD'}
        algP = {'doSampleFUnique', 0, 'doSampleUniqueZ', 1, 'doSplitMerge', 0, 'RJ.birthPropDistr', 'DataDriven'};       

    case 'SM'                
        algP = {'doSampleFUnique', 0, 'doSplitMerge', 1};         
    case 'SM+DD'               
        algP = {'doSampleFUnique', 1, 'doSampleUniqueZ', 0, 'doSplitMerge', 1, 'RJ.birthPropDistr', 'DataDriven'};                      
    case 'SM+DD+AnnealLin'        
        algP = {'doSampleFUnique', 1, 'doSampleUniqueZ', 0, 'doSplitMerge', 1, 'RJ.birthPropDistr', 'DataDriven', 'doAnneal', 'Lin', 'Anneal.T0', T0, 'Anneal.Tf', Tf};               

    case 'SM+zDD'               
        algP = {'doSampleFUnique', 0, 'doSampleUniqueZ', 1, 'doSplitMerge', 1, 'RJ.birthPropDistr', 'DataDriven'};       
    case 'SM+zDD+AnnealExp'        
        algP = {'doSampleFUnique', 0, 'doSampleUniqueZ', 1, 'doSplitMerge', 1, 'RJ.birthPropDistr', 'DataDriven', 'doAnneal', 'Exp', 'Anneal.T0', T0, 'Anneal.Tf', Tf};               
    case 'SM+zDD+AnnealLin'        
        algP = {'doSampleFUnique', 0, 'doSampleUniqueZ', 1, 'doSplitMerge', 1, 'RJ.birthPropDistr', 'DataDriven', 'doAnneal', 'Lin', 'Anneal.T0', T0, 'Anneal.Tf', Tf};               


    % NON VALID OPTIONS (for experiments only)
    case 'SMnoqrev'               
        algP = {'doSampleFUnique', 0, 'doSMNoQRev', 1}; 
    case 'SMnoqrev+DD'               
        algP = {'doSampleFUnique', 1, 'doSMNoQRev', 1, 'theta.birthPropDistr', 'DataDriven'}; 
    otherwise
        error( 'unrecognized inference type' );
end
algP(end+1:end+2) = {'TimeLimit', TimeLimit};


switch initName
    case 'cheat' % use a known set of ground-truth labels to initialize the state sequence
        initP = {'InitFunc', @initBPHMMCheat, 'Cheat.nRepeats', 1};
    case {'seq','seq2'}
        initP = {'InitFunc', @initBPHMMSeq, 'F.nUniquePerObj', 2, 'nSubsetObj', 20};
    case 'one'
        initP = {'F.nTotal', 1};
    case 'unique5'
        initP = {'F.nUniquePerObj', 5};
end

runBPHMM( dataP, modelP, {jobID, taskID}, algP, initP );
