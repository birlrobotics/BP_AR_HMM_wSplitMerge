function model = defaultModelParams_BPHMM( data )

switch data.getObsType()
    case 'Gaussian'                            
        model.obsM.precMu   = 1;            % precMu            : precision parameter for mean
        model.obsM.degFree  = 3;            % degFree           : # deg of freedom
        model.obsM.Scoef    = 1;            % Scoef             : scalar multiplier for Scale Matrix
        model.obsM.doEmpCovScalePrior = 0;  % doEmpCovScalePrior: boolean.
                                            %                     true  = set ScaleMat to Scoef*emp covariance
                                            %                     false = set ScaleMat to Scoef*eye(D)
    case 'AR'
        model.obsM.doEmpCov             = 0;
        model.obsM.doEmpCovFirstDiff    = 1;
        model.obsM.degFree              = data.D+2; % D+2. Start with 6 degrees of freedom. Could go to 12.
        if ~isempty( strmatch( '13_29',data.seqNames )  )
            model.obsM.Scoef            = 0.5; % Default=0.5 Previous work used 5. With better inference, need not be so vague
        else
            model.obsM.Scoef            = 0.5; % Default=0.5 Previous work used 5. With better inference, need not be so vague
        end
end

% ------------------------------- HMM params
model.hmmM.alpha = 1;             % HMM Transition Model: Prior
model.hmmM.kappa = 50;            % Default 25
model.hmmM.prior.a_alpha = 0.01;
model.hmmM.prior.b_alpha = 0.01;
model.hmmM.prior.a_kappa = 0.01;
model.hmmM.prior.b_kappa = 0.01;

% ================================================== BETA PROCESS MODEL
% ------------------------------- GAMMA: Mass param for IBP. See 2014OAOS Supplement paper, Appendix F.1. Sampling IBP Hyperparameters
model.bpM.gamma = 5;           % IBP Mass Parameter
model.bpM.prior.a_mass = 0.01; % Used with the two-parameter version of the IPB
model.bpM.prior.b_mass = 0.01;

% ------------------------------- c0   : Concentration param for IBP
model.bpM.c = 1; 
model.bpM.prior.a_conc = 0.01;
model.bpM.prior.b_conc = 0.01;


