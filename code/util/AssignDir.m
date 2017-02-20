%%*************************** Documentation **************************************
% Function used for Snap Assembly strategies. 
% Assigns appropriate strategy path based on the kind of approach and robot used.
%
% If new categories are created for HIRO/BAXTER, please do a global
% search for HSA/BSA and make sure that any string comparisons include the new
% categories. Need to improve the way we do this. There are several
% instances in the code, where parameters change according to the type of
% experiment that is being run.
%
% I.e. loadData.m, SnapData3.m, InsertStates3.m, WritePrimitivesToFile.m,
% cleanUp.m, WriteCompositesToFile.m, plotMotionCompositions.m,
% GradientClassification.m, CustomizePlotLength.m
%**********************************************************************************
function StratTypeFolder = AssignDir(StrategyType)

    % Assign a directory path based on the StrategyType used. 
    %--------------------- SIMULATION: Pivot Approach------------------------------------------------
    %------------------------ HIRO ------------------------------------------------
    
    %% Simulations: Straight Line Approach w/ HIRO
    if strcmp(StrategyType,'SIM_HIRO_SLA')
        StratTypeFolder = 'PositionControl/StraightLineApproach-NewIKinParams/';            % Straight Line with new IKin params
    
    elseif strcmp(StrategyType,'SIM_HIRO_SLA_NOISE')
        StratTypeFolder = 'PositionControl/StraightLineApproach-NewIkinParams-Noise/';      % Straight Line with new IKin params with noise
    
    elseif strcmp(StrategyType,'SIM_PA10_ONE_SL_SUCCESS')
        StratTypeFolder = 'ForceControl/SIM_StraightLineApproach/';                         % Used with PA10 Simulation

    %% Simulations: Pivot Approach w/ PA10/HIRO         
    elseif strcmp(StrategyType,'SIM_HIRO_ONE_PA_SUCCESS')
        StratTypeFolder = 'PositionControl/PivotApproach-NewIkinParams/';                   % Pivot approach with new IKin Params
    
    elseif strcmp(StrategyType,'SIM_HIRO_ONE_PA_NOISE')
        StratTypeFolder = 'PositionControl/PivotApproach-NewIKin-Noise/';                   % Pivot approach with new IKin Params with noise        
    
    elseif strcmp(StrategyType,'SIM_PA10_ONE_PA_SUCCESS')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         % Used with PA10 PivotApproach Simulation

    %% Simulations: Side Approach w/ HIRO/Baxter     
    elseif strcmp(StrategyType,'SIM_HIRO_ONE_SA_SUCCESS')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         % Used with HIRO SideApproach Simulation and Physical
    
    %% Simulations: Side Approach Error Characterization
    elseif strcmp(StrategyType, 'SIM_HIRO_ONE_SA_ERROR_CHARAC_LoopBack_x')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         % Used with HIRO SideApproach to compute error characteristics
    
    elseif strcmp(StrategyType, 'SIM_HIRO_ONE_SA_ERROR_CHARAC_LoopBack_y')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         
        
    elseif strcmp(StrategyType, 'SIM_HIRO_ONE_SA_ERROR_CHARAC_Prob')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         
    
    elseif strcmp(StrategyType, 'SIM_HIRO_ONE_SA_ERROR_CHARAC_SVM')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');
        
    elseif strcmp(StrategyType, 'SIM_HIRO_ONE_SA_FAILURE')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');
        
    %% Simulations: Dual Arm Side Approach 
    elseif strcmp(StrategyType, 'SIM_HIRO_TWO_SA_SUCCESS')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         % Used with simulations of the HIRO robot performing a side approach strategy with two arms: right and left under a coordination policy of push -hold.           

    %% Simulations: Dual Arm Side FAILURE 
    elseif strcmp(StrategyType, 'SIM_HIRO_TWO_SA_FAILURE')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');        
        
    %------------------------ REAL HIRO---------------------------------
    %% Robot: Side Approach
    elseif strcmp(StrategyType, 'HIRO_SideApproach')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         % Used with real hiro robot performing the side approach strategy.
    
    %% Robot: Side Approach Error Characterization
    elseif strcmp(StrategyType, 'REAL_HIRO_ONE_SA_ERROR_CHARAC')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');                         % Used with real hiro robot performing error deviations in the side approach strategy.
        
    elseif strcmp(StrategyType, 'REAL_HIRO_ONE_SA_FAILURE')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');        

    %------------------------ SIM BAXTER ---------------------------------------------
    elseif strcmp(StrategyType, 'SIM_BAXTER_ONE_SA_SUCCES')
        StratTypeFolder = strcat('PositionControl/',StrategyType,'/');

    elseif strcmp(StrategyType, 'SIM_BAXTER_SA_DUAL')
        StratTypeFolder = strcat('PositionControl/',StrategyType,'/');

    %------------------------ REAL BAXTER ---------------------------------------------
    elseif strcmp(StrategyType, 'REAL_BAXTER_ONE_SA_SUCCESS')
        StratTypeFolder = strcat('PositionControl/',StrategyType,'/');

    elseif strcmp(StrategyType, 'REAL_BAXTER_TWO_SA_SUCCESS')
        StratTypeFolder = strcat('ForceControl/',StrategyType,'/');
    
    else
        StratTypeFolder = '';
%        FolderName='';
    end
end
