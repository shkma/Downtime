%--------------------------------------------------
%---- Compute the rate of exceedance of Sa(T) -----
%--------------------------------------------------
%%
% Time unit: Days

clear, clc

NumIM       =  10;       % Number of intensity measure (return period)
NumGM       =  40;       % Number of ground motions per each RP (Return Period)
NumStory    =  6;        % Number of stories of building
NumFloor    =  7;        % Number of floors of building including ground floor and roof
g           = 386.0;     % Gravity acceleration (in/sec^2)

wh = (8.0 + 7.0) /24.0;  % workday hours (within a day) (8 hrs for day-time crew; 7 hrs for night-time crew) per Mitrani-Resier (2007)
wr = 5.0 / 7.0;          % workday ratio (Mon-Fri: Work days; Sat-Sun: Holidays) per Dong and Frangopol (2016)
cn = 15.0;               % number of crews (per floor) - 15 per floor per Dong and Frangopol (2016)

R_COT_Fast  =  2.0;      % Change-of-trade delay (assumed to be the same for each trade)
R_COT_Slow  = 14.0;      % Change-of-trade delay (assumed to be the same for each trade)

ave_lease_rate = 1.33  /30.0;   % Average lease rate of $1.33/ft2/month, modified to the unit of day (month->day)
E_U_U_m        = ave_lease_rate*18900.0;            % Mean rent per operational unit
E_DTL_TAG      = ave_lease_rate*18900.0 * NumStory; % Mean rent for all operational units
E_DTL_Dem = E_DTL_TAG; E_DTL_Col = E_DTL_TAG;

y_drift_Col = 0.05;    % 0.10; % Story drift ratio that causes collapse
Disp_limit  = 26.9;    % 20.7; % 32.6; % 38.3; % inch; for isolator displacement, D_Ultimate (collapse)

RP          = [43, 144, 289, 475, 949, 1485, 2475, 3899, 7462, 10000]; % Return period

% Specified for different systems%% --------------------------------------------------------------------------------
system        = 4;     % '1'=Nonisolated SCBF;                   '2'=Nonisolated SMRF;
                       % '3'=Isolated SCBF (RI=1) (Lower Bound); '4'=Isolated SCBF (RI=2) (Lower Bound);
                       % '5'=Isolated SMRF (RI=1) (Lower Bound); '6'=Isolated SMRF (RI=2) (Lower Bound);
                       % '7'=Isolated Upper Bound
                       
iso_size      = 1;     % '1'=TFP-1 or DC-1;
                       % '3'=TFP-3;
                       
wall_or_none  = 0;     % '0'=no moat wall;
                       % '1'=moat wall;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     (system == 1) || (system == 2)
        ImpRows_Drift       = 2:7;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
        ImpRows_AccelAll    = 2:8;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
elseif (system == 3) || (system == 4) || (system == 5) || (system == 6) || (system == 7)
        ImpRows_Drift       = 3:8;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
        ImpRows_AccelAll    = 3:9;  % IMPORTANT: PGA is not considered for isolated structures,, elevator is sensitive to PFA at first floor (ie., isolated)
        ImpRows_TFPdisp     = 2:3;  % Rows that contain the seismic response (max response is taken from these rows; change this per structural system & EDP)
end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Values used to calculate collapse probability (Baker JW,2015) ---------
num_gms                   = [NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM,NumGM];
returnPeriod              = [43,144,289,475,949,1485,2475,3899,7462,10000];
num_collapse              = zeros( 1, NumIM );       % Vector for having number of collapse 
num_failure_StoryDrift    = zeros( 1, NumIM );       % Vector for having number of failure of story drift
num_failure_NumericalProb = zeros( 1, NumIM );       % Vector for having number of numerical failure
MaxDriftVector            = zeros( NumGM, NumIM );   % Vector that install maximum drift values to be checked how many collapse happens out of 40 GMs
MaxResDriftVector         = zeros( NumGM, NumIM );   % Vector that install maximum residual drift values to be checked how many collapse happens out of 40 GMs
MaxAccelFloorVector       = zeros( NumGM, NumIM );   % Vector that install maximum floor accel values to be checked how many collapse happens out of 40 GMs
MaxAccelRoofVector        = zeros( NumGM, NumIM );   % Vector that install maximum roof accel values to be checked how many collapse happens out of 40 GMs
MaxAccelAllVector         = zeros( NumGM, NumIM );   % Vector that install maximum roof accel values to be checked how many collapse happens out of 40 GMs

n_NumericalProb_matrix    = zeros( NumGM, NumIM );   % Matrix that install identification of if there is a numerical problem or not
n_Collapse_matrix         = zeros( NumGM, NumIM );   % Matrix that install identification of if there is a collapse or not
n_StoryDrift_matrix       = zeros( NumGM, NumIM );   % Matrix that install identification of if there is a excessive story drift or not

% Store Sa(T1) for different systems ------------------------------------
if system == 1      % Nonisolated SCBF
    Sa_T1     = [0.31,0.63,0.88,1.09,1.39,1.60,1.85,2.19,2.63,2.85];
    Period_T  = 0.524;
    Sa_MCE_T1 = 1.500;        % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;          % Modeling uncertainyy
    beta_TD   = 0.2;          % Test data uncertainyy
    beta_DR   = 0.2;          % Design requirement
    HazardCurve_Name = 'SeismicHazardData_0.524sec.txt';  % Specify name of file containing seismic hazard data
elseif system == 2  % Nonisolated SMRF
    Sa_T1     = [0.16,0.36,0.54,0.68,0.90,1.06,1.25,1.44,1.74,1.88];
    Period_T  = 1.186;
    Sa_MCE_T1 = 0.756302521;  % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;          % Modeling uncertainyy
    beta_TD   = 0.2;          % Test data uncertainyy
    beta_DR   = 0.1;          % Design requirement
    HazardCurve_Name = 'SeismicHazardData_1.186sec.txt';  % Specify name of file containing seismic hazard data
elseif (system == 3) || (system == 4) || (system == 5) || (system == 6)  % Isolated (lower bound)
    Sa_T1     = [0.02,0.06,0.10,0.16,0.24,0.31,0.39,0.44,0.53,0.58];
    Period_T  = 3.660;
    Sa_MCE_T1 = 0.246;  % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;    % Modeling uncertainyy
    beta_TD   = 0.2;    % Test data uncertainyy
    beta_DR   = 0.1;    % Design requirement
    HazardCurve_Name = 'SeismicHazardData_3.660sec.txt';  % Specify name of file containing seismic hazard data
elseif system == 7  % Isolated
    Sa_T1     = [0.04,0.10,0.17,0.24,0.35,0.42,0.50,0.57,0.69,0.75];
    Period_T  = 2.990;
    Sa_MCE_T1 = 0.301;  % Spectral acceleration at T1 or TM
    beta_MDL  = 0.2;    % Modeling uncertainyy
    beta_TD   = 0.2;    % Test data uncertainyy
    beta_DR   = 0.1;    % Design requirement
    HazardCurve_Name = 'SeismicHazardData_2.990sec.txt';  % Specify name of file containing seismic hazard data
end

% Construct collapse fragility curve --------------------------------------
% & store maximum residual drift along the height -------------------------
% & store maximum story drift at each story -------------------------------
% & store maximum PGA, PFA & PRF at each floor ----------------------------
fprintf(strcat('Compute collapse fragility curve - - - - - - - - - - - - - -  \n'));
for i_RP = 1 : NumIM
    dirname  = strcat('ReturnPeriodID_', num2str(i_RP), '_dynamicData');
    
    % Display the progress of processing data to generate collapse fragility curve
    fprintf(strcat('~ Progress(1): \\', num2str(i_RP), 'th RP out of \\', num2str(NumIM), ' RP is currently being processed ~\n'));
    
    % Location of folder that contains the data for analysis (change this depending on location of analysis files)
    if     i_RP == 1  % (RP=   43 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=43');
    elseif i_RP == 2  % (RP=  144 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=144');
    elseif i_RP == 3  % (RP=  289 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=289');
    elseif i_RP == 4  % (RP=  475 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=475');
    elseif i_RP == 5  % (RP=  949 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=949');
    elseif i_RP == 6  % (RP=  1485 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=1485');
    elseif i_RP == 7  % (RP=  2475 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=2475');
    elseif i_RP == 8  % (RP=  3899 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=3899');
    elseif i_RP == 9  % (RP=  7462 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=7462');
    elseif i_RP == 10 % (RP= 10000 yrs)
       DataFolder = strcat('PEER_GroundMotionData_T=', num2str(Period_T,'%5.3f'), 'sec_PR=10000');
    end
        
        for i_gm = 1 : NumGM
            data_Drift      = load( strcat(char(dirname), '/GM_', num2str(i_gm), '_DriftStories.out') ); % Load story drift data (also for check if there's a collapse)
            data_Accel      = load( strcat(char(dirname), '/GM_', num2str(i_gm), '_AccelStories.out') ); % Load floor accel data
            data_TFPdisp    = load( strcat(char(dirname), '/GM_', num2str(i_gm), '_TFPdisp.out') );      % Load triple FP disp data.. Seismic frame
            endTime_Anal    = data_Drift(end,1);  % End (last) time of the time history response data (to be used to check if analysis completes or not; if not=collapse)
            
            % Collapse if defined as exceeding limit of story drift ratio (5% for SCBF) or
            % if there is numerical stability problem (analysis terminates)
                
                % % Load 'NumberGMpoint.txt' to be used in parametric study
                load( strcat('../', strcat(char(DataFolder), '/', 'NumberGMpoint.txt') ) );
                % % Load 'TimeStepGM.txt' to be used in parametric study
                load( strcat('../', strcat(char(DataFolder), '/', 'TimeStepGM.txt') ) );
                % % Calculate max duration of ground motion
                endTime_GM = NumberGMpoint(i_gm) * TimeStepGM(i_gm); % Max duration of ground motion
                
                % Load horizontal displacmeent of isolator (to be checked if it exceeds D_Ultimate)
                MaxTFPdisp = max( max( abs(data_TFPdisp(:, ImpRows_TFPdisp)) ) );    % Obtain max horizontal TFP response (TFPdisp at SF)
                
                % Collapse is defined either by exceeding y_drift_Col OR D_Ultimate in isolator OR Uplift_limit of uplift
                temp_MaxDrift = max( max( abs(data_Drift(:, ImpRows_Drift)) ) );     % Obtain max response (PSDR)
                
                % % For Collapse - - - - -
                    if endTime_Anal < (endTime_GM - 0.0) || MaxTFPdisp > Disp_limit || temp_MaxDrift > y_drift_Col % || MaxTFPdisp > Disp_limit || MaxUpliftdisp > Uplift_limit
                       n_Collapse_matrix(i_gm,i_RP) = 1;                             % identification that there was a collapse
                    else
                       n_Collapse_matrix(i_gm,i_RP) = 0;                             % identification that there wasn't a collapse
                    end
                  
%                 % % For Collapse (excluding collapse at isolator)- - - - -
%                     if temp_MaxDrift > y_drift_Col || MaxTFPdisp > Disp_limit || MaxUpliftdisp > Uplift_limit
%                        n_Collapse_matrix_noIso(i_gm,i_RP) = 1;                                 % identification that there was a collapse
%                     else
%                        n_Collapse_matrix_noIso(i_gm,i_RP) = 0;                                 % identification that there wasn't a collapse
%                     end
                    
                % % For max. residual story drift along the height of the
                % building given IM=im (=return period of i_RP)  - - - - -
                MaxResDriftVector(i_gm,i_RP)= max( abs(data_Drift(end, ImpRows_Drift)) ); % Obtain max response (RSDR)

                % % For peak seismic response at each story or floor - - - - -
                
                for i_story = 1 : length(ImpRows_Drift)
                    % % For story based peak response - - - - - (story drift)
                    if     (system == 1) || (system == 2) % i.e., non-isolated buildings
                        eval(['MaxDriftVector_Story_',num2str(i_story),'(i_gm,i_RP)', '= max( max( abs(data_Drift(:, i_story+1)) ) )']) % Obtain max PSDR at i_story
                    elseif (system == 3) || (system == 4) || (system == 5) || (system == 6) || (system == 7) % i.e., isolated building
                        eval(['MaxDriftVector_Story_',num2str(i_story),'(i_gm,i_RP)', '= max( max( abs(data_Drift(:, i_story+2)) ) )']) % Obtain max PSDR at i_story
                    end
                end
                
                for i_floor = 1 : length(ImpRows_AccelAll)
                    % % For floor (story) based peak response - - - - - (floor & roof acceleration)
                    if     (system == 1) || (system == 2) % i.e., non-isolated buildings
                        eval(['MaxAccelVector_Floor_',num2str(i_floor),'(i_gm,i_RP)', '= max( max( abs(data_Accel(:, i_floor+1)) ) )/g']) % Obtain max PGA, PFA, PRA at i_floor
                    elseif (system == 3) || (system == 4) || (system == 5) || (system == 6) || (system == 7) % i.e., isolated building
                        eval(['MaxAccelVector_Floor_',num2str(i_floor),'(i_gm,i_RP)', '= max( max( abs(data_Accel(:, i_floor+2)) ) )/g']) % Obtain max PGA, PFA, PRA at i_floor
                    end
                end
                
        end
        
    % Count number of GMs that causes collapse (out of 40) in current return period
    num_collapse(1, i_RP) = sum(n_Collapse_matrix(:,i_RP) > 0);  % Count and save in "num_collapse" (used later for maximum likelifood method)
            
end

%%
% Estimate params of collapse fragility function using MLE method (equation 11 in Baker, 2015, Earthquake Spectra) --------------------------
[theta_hat_mle, beta_hat_mle] = fn_mle_pc(Sa_T1, num_gms, num_collapse);
% Compute collapse fragility curve (function) PC|IM using estimated parameters --------------------------
Sa_Range = 0.001:0.001:50;   % IM levels to plot fragility function
PC_Sa    = normcdf((log(Sa_Range/theta_hat_mle))/beta_hat_mle); % compute fragility function using Eq. 1 and estimated parameters
PC_Sa_write = [Sa_Range', PC_Sa'];  % Write out the output in text file
dlmwrite(strcat('PC_Sa.txt'), PC_Sa_write, 'delimiter', '\t', 'precision', 8); % Output the data of fragility curves

PC_IM = zeros(NumIM,2);

for i = 1:NumIM
    PC_IM(i,1)=i; 
    PC_IM(i,2)=normcdf((log(Sa_T1(i)/theta_hat_mle))/beta_hat_mle); 
    
        if num_collapse(1, i) == NumGM % Added 27June2019, to consider Probability of collapse of 100% (all GM caused collapse)
           PC_IM(i,2) = 1.;
        end
    
end

% Compute probability that the building is being considered to be demolished --------------------------
% First, compute probability density function of maximum residual story drift ratio, fRSDR|IM, along the height --------------------------
x_RSDR_pdf = [0.0002:0.0002:0.06]';
    % Assumed fragility curve for decision of demolition of building based
    % on Ramirez and Miranda (2012)...
    theta_hat_Demolish = 0.015; beta_hat_Demolish = 0.3; % From Ramirez and Miranda (2012)
    PD_RSDR = normcdf((log(x_RSDR_pdf/theta_hat_Demolish))/beta_hat_Demolish); % compute fragility function using Eq. 1 and estimated parameters
    vectorPD_RSDR = [x_RSDR_pdf, PD_RSDR];
    
    PD_IM_NC = zeros(NumIM,2);
    vector_mu_sigma_IM = zeros(NumIM,3);
    
for i_RP = 1 : NumIM
    
    MaxResDriftVector_NoCollapse      = MaxResDriftVector( n_Collapse_matrix(:,i_RP) < 1, i_RP ); % Remove any cases that collapse occurs
    MaxResDriftVector_NoCollapse_Log  = log(MaxResDriftVector_NoCollapse);
    mu_RSDR_IM                        = mean(MaxResDriftVector_NoCollapse);     % Mean of RSDR
    sigma_lnRSDR_IM                   = std(MaxResDriftVector_NoCollapse_Log);  % Standard deviation of lnPSDR
                
            if isscalar(MaxResDriftVector_NoCollapse) == 1
               sigma_lnRSDR_IM = 0.009999;
            end

            if isempty(MaxResDriftVector_NoCollapse) == 1
               mu_RSDR_IM      = y_drift_Col; % 9.999; % Corrected: 02.Oct.2019
               sigma_lnRSDR_IM = 0.009999;    % 0.09999; Corrected: 02.Oct.2019
            end
            
    vector_mu_sigma_IM(i_RP,1)                 = i_RP;
    vector_mu_sigma_IM(i_RP,2)                 = mu_RSDR_IM;
    vector_mu_sigma_IM(i_RP,3)                 = sigma_lnRSDR_IM;
    
    vectorPDF_RSDR_IM_col1                     = x_RSDR_pdf;                                          % column=1 of vector of PDF of RSDR at IM=im=i_RP
    vectorPDF_RSDR_IM_col2                     = lognpdf(x_RSDR_pdf,log(mu_RSDR_IM),sigma_lnRSDR_IM); % column=2 of vector of PDF of RSDR at IM=im=i_RP
    vectorPDF_RSDR_IM                          = [vectorPDF_RSDR_IM_col1, vectorPDF_RSDR_IM_col2];
    eval(['vectorPDF_RSDR_IM_',num2str(i_RP), '= [vectorPDF_RSDR_IM_col1, vectorPDF_RSDR_IM_col2]']);
    
    PD_RSDR_x_fRSDR_IM_col1 = x_RSDR_pdf;
    PD_RSDR_x_fRSDR_IM_col2 = vectorPD_RSDR(:,2) .* vectorPDF_RSDR_IM(:,2);
    
    PD_IM_NC(i_RP,1) = i_RP;
    PD_IM_NC(i_RP,2) = trapz(PD_RSDR_x_fRSDR_IM_col1, PD_RSDR_x_fRSDR_IM_col2); % Probability that building is being considered to be demolished
    
end


% Compute mean values of loss for each source of loss ------------------------------------------------------------------------------------------
% Source of losses = "Structural repair time";
%                    "Non-structural repair time (drift)";
%                    "Non-structural repair time (acc)"
%                    "Demolish time";
%                    "Collapse time"


% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
%% COMPUTE: "Structural repair time"  - - - - - - - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %

x_PSDR_pdf = [0.00002:0.00002:3.000]'; % [0.0002:0.0002:3.000]';  0.0002 -> 0.00002 (modified 25Sep2019)

for i_RP = 1 : NumIM
    
%          sigma_m_sigma_n_integral_Eq4_Struct                                      = zeros(NumStory,2);    % Reset the vector for each i_RP
    eval(['sigma_m_sigma_n_integral_Eq4_'                ,num2str(i_RP),'_Struct', '= zeros(NumStory,2)']); % modified 26/Sep/2019, considering "i_RP"                                            = zeros(NumStory,2);    % Reset the vector for each i_RP
    eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_Struct', '= zeros(NumStory,2)']); % Reset the vector for each i_RP
    
    for i_story = 1 : NumStory   % consider each story
        
        [m] = info_num_Components_Structural(i_story, system);  % from outside of main MATLAB file.. from i_story, obtain total # of damageable component IDs
        
        MaxPSDR_Vector_NoCollapse = eval(['MaxDriftVector_Story_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']); % Added 25June2019
        mean_PSDR     = eval(['mean(MaxPSDR_Vector_NoCollapse','(:,1))']); % Added 25June2019
%       mean_PSDR     = eval(['mean(MaxDriftVector_Story_',num2str(i_story),'(:,i_RP))']);
        std_log_PSDR  = eval(['std(log(MaxPSDR_Vector_NoCollapse','(:,1)))']);
%       std_log_PSDR  = eval(['std(log(MaxDriftVector_Story_',num2str(i_story),'(:,i_RP)))']);
            
            if isscalar(MaxPSDR_Vector_NoCollapse) == 1
               std_log_PSDR = 0.009999;
            end
        
            if isempty(MaxPSDR_Vector_NoCollapse) == 1
               mean_PSDR    = y_drift_Col; % 9.999; Corrected: 02.Oct.2019
               std_log_PSDR = 0.009999; % 0.09999; Corrected: 02.Oct.2019
            end
        
        PDF_PSDR_IM   = lognpdf(x_PSDR_pdf,log(mean_PSDR),std_log_PSDR);
        
        sigma_n_integral_Eq4                 = zeros(length(m),2);      % Reset the vector for each i_story
        sigma_n_integral_Eq4_ProbExceed_DSij = zeros(length(m),2);      % Reset the vector for each i_story
        
      for i_m = m                    % consider number of types(?) of damageable components at a story.. note: 1 to m corresponds to component ID; "Sigma na(m)" in MR(2007)
          
          [n] = info_num_DamageStates_Structural(i_m);  % then, from damageable each component ID, obtain number of damage states
          
          integral_Eq4                 = zeros(n+1,2); % Reset the vector for each i_m
          integral_Eq4_ProbExceed_DSij = zeros(n+1,2); % Reset the vector for each i_m
          
         for i_n = 0 : n             % consider number of damage states a component may experience
             
             [PDS_ij_EDP, xm_RepTime, numCompPerStory] = info_Comp_Fragility_Structural_DownTime(i_n, i_m, x_PSDR_pdf, system, i_story); % then, from damageable each component ID, obtain number of damage states
             
             integral_Eq4(i_n+1,1)                   = i_n+1;
             integral_Eq4(i_n+1,2)                   = trapz(x_PSDR_pdf, numCompPerStory*xm_RepTime*(PDS_ij_EDP.*PDF_PSDR_IM))    /(wh*wr*cn); % "numCompPerStory"="Nu_i(m)" in MR(2007)
             integral_Eq4_ProbExceed_DSij(i_n+1,1)   = i_n+1;
             integral_Eq4_ProbExceed_DSij(i_n+1,2)   = trapz(x_PSDR_pdf, (1.-PDS_ij_EDP).*PDF_PSDR_IM) * 100.0; % Prob. exceeding DSij for check of change-of-trade
                                                       % "PDS_ij_EDP" -> "(1.-PDS_ij_EDP)" modified on 26/Sep/2019
         end
         
         sigma_n_integral_Eq4(i_m,1)                 = i_m;
         sigma_n_integral_Eq4(i_m,2)                 = sum(integral_Eq4(:,2));  % Sum of all damageable damage states of a specific component.
         sigma_n_integral_Eq4_ProbExceed_DSij(i_m,1) = i_m;
         sigma_n_integral_Eq4_ProbExceed_DSij(i_m,2) = integral_Eq4_ProbExceed_DSij(1,2); % Take prob of having damage DSi1 for check of change-of-trade.
         
      end
      
%            sigma_m_sigma_n_integral_Eq4_Struct(i_story,1) = i_story;
%            sigma_m_sigma_n_integral_Eq4_Struct(i_story,2) = sum(sigma_n_integral_Eq4(:,2));  % Sum of all damageable components on a specific floor.
      eval(['sigma_m_sigma_n_integral_Eq4_',                num2str(i_RP),'_Struct','(i_story,1)', '= i_story']); % modified 26/Sep/2019
      eval(['sigma_m_sigma_n_integral_Eq4_',                num2str(i_RP),'_Struct','(i_story,2)', '= sum(sigma_n_integral_Eq4(:,2))']);  % Sum of all damageable components on a specific floor. % modified 26/Sep/2019
      eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_Struct','(i_story,1)', '= i_story']);
      eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_Struct','(i_story,2)', '= max(sigma_n_integral_Eq4_ProbExceed_DSij(:,2))>5.0']);  % If max of probability of all damageable components on a specific floor exceeds 5%; for check of change-of-trade..
      
    end
       
       % Added 28.Sep.2019 -----
       dlmwrite(strcat('sigma_m_sigma_n_integral_Eq4_',                 num2str(i_RP),'_Struct','.txt'), eval(['sigma_m_sigma_n_integral_Eq4_',                 num2str(i_RP), '_Struct']), 'delimiter', '\t', 'precision', 8);
       dlmwrite(strcat('sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_', num2str(i_RP),'_Struct','.txt'), eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_', num2str(i_RP), '_Struct']), 'delimiter', '\t', 'precision', 8);

end



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
%% COMPUTE: "Non-structural repair loss (Drift)"  - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %

x_PSDR_pdf = [0.00002:0.00002:3.000]'; % [0.0002:0.0002:3.000]';  0.0002 -> 0.00002 (modified 25Sep2019)

for i_RP = 1 : NumIM
    
%          sigma_m_sigma_n_integral_Eq4_NonStructDrift                                             = zeros(NumStory,2);    % Reset the vector for each i_RP
    eval(['sigma_m_sigma_n_integral_Eq4_'                ,num2str(i_RP),'_NonStructDrift', '= zeros(NumStory,2)']); % modified 26/Sep/2019, considering "i_RP"                                            = zeros(NumStory,2);    % Reset the vector for each i_RP
    eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructDrift', '= zeros(NumStory,2)']); % Reset the vector for each i_RP
    
    for i_story = 1 : NumStory   % consider each story
        
        [m] = info_num_Components_NonStructural_Drift(i_story);  % from outside of main MATLAB file.. from i_story, obtain total # of damageable component IDs
        
        MaxPSDR_Vector_NoCollapse = eval(['MaxDriftVector_Story_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']); % Added 25June2019
        mean_PSDR     = eval(['mean(MaxPSDR_Vector_NoCollapse','(:,1))']); % Added 25June2019
%       mean_PSDR     = eval(['mean(MaxDriftVector_Story_',num2str(i_story),'(:,i_RP))']);
        std_log_PSDR  = eval(['std(log(MaxPSDR_Vector_NoCollapse','(:,1)))']);
%       std_log_PSDR  = eval(['std(log(MaxDriftVector_Story_',num2str(i_story),'(:,i_RP)))']);
            
            if isscalar(MaxPSDR_Vector_NoCollapse) == 1
               std_log_PSDR = 0.009999;
            end

            if isempty(MaxPSDR_Vector_NoCollapse) == 1
               mean_PSDR    = y_drift_Col; % 9.999; % Corrected: 02.Oct.2019
               std_log_PSDR = 0.009999;    % 0.09999; Corrected: 02.Oct.2019
            end

        PDF_PSDR_IM   = lognpdf(x_PSDR_pdf,log(mean_PSDR),std_log_PSDR);
        
        sigma_n_integral_Eq4                 = zeros(length(m),2);  % Reset the vector for each i_story
        sigma_n_integral_Eq4_ProbExceed_DSij = zeros(length(m),2);  % Reset the vector for each i_story
        
      for i_m = m          % consider number of types(?) of damageable components at a story.. note: 1 to m corresponds to component ID
        
          [n] = info_num_DamageStates_NonStructural_Drift(i_m);  % then, from damageable each component ID, obtain number of damage states
              
          integral_Eq4                 = zeros(n+1,2); % Reset the vector for each i_m
          integral_Eq4_ProbExceed_DSij = zeros(n+1,2); % Reset the vector for each i_m
          
         for i_n = 0 : n   % consider number of damage states a component may experience
             
             [PDS_ij_EDP, xm_RepTime, numCompPerStory] = info_Comp_Fragility_NonStructural_Drift_DownTime(i_n, i_m, x_PSDR_pdf); % then, from damageable each component ID, obtain number of damage states
             
             integral_Eq4(i_n+1,1)                   = i_n;
             integral_Eq4(i_n+1,2)                   = trapz(x_PSDR_pdf, numCompPerStory*xm_RepTime*(PDS_ij_EDP.*PDF_PSDR_IM))    /(wh*wr*cn); % "numCompPerStory"="Nu_i(m)" in MR(2007);
             integral_Eq4_ProbExceed_DSij(i_n+1,1)   = i_n+1;
             integral_Eq4_ProbExceed_DSij(i_n+1,2)   = trapz(x_PSDR_pdf, (1.-PDS_ij_EDP).*PDF_PSDR_IM) * 100.0; % Prob. exceeding DSij for check of change-of-trade
                                                       % "PDS_ij_EDP" -> "(1.-PDS_ij_EDP)" modified on 26/Sep/2019
                                                       
         end
         
         sigma_n_integral_Eq4(i_m,1)                 = i_m;
         sigma_n_integral_Eq4(i_m,2)                 = sum(integral_Eq4(:,2));
         sigma_n_integral_Eq4_ProbExceed_DSij(i_m,1) = i_m;
         sigma_n_integral_Eq4_ProbExceed_DSij(i_m,2) = integral_Eq4_ProbExceed_DSij(1,2); % Take prob of having damage DSi1 for check of change-of-trade.
          
      end
      
%            sigma_m_sigma_n_integral_Eq4_NonStructDrift(i_story,1)                                                = i_story;
%            sigma_m_sigma_n_integral_Eq4_NonStructDrift(i_story,2)                                                = sum(sigma_n_integral_Eq4(:,2));
      eval(['sigma_m_sigma_n_integral_Eq4_',                num2str(i_RP),'_NonStructDrift','(i_story,1)', '= i_story']); % modified 26/Sep/2019
      eval(['sigma_m_sigma_n_integral_Eq4_',                num2str(i_RP),'_NonStructDrift','(i_story,2)', '= sum(sigma_n_integral_Eq4(:,2))']);  % Sum of all damageable components on a specific floor. % modified 26/Sep/2019
      eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructDrift','(i_story,1)', '= i_story']);
      eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructDrift','(i_story,2)', '= sum(sigma_n_integral_Eq4_ProbExceed_DSij(:,2)>5.0)']);  % Count number of damageable components whose probability > 5%; for check of change-of-trade..
       
    end
       
       % Added 28.Sep.2019 -----
       dlmwrite(strcat('sigma_m_sigma_n_integral_Eq4_',                 num2str(i_RP),'_NonStructDrift','.txt'), eval(['sigma_m_sigma_n_integral_Eq4_',                 num2str(i_RP), '_NonStructDrift']), 'delimiter', '\t', 'precision', 8);
       dlmwrite(strcat('sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_', num2str(i_RP),'_NonStructDrift','.txt'), eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_', num2str(i_RP), '_NonStructDrift']), 'delimiter', '\t', 'precision', 8);

end



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
%% COMPUTE: "Non-structural repair loss (Accel)" - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %

x_Accel_pdf = [0.0002:0.0002:5.000]';

for i_RP = 1 : NumIM
    
%          sigma_m_sigma_n_integral_Eq4_NonStructAccel                                      = zeros(length(ImpRows_AccelAll),2); % Reset the vector for each i_RP
    eval(['sigma_m_sigma_n_integral_Eq4_'                ,num2str(i_RP),'_NonStructAccel', '= zeros(NumFloor,2)']); % modified 26/Sep/2019, considering "i_RP"                                            = zeros(NumStory,2);    % Reset the vector for each i_RP
    eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructAccel', '= zeros(NumFloor,2)']); % Reset the vector for each i_RP
    
    for i_floor = 1 : NumFloor   % consider each floor
        
        [m] = info_num_Components_NonStructural_Accel(i_floor);  % from outside of main MATLAB file.. from i_floor, obtain total # of damageable component IDs
        
        MaxAccel_Vector_NoCollapse = eval(['MaxAccelVector_Floor_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']); % Added 25June2019
        mean_Accel     = eval(['mean(MaxAccel_Vector_NoCollapse','(:,1))']); % Added 25June2019
%       mean_Accel     = eval(['mean(MaxAccelVector_Floor_',num2str(i_floor),'(:,i_RP))']);
        std_log_Accel  = eval(['std(log(MaxAccel_Vector_NoCollapse','(:,1)))']);
%       std_log_Accel  = eval(['std(log(MaxAccelVector_Floor_',num2str(i_floor),'(:,i_RP)))']);

            if isscalar(MaxAccel_Vector_NoCollapse) == 1
               std_log_Accel = 0.009999;    % 0.09999; Corrected: 02.Oct.2019
            end

            if isempty(MaxAccel_Vector_NoCollapse) == 1
               mean_Accel    = x_Accel_pdf(end); % 9.999; % Corrected: 02.Oct.2019
               std_log_Accel = 0.009999;    % 0.09999; Corrected: 02.Oct.2019
            end
            
        PDF_Accel_IM   = lognpdf(x_Accel_pdf,log(mean_Accel),std_log_Accel);
        
        sigma_n_integral_Eq4                 = zeros(length(m),2);  % Reset the vector for each i_floor
        sigma_n_integral_Eq4_ProbExceed_DSij = zeros(length(m),2);  % Reset the vector for each i_story
        
      for i_m = m            % consider number of types(?) of damageable components at a floor.. note: 1 to m corresponds to component ID
        
          [n] = info_num_DamageStates_NonStructural_Accel(i_m);  % then, from damageable each component ID, obtain number of damage states
              
          integral_Eq4                 = zeros(n+1,2); % Reset the vector for each i_m
          integral_Eq4_ProbExceed_DSij = zeros(n+1,2); % Reset the vector for each i_m

         for i_n = 0 : n         % consider number of damage states a component may experience
             
             [PDS_ij_EDP, xm_RepTime, numCompPerStory] = info_Comp_Fragility_NonStructural_Accel_DownTime(i_n, i_m, x_Accel_pdf); % then, from damageable each component ID, obtain number of damage states
             
             integral_Eq4(i_n+1,1)                   = i_n+1;
             integral_Eq4(i_n+1,2)                   = trapz(x_Accel_pdf, numCompPerStory*xm_RepTime*(PDS_ij_EDP.*PDF_Accel_IM));
             integral_Eq4_ProbExceed_DSij(i_n+1,1)   = i_n+1;
             integral_Eq4_ProbExceed_DSij(i_n+1,2)   = trapz(x_Accel_pdf, (1.-PDS_ij_EDP).*PDF_Accel_IM) * 100.0; % Prob. exceeding DSi1 for check of change-of-trade
                                                       % "PDS_ij_EDP" -> "(1.-PDS_ij_EDP)" modified on 26/Sep/2019
         end
         
         sigma_n_integral_Eq4(i_m,1)                 = i_m;
         sigma_n_integral_Eq4(i_m,2)                 = sum(integral_Eq4(:,2));
         sigma_n_integral_Eq4_ProbExceed_DSij(i_m,1) = i_m;
         sigma_n_integral_Eq4_ProbExceed_DSij(i_m,2) = integral_Eq4_ProbExceed_DSij(1,2); % Take prob of having damage DSi1 for check of change-of-trade.
         
      end
      
%            sigma_m_sigma_n_integral_Eq4_NonStructAccel(i_floor,1)                                         = i_floor;
%            sigma_m_sigma_n_integral_Eq4_NonStructAccel(i_floor,2)                                         = sum(sigma_n_integral_Eq4(:,2));
      eval(['sigma_m_sigma_n_integral_Eq4_',                num2str(i_RP),'_NonStructAccel','(i_floor,1)', '= i_floor']); % modified 26/Sep/2019
      eval(['sigma_m_sigma_n_integral_Eq4_',                num2str(i_RP),'_NonStructAccel','(i_floor,2)', '= sum(sigma_n_integral_Eq4(:,2))']);  % Sum of all damageable components on a specific floor. % modified 26/Sep/2019
      eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructAccel','(i_floor,1)', '= i_floor']);
      eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructAccel','(i_floor,2)', '= sum(sigma_n_integral_Eq4_ProbExceed_DSij(:,2)>5.0)']);  % Count number of damageable components whose probability > 5%; for check of change-of-trade..
      
    end
       
       % Added 28.Sep.2019 -----
       dlmwrite(strcat('sigma_m_sigma_n_integral_Eq4_',                 num2str(i_RP),'_NonStructAccel','.txt'), eval(['sigma_m_sigma_n_integral_Eq4_',                 num2str(i_RP), '_NonStructAccel']), 'delimiter', '\t', 'precision', 8);
       dlmwrite(strcat('sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_', num2str(i_RP),'_NonStructAccel','.txt'), eval(['sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_', num2str(i_RP), '_NonStructAccel']), 'delimiter', '\t', 'precision', 8);

end



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % %
%% COMPUTE: Change-of-Trade Time per Operational Unit (i.e., Story) - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % %
% Calculate: N_T(m) * E[R_COT]

for i_RP = 1 : NumIM
    
    eval(['NT_m_',num2str(i_RP), '= zeros(NumStory,2)']);  % NT(m) is the number of changes of trade in operational unit m
    
%   for i_floor = 1 : NumFloor   % consider each floor
    for i_story = 1 : NumStory   % consider each story
        
        eval(['NT_m_',num2str(i_RP),'(i_story,1)    = i_story']);
        eval(['NT_m_',num2str(i_RP),'(i_story,2)', '= sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_Struct',        '(i_story  ,2)',...
                                                   '+ sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructDrift','(i_story  ,2)',...
                                                   '+ sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructAccel','(i_story+1,2)']); % "+1" is needed to be consistent with operational unit (=story)
                                                   % Number of changes of trade in operational unit m; Number of components exceeding DSi1>1 percent minus 1 to get NT(m)
       
    end
       
        eval(['NT_m_',num2str(i_RP),'(1,      2)', '= NT_m_',num2str(i_RP),'(1,2)',...
                                                   '+ sigma_m_sigma_n_integral_Eq4_ProbExceed_DSij_',num2str(i_RP),'_NonStructAccel','(1,2)']); % Number of components exceeding DSi1>1 percent minus 1 to get NT(m)
            
            for i_story = 1 : NumStory   % consider each story
               if eval(['NT_m_',num2str(i_RP),'(i_story,2)', '> 0.999999999']) 
                  
                  eval(['NT_m_',num2str(i_RP),'(i_story,2)', '= NT_m_',num2str(i_RP),'(i_story,2)',' - 1.0']);  % Num of change-of-trade is num damageable component minus 1
                  
               end
               
               % Compute Change-of-Trade Time per Operational Unit (i.e., Story): N_T(m) * E[R_COT]
                 eval(['Time_COT_', num2str(i_RP), '_Fast', '(i_story,1) = i_story']);                                              % in case of fast-track repair scheme
                 eval(['Time_COT_', num2str(i_RP), '_Fast', '(i_story,2) = NT_m_', num2str(i_RP), '(i_story,2)', '* R_COT_Fast']);  % in case of fast-track repair scheme
                 eval(['Time_COT_', num2str(i_RP), '_Slow', '(i_story,1) = i_story']);                                              % in case of slow-track repair scheme
                 eval(['Time_COT_', num2str(i_RP), '_Slow', '(i_story,2) = NT_m_', num2str(i_RP), '(i_story,2)', '* R_COT_Slow']);  % in case of slow-track repair scheme
               
            end
            
       % Added 28.Sep.2019 -----
       dlmwrite(strcat('Time_COT_', num2str(i_RP), '_Fast','.txt'), eval(['Time_COT_', num2str(i_RP), '_Fast']), 'delimiter', '\t', 'precision', 8);
       dlmwrite(strcat('Time_COT_', num2str(i_RP), '_Slow','.txt'), eval(['Time_COT_', num2str(i_RP), '_Slow']), 'delimiter', '\t', 'precision', 8);

end



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % % % %
%% COMPUTE: Repair time per Operational Unit (i.e., Story) with Time_COT - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % % % %
% Calculate: E[R*_U(m)|NC,im]
% Calculate: (?_(i=1)^na(m) { Nu _i (m)*?_(j=1)^(nds _i) E[R_i DM _ij]*P[D _ij NC, im]})/(wh*wr*cn) + N_T (m)*E[R_COT ]


for i_RP = 1 : NumIM
    
    eval(['TimeStory_RP_',num2str(i_RP), '_Fast', '= zeros(NumStory,2)']);  % NT(m) is the number of changes of trade in operational unit m; Fast repair scheme
    eval(['TimeStory_RP_',num2str(i_RP), '_Slow', '= zeros(NumStory,2)']);  % NT(m) is the number of changes of trade in operational unit m; Slow repair scheme
    
    for i_story = 1 : NumStory   % consider each story
        
        eval(['TimeStory_RP_',num2str(i_RP), '_Fast', '(i_story, 1)', '= i_story']);
%       eval(['TimeStory_RP_',num2str(i_RP), '_Fast', '(i_story, 2)', ' = sigma_m_sigma_n_integral_Eq4_Struct(i_story,2)',...
%                                                                     ' + sigma_m_sigma_n_integral_Eq4_NonStructDrift(i_story  ,2)',...
%                                                                     ' + sigma_m_sigma_n_integral_Eq4_NonStructAccel(i_story+1,2)',...; % "+1" to be consistent with operational unit definition (=story not floor)
%                                                                     ' + Time_COT_', num2str(i_RP), '_Fast', '(i_story,2)']);  % Fast-track repair scheme
        eval(['TimeStory_RP_',num2str(i_RP), '_Fast', '(i_story, 2)', ' = sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_Struct',        '(i_story,2)',...
                                                                      ' + sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_NonStructDrift','(i_story,2)',...
                                                                      ' + sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_NonStructAccel','(i_story,2)',...; % "+1" to be consistent with operational unit definition (=story not floor)
                                                                      ' + Time_COT_', num2str(i_RP), '_Fast', '(i_story,2)']);  % Fast-track repair scheme  % modified on 26Sep2019
        
        
        eval(['TimeStory_RP_',num2str(i_RP), '_Slow', '(i_story, 1)', ' = i_story']);
%       eval(['TimeStory_RP_',num2str(i_RP), '_Slow', '(i_story, 2)', ' = sigma_m_sigma_n_integral_Eq4_Struct(i_story,2)',...
%                                                                     ' + sigma_m_sigma_n_integral_Eq4_NonStructDrift(i_story,2)',...
%                                                                     ' + sigma_m_sigma_n_integral_Eq4_NonStructAccel(i_story+1,2)',...; % "+1" to be consistent with operational unit definition (=story not floor)
%                                                                     ' + Time_COT_', num2str(i_RP), '_Slow', '(i_story,2)']);  % Slow-track repair scheme
        eval(['TimeStory_RP_',num2str(i_RP), '_Slow', '(i_story, 2)', ' = sigma_m_sigma_n_integral_Eq4_',  num2str(i_RP),'_Struct',       '(i_story,2)',...
                                                                      ' + sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_NonStructDrift','(i_story,2)',...
                                                                      ' + sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_NonStructAccel','(i_story,2)',...; % "+1" to be consistent with operational unit definition (=story not floor)
                                                                      ' + Time_COT_', num2str(i_RP), '_Slow', '(i_story,2)']);  % Slow-track repair scheme  % modified on 26Sep2019

    end
    
            % Add nonstructural (accel-sensitive) at first floor to first story to be consistent with the definition of operational unit (=story not floor; elevator only, PGA) 

            eval(['TimeStory_RP_',num2str(i_RP), '_Fast', '(1, 2)',   ' = TimeStory_RP_',num2str(i_RP), '_Fast', '(1, 2)',...
                                                                      ' + sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_NonStructAccel','(1, 2)']); %  ' + sigma_m_sigma_n_integral_Eq4_NonStructAccel(1,2)']);
            
            eval(['TimeStory_RP_',num2str(i_RP), '_Slow', '(1, 2)'    ' = TimeStory_RP_',num2str(i_RP), '_Slow', '(1, 2)',...
                                                                      ' + sigma_m_sigma_n_integral_Eq4_', num2str(i_RP),'_NonStructAccel','(1, 2)']); %  ' + sigma_m_sigma_n_integral_Eq4_NonStructAccel(1,2)']);
     
     % Added 28.Sep.2019 -----
     dlmwrite(strcat('TimeStory_RP_',num2str(i_RP), '_Fast','.txt'), eval(['TimeStory_RP_',num2str(i_RP), '_Fast']), 'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('TimeStory_RP_',num2str(i_RP), '_Slow','.txt'), eval(['TimeStory_RP_',num2str(i_RP), '_Slow']), 'delimiter', '\t', 'precision', 8);

end



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % %
%% COMPUTE: Mean total repair time (rational component of repair time)  - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % %
% Calculate: E[R_T|NC,im] (fast or slow repair schemes)
    
Time_RP_Fast = zeros(NumIM,2);  % NT(m) is the number of changes of trade in operational unit m; Fast repair scheme
Time_RP_Slow = zeros(NumIM,2);  % NT(m) is the number of changes of trade in operational unit m; Slow repair scheme

for i_RP = 1 : NumIM

    Time_RP_Fast(i_RP, 1) = i_RP; 
    Time_RP_Fast(i_RP, 2) = eval(['max(','TimeStory_RP_',num2str(i_RP), '_Fast', '(:, 2)',')']);  % Pick max from all operational units (=stories); Fast repair scheme
    
    Time_RP_Slow(i_RP, 1) = i_RP; 
    Time_RP_Slow(i_RP, 2) = eval(['sum(','TimeStory_RP_',num2str(i_RP), '_Slow', '(:, 2)',')']);  % Pick max from all operational units (=stories); Slow repair scheme
         
     % Added 28.Sep.2019 -----
     dlmwrite(strcat('Time_RP_Fast','.txt'), Time_RP_Fast, 'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('Time_RP_Slow','.txt'), Time_RP_Slow, 'delimiter', '\t', 'precision', 8);

end



% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % %
%% COMPUTE: Mean mobilization time (irrational component of repair time)  - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % % % % % % % % % % % %
% Calculate: R_T0_NC_im = E[R_T0|NC,im] (based on post-earthquake tagging actions)
   
% Rapid evaluation - - - - -
% P[TAG=G|im,RE]=P(Light external structural damage|im,NC)
% P[TAG=Y|im,RE]=P(Moderate external structural damage|im,NC)
% P[TAG=R|im,RE]=P(Severe external structural damage|im,NC)

% Detailed evaluation - - - - -
% P[TAG=G|im,DE]=P(Non-severe internal structural damage|im,NC)
% P[TAG=R|im,DE]=P(Severe internal structural damage|im,NC)

P_TAG_G_NC_im = zeros(NumIM,2);  P_TAG_Y_NC_im = zeros(NumIM,2);  P_TAG_R_NC_im = zeros(NumIM,2); % Modified on 29.Sep.2019
P_TAG_RE_G    = zeros(NumIM,2);  P_TAG_RE_Y    = zeros(NumIM,2);  P_TAG_RE_R    = zeros(NumIM,2);
P_TAG_DE_G    = zeros(NumIM,2);                                   P_TAG_DE_R    = zeros(NumIM,2);

i_n_RE_G = 1; i_n_RE_Y = 2; i_n_RE_R = 3;  % Damage state to pick for fragility curves for Green, Yellow and Red tags

i_n_DE_R = 3;  % Damage state to pick for fragility curves for Green, Yellow and Red tags

if system == 1  || system == 3  || system == 4 % if SCBF (isolated OR non-isolated)
    i_m_RE = 4;  % Round HSS (for SCBFs)
else
    i_m_RE = 5;  % Moment connection; one-sided; <= W27 (for SMRFs)
end

i_m_DE = 8; % Corrugated slab


for i_RP = 1 : NumIM
    
    P_TAG_RE_G_Temp = zeros(NumStory,2);   P_TAG_RE_Y_Temp = zeros(NumStory,2);   P_TAG_RE_R_Temp = zeros(NumStory,2);
    P_TAG_DE_G_Temp = zeros(NumStory,2);   P_TAG_DE_R_Temp = zeros(NumStory,2);
    
    for i_story = 1 : NumStory   % consider each story
        
        % Take PDF of peak story drift ratio at each story of each return period - - - - -
        MaxPSDR_Vector_NoCollapse = eval(['MaxDriftVector_Story_',num2str(i_story),'( n_Collapse_matrix(:,i_RP) < 1, i_RP )']);
        mean_PSDR                 = eval(['mean(MaxPSDR_Vector_NoCollapse','(:,1))']);
        std_log_PSDR              = eval(['std(log(MaxPSDR_Vector_NoCollapse','(:,1)))']);
            
            if isscalar(MaxPSDR_Vector_NoCollapse) == 1
               std_log_PSDR = 0.009999;
            end
            
            if isempty(MaxPSDR_Vector_NoCollapse) == 1
               mean_PSDR    = y_drift_Col; % 9.999;   Corrected: 02.Oct.2019
               std_log_PSDR = 0.009999;    % 0.09999; Corrected: 02.Oct.2019
            end
            
        PDF_PSDR_IM = lognpdf(x_PSDR_pdf,log(mean_PSDR),std_log_PSDR);
        
        P_TAG_RE_G_Temp(i_story, 1) = i_story;  P_TAG_RE_Y_Temp(i_story, 1) = i_story;  P_TAG_RE_R_Temp(i_story, 1) = i_story;
        
        % % Rapid Evaluations - - - - - - - - - -
%         % Green Tagging (Rapid evaluation) - - -
%         [PDS_ij_EDP]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_RE_G, i_m_RE, x_PSDR_pdf, system, i_story);
%         P_TAG_RE_G_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP.*PDF_PSDR_IM));   % Integrate to obtain probability of having Green tag 
%         % Yellow Tagging (Rapid evaluation) - - -
%         [PDS_ij_EDP]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_RE_Y, i_m_RE, x_PSDR_pdf, system, i_story);
%         P_TAG_RE_Y_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP.*PDF_PSDR_IM));   % Integrate to obtain probability of having Yellow
%         % Red Tagging (Rapid evaluation) - - -
%         [PDS_ij_EDP]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_RE_R, i_m_RE, x_PSDR_pdf, system, i_story);
%         P_TAG_RE_R_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP.*PDF_PSDR_IM));   % Integrate to obtain probability of having Red tag 
        
        % Modified on 29.Sep.2019
        % Red Tagging (Rapid evaluation) - - -
        [PDS_ij_EDP_RE_R]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_RE_R, i_m_RE, x_PSDR_pdf, system, i_story);
        P_TAG_RE_R_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP_RE_R.*PDF_PSDR_IM));   % Integrate to obtain probability of having Red tag
        [c, index_RE_R]       = max(abs(P_TAG_RE_R_Temp(:,2)));
        % Yellow Tagging (Rapid evaluation) - - -
        [PDS_ij_EDP_RE_Y]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_RE_Y, i_m_RE, x_PSDR_pdf, system, i_story);
        P_TAG_RE_Y_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP_RE_Y.*PDF_PSDR_IM));   % Integrate to obtain probability of having Yellow
        % Green Tagging (Rapid evaluation) - - -
        [PDS_ij_EDP_RE_G]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_RE_G, i_m_RE, x_PSDR_pdf, system, i_story);
        P_TAG_RE_G_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP_RE_G.*PDF_PSDR_IM));   % Integrate to obtain probability of having Green tag 
        
        % % Detailed Evaluations - - - - - - - - - -
        % Red Tagging (Detailed evaluation) - - -
        [PDS_ij_EDP]     = info_Comp_Fragility_Structural_DownTime_Tagging(i_n_DE_R, i_m_DE, x_PSDR_pdf, system, i_story);
        P_TAG_DE_R_Temp(i_story, 2) = trapz(x_PSDR_pdf, (PDS_ij_EDP.*PDF_PSDR_IM));   % Integrate to obtain probability of having Red tag 
        
    end
    
    P_TAG_G_NC_im(i_RP,1) = RP(i_RP);   P_TAG_Y_NC_im(i_RP,1) = RP(i_RP);   P_TAG_R_NC_im(i_RP,1) = RP(i_RP); 
    
    P_TAG_RE_G(i_RP,1)    = RP(i_RP);   P_TAG_RE_G(i_RP,2)    = P_TAG_RE_G_Temp(index_RE_R,2);   % = max(P_TAG_RE_G_Temp(:,2)); % Modified on 29.Sep.2019
    P_TAG_RE_Y(i_RP,1)    = RP(i_RP);   P_TAG_RE_Y(i_RP,2)    = P_TAG_RE_Y_Temp(index_RE_R,2);   % = max(P_TAG_RE_Y_Temp(:,2)); % Modified on 29.Sep.2019
    P_TAG_RE_R(i_RP,1)    = RP(i_RP);   P_TAG_RE_R(i_RP,2)    = P_TAG_RE_R_Temp(index_RE_R,2);   % = max(P_TAG_RE_R_Temp(:,2)); % Modified on 29.Sep.2019
    
    P_TAG_DE_R(i_RP,1)    = RP(i_RP);   P_TAG_DE_R(i_RP,2)    = max(P_TAG_DE_R_Temp(:,2));
    P_TAG_DE_G(i_RP,1)    = RP(i_RP);   P_TAG_DE_G(i_RP,2)    = 1.0 - P_TAG_DE_R(i_RP,2); % max(P_TAG_DE_G_Temp(:,2));

end

P_TAG_G_NC_im(:,2) = P_TAG_RE_G(:,2) + P_TAG_DE_G(:,2) .* P_TAG_RE_Y(:,2);
P_TAG_Y_NC_im(:,2) = P_TAG_RE_Y(:,2);
P_TAG_R_NC_im(:,2) = P_TAG_RE_R(:,2) + P_TAG_DE_R(:,2) .* P_TAG_RE_Y(:,2);
         
     % Added: 28.Sep.2019; modified: 29.Sep.2019 -----
     dlmwrite(strcat('P_TAG_RE_G',   '.txt'), P_TAG_RE_G,    'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('P_TAG_RE_Y',   '.txt'), P_TAG_RE_Y,    'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('P_TAG_RE_R',   '.txt'), P_TAG_RE_R,    'delimiter', '\t', 'precision', 8);

     dlmwrite(strcat('P_TAG_DE_G',   '.txt'), P_TAG_DE_G,    'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('P_TAG_DE_R',   '.txt'), P_TAG_DE_R,    'delimiter', '\t', 'precision', 8);
     
     dlmwrite(strcat('P_TAG_G_NC_im','.txt'), P_TAG_G_NC_im, 'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('P_TAG_Y_NC_im','.txt'), P_TAG_Y_NC_im, 'delimiter', '\t', 'precision', 8);
     dlmwrite(strcat('P_TAG_R_NC_im','.txt'), P_TAG_R_NC_im, 'delimiter', '\t', 'precision', 8);
     
% Compute mobilization time - - - - -
R_T0_Tag_G =  10.0; % 10 days
R_T0_Tag_Y =  30.0; % 1 month = 30 days
R_T0_Tag_R = 180.0; % 3 month = 180 days

% R_T0_NC_im = R_T0_Tag_G * P_TAG_G_NC_im  +...
%              R_T0_Tag_Y * P_TAG_Y_NC_im  +...
%              R_T0_Tag_R * P_TAG_R_NC_im;
% Modified: 02.Oct.2019 (removed yellow)
R_T0_NC_im = R_T0_Tag_G * P_TAG_G_NC_im  +...
             R_T0_Tag_R * P_TAG_R_NC_im;
             
     % Added: 29.Sep.2019 -----
     dlmwrite(strcat('R_T0_NC_im',   '.txt'), R_T0_NC_im,    'delimiter', '\t', 'precision', 8);

% % % % % % % % % % % % % % % % % % %
%% COMPUTE: Total Downtime - - - - - -
% % % % % % % % % % % % % % % % % % %
% Compute R_TplusT0_NC_im (for fast and slow repair schemes)

% Fast track repair scheme - - -
R_TplusT0_NC_im_Fast       = zeros(NumIM,2);
R_TplusT0_NC_im_Fast(:, 1) = returnPeriod(:, 1);
R_TplusT0_NC_im_Fast(:, 2) = Time_RP_Fast(:, 2) + R_T0_NC_im(:, 2);

% Slow track repair scheme - - -
R_TplusT0_NC_im_Slow       = zeros(NumIM,2);
R_TplusT0_NC_im_Slow(:, 1) = returnPeriod(:, 1);
R_TplusT0_NC_im_Slow(:, 2) = Time_RP_Slow(:, 2) + R_T0_NC_im(:, 2);



% ------------------------------------------------------------------------------------------
% MULTIPLY PROBABILITY OF DEMOLISH & COLLAPSE ----------------------------------------------
% ------------------------------------------------------------------------------------------

% Source of downtimes = "Repair and tagging downtime";
%                       "Demolish downtime";
%                       "Collapse downtime"

Downtime_Repair_RP_Fast  =  zeros(NumIM+1,3);
Downtime_Repair_RP_Slow  =  zeros(NumIM+1,3);
Downtime_Demolish_RP     =  zeros(NumIM+1,3);
Downtime_Collapse_RP     =  zeros(NumIM+1,3);
Downtime_Total_RP_Fast   =  zeros(NumIM+1,3);
Downtime_Total_RP_Slow   =  zeros(NumIM+1,3);

Time_Demolish  =  38.0 * (365.0/12.0);   % 1 year is 12 months=365days. Months is changed to days.
Time_Collapse  =  Time_Demolish;

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Downtime due to repair and tagging"  - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Downtime_Repair_RP_Fast(2:end,2)  = R_TplusT0_NC_im_Fast(:,2) .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));
Downtime_Repair_RP_Slow(2:end,2)  = R_TplusT0_NC_im_Slow(:,2) .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Downtime due to demolish"            - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Downtime_Demolish_RP(2:end,2)    = Time_Demolish * PD_IM_NC(:,2) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Downtime due to collapse"            - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Downtime_Collapse_RP(2:end,2)    = Time_Collapse * PC_IM(:,2);

for i = 1:NumIM
   
    Downtime_Repair_RP_Fast(i+1,1)  =  Sa_T1(i);  Downtime_Repair_RP_Fast(i+1,3)  = RP(i);
    Downtime_Repair_RP_Slow(i+1,1)  =  Sa_T1(i);  Downtime_Repair_RP_Slow(i+1,3)  = RP(i);
    Downtime_Demolish_RP(i+1,1)     =  Sa_T1(i);  Downtime_Demolish_RP(i+1,3)     = RP(i);
    Downtime_Collapse_RP(i+1,1)     =  Sa_T1(i);  Downtime_Collapse_RP(i+1,3)     = RP(i);
    Downtime_Total_RP_Fast(i+1,1)   =  Sa_T1(i);  Downtime_Total_RP_Fast(i+1,3)   = RP(i);
    Downtime_Total_RP_Slow(i+1,1)   =  Sa_T1(i);  Downtime_Total_RP_Slow(i+1,3)   = RP(i);
    
end

% ----------------------------------------------------------------------------
% COMPUTE TOTAL DOWNTIME FOR IM=im  ----------------------------------------------
% ----------------------------------------------------------------------------
% Source of losses = "Structural repair loss";
Downtime_Total_RP_Fast(2:end,2) = Downtime_Repair_RP_Fast(2:end,2) + Downtime_Demolish_RP(2:end,2) + Downtime_Collapse_RP(2:end,2);
Downtime_Total_RP_Slow(2:end,2) = Downtime_Repair_RP_Slow(2:end,2) + Downtime_Demolish_RP(2:end,2) + Downtime_Collapse_RP(2:end,2);



% ----------------------------------------------------------------------------
% OUTPUT VULNERABILITY FUNCTIONS  --------------------------------------------
% ----------------------------------------------------------------------------
% Added 28.Sep.2019 -----
dlmwrite(strcat('Downtime_Repair_RP_Fast','.txt'),  Downtime_Repair_RP_Fast,  'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Downtime_Repair_RP_Slow','.txt'),  Downtime_Repair_RP_Slow,  'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Downtime_Demolish_RP','.txt'),     Downtime_Demolish_RP,     'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Downtime_Collapse_RP','.txt'),     Downtime_Collapse_RP,     'delimiter', '\t', 'precision', 8);

% Normalized by total repair cost - - - - -
dlmwrite(strcat('Downtime_Total_RP_Fast.txt'),  Downtime_Total_RP_Fast,  'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Downtime_Total_RP_Slow.txt'),  Downtime_Total_RP_Slow,  'delimiter', '\t', 'precision', 8);


% ----------------------------------------------------------------------------
%% COMPUTE EXPECTED ANNUAL DOWNTIME (EADT) -----------------------------------------
% ----------------------------------------------------------------------------
fprintf(strcat('Compute Expected Annual Down Time (EADT) - - - - -  \n'));

HazardCurveData  =  load(strcat(char(HazardCurve_Name))); % Load selected seismic hazard data (differs per each period - system)

Sa_1=zeros(length(Sa_T1)+1,1); Sa_2=zeros(length(Sa_T1)+1,1);

Sa_1(2:end,1)    =  Sa_T1(:);
Sa_2(1:end-1,1)  =  Sa_T1(:);
Sa_2(end,1)      =  99.9;

Sa_Ave           = (Sa_1+Sa_2)/2.0;

for choice = 1:2

switch choice
    case 1  % Downtime_Repair_RP
        Loss_Vulnerability  =  Downtime_Total_RP_Fast(:,2);
        type_of_EADT         =  'Downtime_Total_RP_Fast';
    case 2  % Loss_NonStructural_Drift_RP
        Loss_Vulnerability  =  Downtime_Total_RP_Slow(:,2);
        type_of_EADT         =  'Downtime_Total_RP_Slow';

end
    
    MeanLoss_Ave = mean([Loss_Vulnerability(1:end-1)';Loss_Vulnerability(2:end)'])'; % cdf('Lognormal', Sa_Ave(:,1), log(medianSa), beta);  
    
    lamb_Sa_1=zeros(length(Sa_T1),1); lamb_Sa_2=zeros(length(Sa_T1),1);
    
    for s = 1 : NumIM
        [c, index]        = min(abs(HazardCurveData(:,1)-Sa_T1(s)));
        lamb_Sa_1(s, 1)   = HazardCurveData(index,2);
    end
        
        lamb_Sa_2(1:end-1,1) = lamb_Sa_1(2:end,1);
        lamb_Sa_2(  end,  1) = 1.e-9;   % arbitrary very small frequency of occurance

    Delta_lamb_Sa         =  abs(lamb_Sa_2 - lamb_Sa_1);
    lamb_F_i              =  Delta_lamb_Sa .* MeanLoss_Ave;  
    lamb_F                =  sum(lamb_F_i);           % Expected Annual Loss (EADT)
    EADT                   =  lamb_F;                  % Expected Annual Loss (EADT)
    
    dlmwrite(strcat('EADT_',  char(type_of_EADT), '.txt'),    EADT, 'delimiter', '\t', 'precision', 8); % for PSDR
    eval(['EADT_',            char(type_of_EADT), '= EADT']);                  % Obtain EADT

end


% ----------------------------------------------------------------------------
%% COMPUTE INDIRECT LOSS FOR IM=im  ----------------------------------------------
% ----------------------------------------------------------------------------
% Loss from time for repair and tagging -----
% Calculate: E[DTL|NC,im] (fast or slow repair schemes)
% Note: equations are different (modified) from those in Mitrani-Resier (2007) after
% some considerations.

Loss_RP_Fast_temp = zeros(NumIM,2); Loss_RP_Fast = zeros(NumIM,2);  % Fast repair scheme
Loss_RP_Slow_temp = zeros(NumIM,2); Loss_RP_Slow = zeros(NumIM,2);  % Slow repair scheme

for i_RP = 1 : NumIM

    Loss_RP_Fast_temp(i_RP, 1) = RP(i_RP);  Loss_RP_Fast(i_RP, 1) = RP(i_RP);
    Loss_RP_Slow_temp(i_RP, 1) = RP(i_RP);  Loss_RP_Slow(i_RP, 1) = RP(i_RP); 

    eval(['LossStory_RP_',num2str(i_RP), '_Fast', '= zeros(NumStory,2)']); % Fast repair scheme
    eval(['LossStory_RP_',num2str(i_RP), '_Slow', '= zeros(NumStory,2)']); % Slow repair scheme
    
    for i_story = 1 : NumStory   % consider each story
        eval(['LossStory_RP_',num2str(i_RP), '_Fast', '(i_story,2)','=','TimeStory_RP_',num2str(i_RP), '_Fast', '(i_story,2)','* E_U_U_m']);
        eval(['LossStory_RP_',num2str(i_RP), '_Slow', '(i_story,2)','=','TimeStory_RP_',num2str(i_RP), '_Slow', '(i_story,2)','* E_U_U_m']);
    end
    
    Loss_RP_Fast_temp(i_RP, 2) = eval(['max(','LossStory_RP_',num2str(i_RP), '_Fast', '(:, 2)',')']);  % Pick max from all operational units (=stories); Fast repair scheme
    Loss_RP_Slow_temp(i_RP, 2) = eval(['sum(','LossStory_RP_',num2str(i_RP), '_Slow', '(:, 2)',')']);  % Pick max from all operational units (=stories); Slow repair scheme
    
end

Loss_RP_Fast = Loss_RP_Fast_temp  +  E_DTL_TAG*R_T0_NC_im;
Loss_RP_Slow = Loss_RP_Slow_temp  +  E_DTL_TAG*R_T0_NC_im;

% Loss from time for repair, demolish and collapse -----
% Then, calculate: E[DTL|im] (fast or slow repair schemes)

% Source of downtimes = "Repair and tagging downtime";
%                       "Demolish downtime";
%                       "Collapse downtime"

Loss_Repair_RP_Fast  =  zeros(NumIM+1,3);
Loss_Repair_RP_Slow  =  zeros(NumIM+1,3);
Loss_Demolish_RP     =  zeros(NumIM+1,3);
Loss_Collapse_RP     =  zeros(NumIM+1,3);
Loss_Total_RP_Fast   =  zeros(NumIM+1,3);
Loss_Total_RP_Slow   =  zeros(NumIM+1,3);

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Downtime due to repair and tagging"  - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_Repair_RP_Fast(2:end,2)  = Loss_RP_Fast(:,2) .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));
Loss_Repair_RP_Slow(2:end,2)  = Loss_RP_Slow(:,2) .* (1.-PD_IM_NC(:,2)) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Downtime due to demolish"            - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_Demolish_RP(2:end,2)    = E_DTL_Dem * Time_Demolish * PD_IM_NC(:,2) .* (1.-PC_IM(:,2));

% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
% COMPUTE: "Downtime due to collapse"            - - - - - - -
% % % % % % % % % % % % % % % % % % % % %  % % % % % % % % % %
Loss_Collapse_RP(2:end,2)    = E_DTL_Col * Time_Collapse* PC_IM(:,2);

for i = 1:NumIM
   
    Loss_Repair_RP_Fast(i+1,1)  =  Sa_T1(i);  Loss_Repair_RP_Fast(i+1,3)  = RP(i);
    Loss_Repair_RP_Slow(i+1,1)  =  Sa_T1(i);  Loss_Repair_RP_Slow(i+1,3)  = RP(i);
    Loss_Demolish_RP(i+1,1)     =  Sa_T1(i);  Loss_Demolish_RP(i+1,3)     = RP(i);
    Loss_Collapse_RP(i+1,1)     =  Sa_T1(i);  Loss_Collapse_RP(i+1,3)     = RP(i);
    Loss_Total_RP_Fast(i+1,1)   =  Sa_T1(i);  Loss_Total_RP_Fast(i+1,3)   = RP(i);
    Loss_Total_RP_Slow(i+1,1)   =  Sa_T1(i);  Loss_Total_RP_Slow(i+1,3)   = RP(i);
    
end

% ----------------------------------------------------------------------------
% COMPUTE TOTAL DOWNTIME FOR IM=im  ----------------------------------------------
% ----------------------------------------------------------------------------
% Source of losses = "Structural repair loss";
Loss_Total_RP_Fast(2:end,2) = Loss_Repair_RP_Fast(2:end,2) + Loss_Demolish_RP(2:end,2) + Loss_Collapse_RP(2:end,2);
Loss_Total_RP_Slow(2:end,2) = Loss_Repair_RP_Slow(2:end,2) + Loss_Demolish_RP(2:end,2) + Loss_Collapse_RP(2:end,2);


% ----------------------------------------------------------------------------
% OUTPUT VULNERABILITY FUNCTIONS  --------------------------------------------
% ----------------------------------------------------------------------------
% Added 28.Sep.2019 -----
dlmwrite(strcat('Loss_Repair_RP_Fast','.txt'),  Loss_Repair_RP_Fast,  'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Repair_RP_Slow','.txt'),  Loss_Repair_RP_Slow,  'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Demolish_RP','.txt'),     Loss_Demolish_RP,     'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Collapse_RP','.txt'),     Loss_Collapse_RP,     'delimiter', '\t', 'precision', 8);

% Normalized by total repair cost - - - - -
dlmwrite(strcat('Loss_Total_RP_Fast.txt'),      Loss_Total_RP_Fast,   'delimiter', '\t', 'precision', 8);
dlmwrite(strcat('Loss_Total_RP_Slow.txt'),      Loss_Total_RP_Slow,   'delimiter', '\t', 'precision', 8);



% ----------------------------------------------------------------------------
%% COMPUTE EXPECTED ANNUAL LOSS DUE TO DOWNTIME (EADTL) -----------------------------------------
% ----------------------------------------------------------------------------
fprintf(strcat('Compute Expected Annual Loss due to DOWN TIME (EADTL) - - - - -  \n'));

HazardCurveData  =  load(strcat(char(HazardCurve_Name))); % Load selected seismic hazard data (differs per each period - system)

Sa_1=zeros(length(Sa_T1)+1,1); Sa_2=zeros(length(Sa_T1)+1,1);

Sa_1(2:end,1)    =  Sa_T1(:);
Sa_2(1:end-1,1)  =  Sa_T1(:);
Sa_2(end,1)      =  99.9;

Sa_Ave           = (Sa_1+Sa_2)/2.0;

for choice = 1:2

switch choice
    case 1  % Downtime_Repair_RP
        Loss_Vulnerability  =  Loss_Total_RP_Fast(:,2);
        type_of_EADTL         =  'Loss_Total_RP_Fast';
    case 2  % Loss_NonStructural_Drift_RP
        Loss_Vulnerability  =  Loss_Total_RP_Slow(:,2);
        type_of_EADTL         =  'Loss_Total_RP_Slow';

end
    
    MeanLoss_Ave = mean([Loss_Vulnerability(1:end-1)';Loss_Vulnerability(2:end)'])'; % cdf('Lognormal', Sa_Ave(:,1), log(medianSa), beta);  
    
    lamb_Sa_1=zeros(length(Sa_T1),1); lamb_Sa_2=zeros(length(Sa_T1),1);
    
    for s = 1 : NumIM
        [c, index]        = min(abs(HazardCurveData(:,1)-Sa_T1(s)));
        lamb_Sa_1(s, 1)   = HazardCurveData(index,2);
    end
        
        lamb_Sa_2(1:end-1,1) = lamb_Sa_1(2:end,1);
        lamb_Sa_2(  end,  1) = 1.e-9;   % arbitrary very small frequency of occurance

    Delta_lamb_Sa         =  abs(lamb_Sa_2 - lamb_Sa_1);
    lamb_F_i              =  Delta_lamb_Sa .* MeanLoss_Ave;  
    lamb_F                =  sum(lamb_F_i);           % Expected Annual Loss (EADTL)
    EADTL                 =  lamb_F;                  % Expected Annual Loss (EADTL)
    
    dlmwrite(strcat('EADTL_',  char(type_of_EADTL), '.txt'),    EADTL, 'delimiter', '\t', 'precision', 8); % for PSDR
    eval(['EADTL_',            char(type_of_EADTL), '= EADTL']);                  % Obtain EADTL

end

