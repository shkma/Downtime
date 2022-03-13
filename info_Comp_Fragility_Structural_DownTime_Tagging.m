function [PDS_ij_EDP] = info_Comp_Fragility_Structural_DownTime_Tagging(i_n, i_m, x_PSDR_pdf, system, i_story)

% This function file returns from story number IDs (for each story):
% n = number of damage states a component may experience (this is per m)

% i_m = ID of component

    % % % % % % % % % % % % % % % %
    if i_m == 4   % Round HSS brace
    % % % % % % % % % % % % % % % %
             
             % Determining different global slenderness ratio x_GP (=KL/r; Lignos and Karamanci, 2013) 
             % for different systems and different stories (K=1; L=all same; r=varies) - - - - - - - - -
             K=1.0; L = 0.70*sqrt((15*12)^2+(12*12)^2);% Unit=in
             
               if     (i_story == 1 && system == 1) || (i_story == 2 && system == 1)
                  r = 2.89;
               elseif (i_story == 3 && system == 1) || (i_story == 4 && system == 1)
                  r = 2.32;
               elseif (i_story == 5 && system == 1)
                  r = 1.96;
               elseif (i_story == 6 && system == 1)
                  r = 1.69;
                  
               elseif (i_story == 1 && system == 3) || (i_story == 2 && system == 3)
                  r = 2.89;
               elseif (i_story == 3 && system == 3) || (i_story == 4 && system == 3)
                  r = 2.49;
               elseif (i_story == 5 && system == 3)
                  r = 1.96;
               elseif (i_story == 6 && system == 3)
                  r = 1.48;
                  
               elseif (i_story == 1 && system == 4) || (i_story == 2 && system == 4)
                  r = 2.18;
               elseif (i_story == 3 && system == 4) || (i_story == 4 && system == 4)
                  r = 1.96;
               elseif (i_story == 5 && system == 4)
                  r = 1.61;
               elseif (i_story == 6 && system == 4)
                  r = 1.34;
               end
               
               x_GP = K*L/r;
             %   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                   
       if (i_story == 1 && system == 1) || (i_story == 2 && system == 1) || (i_story == 1 && system == 3) || (i_story == 2 && system == 3) %  for first story non-isolated SCBF & isolated SCBF w. Ri=1, thicker HSS is more expensive for DS2 & DS3.
            
           if i_n == 0  % Round HSS (60kg/m < brace weight < 147kg/m); Round HSS

              xm_EDP=0.41 /100.; beta_EDP=0.51;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (KL/r) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_ij_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = 1.0 - F_DS_ij*F_DS_ij_GP;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.41 /100.; beta_EDP=0.51;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.96 /100.; beta_EDP=0.45;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (KL/r) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = 1.0 - F_DS_i2*F_DS_i2_GP;    % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.96 /100.; beta_EDP=0.45;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=2.75 /100.; beta_EDP=0.51;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = F_DS_i1*F_DS_i1_GP - F_DS_i2*F_DS_i2_GP;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=2.75 /100.; beta_EDP=0.51;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_ij_GP    = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = F_DS_ij*F_DS_ij_GP;             % j=n, i.e. biggest damage

           end
       
       else  % for others, the time will be lower.
           
           if i_n == 0  % Round HSS (brace weight < 60kg/m); Round HSS

              xm_EDP=0.41 /100.; beta_EDP=0.51;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_ij_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = 1.0 - F_DS_ij*F_DS_ij_GP;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.41 /100.; beta_EDP=0.51;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.96 /100.; beta_EDP=0.45;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=63.6; beta_GP=0.46;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = 1.0 - F_DS_i2*F_DS_i2_GP;    % 1<=j<=n, i.e. some damage

           elseif i_n == 2

              xm_EDP=0.96 /100.; beta_EDP=0.45;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=2.75 /100.; beta_EDP=0.51;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=66.1; beta_GP=0.45;
                    F_DS_i1_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_i2_GP = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
              
              PDS_ij_EDP = F_DS_i1*F_DS_i1_GP - F_DS_i2*F_DS_i2_GP;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=2.75 /100.; beta_EDP=0.51;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters

                    % Geometric parameters (D/t/lambda_hd) - - - - -
                    xm_GP=68.9; beta_GP=0.40;
                    F_DS_ij_GP    = normcdf((log(x_GP/xm_GP))/beta_GP); % compute fragility function using Eq. 1 and estimated parameters
                    %  - - - - - - - - - - - - - - - - - - - - - - -
                    
              PDS_ij_EDP = F_DS_ij*F_DS_ij_GP;    % j=n, i.e. biggest damage

           end
       
       end
           
           
    % % % % % % % % % % % % % % % % % % %
    elseif i_m == 5  % Moment connections 
    % % % % % % % % % % % % % % % % % % %
        
        if system == 1  || system == 3  || system == 4 % if SCBF (isolated OR non-isolated) 
            
        
           if i_n == 0  % Moment connection; one-sided; <= W27

              xm_EDP=0.03; beta_EDP=0.30;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.03; beta_EDP=0.30;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.04; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0   -   F_DS_i2;   %   F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage % % Modified for mobilization time or tagging..8.Sep.2019

           elseif i_n == 2

              xm_EDP=0.04; beta_EDP=0.30;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.05; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=0.05; beta_EDP=0.30;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage

           end
       
           
        elseif system == 2  || system == 5  || system == 6 % if SMRF (isolated OR non-isolated) 
            
        
           if i_n == 0  % RBS connection; one-sided; <= W27

              xm_EDP=0.01;   beta_EDP=0.17;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage

           elseif i_n == 1

              xm_EDP=0.01;   beta_EDP=0.17;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.0216; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = 1.0   -   F_DS_i2;   %   F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage % % Modified for mobilization time or tagging..8.Sep.2019

           elseif i_n == 2

              xm_EDP=0.0216; beta_EDP=0.30;
              F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              xm_EDP=0.05; beta_EDP=0.30;
              F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage

           elseif i_n == 3

              xm_EDP=0.05;   beta_EDP=0.30;
              F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
              PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage

           end
       
        end


    % % % % % % % % % % % % % % % % % %
    elseif i_m == 8   % Corrugated slab
    % % % % % % % % % % % % % % % % % %
        
       if i_n == 0  % Corrugated slab

          xm_EDP=0.00375; beta_EDP=0.13;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = 1.0 - F_DS_ij;   % j=0, i.e. no damage
          
       elseif i_n == 1
          
          xm_EDP=0.00375; beta_EDP=0.13;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.01; beta_EDP=0.22; 
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 2
               
          xm_EDP=0.01; beta_EDP=0.22;
          F_DS_i1    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          xm_EDP=0.05; beta_EDP=0.35;
          F_DS_i2    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_i1 - F_DS_i2;   % 1<=j<=n, i.e. some damage
          
       elseif i_n == 3
          
          xm_EDP=0.05; beta_EDP=0.35;
          F_DS_ij    = normcdf((log(x_PSDR_pdf/xm_EDP))/beta_EDP); % compute fragility function using Eq. 1 and estimated parameters
          PDS_ij_EDP = F_DS_ij;             % j=n, i.e. biggest damage
          
       end
   
    end
    
end

