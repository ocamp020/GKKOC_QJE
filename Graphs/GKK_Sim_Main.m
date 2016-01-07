% Graphs GKK Wealth Tax
% Sergio Ocampo Diaz


 javaaddpath('jxl.jar');
 javaaddpath('MXL.jar');

 import mymxl.*;
 import jxl.*;   

%%  Parameters

% Grids
    n_a = 201 ; 
    n_z = 7   ;
    n_l = 5   ;
    n_e = 5   ;
    
    Max_Age = 81 ;
    Ret_Age = 45 ;

% Utility and technology
    sigma = 4.00  ;
    gamma = 0.4494; 
    mu    = 0.9   ;
    delta = 0     ;
    
% Taxes
    Threshold_Factor = 0.00 ;
	tauC    = 0.075         ;
    tauK    = 0.25          ;
	
% Set Tax Regime
    % tauPL = 0.185     ;
    % psi   = 0.77      ;
    tauPL = 0.0         ;
    psi   = 0.776       ;
    tauW  = 0.017072675596579098 ;
        
% Age brackets 
    age = [5 , 15, 25, 35, 45, 55, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;
    
% Percentiles
prctl = [10, 25, 50, 75, 90, 99, 99.9];
    
%% Result_Folders

    if ((tauPL==0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Files/') ;
        Bench_Folder  = '../NSU_LT_Results/Bench_Files/' ;
    elseif ((tauPL~=0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Files/') ;
        Bench_Folder  = '../NSU_PT_Results/Bench_Files/' ;
    elseif ((tauPL==0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Files/') ;
        Bench_Folder  = '../SU_LT_Results/Bench_Files/' ;
    elseif ((tauPL~=0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Files/') ;
        Bench_Folder  = '../SU_PT_Results/Bench_Files/' ;
    end
    

%% Grids 
   
% A grid
    eval(['load ',Result_Folder,'agrid']);
    
    A_mat = repmat(agrid,[Max_Age,1,n_z,n_l,n_e]);
    
% Z grid
    eval(['load ',Result_Folder,'zgrid']);
    
    Z_mat = repmat(reshape(zgrid,[1,1,n_z,1,1]),[Max_Age,n_a,1,n_l,n_e]);
    
%% Read files - Benchmark Economy
% Distribution
    eval(['load ',Bench_Folder,'DBN'])
    DBN_bench = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear DBN
    
% Wage
    eval(['load ',Bench_Folder,'wage'])
    W_bench = wage ;
    clear wage
    
% Interst Rate
    eval(['load ',Bench_Folder,'rr'])
    R_bench = rr ;
    clear rr
    
% Interst Rate
    eval(['load ',Bench_Folder,'EBAR'])
    E_bench = EBAR ;
    clear EBAR
    
    
%% Read files - Experimental Economy

% Distribution
    eval(['load ',Result_Folder,'Exp_results_DBN']) ;
    DBN_exp = reshape(Exp_results_DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_DBN

% Wage
    eval(['load ',Result_Folder,'Exp_results_wage'])
    W_exp = Exp_results_wage ;
    clear Exp_results_wage
    
% Interst Rate
    eval(['load ',Result_Folder,'Exp_results_rr'])
    R_exp = Exp_results_rr ;
    clear Exp_results_rr
    
% Wealth Taxes
    Threshold = Threshold_Factor*E_bench ;
    
    
%% Simulate the economy

kubu_folder = '../../kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/more_stats/' ;

kubu_switch = 1 ;

if kubu_switch==1
    eval(['load ',kubu_folder,'panelage_bench']) ;
    
    eval(['load ',kubu_folder,'panela_bench']) ;
    eval(['load ',kubu_folder,'panel_return_bench']) ;
    eval(['load ',kubu_folder,'panel_hours_bench']) ;
    eval(['load ',kubu_folder,'panel_cons_bench']) ;
    eval(['load ',kubu_folder,'panel_aprime_bench']) ;
    eval(['load ',kubu_folder,'panelz_bench']) ;
%     clear panel*
    
    eval(['load ',kubu_folder,'panelage_exp']) ;
    
    eval(['load ',kubu_folder,'panela_exp']) ;
    eval(['load ',kubu_folder,'panel_return_exp']) ;
    eval(['load ',kubu_folder,'panel_hours_exp']) ;
    eval(['load ',kubu_folder,'panel_cons_exp']) ;
    eval(['load ',kubu_folder,'panel_aprime_exp']) ;
    eval(['load ',kubu_folder,'panelz_exp']) ;
    
    Sim_age_ben = panelage_bench; 
    Sim_A_ben   = panela_bench; 
    Sim_R_ben   = 100*panel_return_bench;
    Sim_H_ben   = panel_hours_bench; 
    Sim_C_ben   = panel_cons_bench ; 
    Sim_Ap_ben  = panel_aprime_bench ;
    Sim_Z_ben   = panelz_bench ;
    
    Sim_age_exp = panelage_exp; 
    Sim_A_exp   = panela_exp; 
    Sim_R_exp   = 100*panel_return_exp;
    Sim_H_exp   = panel_hours_exp; 
    Sim_C_exp   = panel_cons_exp ; 
    Sim_Ap_exp  = panel_aprime_exp ;
    Sim_Z_exp   = panelz_exp ;
    
    clear panel*
else 
% Load simulation resutls 
    eval(['load ',Result_Folder,'Simul/Sim_age_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_A_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_C_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_H_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_R_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_Z_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_Ap_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_R_at_ben'])
    eval(['load ',Result_Folder,'Simul/Sim_Yh_ben'])

    eval(['load ',Result_Folder,'Simul/Sim_age_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_A_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_C_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_H_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_R_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_Z_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_Ap_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_R_at_exp'])
    eval(['load ',Result_Folder,'Simul/Sim_Yh_exp'])
    
    Sim_R_ben = 100*(Sim_R_ben-1) ;
    Sim_R_exp = 100*(Sim_R_exp-1) ;
end 
    
    Sim_S_ben = 100*(Sim_Ap_ben./Sim_A_ben-1) ;
    Sim_S_exp = 100*(Sim_Ap_exp./Sim_A_exp-1) ;
    
    
%%  Mean and std by age and age-z
    for i=1:Max_Age
        ind = (Sim_age_ben==i) ;
        n   = sum(ind)         ;
        sim_mean_A_ben(i)  = mean(Sim_A_ben(ind))  ;
        sim_mean_C_ben(i)  = mean(Sim_C_ben(ind))  ;
        sim_mean_H_ben(i)  = mean(Sim_H_ben(ind))  ;
        sim_mean_Ap_ben(i) = mean(Sim_Ap_ben(ind)) ;
        sim_mean_S_ben(i)  = 100*(sim_mean_Ap_ben(i)/sim_mean_A_ben(i)-1) ;
        sim_mean_R_ben(i)     = mean(Sim_R_ben(ind))    ;
        %sim_mean_R_at_ben(i)  = mean(Sim_R_at_ben(ind))  ;
        sim_mean_wR_ben(i)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
        %sim_mean_wR_at_ben(i) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;
        sim_mean_S_ben_aux(i) = mean(Sim_S_ben(ind)) ;
        
        sim_std_A_ben(i)  = std(Sim_A_ben(ind)) ;
        sim_std_C_ben(i)  = std(Sim_C_ben(ind)) ;
        sim_std_H_ben(i)  = std(Sim_H_ben(ind)) ;
        sim_std_Ap_ben(i) = std(Sim_Ap_ben(ind));
        sim_std_R_ben(i)  = std(Sim_R_ben(ind)) ;
        %sim_std_R_at_ben(i)  = std(Sim_R_at_ben(ind)) ;
        sim_std_wR_ben(i) = ( sum( (Sim_R_ben(ind) - sim_mean_wR_ben(i)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;
        %sim_std_wR_at_ben(i) = ( sum( (Sim_R_at_ben(ind) - sim_mean_wR_at_ben(i)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;
        sim_var_lnC_ben(i) = var(log(Sim_C_ben(ind))) ;
        
        ind = (Sim_age_exp==i) ;
        n   = sum(ind)         ;
        sim_mean_A_exp(i)  = mean(Sim_A_exp(ind))  ;
        sim_mean_C_exp(i)  = mean(Sim_C_exp(ind))  ;
        sim_mean_H_exp(i)  = mean(Sim_H_exp(ind))  ;
        sim_mean_Ap_exp(i) = mean(Sim_Ap_exp(ind)) ;
        sim_mean_S_exp(i)  = 100*(sim_mean_Ap_exp(i)/sim_mean_A_exp(i)-1) ;
        sim_mean_R_exp(i)     = mean(Sim_R_exp(ind))     ;
        %sim_mean_R_at_exp(i)  = mean(Sim_R_at_exp(ind))  ;
        sim_mean_wR_exp(i)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
        %sim_mean_wR_at_exp(i) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;
        sim_mean_S_exp_aux(i) = mean(Sim_S_exp(ind)) ;
        
        sim_std_A_exp(i)  = std(Sim_A_exp(ind)) ;
        sim_std_C_exp(i)  = std(Sim_C_exp(ind)) ;
        sim_std_H_exp(i)  = std(Sim_H_exp(ind)) ;
        sim_std_Ap_exp(i) = std(Sim_Ap_exp(ind));
        sim_std_R_exp(i)  = std(Sim_R_exp(ind)) ;
        %sim_std_R_at_exp(i)  = std(Sim_R_at_exp(ind)) ;
        sim_std_wR_exp(i)    = ( sum( (Sim_R_exp(ind) - sim_mean_wR_exp(i)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ;
        %sim_std_wR_at_exp(i) = ( sum( (Sim_R_at_exp(ind) - sim_mean_wR_at_exp(i)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ;
        
    
        for j=1:n_z
            ind = ((Sim_age_ben==i).*(Sim_Z_ben==j))==1 ;
            n   = sum(ind)                         ;
            sim_mean_A_ben_AZ(i,j)  = mean(Sim_A_ben(ind))  ;
            sim_mean_C_ben_AZ(i,j)  = mean(Sim_C_ben(ind))  ;
            sim_mean_H_ben_AZ(i,j)  = mean(Sim_H_ben(ind))  ;
            sim_mean_Ap_ben_AZ(i,j) = mean(Sim_Ap_ben(ind)) ;
            sim_mean_S_ben_AZ(i,j)  = 100*(sim_mean_Ap_ben_AZ(i,j)/sim_mean_A_ben_AZ(i,j)-1) ;
            sim_mean_R_ben_AZ(i,j)  = mean(Sim_R_ben(ind))  ;
            %sim_mean_R_at_ben_AZ(i,j)  = mean(Sim_R_at_ben(ind))  ;
            sim_mean_wR_ben_AZ(i,j)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
            %sim_mean_wR_at_ben_AZ(i,j) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;

            sim_std_A_ben_AZ(i,j)  = std(Sim_A_ben(ind)) ;
            sim_std_C_ben_AZ(i,j)  = std(Sim_C_ben(ind)) ;
            sim_std_H_ben_AZ(i,j)  = std(Sim_H_ben(ind)) ;
            sim_std_Ap_ben_AZ(i,j) = std(Sim_Ap_ben(ind));
            sim_std_R_ben_AZ(i,j)  = std(Sim_R_ben(ind)) ;
            %sim_std_R_at_ben_AZ(i,j)  = std(Sim_R_at_ben(ind) - sim_mean_R_at_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_wR_ben_AZ(i,j)    = ( sum( (Sim_R_ben(ind) - sim_mean_wR_ben_AZ(i,j)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;
            %sim_std_wR_at_ben_AZ(i,j) = ( sum( (Sim_R_at_ben(ind) - sim_mean_wR_at_ben_AZ(i,j)).^2 .*Sim_A_ben(ind) )/sum(Sim_A_ben(ind)) )^0.5 ;

            ind = ((Sim_age_exp==i).*(Sim_Z_exp==j))==1 ;
            n   = sum(ind)                         ;
            sim_mean_A_exp_AZ(i,j)  = mean(Sim_A_exp(ind))  ;
            sim_mean_C_exp_AZ(i,j)  = mean(Sim_C_exp(ind))  ;
            sim_mean_H_exp_AZ(i,j)  = mean(Sim_H_exp(ind))  ;
            sim_mean_Ap_exp_AZ(i,j) = mean(Sim_Ap_exp(ind)) ;
            sim_mean_S_exp_AZ(i,j)  = 100*(sim_mean_Ap_exp_AZ(i,j)/sim_mean_A_exp_AZ(i,j)-1) ;
            sim_mean_R_exp_AZ(i,j)     = mean(Sim_R_exp(ind))     ;
            %sim_mean_R_at_exp_AZ(i,j)  = mean(Sim_R_at_exp(ind))  ;
            sim_mean_wR_exp_AZ(i,j)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
            %sim_mean_wR_at_exp_AZ(i,j) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;

            sim_std_A_exp_AZ(i,j)  = std(Sim_A_exp(ind)) ;
            sim_std_C_exp_AZ(i,j)  = std(Sim_C_exp(ind)) ;
            sim_std_H_exp_AZ(i,j)  = std(Sim_H_exp(ind)) ; 
            sim_std_Ap_exp_AZ(i,j) = std(Sim_Ap_exp(ind));
            sim_std_R_exp_AZ(i,j)  = std(Sim_R_exp(ind)) ;
            %sim_std_R_at_exp_AZ(i,j)  = std(Sim_R_at_exp(ind) - sim_mean_R_at_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
            sim_std_wR_exp_AZ(i,j)    = ( sum( (Sim_R_exp(ind) - sim_mean_wR_exp_AZ(i,j)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ;
            %sim_std_wR_at_exp_AZ(i,j) = ( sum( (Sim_R_at_exp(ind) - sim_mean_wR_at_exp_AZ(i,j)).^2 .*Sim_A_exp(ind) )/sum(Sim_A_exp(ind)) )^0.5 ; 
        end
    end

    
%% 
    for j=1:n_age
        
        
        for i=1:numel(prctl)
            if j==1
                ind = (Sim_age_ben<age(j))  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_age_ben<age(j)))==1  ;
            else
                ind = (Sim_age_ben>=age(j-1))  ;
            end
            prc_A_ben(j,i)    = prctile(Sim_A_ben(ind),prctl(i))  ;
            prc_C_ben(j,i)    = prctile(Sim_C_ben(ind),prctl(i))  ;
            prc_H_ben(j,i)    = prctile(Sim_H_ben(ind),prctl(i))  ;
            prc_Ap_ben(j,i)   = prctile(Sim_Ap_ben(ind),prctl(i)) ;
            prc_R_ben(j,i)    = prctile(Sim_R_ben(ind),prctl(i))    ;
            %prc_R_at_ben(j,i) = prctile(Sim_R_at_ben(ind),prctl(i)) ;
            
            if j==1
                ind = (Sim_age_exp<age(j))  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_age_exp<age(j)))==1  ;
            else
                ind = (Sim_age_exp>=age(j-1))  ;
            end
            prc_A_exp(j,i)    = prctile(Sim_A_exp(ind),prctl(i))  ;
            prc_C_exp(j,i)    = prctile(Sim_C_exp(ind),prctl(i))  ;
            prc_H_exp(j,i)    = prctile(Sim_H_exp(ind),prctl(i))  ;
            prc_Ap_exp(j,i)   = prctile(Sim_Ap_exp(ind),prctl(i)) ;
            prc_R_exp(j,i)    = prctile(Sim_R_exp(ind),prctl(i))    ;
            %prc_R_at_exp(j,i) = prctile(Sim_R_at_exp(ind),prctl(i)) ;
        end
        
        for i=1:n_z
            if j==1
                ind = ((Sim_age_ben<age(j)).*(Sim_Z_ben==i))==1  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_age_ben<age(j)).*(Sim_Z_ben==i))==1  ;
            else
                ind = ((Sim_age_ben>=age(j-1)).*(Sim_Z_ben==i))==1  ;
            end
            n   = sum(ind)                              ;
            A_AZ_ben(j,i)   = mean(Sim_A_ben(ind))  ;
            C_AZ_ben(j,i)   = mean(Sim_C_ben(ind))  ;
            H_AZ_ben(j,i)   = mean(Sim_H_ben(ind))  ;
            Ap_AZ_ben(j,i)  = mean(Sim_Ap_ben(ind)) ;
            S_AZ_ben(j,i)   = 100*(Ap_AZ_ben(j,i)/A_AZ_ben(j,i)-1) ;
            R_AZ_ben(j,i)     = mean(Sim_R_ben(ind))     ;
            %R_at_AZ_ben(j,i)  = mean(Sim_R_at_ben(ind))  ;
            wR_AZ_ben(j,i)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
            %wR_at_AZ_ben(j,i) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;
        
            if j==1
                ind = ((Sim_age_exp<age(j)).*(Sim_Z_exp==i))==1  ;
            elseif j<n_age && j>1
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_age_exp<age(j)).*(Sim_Z_exp==i))==1  ;
            else
                ind = ((Sim_age_exp>=age(j-1)).*(Sim_Z_exp==i))==1  ;
            end
            n   = sum(ind)                              ;
            A_AZ_exp(j,i)   = mean(Sim_A_exp(ind))  ;
            C_AZ_exp(j,i)   = mean(Sim_C_exp(ind))  ;
            H_AZ_exp(j,i)   = mean(Sim_H_exp(ind))  ;
            Ap_AZ_exp(j,i)  = mean(Sim_Ap_exp(ind)) ;
            S_AZ_exp(j,i)   = 100*(Ap_AZ_exp(j,i)/A_AZ_exp(j,i)-1) ;
            R_AZ_exp(j,i)     = mean(Sim_R_exp(ind))     ;
            %R_at_AZ_exp(j,i)  = mean(Sim_R_at_exp(ind))  ;
            wR_AZ_exp(j,i)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
            %wR_at_AZ_exp(j,i) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;
        end
        
        for i=1:n_z
            if j==1
                age_ind = (Sim_age_ben<age(j))  ;
            elseif j<n_age && j>1
                age_ind = ((Sim_age_ben>=age(j-1)).*(Sim_age_ben<age(j)))==1  ;
            else
                age_ind = (Sim_age_ben>=age(j-1))  ;
            end
            
            if i==1
                A_ind = (Sim_A_ben<=prc_A_ben(j,i))  ;
            elseif i<n_z && i>1
                A_ind = ((Sim_A_ben<=prc_A_ben(j,i)).*(Sim_A_ben>prc_A_ben(j,i-1)))==1  ;
            else
                A_ind = (Sim_A_ben>prc_A_ben(j,i))  ;
            end
            
            ind = (age_ind.*A_ind)==1 ;
            
            n   = sum(ind)                              ;
            
            if n~=0
            A_AW_ben(j,i)   = mean(Sim_A_ben(ind))  ;
            C_AW_ben(j,i)   = mean(Sim_C_ben(ind))  ;
            H_AW_ben(j,i)   = mean(Sim_H_ben(ind))  ;
            Ap_AW_ben(j,i)  = mean(Sim_Ap_ben(ind)) ;
            S_AW_ben(j,i)   = 100*(Ap_AW_ben(j,i)/A_AW_ben(j,i)-1) ;
            R_AW_ben(j,i)     = mean(Sim_R_ben(ind))     ;
            %R_at_AW_ben(j,i)  = mean(Sim_R_at_ben(ind))  ;
            wR_AW_ben(j,i)    = sum(Sim_R_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind))    ;
            %wR_at_AW_ben(j,i) = sum(Sim_R_at_ben(ind).*Sim_A_ben(ind))/sum(Sim_A_ben(ind)) ;
            for k=1:n_z
                eval(strcat('Z',int2str(k),'_AW_ben(j,i) = 100*mean(Sim_Z_ben(ind)==',int2str(k),');'))
                
                if eval(strcat('sum(Sim_Z_ben(ind)==',int2str(k),')'))>n
                        beep
                        [j,i,k]
                        error('this is wrong')
                    end 
            end
            else 
                A_AW_ben(j,i)     = A_AW_ben(j,i-1)    ;
                C_AW_ben(j,i)     = C_AW_ben(j,i-1)    ; 
                H_AW_ben(j,i)     = H_AW_ben(j,i-1)    ;
                Ap_AW_ben(j,i)    = Ap_AW_ben(j,i-1)   ;
                S_AW_ben(j,i)     = S_AW_ben(j,i-1)    ;
                R_AW_ben(j,i)     = R_AW_ben(j,i-1)    ;
                %R_at_AW_ben(j,i)  = R_at_AW_ben(j,i-1) ;
                wR_AW_ben(j,i)    = wR_AW_ben(j,i-1)   ;
                %wR_at_AW_ben(j,i) = wR_at_AW_ben(j,i-1);
                for k=1:n_z
                    eval(strcat('Z',int2str(k),'_AW_ben(j,i) = Z',int2str(k),'_AW_ben(j,i-1);'))
                end
            end
            
            if j==1
                age_ind = (Sim_age_exp<age(j))  ;
            elseif j<n_age && j>1
                age_ind = ((Sim_age_exp>=age(j-1)).*(Sim_age_exp<age(j)))==1  ;
            else
                age_ind = (Sim_age_exp>=age(j-1))  ;
            end
            
            if i==1
                A_ind = (Sim_A_exp<=prc_A_exp(j,i))  ;
            elseif i<n_z && i>1
                A_ind = ((Sim_A_exp<=prc_A_exp(j,i)).*(Sim_A_exp>prc_A_exp(j,i-1)))==1  ;
            else
                A_ind = (Sim_A_exp>prc_A_exp(j,i))  ;
            end
            
            ind = (age_ind.*A_ind)==1 ;
            
            n   = sum(ind)                            ;
            if n~=0
                A_AW_exp(j,i)     = mean(Sim_A_exp(ind))  ;
                C_AW_exp(j,i)     = mean(Sim_C_exp(ind))  ;
                H_AW_exp(j,i)     = mean(Sim_H_exp(ind))  ;
                Ap_AW_exp(j,i)    = mean(Sim_Ap_exp(ind)) ;
                S_AW_exp(j,i)     = 100*(Ap_AW_exp(j,i)/A_AW_exp(j,i)-1) ;
                R_AW_exp(j,i)     = mean(Sim_R_exp(ind))     ;
                %R_at_AW_exp(j,i)  = mean(Sim_R_at_exp(ind))  ;
                wR_AW_exp(j,i)    = sum(Sim_R_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind))    ;
                %wR_at_AW_exp(j,i) = sum(Sim_R_at_exp(ind).*Sim_A_exp(ind))/sum(Sim_A_exp(ind)) ;
                for k=1:n_z
                    eval(strcat('Z',int2str(k),'_AW_exp(j,i) = 100*mean(Sim_Z_exp(ind)==',int2str(k),');'))
                    
                    if eval(strcat('sum(Sim_Z_exp(ind)==',int2str(k),')'))>n
                        beep
                        [j,i,k]
                        error('this is wrong')
                    end 
                end
            else
                A_AW_exp(j,i)     = A_AW_exp(j,i-1)    ;
                C_AW_exp(j,i)     = C_AW_exp(j,i-1)    ; 
                H_AW_exp(j,i)     = H_AW_exp(j,i-1)    ;
                Ap_AW_exp(j,i)    = Ap_AW_exp(j,i-1)   ;
                S_AW_exp(j,i)     = S_AW_exp(j,i-1)    ;
                R_AW_exp(j,i)     = R_AW_exp(j,i-1)    ;
                %R_at_AW_exp(j,i)  = R_at_AW_exp(j,i-1) ;
                wR_AW_exp(j,i)    = wR_AW_exp(j,i-1)   ;
                %wR_at_AW_exp(j,i) = wR_at_AW_exp(j,i-1);
                for k=1:n_z
                    eval(strcat('Z',int2str(k),'_AW_exp(j,i) = Z',int2str(k),'_AW_exp(j,i-1);'))
                end
            end 
        end
    end 

%% Wealth Distribution
    n_w = 1000 ;
    wealth_grid = linspace(log(min(Sim_A_ben)),log(max(Sim_A_ben)),n_w)' ;
    Pareto_ben = NaN(n_w,1) ;
    Pareto_exp = NaN(n_w,1) ;
    N = numel(Sim_A_ben) ;
    log_A_ben = log(Sim_A_ben) ;
    log_A_exp = log(Sim_A_exp) ;
    for i=1:n_w
       Pareto_ben(i) = log(sum(log_A_ben>=wealth_grid(i))/N) ;
       Pareto_exp(i) = log(sum(log_A_exp>=wealth_grid(i))/N) ;
    end
    
    figure;
    plot(wealth_grid,[Pareto_ben Pareto_exp]); xlim([wealth_grid(1),wealth_grid(end)])
    xlabel('Log(Wealth)'); legend('Bench','Exp','location','northeast')
    print('-dpdf','./Simulation/Wealth_Pareto.pdf') ;
    
    prctl_vec = linspace(1,99,99) ;
    prctl_W_ben = NaN(99,1) ;
    prctl_W_exp = NaN(99,1) ;
    for i=1:99
        prctl_W_ben(i) = prctile(log_A_ben,prctl_vec(i)) ;
        prctl_W_exp(i) = prctile(log_A_exp,prctl_vec(i)) ;
    end
    
    figure;
    subplot(2,1,1); plot(prctl_vec,[prctl_W_ben prctl_W_exp]);
    legend('Bench','Exp','location','northeast'); title('Percentiles of Wealth Dist.')
    subplot(2,1,2); plot(prctl_vec,(prctl_W_exp-prctl_W_ben));
    xlabel('percentiles'); title('Change in Percentiles of Wealth Dist.')
    print('-dpdf','./Simulation/Wealth_Prctl.pdf') ;

%% Simulation Graphs  
    z_vec = [2 4 6] ;

    figure; 
        subplot(1,4,1); plot(20:100,[sim_mean_A_ben' sim_mean_A_exp']); xlim([19+1,19+Max_Age]); title('Mean Assets')
        subplot(1,4,2); plot(20:80,[sim_mean_S_ben(1:end-20)' sim_mean_S_exp(1:end-20)']); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate')
        subplot(1,4,3); plot(20:100,[sim_mean_C_ben' sim_mean_C_exp']); xlim([19+1,19+Max_Age]); title('Mean Consumption')
        subplot(1,4,4); plot(20:100,[sim_mean_H_ben' sim_mean_H_exp']); xlim([19+1,19+Max_Age]); title('Mean Hours')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Mean_Variables_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,4,1); plot(20:100,sim_mean_A_ben_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Assets - bench')
        subplot(2,4,2); plot(20:80,sim_mean_S_ben_AZ(1:end-20,z_vec)'); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate - bench')
        subplot(2,4,3); plot(20:100,sim_mean_C_ben_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Consumption - bench')
        subplot(2,4,4); plot(20:100,sim_mean_H_ben_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Hours - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,4,5); plot(20:100,sim_mean_A_exp_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Assets - exp')
        subplot(2,4,6); plot(20:80,sim_mean_S_exp_AZ(1:end-20,z_vec)'); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate - exp')
        subplot(2,4,7); plot(20:100,sim_mean_C_exp_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Consumption - exp')
        subplot(2,4,8); plot(20:100,sim_mean_H_exp_AZ(:,z_vec)'); xlim([19+1,19+Max_Age]); title('Mean Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Mean_Variables_Life_Cycle_by_Z.pdf') ;
    
    figure; 
        subplot(1,3,1); plot(20:100,[sim_std_A_ben' sim_std_A_exp']); xlim([19+1,19+Max_Age]); title('Std Assets')
        subplot(1,3,2); plot(20:100,[sim_std_C_ben' sim_std_C_exp']); xlim([19+1,19+Max_Age]); title('Std Consumption')
        subplot(1,3,3); plot(20:100,[sim_std_H_ben' sim_std_H_exp']); xlim([19+1,19+Max_Age]); title('Std Hours')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Std_Variables_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,3,1); plot(20:100,sim_std_A_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Assets - bench')
        subplot(2,3,2); plot(20:100,sim_std_C_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Consumption - bench')
        subplot(2,3,3); plot(20:100,sim_std_H_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Hours - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,3,4); plot(20:100,sim_std_A_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Assets - exp')
        subplot(2,3,5); plot(20:100,sim_std_C_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Consumption - exp')
        subplot(2,3,6); plot(20:100,sim_std_H_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Std_Variables_Life_Cycle_by_Z.pdf') ;
    
    figure; 
        subplot(1,2,1); plot(20:100,[(sim_mean_wR_ben-1)' (sim_mean_wR_exp-1)']); xlim([19+1,19+Max_Age]); title('Mean Return')
        subplot(1,2,2); plot(20:100,[sim_std_wR_ben' sim_std_wR_exp']); xlim([19+1,19+Max_Age]); title('Std Return')
        %subplot(2,2,3); plot(20:100,[(sim_mean_wR_at_ben-1)' (sim_mean_wR_at_exp-1)']); xlim([19+1,19+Max_Age]); title('Mean Return - After tax')
        %subplot(2,2,4); plot(20:100,[sim_std_wR_at_ben' sim_std_wR_at_exp']); xlim([19+1,19+Max_Age]); title('Std Return - After tax')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Return_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,(sim_mean_wR_ben_AZ(:,z_vec)-1)); xlim([19+1,19+Max_Age]); title('Mean Return - bench')
        subplot(2,2,2); plot(20:100,sim_std_wR_ben_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,2,3); plot(20:100,(sim_mean_wR_exp_AZ(:,z_vec)-1)); xlim([19+1,19+Max_Age]); title('Mean Return - exp')
        subplot(2,2,4); plot(20:100,sim_std_wR_exp_AZ(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Return_Life_Cycle_by_Z.pdf') ;

    figure;
        plot(20:(19+Ret_Age-1), [sim_mean_H_ben(1:Ret_Age-1)' sim_mean_H_exp(1:Ret_Age-1)'])
        xlim([19+1,19+Ret_Age-1]); title('Mean Hours')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','./Simulation/Mean_Hours_Life_Cycle.pdf') ;
    
    figure; 
        subplot(1,2,1); plot(20:(19+Ret_Age-1),sim_mean_H_ben_AZ(1:Ret_Age-1,z_vec)); 
        xlim([19+1,19+Ret_Age-1]); title('Mean Hours - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(1,2,2); plot(20:(19+Ret_Age-1),sim_mean_H_exp_AZ(1:Ret_Age-1,z_vec)); 
        xlim([19+1,19+Ret_Age-1]); title('Mean Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','./Simulation/Mean_Hours_Life_Cycle_by_Z.pdf') ;
    
    
%% Simulation Tables
    prc_file = './Simulation/prc_Tables.xls' ;
    AZ_file  = './Simulation/AZ_Tables.xls'  ;
    AW_file  = './Simulation/AW_Tables.xls'  ;

% Titles 
    % Weights by z
        Weights_Z_bench = NaN(1,n_z) ;
        Weights_Z_exp   = NaN(1,n_z) ;
        for i=1:n_z
            Weights_Z_bench(i) = sum(sum(sum(sum(DBN_bench(:,:,i,:,:))))) ;
            Weights_Z_exp(i)   = sum(sum(sum(sum(DBN_exp(:,:,i,:,:)))))   ;
        end 

    % Weights by age
        Weights_Age_bench = zeros(n_age,1) ;
        Weights_Age_exp   = zeros(n_age,1) ;
        j = 1;
        for i=1:Max_Age
            if i<=age(j)
            Weights_Age_bench(j) = Weights_Age_bench(j) + sum(sum(sum(sum(DBN_bench(i,:,:,:,:))))) ;
            Weights_Age_exp(j)   = Weights_Age_exp(j)   + sum(sum(sum(sum(DBN_exp(i,:,:,:,:))))) ;
            end 
            if i==age(j); j=j+1; end 
        end 

    % Weights by age-Z
        Weights_AZ_bench = zeros(n_age,n_z) ;
        Weights_AZ_exp   = zeros(n_age,n_z) ;
        for i=1:n_z
            j=1;
            for l=1:Max_Age
                if l<=age(j)
                Weights_AZ_bench(j,i) = Weights_AZ_bench(j,i) + sum(sum(sum(sum(DBN_bench(l,:,i,:,:))))) ;
                Weights_AZ_exp(j,i)   = Weights_AZ_exp(j,i)   + sum(sum(sum(sum(DBN_exp(l,:,i,:,:))))) ;
                end 
                if l==age(j); j=j+1; end 
            end
        end 

    Weights_Z   = 100*Weights_Z_bench   ; 
    Weights_Age = 100*Weights_Age_bench ; 
    Weights_AZ  = 100*Weights_AZ_bench  ;
    clear Weights_Z_bench Weights_Z_exp Weights_Age_bench Weights_Age_exp Weights_AZ_bench Weights_AZ_exp
    
    z_title{1} = ' ' ;
    for i=1:n_z
        z_title{i+1} = strcat('z',int2str(i),'(',num2str(Weights_Z(i),'%2.1f'),')') ;
    end 
    
    
    age_title{1,1} = strcat('<25 (',num2str(Weights_Age(1),'%2.1f'),')') ;
    age_title{2,1} = strcat('25-34 (',num2str(Weights_Age(2),'%2.1f'),')') ;
    age_title{3,1} = strcat('35-44 (',num2str(Weights_Age(3),'%2.1f'),')') ;
    age_title{4,1} = strcat('45-54 (',num2str(Weights_Age(4),'%2.1f'),')') ;
    age_title{5,1} = strcat('55-64 (',num2str(Weights_Age(5),'%2.1f'),')') ;
    age_title{6,1} = strcat('65-74 (',num2str(Weights_Age(6),'%2.1f'),')') ;
    age_title{7,1} = strcat('>75 (',num2str(Weights_Age(7),'%2.1f'),')') ;
    
% Percentiles 
    p_title{1} = ' ' ;
    for i=1:n_z
        p_title{i+1} = strcat('p',num2str(prctl(i))) ;
    end
    
    Mat_A_ben = [{'prc_A_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_A_ben)] ;
    Mat_C_ben = [{'prc_C_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_C_ben)] ;
    Mat_H_ben = [{'prc_H_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_H_ben)] ;
    Mat_R_ben = [{'prc_R_ben',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_ben)] ;
    %Mat_R_at_ben = [{'prc_R_at_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(prc_R_at_ben)] ;
    
    Mat_A_exp = [{'prc_A_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_A_exp)] ;
    Mat_C_exp = [{'prc_C_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_C_exp)] ;
    Mat_H_exp = [{'prc_H_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_H_exp)] ;
    Mat_R_exp = [{'prc_R_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_exp)] ;
    %Mat_R_at_exp = [{'prc_R_at_exp',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_at_exp)] ;
    
    Mat_A_diff = [{'prc_A_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(100*(prc_A_exp./prc_A_ben-1))] ;
    Mat_C_diff = [{'prc_C_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(100*(prc_C_exp./prc_C_ben-1))] ;
    Mat_H_diff = [{'prc_H_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(100*(prc_H_exp./prc_H_ben-1))] ;
    Mat_R_diff = [{'prc_R_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_exp-prc_R_ben)] ;
    %Mat_R_at_diff = [{'prc_R_at_diff',' ',' ',' ',' ',' ',' ',' '};p_title;age_title,num2cell(prc_R_at_exp-prc_R_at_ben)] ;
    
   
    Mat_7 = cell(1,8*3+3) ;
    Mat_8 = cell(9,1) ;
    
    Mat = [Mat_7;Mat_8 Mat_A_ben Mat_8 Mat_A_exp Mat_8 Mat_A_diff] ;
    status = xlwrite(prc_file,Mat,'Assets') ;
    
    Mat = [Mat_7;Mat_8 Mat_C_ben Mat_8 Mat_C_exp Mat_8 Mat_C_diff] ;
    status = xlwrite(prc_file,Mat,'Cons') ;
    
    Mat = [Mat_7;Mat_8 Mat_H_ben Mat_8 Mat_H_exp Mat_8 Mat_H_diff] ;
    status = xlwrite(prc_file,Mat,'Hours') ;
    
    Mat = [Mat_7;Mat_8 Mat_R_ben Mat_8 Mat_R_exp Mat_8 Mat_R_diff] ;
    status = xlwrite(prc_file,Mat,'Return') ;
    
%     Mat = [Mat_7;Mat_8 Mat_R_at_ben Mat_8 Mat_R_at_exp Mat_8 Mat_R_at_diff] ;
%     status = xlwrite(prc_file,Mat,'Return - AT') ;
    
% AZ Tables
    Mat_A_ben = [{'AZ_A_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_ben)] ;
    Mat_C_ben = [{'AZ_C_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_ben)] ;
    Mat_H_ben = [{'AZ_H_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_ben)] ;
    Mat_S_ben = [{'AZ_S_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_ben)] ;
    Mat_wR_ben = [{'AZ_wR_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell((wR_AZ_ben))] ;
    %Mat_wR_at_ben = [{'AZ_wR_at_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell((wR_at_AZ_ben-1))] ;
    
    Mat_A_exp = [{'AZ_A_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_exp)] ;
    Mat_C_exp = [{'AZ_C_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_exp)] ;
    Mat_H_exp = [{'AZ_H_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_exp)] ;
    Mat_S_exp = [{'AZ_S_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_exp)] ;
    Mat_wR_exp = [{'AZ_wR_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(wR_AZ_exp)] ;
    %Mat_wR_at_exp = [{'AZ_wR_at_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(wR_at_AZ_exp)] ;
    
    Mat_A_diff = [{'AZ_A_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(A_AZ_exp./A_AZ_ben-1))] ;
    Mat_C_diff = [{'AZ_C_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(C_AZ_exp./C_AZ_ben-1))] ;
    Mat_H_diff = [{'AZ_H_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(100*(H_AZ_exp./H_AZ_ben-1))] ;
    Mat_S_diff = [{'AZ_S_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_exp-S_AZ_ben)] ;
    Mat_wR_diff = [{'AZ_wR_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell((wR_AZ_exp-wR_AZ_ben))] ;
    %Mat_wR_at_diff = [{'AZ_wR_at_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell((wR_at_AZ_exp-wR_at_AZ_ben))] ;
    
    Mat = [Mat_7;Mat_8 Mat_A_ben Mat_8 Mat_A_exp Mat_8 Mat_A_diff] ;
    status = xlwrite(AZ_file,Mat,'Assets') ;
    
    Mat = [Mat_7;Mat_8 Mat_S_ben Mat_8 Mat_S_exp Mat_8 Mat_S_diff] ;
    status = xlwrite(AZ_file,Mat,'Saving Rate') ;
    
    Mat = [Mat_7;Mat_8 Mat_C_ben Mat_8 Mat_C_exp Mat_8 Mat_C_diff] ;
    status = xlwrite(AZ_file,Mat,'Cons') ;
    
    Mat = [Mat_7;Mat_8 Mat_H_ben Mat_8 Mat_H_exp Mat_8 Mat_H_diff] ;
    status = xlwrite(AZ_file,Mat,'Hours') ;
    
    Mat = [Mat_7;Mat_8 Mat_wR_ben Mat_8 Mat_wR_exp Mat_8 Mat_wR_diff] ;
    status = xlwrite(AZ_file,Mat,'Return') ;
    
%     Mat = [Mat_7;Mat_8 Mat_wR_at_ben Mat_8 Mat_wR_at_exp Mat_8 Mat_wR_at_diff] ;
%     status = xlwrite(AZ_file,Mat,'Return - AT') ;
    
% AW Tables 
    AW_title{1} = ' ' ;
        AW_title{2} = strcat('A<p',num2str(prctl(1))) ;
    for i=2:n_z-1
        AW_title{i+1} = strcat('p',num2str(prctl(i-1)),'<A<p',num2str(prctl(i))) ;
    end 
        AW_title{n_z+1} = strcat('p',num2str(prctl(i)),'<A') ;
    
    
    Mat_A_ben = [{'AW_A_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(A_AW_ben)] ;
    Mat_C_ben = [{'AW_C_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(C_AW_ben)] ;
    Mat_H_ben = [{'AW_H_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(H_AW_ben)] ;
    Mat_S_ben = [{'AW_S_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(S_AW_ben)] ;
    Mat_wR_ben = [{'AW_wR_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(wR_AW_ben)] ;
    %Mat_wR_at_ben = [{'AW_wR_at_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell((wR_at_AW_ben-1))] ;
    
    Mat_A_exp = [{'AW_A_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(A_AW_exp)] ;
    Mat_C_exp = [{'AW_C_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(C_AW_exp)] ;
    Mat_H_exp = [{'AW_H_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(H_AW_exp)] ;
    Mat_S_exp = [{'AW_S_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(S_AW_exp)] ;
    Mat_wR_exp = [{'AW_wR_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(wR_AW_exp)] ;
    %Mat_wR_at_exp = [{'AW_wR_at_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell((wR_at_AW_exp-1))] ;
    
    Mat_A_diff = [{'AW_A_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(100*(A_AW_exp./A_AW_ben-1))] ;
    Mat_C_diff = [{'AW_C_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(100*(C_AW_exp./C_AW_ben-1))] ;
    Mat_H_diff = [{'AW_H_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(100*(H_AW_exp./H_AW_ben-1))] ;
    Mat_S_diff = [{'AW_S_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(S_AW_exp-S_AW_ben)] ;
    Mat_wR_diff = [{'AW_wR_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(wR_AW_exp-wR_AW_ben)] ;
    %Mat_wR_at_diff = [{'AW_wR_at_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell((wR_at_AW_exp-wR_at_AW_ben))] ;
    
        Mat_Z1_ben = [{'AW_Z1_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z1_AW_ben)] ;
        Mat_Z2_ben = [{'AW_Z2_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z2_AW_ben)] ;
        Mat_Z3_ben = [{'AW_Z3_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z3_AW_ben)] ;
        Mat_Z4_ben = [{'AW_Z4_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z4_AW_ben)] ;
        Mat_Z5_ben = [{'AW_Z5_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z5_AW_ben)] ;
        Mat_Z6_ben = [{'AW_Z6_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z6_AW_ben)] ;
        Mat_Z7_ben = [{'AW_Z7_ben',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z7_AW_ben)] ;
        
        Mat_Z1_exp = [{'AW_Z1_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z1_AW_exp)] ;
        Mat_Z2_exp = [{'AW_Z2_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z2_AW_exp)] ;
        Mat_Z3_exp = [{'AW_Z3_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z3_AW_exp)] ;
        Mat_Z4_exp = [{'AW_Z4_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z4_AW_exp)] ;
        Mat_Z5_exp = [{'AW_Z5_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z5_AW_exp)] ;
        Mat_Z6_exp = [{'AW_Z6_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z6_AW_exp)] ;
        Mat_Z7_exp = [{'AW_Z7_exp',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z7_AW_exp)] ;
 
        Mat_Z1_diff = [{'AW_Z1_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z1_AW_exp-Z1_AW_ben)] ;
        Mat_Z2_diff = [{'AW_Z2_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z2_AW_exp-Z2_AW_ben)] ;
        Mat_Z3_diff = [{'AW_Z3_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z3_AW_exp-Z3_AW_ben)] ;
        Mat_Z4_diff = [{'AW_Z4_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z4_AW_exp-Z4_AW_ben)] ;
        Mat_Z5_diff = [{'AW_Z5_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z5_AW_exp-Z5_AW_ben)] ;
        Mat_Z6_diff = [{'AW_Z6_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z6_AW_exp-Z6_AW_ben)] ;
        Mat_Z7_diff = [{'AW_Z7_diff',' ',' ',' ',' ',' ',' ',' '};AW_title;age_title,num2cell(Z7_AW_exp-Z7_AW_ben)] ;

    Mat = [Mat_7;Mat_8 Mat_A_ben Mat_8 Mat_A_exp Mat_8 Mat_A_diff] ;
    status = xlwrite(AW_file,Mat,'Assets') ;
    
    Mat = [Mat_7;Mat_8 Mat_S_ben Mat_8 Mat_S_exp Mat_8 Mat_S_diff] ;
    status = xlwrite(AW_file,Mat,'Saving Rate') ;
    
    Mat = [Mat_7;Mat_8 Mat_C_ben Mat_8 Mat_C_exp Mat_8 Mat_C_diff] ;
    status = xlwrite(AW_file,Mat,'Cons') ;
    
    Mat = [Mat_7;Mat_8 Mat_H_ben Mat_8 Mat_H_exp Mat_8 Mat_H_diff] ;
    status = xlwrite(AW_file,Mat,'Hours') ;
    
    Mat = [Mat_7;Mat_8 Mat_wR_ben Mat_8 Mat_wR_exp Mat_8 Mat_wR_diff] ;
    status = xlwrite(AW_file,Mat,'Return') ;
    
%     Mat = [Mat_7;Mat_8 Mat_wR_at_ben Mat_8 Mat_wR_at_exp Mat_8 Mat_wR_at_diff] ;
%     status = xlwrite(AW_file,Mat,'Return - AT') ;
    
    Mat = [Mat_7;Mat_8 Mat_Z1_ben Mat_8 Mat_Z1_exp Mat_8 Mat_Z1_diff;
           Mat_7;Mat_8 Mat_Z2_ben Mat_8 Mat_Z2_exp Mat_8 Mat_Z2_diff;
           Mat_7;Mat_8 Mat_Z3_ben Mat_8 Mat_Z3_exp Mat_8 Mat_Z3_diff;
           Mat_7;Mat_8 Mat_Z4_ben Mat_8 Mat_Z4_exp Mat_8 Mat_Z4_diff;
           Mat_7;Mat_8 Mat_Z5_ben Mat_8 Mat_Z5_exp Mat_8 Mat_Z5_diff;
           Mat_7;Mat_8 Mat_Z6_ben Mat_8 Mat_Z6_exp Mat_8 Mat_Z6_diff;
           Mat_7;Mat_8 Mat_Z7_ben Mat_8 Mat_Z7_exp Mat_8 Mat_Z7_diff;] ;
    status = xlwrite(AW_file,Mat,'Z Dist.') ;
    
    for k=1:n_z
    for j=1:n_age
    for i=1:n_z
        eval( strcat('Z_prc',int2str(k),'_ben(j,i)=Z',int2str(i),'_AW_ben(j,k);') ) ;
        eval( strcat('Z_prc',int2str(k),'_exp(j,i)=Z',int2str(i),'_AW_exp(j,k);') ) ;
    end 
    end 
    end
    
        Mat_Z1_ben = [{'AW_Z1_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc1_ben)] ;
        Mat_Z2_ben = [{'AW_Z2_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc2_ben)] ;
        Mat_Z3_ben = [{'AW_Z3_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc3_ben)] ;
        Mat_Z4_ben = [{'AW_Z4_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc4_ben)] ;
        Mat_Z5_ben = [{'AW_Z5_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc5_ben)] ;
        Mat_Z6_ben = [{'AW_Z6_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc6_ben)] ;
        Mat_Z7_ben = [{'AW_Z7_ben',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc7_ben)] ;
        
        Mat_Z1_exp = [{'AW_Z1_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc1_exp)] ;
        Mat_Z2_exp = [{'AW_Z2_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc2_exp)] ;
        Mat_Z3_exp = [{'AW_Z3_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc3_exp)] ;
        Mat_Z4_exp = [{'AW_Z4_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc4_exp)] ;
        Mat_Z5_exp = [{'AW_Z5_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc5_exp)] ;
        Mat_Z6_exp = [{'AW_Z6_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc6_exp)] ;
        Mat_Z7_exp = [{'AW_Z7_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc7_exp)] ;
 
        Mat_Z1_diff = [{'AW_Z1_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc1_exp-Z_prc1_ben)] ;
        Mat_Z2_diff = [{'AW_Z2_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc2_exp-Z_prc2_ben)] ;
        Mat_Z3_diff = [{'AW_Z3_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc3_exp-Z_prc3_ben)] ;
        Mat_Z4_diff = [{'AW_Z4_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc4_exp-Z_prc4_ben)] ;
        Mat_Z5_diff = [{'AW_Z5_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc5_exp-Z_prc5_ben)] ;
        Mat_Z6_diff = [{'AW_Z6_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc6_exp-Z_prc6_ben)] ;
        Mat_Z7_diff = [{'AW_Z7_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Z_prc7_exp-Z_prc7_ben)] ;

    Mat = [Mat_7;Mat_8 Mat_Z1_ben Mat_8 Mat_Z1_exp Mat_8 Mat_Z1_diff;
           Mat_7;Mat_8 Mat_Z2_ben Mat_8 Mat_Z2_exp Mat_8 Mat_Z2_diff;
           Mat_7;Mat_8 Mat_Z3_ben Mat_8 Mat_Z3_exp Mat_8 Mat_Z3_diff;
           Mat_7;Mat_8 Mat_Z4_ben Mat_8 Mat_Z4_exp Mat_8 Mat_Z4_diff;
           Mat_7;Mat_8 Mat_Z5_ben Mat_8 Mat_Z5_exp Mat_8 Mat_Z5_diff;
           Mat_7;Mat_8 Mat_Z6_ben Mat_8 Mat_Z6_exp Mat_8 Mat_Z6_diff;
           Mat_7;Mat_8 Mat_Z7_ben Mat_8 Mat_Z7_exp Mat_8 Mat_Z7_diff;] ;
    status = xlwrite(AW_file,Mat,'Z Dist. by prc') ;

    
%% Wealth stats 

N = numel(Sim_A_ben) ;
% Wealth of people with more tna 10.000 units of assets 
    n = sum(Sim_A_ben>=5000) ;
    Av_Wealth = sum(Sim_A_ben.*(Sim_A_ben>=5000))/n                 ;
    Wealth_Share = sum(Sim_A_ben.*(Sim_A_ben>=5000))/sum(Sim_A_ben) ;
    col={'Wealth','Share of Pop','Av. Wealth','Wealth Share'} ;
    Mat=[col;num2cell([5000,100*n/N,Av_Wealth,100*Wealth_Share])] ;
    disp(Mat)

% Wealth of people with more tna 10.000 units of assets 
    n = sum(Sim_A_ben>=10000) ;
    Av_Wealth = sum(Sim_A_ben.*(Sim_A_ben>=10000))/n                 ;
    Wealth_Share = sum(Sim_A_ben.*(Sim_A_ben>=10000))/sum(Sim_A_ben) ;
    col={'Wealth','Share of Pop','Av. Wealth','Wealth Share'} ;
    Mat=[col;num2cell([10000,100*n/N,Av_Wealth,100*Wealth_Share])] ;
    disp(Mat)
    
% Wealth of people with more tna 50.000 units of assets 
    n = sum(Sim_A_ben>=50000) ;
    Av_Wealth = sum(Sim_A_ben.*(Sim_A_ben>=50000))/n                 ;
    Wealth_Share = sum(Sim_A_ben.*(Sim_A_ben>=50000))/sum(Sim_A_ben) ;
    col={'Wealth','Share of Pop','Av. Wealth','Wealth Share'} ;
    Mat=[col;num2cell([50000,100*n/N,Av_Wealth,100*Wealth_Share])] ;
    disp(Mat)
    
% Wealth of top ten wealth holders
    Wealth_top_10 = sort(Sim_A_ben);
    Wealth_top_10 = sum(Wealth_top_10(end-10:end)) ;
    Wealth_top_10_Share = 100* Wealth_top_10 / sum(Sim_A_ben) ; 
    
% Average Labor income
    Av_labor_income_no_ret = sum(Sim_Yh_ben.*Sim_age_ben<Ret_Age.*Sim_H_ben>0)/sum(Sim_age_ben<Ret_Age.*Sim_H_ben>0) ;
    
% Ratio of wealth of top ten earners to labor income:
    Ratio_top10Wealth_Labor_income = (Wealth_top_10/10)/Av_labor_income_no_ret ;
    
% Percentiles of the asset distribution and distribution of z by them
    for i=1:n_z 
        prc_totA_ben(i) = prctile(Sim_A_ben,prctl(i)) ;
        prc_totA_exp(i) = prctile(Sim_A_exp,prctl(i)) ;
    end 
    
    for i=1:n_z
        if i==1
           ind = (Sim_A_ben<=prc_totA_ben(i))  ;
        elseif i<n_z && i>1
           ind = ((Sim_A_ben<=prc_totA_ben(i)).*(Sim_A_ben>prc_totA_ben(i-1)))==1  ;
        else
           ind = (Sim_A_ben>prc_totA_ben(i))  ;
        end
        n=sum(ind);
        
        for k=1:n_z
            eval(strcat('Z_Aprc_ben(i,k) = 100*sum(Sim_Z_ben(ind)==',int2str(k),')/n;'))

            if eval(strcat('sum(Sim_Z_ben(ind)==',int2str(k),')'))>n
                beep
                [j,i,k]
                error('this is wrong')
            end 
        end
        
        if i==1
           ind = (Sim_A_exp<=prc_totA_exp(i))  ;
        elseif i<n_z && i>1
           ind = ((Sim_A_exp<=prc_totA_exp(i)).*(Sim_A_exp>prc_totA_exp(i-1)))==1  ;
        else
           ind = (Sim_A_exp>prc_totA_exp(i))  ;
        end
        n=sum(ind);
        
        for k=1:n_z
            eval(strcat('Z_Aprc_exp(i,k) = 100*sum(Sim_Z_exp(ind)==',int2str(k),')/n;'))

            if eval(strcat('sum(Sim_Z_exp(ind)==',int2str(k),')'))>n
                beep
                [j,i,k]
                error('this is wrong')
            end 
        end
    end    
    
    z_title = {' ','z1','z2','z3','z4','z5','z6','z7'} ;
    Mat_Z_ben = [{'Z_dist_by_A_prct_ben',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_ben)] ;
    Mat_Z_exp = [{'Z_dist_by_A_prct_exp',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_exp)] ;
    Mat_Z_diff = [{'Z_dist_by_A_prct_diff',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_exp-Z_Aprc_ben)] ;
    Mat = [Mat_7;Mat_8 Mat_Z_ben Mat_8 Mat_Z_exp Mat_8 Mat_Z_diff] ;
    status = xlwrite(AW_file,Mat,'Z Dist. Tot. A') ; 

%% Read files - Benchmark Economy
% A prime
    eval(['load ',Bench_Folder,'aprime'])
    Ap_bench = reshape(aprime,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear aprime
    
    % Savings rate
    S_bench = 100*(Ap_bench./A_mat-1) ; 

% Hours
    eval(['load ',Bench_Folder,'hours'])
    H_bench = reshape(hours,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear hours

% Consumption
    eval(['load ',Bench_Folder,'cons'])
    C_bench = reshape(cons,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear cons
    
% Value
    eval(['load ',Bench_Folder,'value'])
    V_bench = reshape(value,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear value
    
% Distribution
    eval(['load ',Bench_Folder,'DBN'])
    DBN_bench = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear DBN
    
% Wage
    eval(['load ',Bench_Folder,'wage'])
    W_bench = wage ;
    clear wage
    
% Interst Rate
    eval(['load ',Bench_Folder,'rr'])
    R_bench = rr ;
    clear rr
    
% Interst Rate
    eval(['load ',Bench_Folder,'EBAR'])
    E_bench = EBAR ;
    clear EBAR
    
    
%% Read files - Experimental Economy

% A prime
    eval(['load ',Result_Folder,'Exp_results_aprime']) ;
    Ap_exp = reshape(Exp_results_aprime,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_aprime
    
    % Savings rate
    S_exp = 100*(Ap_exp./A_mat-1) ; 

% Hours
    eval(['load ',Result_Folder,'Exp_results_hours']) ;
    H_exp = reshape(Exp_results_hours,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_hours

% Consumption
    eval(['load ',Result_Folder,'Exp_results_cons']) ;
    C_exp = reshape(Exp_results_cons,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_cons
    
% Value
    eval(['load ',Result_Folder,'Exp_results_value']) ;
    V_exp = reshape(Exp_results_value,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_value
    
% Distribution
    eval(['load ',Result_Folder,'Exp_results_DBN']) ;
    DBN_exp = reshape(Exp_results_DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_DBN

% Wage
    eval(['load ',Result_Folder,'Exp_results_wage'])
    W_exp = Exp_results_wage ;
    clear Exp_results_wage
    
% Interst Rate
    eval(['load ',Result_Folder,'Exp_results_rr'])
    R_exp = Exp_results_rr ;
    clear Exp_results_rr
    
% Wealth Taxes
    Threshold = Threshold_Factor*E_bench ; 
    
% Return
    Return_mat =  1 + mu * R_bench * Z_mat.^mu .* A_mat.^(mu-1) - delta  ;
    
a%% Moments 

    mean_A_dbn  = sum(sum(sum(sum(sum(A_mat.*DBN_bench)))))      ;
    mean_A_sim  = mean(Sim_A_ben)            ; 
    
    mean_C_dbn  = sum(sum(sum(sum(sum(C_bench.*DBN_bench)))))    ;
    mean_C_sim  = mean(Sim_C_ben)            ;
    
    mean_H_dbn  = sum(sum(sum(sum(sum(H_bench.*DBN_bench)))))    ;
    mean_H_sim  = mean(Sim_H_ben)            ;
    
    mean_Ap_dbn = sum(sum(sum(sum(sum(Ap_bench.*DBN_bench)))))   ;
    mean_Ap_sim = mean(Sim_Ap_ben)           ;
    
    mean_R_dbn  = sum(sum(sum(sum(sum(Return_mat.*DBN_bench))))) ;
    mean_R_sim  = mean(Sim_R_ben)            ;
    
    std_A_dbn   = ( sum(sum(sum(sum(sum( ((A_mat-mean_A_dbn).^2    ).*DBN_bench))))) ).^0.5 ;
    std_C_dbn   = ( sum(sum(sum(sum(sum( ((C_bench-mean_C_dbn).^2  ).*DBN_bench))))) ).^0.5 ;
    std_H_dbn   = ( sum(sum(sum(sum(sum( ((H_bench-mean_H_dbn).^2  ).*DBN_bench))))) ).^0.5 ;
    std_Ap_dbn  = ( sum(sum(sum(sum(sum( ((Ap_bench-mean_Ap_dbn).^2).*DBN_bench))))) ).^0.5 ;
    std_R_dbn   = ( sum(sum(sum(sum(sum( ((Return_mat-mean_R_dbn).^2  ).*DBN_bench))))) ).^0.5 ;
    
    std_A_sim   = std(Sim_A_ben)  ;
    std_C_sim   = std(Sim_C_ben)  ;
    std_H_sim   = std(Sim_H_ben)  ;
    std_Ap_sim  = std(Sim_Ap_ben) ;
    std_R_sim   = std(Sim_R_ben)  ;
    
    cols = {' ','mean_dbn','mean_sim',' ','std_dbn','std_sim'};
    rows = {'A';'R';'Ap';'C';'H'} ;
    Mat = [mean_A_dbn  mean_A_sim  NaN std_A_dbn  std_A_sim ;
           mean_R_dbn  mean_R_sim  NaN std_R_dbn  std_R_sim ;
           mean_Ap_dbn mean_Ap_sim NaN std_Ap_dbn std_Ap_sim;
           mean_C_dbn  mean_C_sim  NaN std_C_dbn  std_C_sim ;
           mean_H_dbn  mean_H_sim  NaN std_H_dbn  std_H_sim ];
    Mat = [cols;rows num2cell(Mat)]   
    
    
    N = numel(Sim_age_ben) ;
    for i=1:Max_Age 
        Share_age_dbn(i,1) = sum(sum(sum(sum(DBN_bench(i,:,:,:,:))))) ;
        Share_age_sim(i,1) = sum(Sim_age_ben==i)/N  ;
    end 
    
    for i=1:n_z
        Share_Z_dbn(i,1) = sum(sum(sum(sum(DBN_bench(:,:,i,:,:))))) ;
        Share_Z_sim(i,1) = sum(Sim_Z_ben==i)/N   ;
    end 
    %[(1:n_z)' 100*(Share_Z_dbn) 100*(Share_Z_sim) 100*(Share_Z_dbn-Share_Z_sim)]
    
    for i=1:n_a
        Share_A_dbn(i,1) = sum(sum(sum(sum(DBN_bench(:,i,:,:,:))))) ;
        if i<n_a
        Share_A_sim(i,1) = sum((Sim_A_ben>=agrid(i)).*(Sim_A_ben<agrid(i+1)))/N   ;
        else
        Share_A_sim(i,1) = sum(Sim_A_ben>=agrid(i))/N   ;
        end
    end
    %[(1:n_a)' 100*(Share_A_dbn) 100*(Share_A_sim) 100*(Share_A_dbn-Share_A_sim)]

    
    
%% Gini Coefficient with Burhan's data 
% I use Deaton's formulas for the gini. 
% The formulas are found in Deaton (1997) "The analysis of household surveys"
% The formulas are in page 139.



kubu_folder = '../../kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/more_stats/' ;

eval(['load ',kubu_folder,'panela_bench']) ;
    
N  = numel(panela_bench) ;
mu = mean(panela_bench)  ;
panela_bench_sort = sort(panela_bench,'descend') ;
index = 1:N ;

G_bench = (N+1)/(N-1) - 2*sum(panela_bench_sort.*index)/(mu*N*(N-1)) ;


eval(['load ',kubu_folder,'panela_exp']) ;
    
N  = numel(panela_exp) ;
mu = mean(panela_exp)  ;
panela_exp_sort = sort(panela_exp,'descend') ;
index = 1:N ;

G_exp = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu*N*(N-1)) ;


cols = {' ','Bench','Exp'};
rows = {'Gini'} ;
Mat = [cols; rows num2cell([G_bench G_exp])]

eval(['load ',kubu_folder,'panela_bench']) ;
eval(['load ',kubu_folder,'panel_return_bench']) ;
eval(['load ',kubu_folder,'panel_hours_bench']) ;
eval(['load ',kubu_folder,'panel_cons_bench']) ;
eval(['load ',kubu_folder,'panel_aprime_bench']) ;
eval(['load ',kubu_folder,'panelz_bench']) ;
eval(['load ',kubu_folder,'panelz_exp']) ;


    mean_A_dbn  = sum(sum(sum(sum(sum(A_mat.*DBN_bench)))))      ;
    mean_A_sim  = mean(panela_bench)            ; 
    
    mean_C_dbn  = sum(sum(sum(sum(sum(C_bench.*DBN_bench)))))    ;
    mean_C_sim  = mean(panel_cons_bench)            ;
    
    mean_H_dbn  = sum(sum(sum(sum(sum(H_bench.*DBN_bench)))))    ;
    mean_H_sim  = mean(panel_hours_bench)            ;
    
    mean_Ap_dbn = sum(sum(sum(sum(sum(Ap_bench.*DBN_bench)))))   ;
    mean_Ap_sim = mean(panel_aprime_bench)           ;
    
    mean_R_dbn  = sum(sum(sum(sum(sum(Return_mat.*DBN_bench))))) ;
    mean_R_sim  = mean(panel_return_bench)            ;
    
    std_A_dbn   = ( sum(sum(sum(sum(sum( ((A_mat-mean_A_dbn).^2    ).*DBN_bench))))) ).^0.5 ;
    std_C_dbn   = ( sum(sum(sum(sum(sum( ((C_bench-mean_C_dbn).^2  ).*DBN_bench))))) ).^0.5 ;
    std_H_dbn   = ( sum(sum(sum(sum(sum( ((H_bench-mean_H_dbn).^2  ).*DBN_bench))))) ).^0.5 ;
    std_Ap_dbn  = ( sum(sum(sum(sum(sum( ((Ap_bench-mean_Ap_dbn).^2).*DBN_bench))))) ).^0.5 ;
    std_R_dbn   = ( sum(sum(sum(sum(sum( ((Return_mat-mean_R_dbn).^2  ).*DBN_bench))))) ).^0.5 ;
    
    std_A_sim   = std(panela_bench)  ;
    std_C_sim   = std(panel_cons_bench)  ;
    std_H_sim   = std(panel_hours_bench)  ;
    std_Ap_sim  = std(panel_aprime_bench) ;
    std_R_sim   = std(panel_return_bench)  ;
    
    cols = {' ','mean_dbn','mean_sim',' ','std_dbn','std_sim'};
    rows = {'A';'R';'Ap';'C';'H'} ;
    Mat = [mean_A_dbn  mean_A_sim  NaN std_A_dbn  std_A_sim ;
           mean_R_dbn  mean_R_sim  NaN std_R_dbn  std_R_sim ;
           mean_Ap_dbn mean_Ap_sim NaN std_Ap_dbn std_Ap_sim;
           mean_C_dbn  mean_C_sim  NaN std_C_dbn  std_C_sim ;
           mean_H_dbn  mean_H_sim  NaN std_H_dbn  std_H_sim ];
    Mat = [cols;rows num2cell(Mat)]  

    
    
% Percentiles of the asset distribution and distribution of z by them
    for i=1:n_z 
        prc_totA_ben(i) = prctile(panela_bench,prctl(i)) ;
        prc_totA_exp(i) = prctile(panela_exp,prctl(i)) ;
    end 
    
    for i=1:n_z
        if i==1
           ind = (panela_bench<=prc_totA_ben(i))  ;
        elseif i<n_z && i>1
           ind = ((panela_bench<=prc_totA_ben(i)).*(panela_bench>prc_totA_ben(i-1)))==1  ;
        else
           ind = (panela_bench>prc_totA_ben(i))  ;
        end
        n=sum(ind);
        
        for k=1:n_z
            eval(strcat('Z_Aprc_ben(i,k) = 100*sum(panelz_bench(ind)==',int2str(k),')/n;'))

            if eval(strcat('sum(panelz_bench(ind)==',int2str(k),')'))>n
                beep
                [j,i,k]
                error('this is wrong')
            end 
        end
        
        if i==1
           ind = (panela_exp<=prc_totA_exp(i))  ;
        elseif i<n_z && i>1
           ind = ((panela_exp<=prc_totA_exp(i)).*(panela_exp>prc_totA_exp(i-1)))==1  ;
        else
           ind = (panela_exp>prc_totA_exp(i))  ;
        end
        n=sum(ind);
        
        for k=1:n_z
            eval(strcat('Z_Aprc_exp(i,k) = 100*sum(panelz_exp(ind)==',int2str(k),')/n;'))

            if eval(strcat('sum(panelz_exp(ind)==',int2str(k),')'))>n
                beep
                [j,i,k]
                error('this is wrong')
            end 
        end
    end   
    
    z_title = {' ','z1','z2','z3','z4','z5','z6','z7'} ;
    Mat_Z_ben = [{'Z_dist_by_A_prct_ben',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_ben)] ;
    Mat_Z_exp = [{'Z_dist_by_A_prct_exp',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_exp)] ;
    Mat_Z_diff = [{'Z_dist_by_A_prct_diff',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_exp-Z_Aprc_ben)] ;
    Mat = [Mat_7;Mat_8 Mat_Z_ben Mat_8 Mat_Z_exp Mat_8 Mat_Z_diff] ;
    status = xlwrite(AW_file,Mat,'Z Dist. Tot. A') ; 
    
    
    