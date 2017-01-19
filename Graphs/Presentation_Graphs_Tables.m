% Simulation Graphs and Tables
% Sergio Ocampo Diaz

%% Clear and load required files
   % Clear and close
        close all; clear all; clc

   % Load files for excel in mac
        javaaddpath('poi_library/poi-3.8-20120326.jar');
        javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
        javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
        javaaddpath('poi_library/xmlbeans-2.3.0.jar');
        javaaddpath('poi_library/dom4j-1.6.1.jar');
        javaaddpath('poi_library/stax-api-1.0.1.jar');
    
    % Line style
    set(groot,'defaultAxesLineStyleOrder','-|--|-.|:')

        
%% Fixed Parameters 

% Number of agents for plots
    N_plots = 20000000 ;

% Grids
    n_a = 201 ; 
    n_z = 9   ;
    n_l = 5   ;
    n_e = 5   ;
    
    Max_Age = 81 ;
    Ret_Age = 45 ;

% Age brackets 
    age = [5 , 15, 25, 35, 45, 55, Max_Age ];
    age_limit = [0, 5, 15, 25, 35, 45, 55, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;
    
    prc         = [10 25 50 60 75 80 90 95 99 99.9] ;

%% Fixed Parameters
    X_Switch  = 1.2;

    EBAR_data = 8.8891*10^12/(122.46*10^6) ; % 2013 total compensation of employees' devided by number of HH's in 2013

    theta_folder = 1.5 ;
    R_gap        = 0.00 ;
    Threshold_Factor = 0.0;
    
% Fixed Parameters
    if X_Switch==1
        n_x = 2   ;
        x_grid(1) = 5.0 ; x_grid(2) = 1 ;
        theta = 1 + (2.5-1)*(0:8)/8 ;
    elseif X_Switch==0
        n_x = 1   ;
        x_grid(1) = 1 ;
        theta = 1.5*ones(1,n_z) ;
    elseif X_Switch==2
        n_x = 2   ;
        x_grid(1) = 5.0 ; x_grid(2) = 1 ;
        theta = 1.5*ones(1,n_z) ;
    elseif X_Switch==1.1 || X_Switch==1.2
        n_x = 3   ;
        x_grid(1) = 5.0 ; x_grid(2) = 1 ; x_grid(3) = 0 ;
        theta = 1 + (2.5-1)*(0:8)/8 ;
    elseif X_Switch==2.1
        n_x = 3   ;
        x_grid(1) = 5.0 ; x_grid(2) = 1 ; x_grid(3) = 0 ;
        theta = 1.5*ones(1,n_z) ;
    elseif X_Switch==3.1
        n_x = 3   ;
        x_grid(1) = 1.0 ; x_grid(2) = 5 ; x_grid(3) = 0 ;
        theta = 1 + (10-1)*(0:8)/8 ;
    end
    
    mu = 0.9 ;
    Dep_Rate = 0.05 ;
    


% Utility and technology
    sigma = 4.00  ;
    gamma = 0.470;
    
% Taxes
	tauC    = 0.075         ;
    tauK    = 0.25          ;
	
% Set Tax Regime
    % tauPL = 0.185     ;
    % psi   = 0.77      ;
    tauPL = 0.0         ;
    psi   = 0.776       ;
    tauW  = 0.0126127644280288 ;
        
    
%% Result_Folders

if X_Switch==0
    mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents/'
    Result_Folder = strcat('../../../../NSU_F_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Markov_Cut/') ;
    Simul_Folder  = strcat('../../../../NSU_F_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Markov_Cut/Simul/') ;
    Bench_Folder  = strcat('../../../../NSU_F_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Markov_Cut/Bench_Files/') ;
    Exp_Folder  = strcat('../../../../NSU_F_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Markov_Cut/Exp_Files/') ;
    Top_Folder    = strcat('../../../../NSU_F_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Markov_Cut/Simul/Top_A/') ;
    mkdir('markov_cut/presentation')
    mkdir('markov_cut/presentation/png')
    cd 'markov_cut/presentation'
    Tables_file     = 'Tables_Presentation_Model_0.xls' ;
    
elseif X_Switch==1
    
        cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_zs/'
        Result_Folder = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_mu90/') ;
        Simul_Folder  = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_mu90/Simul/') ;
        Bench_Folder  = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_mu90/Bench_Files/') ;
        Exp_Folder  = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_mu90/Exp_Files/') ;
        Top_Folder    = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_mu90/Simul/Top_A/') ;
        mkdir('top_mu90/presentation')
        mkdir('top_mu90/presentation/png')
        cd 'top_mu90/presentation'
        Tables_file     = 'Tables_Presentation_Model_1.xls' ;
    
elseif X_Switch==2
    mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_F2/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_F2/'
        Result_Folder = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_rho1/') ;
        Simul_Folder  = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_rho1/Simul/') ;
        Bench_Folder  = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_rho1/Bench_Files/') ;
        Exp_Folder  = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_rho1/Exp_Files/') ;
        Top_Folder    = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Exp_Shock_top_rho1/Simul/Top_A/') ;
        mkdir('r_gap_0/presentation')
        mkdir('r_gap_0/presentation/png')
        cd 'r_gap_0/presentation'
        Tables_file     = 'Tables_Presentation_Model_2.xls' ;
    
elseif X_Switch==1.1
   mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_M11/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_M11/'
        Result_Folder = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock_theta25/') ;
        Simul_Folder  = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock_theta25/Simul/') ;
        Bench_Folder  = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock_theta25/Bench_Files/') ;
        Exp_Folder  = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock_theta25/Exp_Files/') ;
        Top_Folder    = strcat('../../../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock_theta25/Simul/Top_A/') ;
        mkdir('theta25/presentation')
        mkdir('theta25/presentation/png')
        cd 'theta25/presentation'
        Tables_file     = 'Tables_Presentation_Model_11.xls' ;
    
elseif X_Switch==1.2
   mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Model_1.2_bv/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Model_1.2_bv/'
        Result_Folder = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/') ;
        Simul_Folder  = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Simul/') ;
        Bench_Folder  = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Bench_Files/') ;
        Exp_Folder    = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Exp_Files/') ;
        Top_Folder    = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Simul/Top_A/') ;
        Tables_file   = 'Tables_Presentation_Model_12_bv.xls' ;
        
   mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Model_1.2_bv/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Model_1.2_bv/'
        Result_Folder   = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/') ;
        Simul_Folder    = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Simul/') ;
        Simul_Folder_Exp= strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Tau_C_Experiment/Simul/') ;
        Bench_Folder    = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Bench_Files/') ;
        Exp_Folder      = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Tau_C_Experiment/Exp_Files/') ;
        Top_Folder      = strcat('../../NSU_ZS_LT_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_1.2_bv/Simul/Top_A/') ;
        Tables_file     = 'Tables_Presentation_Model_12_bv_tauC.xls' ;
    
elseif X_Switch==2.1
    mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_F2/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_F2/'
        Result_Folder = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock/') ;
        Simul_Folder  = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock/Simul/') ;
        Bench_Folder  = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock/Bench_Files/') ;
        Exp_Folder  = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock/Exp_Files/') ;
        Top_Folder    = strcat('../../../../NSU_F2_LT_Results/R_gap_',num2str(R_gap,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Ret_Shock/Simul/Top_A/') ;
        mkdir('r_gap_0_ret_shock/presentation')
        mkdir('r_gap_0_ret_shock/presentation/png')
        cd 'r_gap_0_ret_shock/presentation'
        Tables_file     = 'Tables_Presentation_Model_21.xls' ;
elseif X_Switch==3.1
    mkdir('/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_M3/')
    cd '/Users/ocamp020/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Top_Agents_M3/'
        Result_Folder = strcat('../../../../NSU_M3_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_3_exp/') ;
        Simul_Folder  = strcat('../../../../NSU_M3_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_3_exp/Simul/') ;
        Bench_Folder  = strcat('../../../../NSU_M3_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_3_exp/Bench_Files/') ;
        Exp_Folder  = strcat('../../../../NSU_M3_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_3_exp/Exp_Files/') ;
        Top_Folder    = strcat('../../../../NSU_M3_Results/Theta_',num2str(theta_folder,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Model_3_exp/Simul/Top_A/') ;
        mkdir('exp_shock/presentation')
        mkdir('exp_shock/presentation/png')
        cd 'exp_shock/presentation'
        Tables_file     = 'Tables_Presentation_Model_13.xls' ;
    
end
      


%% Return percentiles table

    % Bench return 
        eval(['load ',Result_Folder,'zgrid']); eval(['load ',Bench_Folder,'P']); eval(['load ',Bench_Folder,'R']); eval(['load ',Bench_Folder,'EBAR'])
        EBAR_bench = EBAR*0.727853584919652;

        if X_Switch~=0
        eval(['load ',Simul_Folder,'panelx_bench'])      ; panel_x   = panelx_bench(1:N_plots)      ; clear panelx_bench       ;
        else
        panel_x = ones(size(panel_age)) ;
        end
        eval(['load ',Simul_Folder,'panelz_bench'])      ; panel_z   = panelz_bench(1:N_plots)      ; clear panelz_bench       ;
        eval(['load ',Simul_Folder,'panela_bench'])      ; panel_a   = panela_bench(1:N_plots)      ; clear panela_bench       ;    
        eval(['load ',Simul_Folder,'panelK_bench'])      ; panel_K   = panelK_bench(1:N_plots)      ; clear panelK_bench       ;
        eval(['load ',Simul_Folder,'panelage_bench'])    

        %  X_grid
        if X_Switch==1 || X_Switch==2
            % panel_xz  = exp(log(zgrid(panel_z)).*x_grid(panel_x)) ;
            panel_xz = NaN(size(panel_a)) ;
            panel_xz(panel_z<=4)  = exp(log(zgrid(panel_z(panel_z<=4))).*x_grid(1)) ;
            panel_xz(panel_z>4)   = exp(log(zgrid(panel_z(panel_z>4))).*x_grid(panel_x(panel_z>4))) ;
        elseif X_Switch==0
            panel_xz = zgrid(panel_z) ;
        elseif X_Switch==1.1 || X_Switch==1.2 || X_Switch==2.1
            panel_xz = NaN(size(panel_a)) ;
            panel_xz(panel_z<=4)  = exp(log(zgrid(panel_z(panel_z<=4))).*x_grid(1)) ;
            panel_xz(panel_z>4)   = exp(log(zgrid(panel_z(panel_z>4))).*x_grid(panel_x(panel_z>4))) ;
            pabe_xz(panel_x==3)   = 0 ;
        elseif X_Switch==3.1
            panel_xz = NaN(size(panel_a)) ;
            panel_xz(panel_x==1)  = zgrid(panel_z(panel_x==1)) ;
            panel_xz(panel_x==2)  = exp(log(zgrid(panel_z(panel_x==2))).*x_grid(2)) ;
            panel_xz(panel_x==3)  = 0 ;
        end

        % Profits
        if X_Switch==1 || X_Switch==0 || X_Switch==1.1 || X_Switch==1.2 || X_Switch==3.1
            panel_Pr = P*(panel_xz.*panel_K).^mu - (R+Dep_Rate)*panel_K ;
        elseif X_Switch==2 || X_Switch==2.1
            panel_Pr = P*(panel_xz.*panel_K).^mu - (R+R_gap+Dep_Rate)*max(panel_K-panel_a,0) - (R+Dep_Rate)*min(panel_a,panel_K);
        end        

        % Return 
            panel_Ret_K = (R*panel_a + panel_Pr)./panel_a ;

        clear panel_xz panel_K panel_a panel_z panel_x
        
        
    % Bench tax reform 
        eval(['load ',Result_Folder,'zgrid']); eval(['load ',Exp_Folder,'Exp_results_P']); eval(['load ',Exp_Folder,'Exp_results_R']); eval(['load ',Exp_Folder,'Exp_results_EBAR'])
        EBAR_bench = EBAR*0.727853584919652;
        P = Exp_results_P; R = Exp_results_R ;

        if X_Switch~=0
        eval(['load ',Simul_Folder_Exp,'panelx_exp'])      ; panel_x   = panelx_exp(1:N_plots)      ; clear panelx_exp       ;
        else
        panel_x = ones(size(panel_age)) ;
        end
        eval(['load ',Simul_Folder_Exp,'panelz_exp'])      ; panel_z   = panelz_exp(1:N_plots)      ; clear panelz_exp       ;
        eval(['load ',Simul_Folder_Exp,'panela_exp'])      ; panel_a   = panela_exp(1:N_plots)      ; clear panela_exp       ;    
        eval(['load ',Simul_Folder_Exp,'panelK_exp'])      ; panel_K   = panelK_exp(1:N_plots)      ; clear panelK_exp       ;

        %  X_grid
        if X_Switch==1 || X_Switch==2
            % panel_xz  = exp(log(zgrid(panel_z)).*x_grid(panel_x)) ;
            panel_xz = NaN(size(panel_a)) ;
            panel_xz(panel_z<=4)  = exp(log(zgrid(panel_z(panel_z<=4))).*x_grid(1)) ;
            panel_xz(panel_z>4)   = exp(log(zgrid(panel_z(panel_z>4))).*x_grid(panel_x(panel_z>4))) ;
        elseif X_Switch==0
            panel_xz = zgrid(panel_z) ;
        elseif X_Switch==1.1 || X_Switch==1.2 || X_Switch==2.1
            panel_xz = NaN(size(panel_a)) ;
            panel_xz(panel_z<=4)  = exp(log(zgrid(panel_z(panel_z<=4))).*x_grid(1)) ;
            panel_xz(panel_z>4)   = exp(log(zgrid(panel_z(panel_z>4))).*x_grid(panel_x(panel_z>4))) ;
            pabe_xz(panel_x==3)   = 0 ;
        elseif X_Switch==3.1
            panel_xz = NaN(size(panel_a)) ;
            panel_xz(panel_x==1)  = zgrid(panel_z(panel_x==1)) ;
            panel_xz(panel_x==2)  = exp(log(zgrid(panel_z(panel_x==2))).*x_grid(2)) ;
            panel_xz(panel_x==3)  = 0 ;
        end

        % Profits
        if X_Switch==1 || X_Switch==0 || X_Switch==1.1 || X_Switch==1.2 || X_Switch==3.1
            panel_Pr = P*(panel_xz.*panel_K).^mu - (R+Dep_Rate)*panel_K ;
        elseif X_Switch==2 || X_Switch==2.1
            panel_Pr = P*(panel_xz.*panel_K).^mu - (R+R_gap+Dep_Rate)*max(panel_K-panel_a,0) - (R+Dep_Rate)*min(panel_a,panel_K);
        end        

        % Return 
            panel_Ret_W    = 100*(R*panel_a + panel_Pr)./panel_a ;
            panel_Ret_at_W = 100*(((1+R)*panel_a+panel_Pr)*(1-tauW)./panel_a-1) ;
            

        clear panel_xz panel_K panel_a panel_z panel_x
        
        
    % Percentiles
        prc_Ret_K    = prctile(panel_Ret_K,[10,50,90,95,99]);
        prc_Ret_at_K = prc_Ret_K*(1-tauK) ;
        prc_Ret_W    = prctile(panel_Ret_W,[10,50,90,95,99]);
        prc_Ret_at_W = prctile(panel_Ret_at_W,[10,50,90,95,99]);
        
    % Table 
    col_title = {'Return','p10','p50','p90','p95','p99'};
    row_title = {'Benchmark';'Tax Reform';'Benchmark AT';'Tax Reform AT'};
    Mat = [col_title; row_title num2cell([prc_Ret_K;prc_Ret_W;prc_Ret_at_K;prc_Ret_at_W])]
    status = xlwrite(Tables_file,Mat,'Return_Percentile') ;
    
    
    % Table by agregroups with extra percentiles - Benchmark
    Mat = NaN(8,7);
    Mat(1,:) = prctile(panel_Ret_K,[10 25 50 75 90 95 99]);    
    for age_group_counter=1:7
        indx1 = ((panelage_bench>age_limit(age_group_counter)).*(panelage_bench<=age_limit(age_group_counter+1)));
        indx1 = indx1>0;
        Mat(1+age_group_counter,:) = prctile(panel_Ret_K(indx1),[10 25 50 75 90 95 99]);
        clear indx1
    end
    
    col_title = {'age_group','p10','p25','p50','p75','p90','p95','p99'};
    row_title = {'Pop';'Age1';'Age2';'Age3';'Age4';'Age5';'Age6';'Age7'};
    Mat = [col_title;row_title num2cell(Mat)]
    status = xlwrite(Tables_file,Mat,'Ret_Prc_age_group') ;
    
%% Change in composition of top X% by Z
    Switch_PV = 0;

    % Load PV wealth and Z
        eval(['load ',Simul_Folder,'panelz_bench']) ; eval(['load ',Simul_Folder,'panelx_bench']) ; 
        eval(['load ',Simul_Folder_Exp,'panelz_exp'])   ; eval(['load ',Simul_Folder_Exp,'panelx_exp'])   ;
        if Switch_PV==1
        eval(['load ',Simul_Folder,'panelPV_a_bench']) ; eval(['load ',Simul_Folder_Exp,'panelPV_a_exp'])   ;
        Wealth_K = panelPV_a_bench; Wealth_W = panelPV_a_exp;
        else
        eval(['load ',Simul_Folder,'panela_bench']) ; eval(['load ',Simul_Folder_Exp,'panela_exp'])   ;
        Wealth_K = panela_bench; Wealth_W = panela_exp;
        end
        
	% Get Percentiles
        prc_K = prctile(Wealth_K,[50 90 95 99]);
        prc_W = prctile(Wealth_W,[50 90 95 99]);
    % Get composition by top x%
        ii = 1;
        for i=numel(prc_K):-1:1
            ind_K = Wealth_K>=prc_K(i) ;
            ind_W = Wealth_W>=prc_W(i) ;
            for z=1:n_z
                z_share_top_x_K(ii,z) = 100*sum(panelz_bench(ind_K)==z)/sum(ind_K) ;
                z_share_top_x_W(ii,z) = 100*sum(panelz_exp(ind_W)==z)/sum(ind_W) ;
            end
            ii= ii+1;
        end 
	% Tables
        col_title = {'Top X%','z1','z2','z3','z4','z5','z6','z7','z8','z9'};
        row_title = {'Top 1%';'Top 5%';'Top 10%';'Top 50%'};
        Mat = [{'Benchmark'} cell(1,n_z); col_title; row_title num2cell(z_share_top_x_K);
               {'Tax_Reform'} cell(1,n_z); col_title; row_title num2cell(z_share_top_x_W);
               {'Difference'} cell(1,n_z); col_title; row_title num2cell(z_share_top_x_W-z_share_top_x_K);
               {'% Difference'} cell(1,n_z); col_title; row_title num2cell(100*(z_share_top_x_W./z_share_top_x_K-1))]
        if Switch_PV==1 
            status = xlwrite(Tables_file,Mat,'Z shares in Top x - Present Value') ;
        else
            status = xlwrite(Tables_file,Mat,'Z shares in Top x - Book Value') ;
        end 
        
   % Composition by xz
        ii = 1;
        for i=numel(prc_K):-1:1
            ind_K = Wealth_K>=prc_K(i) ;
            ind_W = Wealth_W>=prc_W(i) ;
            xz = 1 ;
            for z=1:n_z
                for x=1:n_x
                xz_share_top_x_K(ii,xz) = 100*sum(panelz_bench(ind_K)==z & panelx_bench(ind_K)==x)/sum(ind_K) ;
                xz_share_top_x_W(ii,xz) = 100*sum(panelz_exp(ind_W)==z & panelx_exp(ind_W)==x)/sum(ind_W) ;
                xz = xz+1;
                end
            end
            ii= ii+1;
        end 
        
    % Tables
        col_title = {'Top X%','z1 x1','z1 x2','z1 x3','z2 x1','z2 x2','z2 x3','z3 x1','z3 x2','z3 x3',...
                              'z4 x1','z4 x2','z4 x3','z5 x1','z5 x2','z5 x3','z6 x1','z6 x2','z6 x3',...
                              'z7 x1','z7 x2','z7 x3','z8 x1','z8 x2','z8 x3','z9 x1','z9 x2','z9 x3'};
        row_title = {'Top 1%';'Top 5%';'Top 10%';'Top 50%'};
        Mat = [{'Benchmark'} cell(1,n_z*n_x); col_title; row_title num2cell(xz_share_top_x_K);
               {'Tax_Reform'} cell(1,n_z*n_x); col_title; row_title num2cell(xz_share_top_x_W);
               {'Difference'} cell(1,n_z*n_x); col_title; row_title num2cell(xz_share_top_x_W-xz_share_top_x_K);
               {'% Difference'} cell(1,n_z*n_x); col_title; row_title num2cell(100*(xz_share_top_x_W./xz_share_top_x_K-1))]
        if Switch_PV==1 
            status = xlwrite(Tables_file,Mat,'XZ shares in Top x - Present Value') ;
        else
            status = xlwrite(Tables_file,Mat,'XZ shares in Top x - Book Value') ;
        end 

%% Tax Reform welfare by age group and productivity

    % Load Distribution
        eval(['load ',Bench_Folder,'DBN'])
        DBN = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e,n_x]) ;

    % Consumption equivalent welfare CE1
        eval(['load ',Bench_Folder,'value']); ValueFunction_bench = reshape(value,[Max_Age,n_a,n_z,n_l,n_e,n_x]) ; clear value
        eval(['load ',Exp_Folder,'Exp_results_value']); ValueFunction_exp = reshape(Exp_results_value,[Max_Age,n_a,n_z,n_l,n_e,n_x]) ; clear Exp_results_value
    Cons_Eq_Welfare=(ValueFunction_exp./ValueFunction_bench).^ ( 1 / ( gamma* (1-sigma)) )-1 ;
    
                                 
    CE_by_agegroup_z  = NaN(n_age,n_z);
    CE_by_agegroup_xz = NaN(n_age,n_z*n_x);
    for age_group_counter=1:n_age
        ii=1 ;
    for zi=1:n_z
        CE_by_agegroup_z(age_group_counter,zi)= ...
             100*sum(sum(sum(sum(sum(Cons_Eq_Welfare(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:).* ...
                    DBN(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:))))))/ ...
                 sum(sum(sum(sum(sum( DBN(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,:)))))) ;
        for xi=1:n_x
            CE_by_agegroup_xz(age_group_counter,ii) = ...
             100*sum(sum(sum(sum(sum(Cons_Eq_Welfare(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,xi).* ...
                    DBN(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,xi))))))/ ...
                 sum(sum(sum(sum(sum( DBN(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,zi,:,:,xi)))))) ;
            ii=ii+1;
        end 
    end
    end
    col_title = {'Age','z1','z2','z3','z4','z5','z6','z7','z8','z9'};
    row_title = {'Age1';'Age2';'Age3';'Age4';'Age5';'Age6';'Age7'};
    Mat = [col_title;row_title num2cell(CE_by_agegroup_z)]
    status = xlwrite(Tables_file,Mat,'CE by Z') ;
    col_title = {'Top X%','z1 x1','z1 x2','z1 x3','z2 x1','z2 x2','z2 x3','z3 x1','z3 x2','z3 x3',...
                              'z4 x1','z4 x2','z4 x3','z5 x1','z5 x2','z5 x3','z6 x1','z6 x2','z6 x3',...
                              'z7 x1','z7 x2','z7 x3','z8 x1','z8 x2','z8 x3','z9 x1','z9 x2','z9 x3'};
    row_title = {'Age1';'Age2';'Age3';'Age4';'Age5';'Age6';'Age7'};
    Mat = [col_title;row_title num2cell(CE_by_agegroup_xz)]
    status = xlwrite(Tables_file,CE_by_agegroup_xz,'CE by XZ') ;
    
    CE_mat = NaN(n_a,n_z,5) ;
    for age=1:5
        for z=1:n_z
            for i=1:n_a
                CE_mat(i,z,age) = 100*sum(sum(sum(Cons_Eq_Welfare(age,i,z,:,:,:).*DBN(age,i,z,:,:,:))))/sum(sum(sum(DBN(age,i,z,:,:,:)))) ; 
            end
        end
    end
    Mat=[];
    for i=1:5
    Mat = [Mat NaN(n_a,1) CE_mat(:,:,i)];
    end 
    status = xlwrite(Tables_file,Mat,'CE by asset') ;
    

    
%% Mean wealth by age

    % Load Distribution
        eval(['load ',Bench_Folder,'DBN'])
        DBN = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e,n_x]) ;
        
    % A grid
        eval(['load ',Result_Folder,'agrid']);
        A_mat = repmat(agrid,[Max_Age,1,n_z,n_l,n_e,n_x]);
        
    % Mean wealth by age
    Wealth_by_age= NaN(Max_Age,n_z) ;
    
    for i=1:Max_Age
        for z=1:n_z
        Wealth_by_age(i,z) = sum(sum(sum(sum(sum(A_mat(i,:,z,:,:,:).*DBN(i,:,z,:,:,:))))))/sum(sum(sum(sum(sum(DBN(i,:,z,:,:,:))))));
        end
    end 
    col_title = {'Age','z1','z2','z3','z4','z5','z6','z7','z8','z9'};
    row_title = num2cell([20:100]');
    Mat= [col_title;row_title num2cell(Wealth_by_age)]
    status = xlwrite(Tables_file,Mat,'Asset by age-Z') ;
    
    

        
%% Gini Coefficient with Burhan's data 
% I use Deaton's formulas for the gini. 
% The formulas are found in Deaton (1997) "The analysis of household surveys"
% The formulas are in page 139.

    % Benchmark
        eval(['load ',Simul_Folder,'panela_bench']) ; 
        N  = numel(panela_bench) ;
        mu_g = mean(panela_bench)  ;
        panela_bench_sort = sort(panela_bench,'descend') ;
        index = 1:N ;
        G_bench = (N+1)/(N-1) - 2*sum(panela_bench_sort.*index)/(mu_g*N*(N-1)) ;
    % Experiment 
        eval(['load ',Simul_Folder_Exp,'panela_exp']) ; 
        N  = numel(panela_exp) ;
        mu_g = mean(panela_exp)  ;
        panela_exp_sort = sort(panela_exp,'descend') ;
        index = 1:N ;
        G_exp = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu_g*N*(N-1)) ;
    % Optimal Tau K
        eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panela_exp']) ; 
        N  = numel(panela_exp) ;
        mu_g = mean(panela_exp)  ;
        panela_exp_sort = sort(panela_exp,'descend') ;
        index = 1:N ;
        G_opt_K = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu_g*N*(N-1)) ;
    % Optimal Tau W
        eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panela_exp']) ; 
        N  = numel(panela_exp) ;
        mu_g = mean(panela_exp)  ;
        panela_exp_sort = sort(panela_exp,'descend') ;
        index = 1:N ;
        G_opt_W = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu_g*N*(N-1)) ;

    % Save results
        cols = {' ','Bench','Exp','Opt_tau_K','Opt_tau_W'};
        rows = {'Gini'} ;
        Mat = [cols; rows num2cell([G_bench G_exp G_opt_K G_opt_W])]

        status = xlwrite(Tables_file,Mat,'Gini') ;

       
%% Pareto Tail

    % Vermuelen Data
    load /Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Graphs/Vermeulen_Data/USA_IM1.dat
    V_Data = USA_IM1;  clear USA_IM1;
    V_Data = V_Data(:,2:4) ;
    [panel_a_V,ind_V] = sort(V_Data(:,2)); weights_a_V = V_Data(ind_V,3); type_a_V = V_Data(ind_V,1);
    panel_aux = panel_a_V; type_aux = type_a_V ;
    i = 1;
    while numel(panel_aux)>=1
        panel_V(i)   = panel_aux(1)                              ;
        weights_V(i) = sum(weights_a_V(panel_a_V==panel_aux(1))) ;
        type_V(i)    = type_aux(1)                               ;
        type_aux     = type_aux(panel_aux~=panel_aux(1))         ;
        panel_aux    = panel_aux(panel_aux~=panel_aux(1))        ;
        i=i+1;
    end 
        % Select only observations with positive wealth
        weights_V = weights_V(panel_V>0); type_V = type_V(panel_V>0);
        panel_V = 1.3572*panel_V(panel_V>0); 
    
    % Bench Files
        eval(['load ',Bench_Folder,'EBAR'])         ; EBAR_bench = EBAR*0.727853584919652 ;
        
            eval(['load ',Simul_Folder,'panela_bench']) ; 
            wealth_bench = EBAR_data/EBAR_bench * sort(panela_bench) ;
        
    % Pareto Tail
        w_min = [1 1000000];
        for ii=1:numel(w_min)
            % Bench
                %x_bench   = -log(wealth_bench(wealth_bench>w_min(ii)))     ;
                %y_bench   =  log(((numel(x_bench)):-1:1)/(numel(x_bench))) ;
                [y_bench,x_bench] = ecdf(wealth_bench(wealth_bench>w_min(ii))) ;
                x_bench = -log(x_bench(2:end)) ; y_bench = log(1-y_bench(1:end-1)) ;
                %x_bench = log(wealth_bench(wealth_bench>w_min(ii)));
                %y_bench = log(1-[0:numel(x_bench)-1]/numel(x_bench));
            % Vermuelen
                x_V_s     = -log(panel_V((panel_V>w_min(ii))&(type_V==1))) ;
                y_V_s     = log(cumsum(weights_V((panel_V>w_min(ii))&(type_V==1)),'reverse')/sum(weights_V(panel_V>w_min(ii)))) ;
                x_V_f     = -log(panel_V((panel_V>w_min(ii))&(type_V==0))) ;
                y_V_f     = log(cumsum(weights_V((panel_V>w_min(ii))&(type_V==0)),'reverse')/sum(weights_V(panel_V>w_min(ii)))) ;
                x_V       = -log(panel_V(panel_V>w_min(ii))) ;
                y_V       = log(cumsum(weights_V((panel_V>w_min(ii))),'reverse')/sum(weights_V(panel_V>w_min(ii)))) ;
            % Maximum Likelihood
            a_ml(ii)  = sum(weights_V(panel_V>w_min(ii))) / sum(weights_V(panel_V>w_min(ii)).*log(panel_V(panel_V>w_min(ii))/w_min(ii))) ;
            % Regression
            mdl_b     = fitlm(x_bench',y_bench')         ;
            C_b       = mdl_b.Coefficients.Estimate(1)   ;
            a_b(ii)   = mdl_b.Coefficients.Estimate(2)   ; 
            mdl_V_s   = fitlm(x_V_s',y_V_s')             ;      
            C_V_s     = mdl_V_s.Coefficients.Estimate(1) ;
            a_V_s(ii) = mdl_V_s.Coefficients.Estimate(2) ;
            mdl_V     = fitlm(x_V',y_V')                 ;
            C_V       = mdl_V.Coefficients.Estimate(1)   ;
            a_V(ii)   = mdl_V.Coefficients.Estimate(2)   ;
            xx        = linspace(min(-x_V),max(-x_V),3)  ;
            % Select x_bench and y_bench
            x_bench_2 = [x_bench(1:200:end-600) ; x_bench(end-599:10:end-50) ; x_bench(end-50:end)];
            y_bench_2 = [y_bench(1:200:end-600) ; y_bench(end-599:10:end-50) ; y_bench(end-50:end)];
            % Figure
            fig_title = ['Pareto Tail Above $',num2str(w_min(ii))] ;
            figure; hold on; 
            % scatter(-x_V_s,y_V_s); scatter(-x_V_f,y_V_f,'*'); 
%             plot(-x_V,y_V,'o','color',[0    0.4470    0.7410]);              plot(xx,C_V-a_V(ii)*xx,'color',[0    0.4470    0.7410],'linewidth',1.5); 
%             plot(-x_bench_2,y_bench_2,'d','color',[0.8500    0.3250    0.0980]); plot(xx,C_b-a_b(ii)*xx,'color',[0.8500    0.3250    0.0980],'linewidth',1.5); hold off;
            plot(-x_V,y_V,'or');             plot(xx,C_V-a_V(ii)*xx,'r','linewidth',1.5); 
            plot(-x_bench_2,y_bench_2,'db'); plot(xx,C_b-a_b(ii)*xx,'b','linewidth',1.5); hold off;
            title(fig_title); xlabel('Wealth (log scale)'); ylabel('Log Counter-CDF'); xlim([log(w_min(ii)),max(-x_V_f)]);
            if w_min(ii)~=1000000
                ww = linspace(log(w_min(ii)),log(panel_V(end)),5) ;
            else
                ww = log([1e6 10e6 100e6 1e9 10e9 50e9]);
            end 
            set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)); ylim([-16 0])
            legend('US Data','Regression Line','Model','Regression Line','location','NorthEast')
            file_name_eps = ['Pareto_Tail_$',num2str(w_min(ii)),'_Presentation.eps'] ;
            file_name_png = ['png/Pareto_Tail_$',num2str(w_min(ii)),'_Presentation.png'] ;
            file_name_fig = ['fig/Pareto_Tail_$',num2str(w_min(ii)),'_Presentation.fig'] ;
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);
            % Table
            AA(:,ii) = [w_min(ii) ; a_b(ii) ; a_V_s(ii) ; a_V(ii) ; a_ml(ii) ] ; 
        end
        
        bench_name = {'Bench',' ',' '} ;
        blank_row  = cell(1,3) ;
        
        row_name = {'w_min';'alpha_bench';'alpha_survey';'alpha_survey_forbes';'alpha_ml'} ;
        Mat = [blank_row ; row_name num2cell(AA)]
        status = xlwrite(Tables_file,Mat,'Pareto_Tail_PV') ;
        
    
clear panela_bench panelPV_a_bench wealth_bench y_bench x_bench
        

%% Optimal Tax Graphs

% Load stats for tau K tau W
    eval(['load ',Result_Folder,'Opt_Tax_K/','Stats_by_tau_k.txt']) ; 
    eval(['load ',Result_Folder,'Opt_Tax_W/','Stats_by_tau_w.txt']) ; 
    % ORDER OF THINGS
        % 1. tauk, 2. tauw, 3. taul, 4. GBAR_K/(GBAR+SSC), 5. KBAR, 6. QBAR, 7. NBAR, 8. YBAR
        % 9. %Y , 10. wage, 11. VBAR, 12. CE2_NB, 13. CE2, 14. KBAR/Y, 
        % 15. Top1%, 16. Top10%, 17. std(log(E)), 18. Av. Hours
    aa = Stats_by_tau_k(:,4) ;
    aa = aa*0.16316577543219390/(0.16316577543219390+0.12869127919744788);
    Stats_by_tau_k(:,4) = aa ;
    aa = Stats_by_tau_w(:,4) ;
    aa = aa*0.16316577543219390/(0.16316577543219390+0.12869127919744788);
    Stats_by_tau_w(:,4) = aa ;
    
  
    G_K_Frac_k = Stats_by_tau_k(:,4) ;
    G_K_Frac_w = Stats_by_tau_w(:,4) ;
    

    % KBAR & QBAR graph
        figure;
        hold on
        plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,5)/Stats_by_tau_k(41,5)-1),'r', 'linewidth',2)
        plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,5)/Stats_by_tau_w(41,5)-1),'b', 'linewidth',2)
        xlabel('Tax Revenue from K / Total Tax Revenue')
        ylabel('Percent Change')
        h = legend('$\bar k, \tau_k$' , '$\bar k, \tau_a$','Location','SouthWest');
        set(h,'Interpreter','latex','FontSize',20)
        grid on
        axis([min(G_K_Frac_k) max(G_K_Frac_k)  -40 40 ])

        hgsave('fig/1.1fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
        print -depsc 1.1fig_KBAR_QBAR_by_CAP_TAX_REV.eps

        plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,6)/Stats_by_tau_k(41,6)-1),'-dr', 'linewidth',1.5)
        plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,6)/Stats_by_tau_w(41,6)-1),'-db', 'linewidth',1.5)
        h2 = legend('$\bar k, \tau_k$' , '$\bar k, \tau_a$','$\bar Q, \tau_k$' , '$\bar Q, \tau_a$','Location','SouthWest');
        set(h2,'Interpreter','latex','FontSize',20)

        hgsave('fig/1.2fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
        print -depsc 1.2fig_KBAR_QBAR_by_CAP_TAX_REV.eps

CE2_NB_k =100*( (Stats_by_tau_k(:,11)./Stats_by_tau_k(41+25,11)).^(1/(gamma*(1-sigma)))-1);
CE2_NB_w =100*( (Stats_by_tau_w(:,11)./Stats_by_tau_k(41+25,11)).^(1/(gamma*(1-sigma)))-1);
% Burhan's graphs on welfare by tau
    [x,tauindx ] = min(abs(Stats_by_tau_k(:,1) - 0.25)) ;

    figure
    axes1 = axes(...
      'FontName','Times New Roman',...
      'FontSize',14);
    hold on
    plot(G_K_Frac_k, (CE2_NB_k-CE2_NB_k(tauindx)),'r', 'linewidth',2)
    plot(G_K_Frac_w, zeros(size(G_K_Frac_k)),'k--', 'linewidth',2)
    grid on
    xlabel('Tax Revenue from K / Total Tax Revenue')
    ylabel('$CE_2$ Welfare Change from \textbf{Benchmark}','Interpreter','latex')
    plot(G_K_Frac_k(tauindx)*ones(2),[min(CE2_NB_k-CE2_NB_k(tauindx))-1 0],'k--', 'linewidth',2)
    axis([min(G_K_Frac_k) max(G_K_Frac_w)  min(CE2_NB_k-CE2_NB_k(tauindx))-1 max(CE2_NB_w-CE2_NB_k(tauindx))+1 ])
    
    annotation('textarrow',[0.52 0.57],[0.54 0.56], 'string','Cap. Income Tax Economy','linewidth',2,'FontSize',12)
    annotation('textarrow',[0.63 0.68],[0.36 0.41], 'string','Benchmark, \tau_k = 25%','linewidth',2,'FontSize',12)
    annotation('textarrow',[0.51+0.127 0.56+0.127],[0.17 0.12], 'string','0.25','linewidth',2,'FontSize',12)


    hgsave('fig/1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -depsc 1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps

%     plot(-0.3416*ones(2),[min(CE2_NB_k-CE2_NB_k(tauindx))-1 8],'k--', 'linewidth',2)
    annotation('textarrow',[0.21 0.18],[0.76 0.72], 'string','Opt. \tau_k = -34.4%','linewidth',2,'FontSize',12)
    hgsave('fig/1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -depsc 1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps

    plot(G_K_Frac_w, (CE2_NB_w-CE2_NB_k(tauindx)),'b', 'linewidth',2)
%     plot(0.4107*ones(2),[min(CE2_NB_k-CE2_NB_k(tauindx))-1 10],'k--', 'linewidth',2)
    annotation('textarrow',[0.80 0.83],[0.82 0.87], 'string','Opt. \tau_a = 3.06%','linewidth',2,'FontSize',12)


    hgsave('fig/1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -depsc 1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps
    
    
    
    return 




%% ==============================================================================
    
%% Composition of Z by welth percentile.
    % Each table shows the composition of the top x% into Z categories

    eval(['load ',Simul_Folder,'panela_bench']) ;    eval(['load ',Simul_Folder,'panelz_bench']);
    eval(['load ',Simul_Folder_Exp,'panela_exp'])   ;    eval(['load ',Simul_Folder_Exp,'panelz_exp'])  ;
    
    table_bench = NaN(10,n_z+1) ; table_exp = table_bench ;
       
    x=1;
    for i=100-[1, 5, 10, 25,50,75, 90,95,99, 100]
       tempa = panela_bench(panela_bench>=prctile(panela_bench,i));
       tempz = panelz_bench(panela_bench>=prctile(panela_bench,i));
       
      table_bench(x,:)=[100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
              sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
              sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
              sum(tempz==7)/sum(tempz>-1) ] ;   
          x=x+1;
    end

    x=1;
    for i=100-[1, 5, 10, 25,50,75, 90,95,99, 100]
       tempa=panela_exp(panela_exp>=prctile(panela_exp,i));
       tempz=panelz_exp(panela_exp>=prctile(panela_exp,i));
          table_exp(x,:)=[100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
              sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
              sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
              sum(tempz==7)/sum(tempz>-1) ] ;   
          x=x+1;

    end

    
    col_title = {'prct','z1','z2','z3','z4','z5','z6','z7',' '} ;
    Mat_1     = cell(10,1) ;
    
    Mat_A = [{'Benchmark',' ',' ',' ',' ',' ',' ',' ',' '},{'Experiment',' ',' ',' ',' ',' ',' ',' ',' '},{'Difference',' ',' ',' ',' ',' ',' ',' ',' '};
              col_title,col_title,col_title;
              num2cell(100*table_bench),Mat_1,num2cell(100*table_exp),Mat_1,num2cell(100*table_exp-100*table_bench),Mat_1] ;
	status = xlwrite(Tables_file,Mat_A,'Reallocation_Wealth_all') ;
          
    disp(Mat_A)
    
% Small table
    table_bench = NaN(4,n_z+1) ; table_exp = table_bench ;
    x=1;
    for i=100-[1, 5, 10, 50]
       tempa = panela_bench(panela_bench>=prctile(panela_bench,i));
       tempz = panelz_bench(panela_bench>=prctile(panela_bench,i));
       
      table_bench(x,:)=[100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
              sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
              sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
              sum(tempz==7)/sum(tempz>-1) ] ;   
          x=x+1;
    end

    x=1;
    for i=100-[1, 5, 10, 50]
       tempa=panela_exp(panela_exp>=prctile(panela_exp,i));
       tempz=panelz_exp(panela_exp>=prctile(panela_exp,i));
          table_exp(x,:)=[100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
              sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
              sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
              sum(tempz==7)/sum(tempz>-1) ] ;   
          x=x+1;
    end

    
    col_title = {'prct','z1','z2','z3','z4','z5','z6','z7',' '} ;
    Mat_1     = cell(4,1) ;
    
    Mat_A = [{'Benchmark',' ',' ',' ',' ',' ',' ',' ',' '},{'Experiment',' ',' ',' ',' ',' ',' ',' ',' '},{'Difference',' ',' ',' ',' ',' ',' ',' ',' '};
              col_title,col_title,col_title;
              num2cell(100*table_bench),Mat_1,num2cell(100*table_exp),Mat_1,num2cell(100*table_exp-100*table_bench),Mat_1] ;
	status = xlwrite(Tables_file,Mat_A,'Reallocation_Wealth') ;
          
    disp(Mat_A)
    
%% Opt Tax 
    % Compute consumption and hours increase for different percentiles
    
    eval(['load ',Simul_Folder,'panel_cons_bench']) ;    eval(['load ',Simul_Folder,'panel_hours_bench']);
    eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panel_cons_exp'])   ;    eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panel_hours_exp'])  ;
    %eval(['load ',Simul_Folder,'panel_cons_exp'])   ;    eval(['load ',Simul_Folder,'panel_hours_exp'])  ;
    
    panel_cons_exp_K  = panel_cons_exp  ;
    panel_hours_exp_K = panel_hours_exp ;
    
    eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panel_cons_exp'])   ;    eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panel_hours_exp'])  ;
    %eval(['load ',Simul_Folder,'panel_cons_exp'])   ;    eval(['load ',Simul_Folder,'panel_hours_exp'])  ;
    
    panel_cons_exp_W  = panel_cons_exp  ;
    panel_hours_exp_W = panel_hours_exp ;

    % Percentile of consumption in benchmark, opt tax k, opt tax w and
    % growth in mean consumption below the prct.
    Cons_prc_change = NaN(9,6) ;
    Hours_prc_change = NaN(9,6) ;
    Hours_by_cons_prc_change = NaN(9,6) ;
    x = 1 ;
    for i=[1, 10, 25,50,75, 90,95,99,100]
      Cons_prc_change(x,:) = [i , prctile(panel_cons_bench,i) , prctile(panel_cons_exp_K,i) , prctile(panel_cons_exp_W,i) , ...
          100*(mean(panel_cons_exp_K(panel_cons_exp_K<=prctile(panel_cons_exp_K,i)))/ mean(panel_cons_bench(panel_cons_bench<=prctile(panel_cons_bench,i)))-1) ... 
          100*(mean(panel_cons_exp_W(panel_cons_exp_W<=prctile(panel_cons_exp_W,i)))/ mean(panel_cons_bench(panel_cons_bench<=prctile(panel_cons_bench,i)))-1)] ;
      Hours_prc_change(x,:) = [i , prctile(panel_hours_bench,i) , prctile(panel_hours_exp_K,i) , prctile(panel_hours_exp_W,i) , ...
          100*(mean(panel_hours_exp_K(panel_hours_exp_K<=prctile(panel_hours_exp_K,i)))/ mean(panel_hours_bench(panel_hours_bench<=prctile(panel_hours_bench,i)))-1) ... 
          100*(mean(panel_hours_exp_W(panel_hours_exp_W<=prctile(panel_hours_exp_W,i)))/ mean(panel_hours_bench(panel_hours_bench<=prctile(panel_hours_bench,i)))-1)] ;
      Hours_by_cons_prc_change(x,:) = [i , prctile(panel_cons_bench,i) , prctile(panel_cons_exp_K,i) , prctile(panel_cons_exp_W,i) , ...
          100*(mean(panel_hours_exp_K(panel_cons_exp_K<=prctile(panel_cons_exp_K,i)))/ mean(panel_hours_bench(panel_cons_bench<=prctile(panel_cons_bench,i)))-1) ... 
          100*(mean(panel_hours_exp_W(panel_cons_exp_W<=prctile(panel_cons_exp_W,i)))/ mean(panel_hours_bench(panel_cons_bench<=prctile(panel_cons_bench,i)))-1)] ;
      x = x + 1 ;
    end

    
    Mat_1     = cell(9,1) ;
    col_title = {'prct','Cprc_bench','Cprc_opt_tk','Cprc_opt_tw','E[C_k<prc]/E[C_b<prc]','E[H_w<prc]/E[H_b<prc]',' '} ;
    Mat_2     = [col_title;num2cell(Cons_prc_change),Mat_1] ;
    col_title = {'prct','Hprc_bench','Hprc_opt_tk','Hprc_opt_tw','E[H_k<prc]/E[H_b<prc]','E[H_w<prc]/E[H_b<prc]',' '} ;
    Mat_3     = [col_title;num2cell(Hours_prc_change),Mat_1] ;
    col_title = {'prct','Cprc_bench','Cprc_opt_tk','Cprc_opt_tw','E[H_k<prc]/E[H_b<prc]','E[H_w<prc]/E[H_b<prc]',' '} ;
    Mat_4     = [col_title;num2cell(Hours_by_cons_prc_change),Mat_1] ;
    
    Mat_A = [Mat_2,Mat_3,Mat_4];
	status = xlwrite(Tables_file,Mat_A,'C&H change(Opt_Tax)') ;
          
    disp(Mat_A)
    
%% Return

Mat_0     = cell(1,20) ;
row_names = {'mean';'std';'mean_at';'std_at';'Weight_mean';'Weight_mean_at';' '} ;
prc_names = {'p10';'p25';'p50';'p75';'p90';'p95';'p99'} ;
col_names = {' ','Moments',' ','prc_return','age1','age2','age3','age4','age5','age6','age7',' ','prc_at_return','age1','age2','age3','age4','age5','age6','age7'};

% load bench 
    eval(['load ',Simul_Folder,'panelage_bench']) ; eval(['load ',Simul_Folder,'panela_bench']) ; 
    eval(['load ',Simul_Folder,'panel_return_bench']) ; eval(['load ',Simul_Folder,'panel_at_return_bench']) ;
    
    Mat_2 = NaN(8,7) ; Mat_3 = NaN(8,7) ;
    
    Mat_2(1,:) = [prctile(panel_return_bench./panela_bench,10) , prctile(panel_return_bench./panela_bench,25) , ...
                  prctile(panel_return_bench./panela_bench,50) , prctile(panel_return_bench./panela_bench,75) , ... 
                  prctile(panel_return_bench./panela_bench,90) , prctile(panel_return_bench./panela_bench,95) , ... 
                  prctile(panel_return_bench./panela_bench,99) ] ;
	Mat_3(1,:) = [prctile(panel_at_return_bench./panela_bench,10) , prctile(panel_at_return_bench./panela_bench,25) , ...
                  prctile(panel_at_return_bench./panela_bench,50) , prctile(panel_at_return_bench./panela_bench,75) , ... 
                  prctile(panel_at_return_bench./panela_bench,90) , prctile(panel_at_return_bench./panela_bench,95) , ... 
                  prctile(panel_at_return_bench./panela_bench,99) ] ;
              
    for age_group_counter=1:7
        indx1 = ((panelage_bench>age_limit(age_group_counter)).*(panelage_bench<=age_limit(age_group_counter+1)));
        indx1=indx1>0;
        temp_return=panel_return_bench(indx1)./panela_bench(indx1);
        Mat_2(1+age_group_counter,:) = [prctile(temp_return,10) , prctile(temp_return,25) , prctile(temp_return,50) , ...
                                        prctile(temp_return,75) , prctile(temp_return,90) , prctile(temp_return,95) , prctile(temp_return,99)] ;
        clear temp_return
        temp_return=panel_at_return_bench(indx1)./panela_bench(indx1);
        Mat_3(1+age_group_counter,:) = [prctile(temp_return,10) , prctile(temp_return,25) , prctile(temp_return,50) , ...
                                        prctile(temp_return,75) , prctile(temp_return,90) , prctile(temp_return,95) , prctile(temp_return,99)] ;
        clear indx1 temp_return
    end
    
    Mat_row_head = [{' ';'Bench'};cell(8,1);{'Exp'};cell(8,1);{'Opt_K'};cell(8,1);{'Opt_W'};cell(8,1)] ;
    Mat_col_A  = {' ','p10','p25','p50','p75','p90','p95','p99'} ; Mat_col_A = [Mat_col_A ,{' '}, Mat_col_A] ;
    Mat_row_A  = {'Pop';'Age1';'Age2';'Age3';'Age4';'Age5';'Age6';'Age7'};
    Mat_00     = cell(size(Mat_col_A)) ;
    Mat_A      = [Mat_row_A num2cell(100*Mat_2)] ;
    Mat_at_A   = [Mat_row_A num2cell(100*Mat_3)] ;
    
    % Small Table
    mean_R_ben = sum(panel_return_bench)/sum(panela_bench) ;
    std_R_ben  = sqrt(sum( ((panel_return_bench./panela_bench - mean_R_ben).^2).*panela_bench )/sum(panela_bench)) ;
    p10_R_ben  = prctile(panel_return_bench./panela_bench,10) ;
    p50_R_ben  = prctile(panel_return_bench./panela_bench,50) ;
    p90_R_ben  = prctile(panel_return_bench./panela_bench,90) ;
    p95_R_ben  = prctile(panel_return_bench./panela_bench,95) ;
    p99_R_ben  = prctile(panel_return_bench./panela_bench,99) ;
    
    mean_at_R_ben = sum(panel_at_return_bench)/sum(panela_bench) ;
    std_at_R_ben  = sqrt(sum( ((panel_at_return_bench./panela_bench - mean_at_R_ben).^2).*panela_bench )/sum(panela_bench)) ;
    p10_at_R_ben  = prctile(panel_at_return_bench./panela_bench,10) ;
    p50_at_R_ben  = prctile(panel_at_return_bench./panela_bench,50) ;
    p90_at_R_ben  = prctile(panel_at_return_bench./panela_bench,90) ;
    p95_at_R_ben  = prctile(panel_at_return_bench./panela_bench,95) ;
    p99_at_R_ben  = prctile(panel_at_return_bench./panela_bench,99) ;
    
    
    Mat_head   = [{'Before_Tax'} cell(1,8) {'After_Tax'} cell(1,7)] ;
    Mat_col    = {' ','Mean','Std','p10','p50','p90','p95','p99'} ; Mat_col = [Mat_col ,{' '}, Mat_col] ;
    Mat_row    = {'Bench';'Exp';'Opt_K';'Opt_W'};
    Mat_ben    = [mean_R_ben std_R_ben p10_R_ben p50_R_ben p90_R_ben p95_R_ben p99_R_ben];
    Mat_at_ben = [mean_at_R_ben std_at_R_ben p10_at_R_ben p50_at_R_ben p90_at_R_ben p95_at_R_ben p99_at_R_ben];
    Mat_B      = [100*Mat_ben] ;
    Mat_at_B   = [100*Mat_at_ben] ;
    
    
for j = 1:3
    if j==1
        % Experiment
        eval(['load ',Simul_Folder_Exp,'panelage_exp']) ; eval(['load ',Simul_Folder_Exp,'panela_exp']) ; 
        eval(['load ',Simul_Folder_Exp,'panel_return_exp']) ; eval(['load ',Simul_Folder_Exp,'panel_at_return_exp']) ;
    elseif j==2
        % Opt Tax K
        eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panelage_exp']) ; eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panela_exp']) ; 
        eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panel_return_exp']) ; eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panel_at_return_exp']) ;
    elseif j==3
        % Opt Tax W
        eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panelage_exp']) ; eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panela_exp']) ; 
        eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panel_return_exp']) ; eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panel_at_return_exp']) ;
    end
    
    
    Mat_2 = NaN(8,7) ; Mat_3 = NaN(8,7) ;
    
    Mat_2(1,:) = [prctile(panel_return_exp./panela_exp,10) , prctile(panel_return_exp./panela_exp,25) , ...
                  prctile(panel_return_exp./panela_exp,50) , prctile(panel_return_exp./panela_exp,75) , ... 
                  prctile(panel_return_exp./panela_exp,90) , prctile(panel_return_exp./panela_exp,95) , ... 
                  prctile(panel_return_exp./panela_exp,99) ] ;
	Mat_3(1,:) = [prctile(panel_at_return_exp./panela_exp,10) , prctile(panel_at_return_exp./panela_exp,25) , ...
                  prctile(panel_at_return_exp./panela_exp,50) , prctile(panel_at_return_exp./panela_exp,75) , ... 
                  prctile(panel_at_return_exp./panela_exp,90) , prctile(panel_at_return_exp./panela_exp,95) , ... 
                  prctile(panel_at_return_exp./panela_exp,99) ] ;
              
    for age_group_counter=1:7
        indx1 = ((panelage_exp>age_limit(age_group_counter)).*(panelage_exp<=age_limit(age_group_counter+1)));
        indx1=indx1>0;
        temp_return=panel_return_exp(indx1)./panela_exp(indx1);
        Mat_2(1+age_group_counter,:) = [prctile(temp_return,10) , prctile(temp_return,25) , prctile(temp_return,50) , ...
                                        prctile(temp_return,75) , prctile(temp_return,90) , prctile(temp_return,95) , prctile(temp_return,99)] ;
        clear temp_return
        temp_return=panel_at_return_exp(indx1)./panela_exp(indx1);
        Mat_3(1+age_group_counter,:) = [prctile(temp_return,10) , prctile(temp_return,25) , prctile(temp_return,50) , ...
                                        prctile(temp_return,75) , prctile(temp_return,90) , prctile(temp_return,95) , prctile(temp_return,99)] ;
        clear indx1 temp_return
    end
    
    Mat_A      = [Mat_A ; cell(1,8) ; Mat_row_A num2cell(100*Mat_2) ] ;
    Mat_at_A   = [Mat_at_A ; cell(1,8) ; Mat_row_A num2cell(100*Mat_3) ] ;
    
    % Small Table
    mean_R_exp = sum(panel_return_exp)/sum(panela_exp) ;
    std_R_exp  = sqrt(sum( ((panel_return_exp./panela_exp - mean_R_exp).^2).*panela_exp )/sum(panela_exp)) ;
    p10_R_exp  = prctile(panel_return_exp./panela_exp,10) ;
    p50_R_exp  = prctile(panel_return_exp./panela_exp,50) ;
    p90_R_exp  = prctile(panel_return_exp./panela_exp,90) ;
    p95_R_exp  = prctile(panel_return_exp./panela_exp,95) ;
    p99_R_exp  = prctile(panel_return_exp./panela_exp,99) ;
    
    mean_at_R_exp = sum(panel_at_return_exp)/sum(panela_exp) ;
    std_at_R_exp  = sqrt(sum( ((panel_at_return_exp./panela_exp - mean_at_R_exp).^2).*panela_exp )/sum(panela_exp)) ;
    p10_at_R_exp  = prctile(panel_at_return_exp./panela_exp,10) ;
    p50_at_R_exp  = prctile(panel_at_return_exp./panela_exp,50) ;
    p90_at_R_exp  = prctile(panel_at_return_exp./panela_exp,90) ;
    p95_at_R_exp  = prctile(panel_at_return_exp./panela_exp,95) ;
    p99_at_R_exp  = prctile(panel_at_return_exp./panela_exp,99) ;
    
    Mat_exp    = [mean_R_exp std_R_exp p10_R_exp p50_R_exp p90_R_exp p95_R_exp p99_R_exp];
    Mat_at_exp = [mean_at_R_exp std_at_R_exp p10_at_R_exp p50_at_R_exp p90_at_R_exp p95_at_R_exp p99_at_R_exp];
    Mat_B      = [Mat_B    ; 100*Mat_exp] ;
    Mat_at_B   = [Mat_at_B ; 100*Mat_at_exp] ;
     
end 

    Mat_A = [Mat_head; Mat_col_A; Mat_A cell(size(Mat_A,1),1) Mat_at_A] ;
    Mat_A = [Mat_row_head Mat_A] ;
    Mat_B = [Mat_head; Mat_col; Mat_row num2cell(Mat_B) cell(4,1) Mat_row num2cell(Mat_at_B)] ;
    
	status = xlwrite(Tables_file,Mat_A,'Return_Stats_1') ;
    status = xlwrite(Tables_file,Mat_B,'Return_Table_2') ;
          
    disp(Mat_A)
    disp(Mat_B)
    
    
%% Gini Coefficient with Burhan's data 
% I use Deaton's formulas for the gini. 
% The formulas are found in Deaton (1997) "The analysis of household surveys"
% The formulas are in page 139.

% Benchmark
    eval(['load ',Simul_Folder,'panela_bench']) ; 
    N  = numel(panela_bench) ;
    mu = mean(panela_bench)  ;
    panela_bench_sort = sort(panela_bench,'descend') ;
    index = 1:N ;
    G_bench = (N+1)/(N-1) - 2*sum(panela_bench_sort.*index)/(mu*N*(N-1)) ;
% Experiment 
    eval(['load ',Simul_Folder,'panela_exp']) ; 
    N  = numel(panela_exp) ;
    mu = mean(panela_exp)  ;
    panela_exp_sort = sort(panela_exp,'descend') ;
    index = 1:N ;
    G_exp = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu*N*(N-1)) ;
% Optimal Tau K
    eval(['load ',Result_Folder,'Opt_Tax_K/Simul/','panela_exp']) ; 
    N  = numel(panela_exp) ;
    mu = mean(panela_exp)  ;
    panela_exp_sort = sort(panela_exp,'descend') ;
    index = 1:N ;
    G_opt_K = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu*N*(N-1)) ;
% Optimal Tau W
    eval(['load ',Result_Folder,'Opt_Tax_W/Simul/','panela_exp']) ; 
    N  = numel(panela_exp) ;
    mu = mean(panela_exp)  ;
    panela_exp_sort = sort(panela_exp,'descend') ;
    index = 1:N ;
    G_opt_W = (N+1)/(N-1) - 2*sum(panela_exp_sort.*index)/(mu*N*(N-1)) ;

% Save results
    cols = {' ','Bench','Exp','Opt_tau_K','Opt_tau_W'};
    rows = {'Gini'} ;
    Mat = [cols; rows num2cell([G_bench G_exp G_opt_K G_opt_W])]

    status = xlwrite(Tables_file,Mat,'Gini') ;

%% Vote for reform
% Distribution
    eval(['load ',Bench_Folder,'DBN'])
    DBN_bench = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear DBN

    eval(['load ',Result_Folder,'Exp_Files/Exp_results_DBN']) ;
    DBN_exp = reshape(Exp_results_DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_DBN
% Vaule Function
    eval(['load ',Bench_Folder,'value'])
    V_bench = reshape(value,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear value

    eval(['load ',Result_Folder,'Exp_Files/Exp_results_value']) ;
    V_exp = reshape(Exp_results_value,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_value
    
% Vote
    yes_votes = 100*sum(sum(sum(sum(sum( (V_exp>=V_bench).*DBN_bench ))))) ;
    
    vote_mat = NaN(7,7) ;
    
    for z=1:n_z
    for age_group_counter=1:7
        V_bench_aux   = V_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,z,:,:)     ;
        V_exp_aux     = V_exp(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,z,:,:)       ;
        DBN_bench_aux = DBN_bench(age_limit(age_group_counter)+1:age_limit(age_group_counter+1),:,z,:,:)   ;
        
        size_group    = sum(sum(sum(sum(sum( DBN_bench_aux ))))) ;
        
        vote_mat(age_group_counter,z) = 100*sum(sum(sum(sum(sum( (V_exp_aux>=V_bench_aux).*DBN_bench_aux )))))/size_group ;
    end
    end
    
% Save Results
    Total = [{'Votes'},num2cell(yes_votes),{' ',' ',' ',' ',' ',' '}] ;
    Mat_0 = {' ',' ',' ',' ',' ',' ',' ',' '} ;
    col_names = {' ','z1','z2','z3','z4','z5','z6','z7'} ;
    row_names = {'Age1';'Age2';'Age3';'Age4';'Age5';'Age6';'Age7'};
    Mat = [Total;Mat_0;col_names; row_names , num2cell(vote_mat)]
    status = xlwrite(Tables_file,Mat,'Votes_Reform') ;
    
%% Pareto Tail
    eval(['load ',Simul_Folder,'panela_bench']) ; 
    eval(['load ',Simul_Folder,'panela_exp'])   ;
    eval(['load ',Bench_Folder,'YBAR'])                        ;
    eval(['load ',Result_Folder,'Exp_Files/Exp_results_YBAR']) ;
    n_w = 1000 ;
    wealth_grid = linspace(log(min(panela_bench)),log(max(panela_bench)),n_w)' ;
    Pareto_ben = NaN(n_w,1) ;
    Pareto_exp = NaN(n_w,1) ;
    N = numel(panela_bench) ;
    log_A_ben = log(panela_bench) ;
    log_A_exp = log(panela_exp) ;
    for i=1:n_w
       Pareto_ben(i) = log(sum(log_A_ben>=wealth_grid(i))/N) ;
       Pareto_exp(i) = log(sum(log_A_exp>=wealth_grid(i))/N) ;
    end
    
    figure;
    plot(wealth_grid,[Pareto_ben Pareto_exp]); xlim([wealth_grid(1),wealth_grid(end)])
    xlabel('Log(Wealth)'); legend('Bench','Exp','location','northeast')
    print('-depsc','Wealth_Pareto.eps') ;
    
    GDP_pc = 54629 ;
    wealth_bench = sort(panela_bench) ; wealth_bench = GDP_pc/YBAR * wealth_bench ;
    wealth_exp   = sort(panela_exp)   ; wealth_exp   = GDP_pc/YBAR * wealth_exp   ;
    Exp_results_YBAR = GDP_pc/YBAR * Exp_results_YBAR ; YBAR =  GDP_pc ;
    nn=[1,floor(0.50*N),floor(0.75*N),floor(0.90*N),floor(0.95*N),floor(0.99*N)];
    ind = [100,50,25,10,5,1] ;
    for ii=1:numel(nn)
        % Bench
        w_min =  wealth_bench(nn(ii))                  ;
        x     = -log(wealth_bench(nn(ii):end)/w_min)   ;
        y     =  log(((N-nn(ii)+1):-1:1)/(N-nn(ii)+1)) ;
        xx    =  linspace(min(-x),max(-x),3)           ;
        mdl   =  fitlm(x',y')                          ;
        cons_bench(ii)     = mdl.Coefficients.Estimate(1)                          ;
        alpha_bench(ii)    = mdl.Coefficients.Estimate(2)                          ;
        alpha_bench_ml(ii) = (N-nn(ii)+1)/sum(log(wealth_bench(nn(ii):end)/w_min)) ;
        SE_bench(ii)       = mdl.Coefficients.SE(2)                                ;
        R2_bench(ii)       = mdl.Rsquared.Ordinary                                 ;
        fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% bench'] ;
        figure; hold on; plot(-x,y,'linewidth',2); plot(xx,cons_bench(ii)-alpha_bench(ii)*xx,'-k'); hold off;
        title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
        file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_bench.fig'] ;
        savefig(file_name); %print('-depsc',file_name) ;
        % Tax Reform
        w_min =  wealth_exp(nn(ii))                    ;
        x     = -log(wealth_exp(nn(ii):end)/w_min)     ;
        y     =  log(((N-nn(ii)+1):-1:1)/(N-nn(ii)+1)) ;
        mdl   =  fitlm(x',y')                          ;
        alpha_exp(ii)    = mdl.Coefficients.Estimate(2)                          ;
        alpha_exp_ml(ii) = (N-nn(ii)+1)/sum(log(wealth_exp(nn(ii):end)/w_min))   ;
        SE_exp(ii)       = mdl.Coefficients.SE(2)                                ;
        R2_exp(ii)       = mdl.Rsquared.Ordinary                                 ;
    end
    
    AA_bench = [100 50 25 10 05 01;
                alpha_bench;
                SE_bench;
                R2_bench;
                alpha_bench_ml;
                wealth_bench(nn);
                wealth_bench(end),NaN,NaN,NaN,NaN,NaN;
                YBAR,NaN,NaN,NaN,NaN,NaN];
	AA_exp   = [100 50 25 10 05 01;
                alpha_exp;
                SE_exp;
                R2_exp;
                alpha_exp_ml;
                wealth_exp(nn);
                wealth_exp(end),NaN,NaN,NaN,NaN,NaN;
                Exp_results_YBAR,NaN,NaN,NaN,NaN,NaN];
          
    bench_name = {'Bench',' ',' ',' ',' ',' ',' '} ;
    Exp_name = {'Exp',' ',' ',' ',' ',' ',' '} ;
    blank_row = cell(1,7) ;
    row_name = {'top_x%';'Pareto_ind_reg';'SE';'R2';'Pareto_ind_ml';'Min_Wealth';'Max_Wealth';'YBAR'} ;
    Mat = [bench_name ; row_name num2cell(AA_bench) ; blank_row ; Exp_name; row_name num2cell(AA_exp) ]
    status = xlwrite(Tables_file,Mat,'Pareto_Tail') ;
    
    
    
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
    print('-depsc','Wealth_Prctl.eps') ;
    

    
% % %% Percentiles of the asset distribution and distribution of z by them
% %     for i=1:n_z 
% %         prc_totA_ben(i) = prctile(panela_bench,prctl(i)) ;
% %         prc_totA_exp(i) = prctile(panela_exp,prctl(i)) ;
% %     end 
% %     
% %     for i=1:n_z
% %         if i==1
% %            ind = (panela_bench<=prc_totA_ben(i))  ;
% %         elseif i<n_z && i>1
% %            ind = ((panela_bench<=prc_totA_ben(i)).*(panela_bench>prc_totA_ben(i-1)))==1  ;
% %         else
% %            ind = (panela_bench>prc_totA_ben(i))  ;
% %         end
% %         n=sum(ind);
% %         
% %         for k=1:n_z
% %             eval(strcat('Z_Aprc_ben(i,k) = 100*sum(panelz_bench(ind)==',int2str(k),')/n;'))
% % 
% %             if eval(strcat('sum(panelz_bench(ind)==',int2str(k),')'))>n
% %                 beep
% %                 [j,i,k]
% %                 error('this is wrong')
% %             end 
% %         end
% %         
% %         if i==1
% %            ind = (panela_exp<=prc_totA_exp(i))  ;
% %         elseif i<n_z && i>1
% %            ind = ((panela_exp<=prc_totA_exp(i)).*(panela_exp>prc_totA_exp(i-1)))==1  ;
% %         else
% %            ind = (panela_exp>prc_totA_exp(i))  ;
% %         end
% %         n=sum(ind);
% %         
% %         for k=1:n_z
% %             eval(strcat('Z_Aprc_exp(i,k) = 100*sum(panelz_exp(ind)==',int2str(k),')/n;'))
% % 
% %             if eval(strcat('sum(panelz_exp(ind)==',int2str(k),')'))>n
% %                 beep
% %                 [j,i,k]
% %                 error('this is wrong')
% %             end 
% %         end
% %     end   
% %     
% %     z_title = {' ','z1','z2','z3','z4','z5','z6','z7'} ;
% %     Mat_Z_ben = [{'Z_dist_by_A_prct_ben',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_ben)] ;
% %     Mat_Z_exp = [{'Z_dist_by_A_prct_exp',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_exp)] ;
% %     Mat_Z_diff = [{'Z_dist_by_A_prct_diff',' ',' ',' ',' ',' ',' ',' '};z_title;AW_title(2:end)',num2cell(Z_Aprc_exp-Z_Aprc_ben)] ;
% %     Mat = [Mat_7;Mat_8 Mat_Z_ben Mat_8 Mat_Z_exp Mat_8 Mat_Z_diff] ;
% %     status = xlwrite(AW_file,Mat,'Z Dist. Tot. A') ; 
% %     

    
    
%% Optimal Tax Graphs

% Load stats for tau K tau W
    eval(['load ',Result_Folder,'Opt_Tax_K/','Stats_by_tau_k.txt']) ; eval(['load ',Result_Folder,'Opt_Tax_W/','Stats_by_tau_w.txt']) ; 
    % ORDER OF THINGS
        % 1. tauk, 2. tauw, 3. taul, 4. GBAR_K/(GBAR+SSC), 5. KBAR, 6. QBAR, 7. NBAR, 8. YBAR
        % 9. %Y , 10. wage, 11. VBAR, 12. CE2_NB, 13. CE2, 14. KBAR/Y, 
        % 15. Top1%, 16. Top10%, 17. std(log(E)), 18. Av. Hours
    
  
    aa = Stats_by_tau_k(:,4) ;
    aa = aa*0.16004634619224012/(0.16004634619224012+0.12643367375356479);
    Stats_by_tau_k(:,4) = aa ;
    aa = Stats_by_tau_w(:,4) ;
    aa = aa*0.16004634619224012/(0.16004634619224012+0.12643367375356479);
    Stats_by_tau_w(:,4) = aa ;
    
    G_K_Frac_k = Stats_by_tau_k(:,4) ;
    G_K_Frac_w = Stats_by_tau_w(:,4) ;
    

% Aggregate variabls as a function of Borrowing/wealth tax revenue
    clf
    hold on
    plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,5)/Stats_by_tau_k(1,5)-1),'r')
    plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,5)/Stats_by_tau_w(1,5)-1),'b')
    plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,6)/Stats_by_tau_k(1,6)-1),'r-d')
    plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,6)/Stats_by_tau_w(1,6)-1),'b-d')
    plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,7)/Stats_by_tau_k(1,7)-1),'ro')
    plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,7)/Stats_by_tau_w(1,7)-1),'bo')
    plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,8)/Stats_by_tau_k(1,8)-1),'rd')
    plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,8)/Stats_by_tau_w(1,8)-1),'bd')
    grid on
    hold off
    legend(' KBAR, \tau_K', 'KBAR, \tau_W', 'QBAR, \tau_K', 'QBAR, \tau_W', ...
        'NBAR, \tau_K', 'NBAR, \tau_W', 'YBAR, \tau_K', 'YBAR, \tau_W','Location','SouthWest')
    xlabel('GBAR_K')

    hgsave('fig/1fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
    print -depsc 1fig_KBAR_QBAR_by_CAP_TAX_REV.eps
    print -dps  1fig_KBAR_QBAR_by_CAP_TAX_REV.eps
    print -dpng 1fig_KBAR_QBAR_by_CAP_TAX_REV.png
    
    % KBAR & QBAR graph
        clf
        hold on
        plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,5)/Stats_by_tau_k(1,5)-1))
        plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,5)/Stats_by_tau_w(1,5)-1))
        xlabel('Tax Revenue from K / Total Tax Revenue')
        ylabel('Percent Change')
        h = legend('$\bar k, \tau_k$' , '$\bar k, \tau_a$','Location','SouthWest');
        set(h,'Interpreter','latex','FontSize',20)
        grid on
        axis([min(G_K_Frac_w) max(G_K_Frac_w)  -50 0 ])


        hgsave('fig/1.1fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
        print -depsc 1.1fig_KBAR_QBAR_by_CAP_TAX_REV.eps
        print -dpng png/1.1fig_KBAR_QBAR_by_CAP_TAX_REV.png

        plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,6)/Stats_by_tau_k(1,6)-1))
        plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,6)/Stats_by_tau_w(1,6)-1))
        h2 = legend('$\bar k, \tau_k$' , '$\bar k, \tau_a$','$\bar Q, \tau_k$' , '$\bar Q, \tau_a$','Location','SouthWest');
        set(h2,'Interpreter','latex','FontSize',20)

        hgsave('fig/1.2fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
        print -depsc 1.2fig_KBAR_QBAR_by_CAP_TAX_REV.eps
        print -dpng png/1.2fig_KBAR_QBAR_by_CAP_TAX_REV.png


% Wage by capital/Wealth tax revenue
    clf
    hold on
    plot(G_K_Frac_k, (Stats_by_tau_k(:,10)),'r')
    plot(G_K_Frac_w, (Stats_by_tau_w(:,10)),'b')
    grid on
    hold off
    legend('\tau_K', '\tau_W')
    xlabel('GBAR_K')
    title('Wage')
    [maxvalK indxK] = max((Stats_by_tau_k(:,11)-Stats_by_tau_k(1,11)));
    text(Stats_by_tau_k(indxK,4)+0.0, maxvalK -0.15,['Opt. \tau_K = ' ,num2str(100*Stats_by_tau_k(indxK,1)),'%'])

    [maxvalW indxW] = max((Stats_by_tau_w(:,10)-Stats_by_tau_w(1,11)));
    text(Stats_by_tau_w(indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W = ' ,num2str(100*Stats_by_tau_w(indxW,2)),'%'])


    hgsave('fig/1fig_Wage_by_CAP_TAX_REV.fig')
    print -depsc 1fig_Wage_by_CAP_TAX_REV.eps
    print -dpng 1fig_Wage_by_CAP_TAX_REV.png
    
    % After-tax wage
        clf
        hold on
        plot(G_K_Frac_k, (Stats_by_tau_k(:,10).*(1-Stats_by_tau_k(:,3))))
        plot(G_K_Frac_w, (Stats_by_tau_w(:,10).*(1-Stats_by_tau_w(:,3))))
        
        [x,wagemax_indx_tauK] = max(Stats_by_tau_k(:,10).*(1-Stats_by_tau_k(:,3)));
        [x,wagemax_indx_tauW] = max(Stats_by_tau_w(:,10).*(1-Stats_by_tau_w(:,3)));
        display('    tauK      tauW   that maximizes after-tax wage')
        disp([Stats_by_tau_k(wagemax_indx_tauK,1) Stats_by_tau_w(wagemax_indx_tauW,2)])
        hgsave('1fig_Wage_by_CAP_TAX_REV.fig')
        print -depsc 1fig_AT_Wage_by_CAP_TAX_REV.eps
        print -dps  1fig_AT_Wage_by_CAP_TAX_REV.eps
        print -dpng 1fig_AT_Wage_by_CAP_TAX_REV.png

% Average utilities (relavite to no taxes)
    clf
    hold on
    plot(G_K_Frac_k, (Stats_by_tau_k(:,11)-Stats_by_tau_k(1,11)),'r')
    plot(G_K_Frac_w, (Stats_by_tau_w(:,11)-Stats_by_tau_w(1,11)),'b')
    grid on; hold off
    xlabel('GBAR_K')
    title('Average Utility')

    [maxvalK indxK] = max((Stats_by_tau_k(:,11)-Stats_by_tau_k(1,11)));
    text(Stats_by_tau_k(indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K = ' ,num2str(100*Stats_by_tau_k(indxK,1)),'%'])

    [maxvalW indxW] = max((Stats_by_tau_w(:,11)-Stats_by_tau_w(1,11)));
    text(Stats_by_tau_w(indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W = ' ,num2str(100*Stats_by_tau_w(indxW,2)),'%'])

    hgsave('1fig_Average_Utility_by_CAP_TAX_REV.fig')
    print -depsc 1fig_Average_Utility_by_CAP_TAX_REV.eps
    print -dps  1fig_Average_Utility_by_CAP_TAX_REV.eps
    print -dpng 1fig_Average_Utility_by_CAP_TAX_REV.png

CE2_NB_k =100*( (Stats_by_tau_k(:,11)./Stats_by_tau_k(26,11)).^(1/(gamma*(1-sigma)))-1);
CE2_NB_w =100*( (Stats_by_tau_w(:,11)./Stats_by_tau_k(26,11)).^(1/(gamma*(1-sigma)))-1);

% CE Newborn computed using average utilities (relavite to no taxes)
    clf
    hold on
    plot(G_K_Frac_k, (CE2_NB_k-CE2_NB_k(1)),'r')
    plot(G_K_Frac_w, (CE2_NB_w-CE2_NB_w(1)),'b')
    grid on; hold off;
    xlabel('GBAR_K')
    title('CE Newborn computed using average utilities')

    [maxvalK indxK] = max((CE2_NB_k-CE2_NB_k(1)));
    text(Stats_by_tau_k(indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K = ' ,num2str(100*Stats_by_tau_k(indxK,1)),'%'])

    [maxvalW indxW] = max((CE2_NB_w-CE2_NB_w(1)));
    text(Stats_by_tau_w(indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W = ' ,num2str(100*Stats_by_tau_w(indxW,2)),'%'])

    hgsave('2fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
    print -depsc 2fig_CE_NEWBORN_by_CAP_TAX_REV.eps
    print -dps  2fig_CE_NEWBORN_by_CAP_TAX_REV.eps
    print -dpng 2fig_CE_NEWBORN_by_CAP_TAX_REV.png

% CE Newborn computed using average utilities (in levels)
    clf
    hold on
    plot(G_K_Frac_k, (CE2_NB_k ),'r')
    plot(G_K_Frac_w, (CE2_NB_w ),'b')
    grid on; hold off;
    xlabel('GBAR_K')
    title('CE Newborn computed using average utilities')

    [maxvalK indxK] = max((CE2_NB_k ));
    text(Stats_by_tau_k(indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K = ' ,num2str(100*Stats_by_tau_k(indxK,1)),'%'])

    [maxvalW indxW] = max((CE2_NB_w ));
    text(Stats_by_tau_w(indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W = ' ,num2str(100*Stats_by_tau_w(indxW,2)),'%'])


    hgsave('3fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
    print -depsc 3fig_CE_NEWBORN_by_CAP_TAX_REV.eps
    print -dps  3fig_CE_NEWBORN_by_CAP_TAX_REV.eps
    print -dpng 3fig_CE_NEWBORN_by_CAP_TAX_REV.png


% Burhan's graphs on welfare by tau
    [x,tauindx ] = min(abs(Stats_by_tau_k(:,1) - 0.25)) ;

    clf
    axes1 = axes(...
      'FontName','Times New Roman',...
      'FontSize',14);
    hold on
    plot(G_K_Frac_k, (CE2_NB_k-CE2_NB_k(tauindx)),'r', 'linewidth',2)
    plot(G_K_Frac_k, zeros(size(G_K_Frac_k)),'k--', 'linewidth',2)
    grid on
    xlabel('Tax Revenue from K / Total Tax Revenue')
    ylabel('$CE_2$ Welfare Change from \textbf{Benchmark}','Interpreter','latex')
    plot(G_K_Frac_k(tauindx)*ones(2),[min(CE2_NB_k-CE2_NB_k(tauindx))-1 0],'k--', 'linewidth',2)
    axis([min(G_K_Frac_k) max(G_K_Frac_k)  min(CE2_NB_k-CE2_NB_k(tauindx))-1 max(CE2_NB_w-CE2_NB_k(tauindx))+1 ])

    annotation('textarrow',[0.71 0.68],[0.55 0.52], 'string','Cap. Income Tax Economy','linewidth',2,'FontSize',12)
    annotation('textarrow',[0.68 0.73],[0.39 0.44], 'string','Benchmark, \tau_k = 25%','linewidth',2,'FontSize',12)
    annotation('textarrow',[0.60 0.65],[0.17 0.12], 'string','0.28','linewidth',2,'FontSize',12)


    hgsave('1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -depsc 1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps
    print -dps 1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps

    annotation('textarrow',[0.2 0.15],[0.71 0.76], 'string','Opt. \tau_k = 1.62%','linewidth',2,'FontSize',12)
    hgsave('1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -depsc 1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps
    print -dps 1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps

    plot(G_K_Frac_w, (CE2_NB_w-CE2_NB_k(tauindx)),'b', 'linewidth',2)
    annotation('textarrow',[0.47 0.47],[0.81 0.86], 'string','Opt. \tau_a = 1.54%','linewidth',2,'FontSize',12)


    hgsave('1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -depsc 1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps
    print -dps 1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps


%% Tables from experiment
% Load
    eval(['load ',Result_Folder,'mean_wealth_by_agegroup_z_bench']) 
    eval(['load ',Result_Folder,'agrid']) 
    eval(['load ',Result_Folder,'aprime_age1_bench']) 
    eval(['load ',Result_Folder,'aprime_age1_exp'])
    eval(['load ',Result_Folder,'aprime_age16_bench'])
    eval(['load ',Result_Folder,'aprime_age16_exp'])
    eval(['load ',Result_Folder,'aprime_age31_bench'])
    eval(['load ',Result_Folder,'aprime_age31_exp'])
    eval(['load ',Result_Folder,'MeanAsset_by_z_med_E_Lambda_age1'])
    eval(['load ',Result_Folder,'CE_by_Asset_z_med_E_Lambda_age1'])
    eval(['load ',Result_Folder,'MeanAsset_by_z_med_E_Lambda_age16'])
    eval(['load ',Result_Folder,'CE_by_Asset_z_med_E_Lambda_age16'])
    eval(['load ',Result_Folder,'MeanAsset_by_z_med_E_Lambda_age31'])
    eval(['load ',Result_Folder,'CE_by_Asset_z_med_E_Lambda_age31'])
    
    axes1 = axes(...
      'FontName','Times New Roman',...
      'FontSize',18);
    
    
%     
    clf
    alow=8;
    ahigh=20;
    hold on
    plot(agrid(alow:ahigh),aprime_age1_exp(1,alow:ahigh)./agrid(alow:ahigh)-aprime_age1_bench(1,alow:ahigh)./agrid(alow:ahigh),'b','LineWidth',2)
    plot(agrid(alow:ahigh),aprime_age1_exp(4,alow:ahigh)./agrid(alow:ahigh)-aprime_age1_bench(4,alow:ahigh)./agrid(alow:ahigh),'k','LineWidth',2)
    plot(agrid(alow:ahigh),aprime_age1_exp(7,alow:ahigh)./agrid(alow:ahigh)-aprime_age1_bench(7,alow:ahigh)./agrid(alow:ahigh),'r','LineWidth',2)
    plot(agrid(alow:ahigh),zeros(size(agrid(alow:ahigh))),'k--')
    axis([min(agrid(alow:ahigh)) max(agrid(alow:ahigh)) ...
        min(aprime_age1_exp(1,alow:ahigh)./agrid(alow:ahigh)-aprime_age1_bench(1,alow:ahigh)./agrid(alow:ahigh))-0.1...
        max(aprime_age1_exp(7,alow:ahigh)./agrid(alow:ahigh)-aprime_age1_bench(7,alow:ahigh)./agrid(alow:ahigh))+0.1 ])
    xlabel('wealth')
    legend('z(1)','z(4)','z(7)')
    grid on

    hgsave('1fig_diff_savings_rate_age1.fig')
    print -depsc 1fig_diff_savings_rate_age1.eps
    print -dps  1fig_diff_savings_rate_age1.eps
    print -dpng 1fig_diff_savings_rate_age1.png
naa=40;
nzz=1;

clf
hold on
plot(agrid(1:naa),CE_by_Asset_z_med_E_Lambda_age1(nzz,1:naa)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age1(nzz),0,'d','LineWidth',4)
plot(agrid(1:naa),zeros(size(agrid(1:naa))),'k--')
xlabel('wealth')
title('Welfare Gain for age=20, z=1')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age1_z1.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age1_z1.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age1_z1.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age1_z1.png


naa=40;
nzz=4;

clf
hold on
plot(agrid(1:naa),CE_by_Asset_z_med_E_Lambda_age1(nzz,1:naa)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age1(nzz),0,'d','LineWidth',4)
plot(agrid(1:naa),zeros(size(agrid(1:naa))),'k--')
xlabel('wealth')
title('Welfare Gain for age=20, z=4')
axis([0 100 -15 25])
grid on

print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age1_z4.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age1_z4.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age1_z4.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age1_z4.png


naa=40;
nzz=6;

clf
hold on
plot(agrid(1:naa),CE_by_Asset_z_med_E_Lambda_age1(nzz,1:naa)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age1(nzz),0,'d','LineWidth',4)
plot(agrid(1:naa),zeros(size(agrid(1:naa))),'k--')
xlabel('wealth')
title('Welfare Gain for age=20, z=7')
axis([0 100 -15 25])
grid on

print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age1_z7.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age1_z7.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age1_z7.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age1_z7.png



clf
hold on
plot(agrid(1:40),CE_by_Asset_z_med_E_Lambda_age16(1,1:40)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age16(1),0,'d','LineWidth',4)
plot(agrid(1:40),zeros(size(agrid(1:40))),'k--')
xlabel('wealth')
title('Welfare Gain for age=35, z=1')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age16_z1.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age16_z1.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age16_z1.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age16_z1.png


clf
hold on
plot(agrid(1:40),CE_by_Asset_z_med_E_Lambda_age16(4,1:40)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age16(4),0,'d','LineWidth',4)
plot(agrid(1:40),zeros(size(agrid(1:40))),'k--')
xlabel('wealth')
title('Welfare Gain for age=35, z=4')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age16_z4.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age16_z4.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age16_z4.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age16_z4.png


clf
hold on
plot(agrid(1:40),CE_by_Asset_z_med_E_Lambda_age16(7,1:40)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age16(7),0,'d','LineWidth',4)
plot(agrid(1:40),zeros(size(agrid(1:40))),'k--')
xlabel('wealth')
title('Welfare Gain for age=35, z=7')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age16_z7.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age16_z7.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age16_z7.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age16_z7.png


clf
hold on
plot(agrid(1:40),CE_by_Asset_z_med_E_Lambda_age31(1,1:40)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age31(1),0,'d','LineWidth',4)
plot(agrid(1:40),zeros(size(agrid(1:40))),'k--')
xlabel('wealth')
title('Welfare Gain for age=50, z=1')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age31_z1.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age31_z1.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age31_z1.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age31_z1.png


clf
hold on
plot(agrid(1:40),CE_by_Asset_z_med_E_Lambda_age31(4,1:40)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age31(4),0,'d','LineWidth',4)
plot(agrid(1:40),zeros(size(agrid(1:40))),'k--')
xlabel('wealth')
title('Welfare Gain for age=50, z=4')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age31_z4.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age31_z4.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age31_z4.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age31_z4.png


clf
hold on
plot(agrid(1:40),CE_by_Asset_z_med_E_Lambda_age31(7,1:40)*100,'LineWidth',2)
plot(MeanAsset_by_z_med_E_Lambda_age31(7),0,'d','LineWidth',4)
plot(agrid(1:40),zeros(size(agrid(1:40))),'k--')
xlabel('wealth')
title('Welfare Gain for age=50, z=7')
axis([0 100 -15 25])
grid on
print -depsc 1fig_Welfare_Gain_by_Wealth_for_Age31_z7.eps
hgsave('1fig_Welfare_Gain_by_Wealth_for_Age31_z7.fig')
print -dps  1fig_Welfare_Gain_by_Wealth_for_Age31_z7.eps
print -dpng 1fig_Welfare_Gain_by_Wealth_for_Age31_z7.png


%------------------------

alow=9;
ahigh=20;
clf
hold on
plot(agrid(alow:ahigh),aprime_age16_exp(1,alow:ahigh)./agrid(alow:ahigh)-aprime_age16_bench(1,alow:ahigh)./agrid(alow:ahigh),'b','LineWidth',2)
plot(agrid(alow:ahigh),aprime_age16_exp(4,alow:ahigh)./agrid(alow:ahigh)-aprime_age16_bench(4,alow:ahigh)./agrid(alow:ahigh),'k','LineWidth',2)
plot(agrid(alow:ahigh),aprime_age16_exp(7,alow:ahigh)./agrid(alow:ahigh)-aprime_age16_bench(7,alow:ahigh)./agrid(alow:ahigh),'r','LineWidth',2)
plot(agrid(alow:ahigh),zeros(size(agrid(alow:ahigh))),'k--')
axis([min(agrid(alow:ahigh)) max(agrid(alow:ahigh)) ...
    min(aprime_age16_exp(1,alow:ahigh)./agrid(alow:ahigh)-aprime_age16_bench(1,alow:ahigh)./agrid(alow:ahigh))-0.1...
    max(aprime_age16_exp(7,alow:ahigh)./agrid(alow:ahigh)-aprime_age16_bench(7,alow:ahigh)./agrid(alow:ahigh))+0.1 ])
xlabel('wealth')
legend('z(1)','z(4)','z(7)')
grid on


hgsave('1fig_diff_savings_rate_age16.fig')
print -depsc 1fig_diff_savings_rate_age16.eps
print -dps  1fig_diff_savings_rate_age16.eps
print -dpng 1fig_diff_savings_rate_age16.png


alow=9;
ahigh=20;
clf
hold on
plot(agrid(alow:ahigh),aprime_age31_exp(1,alow:ahigh)./agrid(alow:ahigh)-aprime_age31_bench(1,alow:ahigh)./agrid(alow:ahigh),'b','LineWidth',2)
plot(agrid(alow:ahigh),aprime_age31_exp(4,alow:ahigh)./agrid(alow:ahigh)-aprime_age31_bench(4,alow:ahigh)./agrid(alow:ahigh),'k','LineWidth',2)
plot(agrid(alow:ahigh),aprime_age31_exp(7,alow:ahigh)./agrid(alow:ahigh)-aprime_age31_bench(7,alow:ahigh)./agrid(alow:ahigh),'r','LineWidth',2)
plot(agrid(alow:ahigh),zeros(size(agrid(alow:ahigh))),'k--')
axis([min(agrid(alow:ahigh)) max(agrid(alow:ahigh)) ...
    min(aprime_age31_exp(1,alow:ahigh)./agrid(alow:ahigh)-aprime_age31_bench(1,alow:ahigh)./agrid(alow:ahigh))-0.1...
    max(aprime_age31_exp(7,alow:ahigh)./agrid(alow:ahigh)-aprime_age31_bench(7,alow:ahigh)./agrid(alow:ahigh))+0.1 ])
xlabel('wealth')
legend('z(1)','z(4)','z(7)')
grid on; hold off;


hgsave('1fig_diff_savings_rate_age31.fig')
print -depsc 1fig_diff_savings_rate_age31.eps
print -dps  1fig_diff_savings_rate_age31.eps
print -dpng 1fig_diff_savings_rate_age31.png


%% Wealth measures

    eval(['load ',Simul_Folder,'panela_bench'])             ; 
    eval(['load ',Simul_Folder,'panel_firm_wealth_bench'])  ; 
    eval(['load ',Simul_Folder,'panelage_bench'])           ; 
    eval(['load ',Bench_Folder,'EBAR'])                     ;
    
    Wage_earnings = 47000; 
    panel_a  = Wage_earnings/EBAR * panela_bench ;
    panel_PV = Wage_earnings/EBAR * panel_firm_wealth_bench ;
    
    mean_w = [mean(panel_a) mean(panel_PV)];
    max_w  = [max(panel_a)  max(panel_PV) ];
    tot_w  = [sum(panel_a)  sum(panel_PV) ];
    
    i=1;
    for prc=100-[0.01 0.10 1.00 10.00 20.00 40.00 50.00 90.00 99.00]
        a_prc  = prctile(panel_a ,prc);
        PV_prc = prctile(panel_PV,prc);
        W_a    = sum(panel_a(panel_a>=a_prc))/tot_w(1)   ;
        W_PV   = sum(panel_PV(panel_PV>=PV_prc))/tot_w(2);
        AA(i,:)= [prc W_a W_PV]                        ;
        i=i+1;
    end
    
    % Age brackets 
    age = [6  11, 16, 21, 26, 31, 36, 41, 46, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;
    
    
    %  Mean and std by age
    for j=1:n_age
        if j==1
            ind = (panelage_bench<age(j))  ;
        elseif j<n_age && j>1
            ind = ((panelage_bench>=age(j-1)).*(panelage_bench<age(j)))==1  ;
        else
            ind = (panelage_bench>=age(j-1))  ;
        end
        mean_wealth_age(j,1) = mean(panel_a(ind))   ;
        mean_wealth_age(j,2) = mean(panel_PV(ind))  ;
    end
    
    for j=1:Max_Age
        ind = (panelage_bench==j)  ;
        
        mean_wealth_age_all(j,1) = mean(panel_a(ind))   ;
        mean_wealth_age_all(j,2) = mean(panel_PV(ind))  ;
    end
    
    BB = [ 19+age' mean_wealth_age ];
    
    CC = [ 19+(1:Max_Age)' mean_wealth_age_all ];
    
    col_name  = {' ','assets','present_value'};
    row_name  = {'Mean';'Max'};
    col_name1 = {'pcrt','assets','present_value'};
    col_name2 = {'Age-Group','assets','present_value'};
    col_name3 = {'Age','assets','present_value'};
    Mat = [col_name;row_name num2cell([mean_w;max_w]);col_name1;num2cell(AA);col_name2;num2cell(BB);col_name2;num2cell(CC)]
    status = xlwrite(Tables_file,Mat,'Wealth_Stats') ;
    
%% Constrained firms
    % Grids    
        % A grid
        eval(['load ',Result_Folder,'agrid']);
        A_mat = repmat(agrid,[Max_Age,1,n_z,n_l,n_e]);

        % Z grid
        eval(['load ',Result_Folder,'zgrid']);
        Z_mat = repmat(reshape(zgrid,[1,1,n_z,1,1]),[Max_Age,n_a,1,n_l,n_e]);

    % Distribution
    eval(['load ',Bench_Folder,'DBN'])
    DBN_bench = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear DBN

    eval(['load ',Result_Folder,'Exp_Files/Exp_results_DBN']) ;
    DBN_exp = reshape(Exp_results_DBN,[Max_Age,n_a,n_z,n_l,n_e]) ;
    clear Exp_results_DBN
    
    % Constrained Firms output and profits
    eval(['load ',Result_Folder,'Exp_Files/constrained_ind_bench']) ;
    constrained_ind_bench = reshape(constrained_ind_bench,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    eval(['load ',Result_Folder,'Exp_Files/constrained_ind_exp']) ;
    constrained_ind_exp = reshape(constrained_ind_exp,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    eval(['load ',Result_Folder,'Exp_Files/firm_output_bench']) ;
    firm_output_bench = reshape(firm_output_bench,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    eval(['load ',Result_Folder,'Exp_Files/firm_output_exp']) ;
    firm_output_exp = reshape(firm_output_exp,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    eval(['load ',Result_Folder,'Exp_Files/firm_profit_bench']) ;
    firm_profit_bench = reshape(firm_profit_bench,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    eval(['load ',Result_Folder,'Exp_Files/firm_profit_exp']) ;
    firm_profit_exp = reshape(firm_profit_exp,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    Capital_bench   = firm_output_bench./Z_mat ;
    Capital_exp     = firm_output_exp./Z_mat   ;
    
    Borrowing_bench = Capital_bench - A_mat    ;
    Borrowing_exp   = Capital_exp   - A_mat    ;
    
    
    % Aggregate measures
    Constrained_firms_tot(1) = sum(sum(sum(sum(sum(constrained_ind_bench.*DBN_bench)))));
    Constrained_firms_tot(2) = sum(sum(sum(sum(sum(constrained_ind_exp.*DBN_exp)))))    ;
    Capital_tot(1)           = sum(sum(sum(sum(sum(Capital_bench.*DBN_bench)))))        ;
    Capital_tot(2)           = sum(sum(sum(sum(sum(Capital_exp.*DBN_exp)))))            ;
    Borrowing_tot(1)         = sum(sum(sum(sum(sum(Borrowing_bench.*DBN_bench)))))      ;
    Borrowing_tot(2)         = sum(sum(sum(sum(sum(Borrowing_exp.*DBN_exp)))))          ;
    Assets_tot(1)            = sum(sum(sum(sum(sum(A_mat.*DBN_bench)))))                ;
    Assets_tot(2)            = sum(sum(sum(sum(sum(A_mat.*DBN_exp)))))                  ;
    for age = 1:Max_Age
        Constrained_firms_age(age,1) = sum(sum(sum(sum(constrained_ind_bench(age,:,:,:,:).*DBN_bench(age,:,:,:,:)))))/sum(sum(sum(sum(DBN_bench(age,:,:,:,:)))));
        Constrained_firms_age(age,2) = sum(sum(sum(sum(constrained_ind_exp(age,:,:,:,:).*DBN_exp(age,:,:,:,:)))))/sum(sum(sum(sum(DBN_exp(age,:,:,:,:)))));
        Capital_age(age,1)           = sum(sum(sum(sum(Capital_bench(age,:,:,:,:).*DBN_bench(age,:,:,:,:)))))/sum(sum(sum(sum(DBN_bench(age,:,:,:,:)))));
        Capital_age(age,2)           = sum(sum(sum(sum(Capital_exp(age,:,:,:,:).*DBN_exp(age,:,:,:,:)))))/sum(sum(sum(sum(DBN_exp(age,:,:,:,:)))));
        Borrowing_age(age,1)         = sum(sum(sum(sum(Borrowing_bench(age,:,:,:,:).*DBN_bench(age,:,:,:,:)))))/sum(sum(sum(sum(DBN_bench(age,:,:,:,:)))));
        Borrowing_age(age,2)         = sum(sum(sum(sum(Borrowing_exp(age,:,:,:,:).*DBN_exp(age,:,:,:,:)))))/sum(sum(sum(sum(DBN_exp(age,:,:,:,:)))));
        for z = 1:n_z
            Constrained_firms_AZ(age,z,1) = sum(sum(sum(constrained_ind_bench(age,:,z,:,:).*DBN_bench(age,:,z,:,:))))/sum(sum(sum(DBN_bench(age,:,z,:,:))));
            Constrained_firms_AZ(age,z,2) = sum(sum(sum(constrained_ind_exp(age,:,z,:,:).*DBN_exp(age,:,z,:,:))))/sum(sum(sum(DBN_exp(age,:,z,:,:))));
            Capital_AZ(age,z,1)           = sum(sum(sum(Capital_bench(age,:,z,:,:).*DBN_bench(age,:,z,:,:))))/sum(sum(sum(DBN_bench(age,:,z,:,:))));
            Capital_AZ(age,z,2)           = sum(sum(sum(Capital_exp(age,:,z,:,:).*DBN_exp(age,:,z,:,:))))/sum(sum(sum(DBN_exp(age,:,z,:,:))));
            Borrowing_AZ(age,z,1)         = sum(sum(sum(Borrowing_bench(age,:,z,:,:).*DBN_bench(age,:,z,:,:))))/sum(sum(sum(DBN_bench(age,:,z,:,:))));
            Borrowing_AZ(age,z,2)         = sum(sum(sum(Borrowing_exp(age,:,z,:,:).*DBN_exp(age,:,z,:,:))))/sum(sum(sum(DBN_exp(age,:,z,:,:))));
        end 
    end 
    
    
    for age=1:7
        DBN_bench_aux       = DBN_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)             ;
        DBN_exp_aux         = DBN_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)               ;
        cons_ind_bench_aux  = constrained_ind_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:) ;
        cons_ind_exp_aux    = constrained_ind_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)   ;
        Capital_bench_aux   = Capital_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)         ;
        Capital_exp_aux     = Capital_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)           ;
        Borrowing_bench_aux = Borrowing_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)       ;
        Borrowing_exp_aux   = Borrowing_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)         ;
        Assets_bench_aux    = A_mat(age_limit(age)+1:age_limit(age+1),:,:,:,:)                 ;
        Assets_exp_aux      = A_mat(age_limit(age)+1:age_limit(age+1),:,:,:,:)                 ;
        size_bench          = 100*sum(sum(sum(sum(sum(sum( DBN_bench_aux )))))) ;
        size_exp            = 100*sum(sum(sum(sum(sum(sum( DBN_exp_aux ))))))   ;
         
        Constrained_firms_A_group(age,1) = sum(sum(sum(sum(sum(cons_ind_bench_aux.*DBN_bench_aux)))))/size_bench ;
        Constrained_firms_A_group(age,2) = sum(sum(sum(sum(sum(cons_ind_exp_aux.*DBN_exp_aux)))))/size_exp       ;
        Capital_A_group(age,1)           = sum(sum(sum(sum(sum(Capital_bench_aux.*DBN_bench_aux)))))/size_bench  ;
        Capital_A_group(age,2)           = sum(sum(sum(sum(sum(Capital_exp_aux.*DBN_exp_aux)))))/size_exp        ;
        Borrowing_A_group(age,1)         = sum(sum(sum(sum(sum(Borrowing_bench_aux.*DBN_bench_aux)))))/size_bench;
        Borrowing_A_group(age,2)         = sum(sum(sum(sum(sum(Borrowing_exp_aux.*DBN_exp_aux)))))/size_exp      ;
        Assets_A_group(age,1)            = sum(sum(sum(sum(sum(Assets_bench_aux.*DBN_bench_aux)))))/size_bench   ;
        Assets_A_group(age,2)            = sum(sum(sum(sum(sum(Assets_exp_aux.*DBN_exp_aux)))))/size_exp         ;
        A_share_A_group(age,1)           = 100*sum(sum(sum(sum(sum(Assets_exp_aux.*DBN_bench_aux)))))/Assets_tot(1)  ;
        A_share_A_group(age,2)           = 100*sum(sum(sum(sum(sum(Assets_exp_aux.*DBN_exp_aux)))))/Assets_tot(2)    ;
        size_A_group(age,1)              = 100*size_bench                                                            ;
        size_A_group(age,2)              = 100*size_exp                                                              ;
        for z=1:n_z
            DBN_bench_aux       = DBN_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)             ;
            DBN_exp_aux         = DBN_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)               ;
            cons_ind_bench_aux  = constrained_ind_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:) ;
            cons_ind_exp_aux    = constrained_ind_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)   ;
            Capital_bench_aux   = Capital_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)         ;
            Capital_exp_aux     = Capital_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)           ;
            Borrowing_bench_aux = Borrowing_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)       ;
            Borrowing_exp_aux   = Borrowing_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)         ;
            Assets_bench_aux    = A_mat(age_limit(age)+1:age_limit(age+1),:,z,:,:)                 ;
            Assets_exp_aux      = A_mat(age_limit(age)+1:age_limit(age+1),:,z,:,:)                 ;
            size_bench          = sum(sum(sum(sum(sum( DBN_bench_aux ))))) ;
            size_exp            = sum(sum(sum(sum(sum( DBN_exp_aux )))))   ;

            Constrained_firms_AZ_group(age,z,1) = sum(sum(sum(sum(cons_ind_bench_aux.*DBN_bench_aux))))/size_bench  ;
            Constrained_firms_AZ_group(age,z,2) = sum(sum(sum(sum(cons_ind_exp_aux.*DBN_exp_aux))))/size_exp        ;
            Capital_AZ_group(age,z,1)           = sum(sum(sum(sum(Capital_bench_aux.*DBN_bench_aux))))/size_bench   ;
            Capital_AZ_group(age,z,2)           = sum(sum(sum(sum(Capital_exp_aux.*DBN_exp_aux))))/size_exp         ;
            Borrowing_AZ_group(age,z,1)         = sum(sum(sum(sum(Borrowing_bench_aux.*DBN_bench_aux))))/size_bench ;
            Borrowing_AZ_group(age,z,2)         = sum(sum(sum(sum(Borrowing_exp_aux.*DBN_exp_aux))))/size_exp       ;
            Assets_AZ_group(age,z,1)            = sum(sum(sum(sum(Assets_bench_aux.*DBN_bench_aux))))/size_bench    ;
            Assets_AZ_group(age,z,2)            = sum(sum(sum(sum(Assets_exp_aux.*DBN_exp_aux))))/size_exp          ;
            A_share_AZ_group(age,z,1)           = 100*sum(sum(sum(sum(Assets_exp_aux.*DBN_bench_aux))))/Assets_tot(1)   ;
            A_share_AZ_group(age,z,2)           = 100*sum(sum(sum(sum(Assets_exp_aux.*DBN_exp_aux))))/Assets_tot(2)     ;
            size_AZ_group(age,z,1)              = 100*size_bench                                                        ;
            size_AZ_group(age,z,2)              = 100*size_exp                                                          ;
            
            % Wealth Shares by Z
            size_Z(z,1)    = 100*sum(sum(sum(sum(DBN_bench(:,:,z,:,:))))) ;
            size_Z(z,2)    = 100*sum(sum(sum(sum(DBN_exp(:,:,z,:,:)))))   ;
            A_share_Z(z,1) = 100*sum(sum(sum(sum(A_mat(:,:,z,:,:).*DBN_bench(:,:,z,:,:)))))/Assets_tot(1) ;
            A_share_Z(z,2) = 100*sum(sum(sum(sum(A_mat(:,:,z,:,:).*DBN_exp(:,:,z,:,:)))))/Assets_tot(2)   ;
        end 
    end

    
    
    % Firms that were constrained and now are not
    ind  = (constrained_ind_bench==1)&(constrained_ind_exp==0) ;
    size = sum(sum(sum(sum(sum(ind.*DBN_bench)))));
    output_ba(1)    = sum(sum(sum(sum(sum(firm_output_bench.*ind.*DBN_bench)))))/size ;
    output_ba(2)    = sum(sum(sum(sum(sum(firm_output_exp.*ind.*DBN_bench)))))/size   ;
    profit_ba(1)    = sum(sum(sum(sum(sum(firm_profit_bench.*ind.*DBN_bench)))))/size ;
    profit_ba(2)    = sum(sum(sum(sum(sum(firm_profit_exp.*ind.*DBN_bench)))))/size   ;
    Capital_ba(1)   = sum(sum(sum(sum(sum(Capital_bench.*ind.*DBN_bench)))))/size ;
    Capital_ba(2)   = sum(sum(sum(sum(sum(Capital_exp.*ind.*DBN_bench)))))/size   ;
    Borrowing_ba(1) = sum(sum(sum(sum(sum(Borrowing_bench.*ind.*DBN_bench)))))/size ;
    Borrowing_ba(2) = sum(sum(sum(sum(sum(Borrowing_exp.*ind.*DBN_bench)))))/size   ;
              
    for age = 1:Max_Age
        size = sum(sum(sum(sum(ind(age,:,:,:,:).*DBN_bench(age,:,:,:,:))))) ;
        output_ba_age(age,1) = sum(sum(sum(sum(sum(firm_output_bench(age,:,:,:,:).*ind(age,:,:,:,:).*DBN_bench(age,:,:,:,:))))))/size ;
        output_ba_age(age,2) = sum(sum(sum(sum(sum(firm_output_exp(age,:,:,:,:).*ind(age,:,:,:,:).*DBN_bench(age,:,:,:,:))))))/size ;
        profit_ba_age(age,1) = sum(sum(sum(sum(sum(firm_profit_bench(age,:,:,:,:).*ind(age,:,:,:,:).*DBN_bench(age,:,:,:,:))))))/size ;
        profit_ba_age(age,2) = sum(sum(sum(sum(sum(firm_profit_exp(age,:,:,:,:).*ind(age,:,:,:,:).*DBN_bench(age,:,:,:,:))))))/size ;
        for z = 1:n_z
            size = sum(sum(sum(ind(age,:,z,:,:).*DBN_bench(age,:,z,:,:)))) ;
            output_ba_AZ(age,z,1) = sum(sum(sum(sum(firm_output_bench(age,:,z,:,:).*ind(age,:,z,:,:).*DBN_bench(age,:,z,:,:)))))/size ;
            output_ba_AZ(age,z,2) = sum(sum(sum(sum(firm_output_exp(age,:,z,:,:).*ind(age,:,z,:,:).*DBN_bench(age,:,z,:,:)))))/size   ;
            profit_ba_AZ(age,z,1) = sum(sum(sum(sum(firm_profit_bench(age,:,z,:,:).*ind(age,:,z,:,:).*DBN_bench(age,:,z,:,:)))))/size ;
            profit_ba_AZ(age,z,2) = sum(sum(sum(sum(firm_profit_exp(age,:,z,:,:).*ind(age,:,z,:,:).*DBN_bench(age,:,z,:,:)))))/size   ;
        end 
    end 

    for age=1:7
        ind_aux             = ind(age_limit(age)+1:age_limit(age+1),:,:,:,:)                   ;
        DBN_bench_aux       = DBN_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)             ;
        DBN_exp_aux         = DBN_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)               ;
        profit_bench_aux    = firm_profit_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)     ;
        profit_exp_aux      = firm_profit_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)       ;
        output_bench_aux    = firm_output_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)     ;
        output_exp_aux      = firm_output_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)       ;
        Capital_bench_aux   = Capital_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)         ;
        Capital_exp_aux     = Capital_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)           ;
        Borrowing_bench_aux = Borrowing_bench(age_limit(age)+1:age_limit(age+1),:,:,:,:)       ;
        Borrowing_exp_aux   = Borrowing_exp(age_limit(age)+1:age_limit(age+1),:,:,:,:)         ;
        size_bench          = sum(sum(sum(sum(sum(sum( ind_aux.*DBN_bench_aux )))))) ;
        size_exp            = sum(sum(sum(sum(sum(sum( ind_aux.*DBN_exp_aux ))))))   ;
        
        output_A_group(age,1)    = sum(sum(sum(sum(sum(sum(output_bench_aux.*ind_aux.*DBN_bench_aux))))))/size_bench ;
        output_A_group(age,2)    = sum(sum(sum(sum(sum(sum(output_exp_aux.*ind_aux.*DBN_bench_aux))))))/size_exp     ;
        profit_A_group(age,1)    = sum(sum(sum(sum(sum(sum(profit_bench_aux.*ind_aux.*DBN_bench_aux))))))/size_bench ;
        profit_A_group(age,2)    = sum(sum(sum(sum(sum(sum(profit_exp_aux.*ind_aux.*DBN_bench_aux))))))/size_exp     ;
        Capital_ba_A_group(age,1)   = sum(sum(sum(sum(sum(sum(Capital_bench_aux.*ind_aux.*DBN_bench_aux))))))/size_bench ;
        Capital_ba_A_group(age,2)   = sum(sum(sum(sum(sum(sum(Capital_exp_aux.*ind_aux.*DBN_bench_aux))))))/size_exp     ;
        Borrowing_ba_A_group(age,1) = sum(sum(sum(sum(sum(sum(Borrowing_bench_aux.*ind_aux.*DBN_bench_aux))))))/size_bench ;
        Borrowing_ba_A_group(age,2) = sum(sum(sum(sum(sum(sum(Borrowing_exp_aux.*ind_aux.*DBN_bench_aux))))))/size_exp     ;
        
        for z=1:n_z
            ind_aux            = ind(age_limit(age)+1:age_limit(age+1),:,z,:,:)                   ;
            DBN_bench_aux      = DBN_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)             ;
            DBN_exp_aux        = DBN_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)               ;
            profit_bench_aux   = firm_profit_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)     ;
            profit_exp_aux     = firm_profit_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)       ;
            output_bench_aux   = firm_output_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)     ;
            output_exp_aux     = firm_output_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)       ;
            Capital_bench_aux  = Capital_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)         ;
            Capital_exp_aux    = Capital_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)           ;
            Borrowing_bench_aux  = Borrowing_bench(age_limit(age)+1:age_limit(age+1),:,z,:,:)     ;
            Borrowing_exp_aux    = Borrowing_exp(age_limit(age)+1:age_limit(age+1),:,z,:,:)       ;
            size_bench           = sum(sum(sum(sum(sum(sum( ind_aux.*DBN_bench_aux )))))) ;
            size_exp             = sum(sum(sum(sum(sum(sum( ind_aux.*DBN_exp_aux ))))))   ;

            output_AZ_group(age,z,1) = sum(sum(sum(sum(sum(output_bench_aux.*ind_aux.*DBN_bench_aux)))))/size_bench ;
            output_AZ_group(age,z,2) = sum(sum(sum(sum(sum(output_exp_aux.*ind_aux.*DBN_bench_aux)))))/size_exp     ;
            profit_AZ_group(age,z,1) = sum(sum(sum(sum(sum(profit_bench_aux.*ind_aux.*DBN_bench_aux)))))/size_bench ;
            profit_AZ_group(age,z,2) = sum(sum(sum(sum(sum(profit_exp_aux.*ind_aux.*DBN_bench_aux)))))/size_exp     ;
            Capital_ba_AZ_group(age,z,1) = sum(sum(sum(sum(sum(Capital_bench_aux.*ind_aux.*DBN_bench_aux)))))/size_bench ;
            Capital_ba_AZ_group(age,z,2) = sum(sum(sum(sum(sum(Capital_exp_aux.*ind_aux.*DBN_bench_aux)))))/size_exp     ;
            Borrowing_ba_AZ_group(age,z,1) = sum(sum(sum(sum(sum(Borrowing_bench_aux.*ind_aux.*DBN_bench_aux)))))/size_bench ;
            Borrowing_ba_AZ_group(age,z,2) = sum(sum(sum(sum(sum(Borrowing_exp_aux.*ind_aux.*DBN_bench_aux)))))/size_exp     ;
        end
    end

    
    % Excel file
        Blank = {' ',' ',' ',' ',' ',' ',' ',' ',' '} ;
        % Constrained firms
        Mat_bench   = [Constrained_firms_A_group(:,1) Constrained_firms_AZ_group(:,:,1)] ;
        Mat_exp     = [Constrained_firms_A_group(:,2) Constrained_firms_AZ_group(:,:,2)] ;
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over all population'},num2cell(Constrained_firms_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over all population'},num2cell(Constrained_firms_tot(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over all population'},num2cell(Constrained_firms_tot(2)-Constrained_firms_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Percentage of firms constrained by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank; Title;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank;
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank;
               diff_title  ; col_title ; row_title , num2cell(Mat_exp-Mat_bench) ; diff_end ] 
        status = xlwrite(Tables_file,Mat,'Constrained Firms (all)') ;
        
        Mat_bench   = [Capital_A_group(:,1) Capital_AZ_group(:,:,1)] ;
        Mat_exp     = [Capital_A_group(:,2) Capital_AZ_group(:,:,2)] ;
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over all population'},num2cell(Capital_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over all population'},num2cell(Capital_tot(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over all population'},num2cell(Capital_tot(2)-Capital_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over all population'},num2cell(100*(Capital_tot(2)./Capital_tot(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average Capital Demand by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank;
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank;
               diff_title  ; col_title ; row_title , num2cell(Mat_exp-Mat_bench) ; diff_end ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(100*(Mat_exp./Mat_bench-1)) ; diff_pc_end] 
        status = xlwrite(Tables_file,Mat,'Capital Firms (all)') ;
        
        Mat_bench   = [Borrowing_A_group(:,1) Borrowing_AZ_group(:,:,1)] ;
        Mat_exp     = [Borrowing_A_group(:,2) Borrowing_AZ_group(:,:,2)] ;
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over all population'},num2cell(Borrowing_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over all population'},num2cell(Borrowing_tot(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over all population'},num2cell(Borrowing_tot(2)-Borrowing_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over all population'},num2cell(100*(Borrowing_tot(2)./Borrowing_tot(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average Borrowing by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank;
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank;
               diff_title  ; col_title ; row_title , num2cell(Mat_exp-Mat_bench) ; diff_end ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(100*(Mat_exp./Mat_bench-1)) ; diff_pc_end] 
        status = xlwrite(Tables_file,Mat,'Borrowing Firms (all)') ;
        
        Mat_bench   = [Assets_A_group(:,1) Assets_AZ_group(:,:,1)] ;
        Mat_exp     = [Assets_A_group(:,2) Assets_AZ_group(:,:,2)] ;
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over all population'},num2cell(Assets_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over all population'},num2cell(Assets_tot(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over all population'},num2cell(Assets_tot(2)-Assets_tot(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over all population'},num2cell(100*(Assets_tot(2)./Assets_tot(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average Assets by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank;
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank;
               diff_title  ; col_title ; row_title , num2cell(Mat_exp-Mat_bench) ; diff_end ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(100*(Mat_exp./Mat_bench-1)) ; diff_pc_end] 
        status = xlwrite(Tables_file,Mat,'Assets Firms (all)') ;
        
        
        Mat_bench   = [A_share_A_group(:,1) A_share_AZ_group(:,:,1)] ;
        Mat_exp     = [A_share_A_group(:,2) A_share_AZ_group(:,:,2)] ;
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over all population'},num2cell([100 A_share_Z(:,1)'])] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end   = [{'Average over all population'},num2cell([100 A_share_Z(:,2)'])] ;
        diff_title  = {'Difference (Wealth-Capital',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over all population'},num2cell([0, A_share_Z(:,2)'-A_share_Z(:,1)'])] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Wealth Shares by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank;
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank;
               diff_title  ; col_title ; row_title , num2cell(Mat_exp-Mat_bench) ; diff_end ] 
        status = xlwrite(Tables_file,Mat,'Wealth Shares AZ') ;
        
        Mat_bench   = [size_A_group(:,1) size_AZ_group(:,:,1)] ;
        Mat_exp     = [size_A_group(:,2) size_AZ_group(:,:,2)] ;
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over all population'},num2cell([100 size_Z(:,1)'])] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end   = [{'Average over all population'},num2cell([100 size_Z(:,2)'])] ;
        diff_title  = {'Difference (Wealth-Capital',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over all population'},num2cell([0, size_Z(:,2)'-size_Z(:,1)'])] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Population Shares by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank;
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank;
               diff_title  ; col_title ; row_title , num2cell(Mat_exp-Mat_bench) ; diff_end ] 
        status = xlwrite(Tables_file,Mat,'Population Shares AZ') ;
        
        % Output and profit Change
        Mat_bench   = [output_A_group(:,1) output_AZ_group(:,:,1)] ;
        Mat_exp     = [output_A_group(:,2) output_AZ_group(:,:,2)] ;
        Mat_diff    = 100*(Mat_exp./Mat_bench-1)                   ; 
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over sample'},num2cell(output_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over sample'},num2cell(output_ba(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over sample'},num2cell(output_ba(2)-output_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over sample'},num2cell(100*(output_ba(2)/output_ba(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average output by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Title_aux   = {'Only firms that were constrained under Capital_Tax that are not constrained under Wealth_Tax',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title; Title_aux ;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank ; 
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank ; 
               diff_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_end  ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_pc_end  ] 
        status = xlwrite(Tables_file,Mat,'Output (Cons-Unc)') ;
        
        Mat_bench   = [profit_A_group(:,1) profit_AZ_group(:,:,1)] ;
        Mat_exp     = [profit_A_group(:,2) profit_AZ_group(:,:,2)] ;
        Mat_diff    = 100*(Mat_exp./Mat_bench-1)                   ; 
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over sample'},num2cell(profit_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over sample'},num2cell(profit_ba(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over sample'},num2cell(profit_ba(2)-profit_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over sample'},num2cell(100*(profit_ba(2)/profit_ba(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average profit by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Title_aux   = {'Only firms that were constrained under Capital_Tax that are not constrained under Wealth_Tax',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title; Title_aux ;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank ; 
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank ; 
               diff_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_end  ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_pc_end  ] 
        status = xlwrite(Tables_file,Mat,'Profits (Cons-Unc)') ;
        
        Mat_bench   = [Capital_ba_A_group(:,1) Capital_ba_AZ_group(:,:,1)] ;
        Mat_exp     = [Capital_ba_A_group(:,2) Capital_ba_AZ_group(:,:,2)] ;
        Mat_diff    = 100*(Mat_exp./Mat_bench-1)                   ; 
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over sample'},num2cell(Capital_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over sample'},num2cell(Capital_ba(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over sample'},num2cell(Capital_ba(2)-Capital_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Capital)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over sample'},num2cell(100*(Capital_ba(2)/Capital_ba(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average Capital_ba by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Title_aux   = {'Only firms that were constrained under Capital_Tax that are not constrained under Wealth_Tax',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title; Title_aux ;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank ; 
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank ; 
               diff_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_end  ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_pc_end  ] 
        status = xlwrite(Tables_file,Mat,'Capital (Cons-Unc)') ;
        
        Mat_bench   = [Borrowing_ba_A_group(:,1) Borrowing_ba_AZ_group(:,:,1)] ;
        Mat_exp     = [Borrowing_ba_A_group(:,2) Borrowing_ba_AZ_group(:,:,2)] ;
        Mat_diff    = 100*(Mat_exp./Mat_bench-1)                   ; 
        bench_title = {'Capital_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        bench_end   = [{'Average over sample'},num2cell(Borrowing_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        exp_title   = {'Wealth_Tax Economy',' ',' ',' ',' ',' ',' ',' ',' '} ;
        exp_end     = [{'Average over sample'},num2cell(Borrowing_ba(2)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_title  = {'Difference (Wealth-Borrowing)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_end    = [{'Average over sample'},num2cell(Borrowing_ba(2)-Borrowing_ba(1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        diff_pc_title  = {'% Difference (Wealth/Borrowing)',' ',' ',' ',' ',' ',' ',' ',' '} ;
        diff_pc_end    = [{'Average over sample'},num2cell(100*(Borrowing_ba(2)/Borrowing_ba(1)-1)),{' ',' ',' ',' ',' ',' ',' '}] ;
        col_title   = {'Age_Group','Total','z1','z2','z3','z4','z5','z6','z7'} ;
        row_title   = {'20-24';'25-34';'35-44';'45-54';'55-64';'65-74';'75-80'};
        Title       = {'Average Borrowing_ba by age-group and Z',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Title_aux   = {'Only firms that were constrained under Borrowing_Tax that are not constrained under Wealth_Tax',' ',' ',' ',' ',' ',' ',' ',' '} ;
        Mat = [Blank ; Title; Title_aux ;
               bench_title ; col_title ; row_title , num2cell(Mat_bench) ; bench_end ; Blank ; 
               exp_title   ; col_title ; row_title , num2cell(Mat_exp)   ; exp_end   ; Blank ; 
               diff_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_end  ; Blank ;
               diff_pc_title  ; col_title ; row_title , num2cell(Mat_diff)  ; diff_pc_end  ] 
        status = xlwrite(Tables_file,Mat,'Borrowing (Cons-Unc)') ;
        
    % Figures
        % Constrained firms
        figure;
        plot([20:100]',Constrained_firms_age,'linewidth',2);
        legend('Capital Tax','Wealth_Tax','location','southeast'); 
        xlabel('age'); xlim([20,100]); ylim([0 1]); title('Percentage of Constrained Firms')
        print('-depsc','Constrained_Firms_Age.eps') ;
        
        figure;
        for z=1:n_z
            subplot(2,4,z); plot([20:100]',[Constrained_firms_AZ(:,z,1),Constrained_firms_AZ(:,z,2)],'linewidth',2);
            aa = ['Cons. Firms (Z',num2str(z),')'] ; title(aa); xlabel('age'); xlim([20,100]); ylim([0 1]);
        end 
            legend('K Tax','W Tax','location','south');
        print('-depsc','Constrained_Firms_Age_Z.eps') ;



