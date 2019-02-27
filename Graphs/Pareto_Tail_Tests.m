% Pareto Tail tests
% Sergio Ocampo Diaz

%% Clear and load required files
   % Clear and close
        close all;
        clear; clc;

   % Load files for excel in mac
        javaaddpath('poi_library/poi-3.8-20120326.jar');
        javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
        javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
        javaaddpath('poi_library/xmlbeans-2.3.0.jar');
        javaaddpath('poi_library/dom4j-1.6.1.jar');
        javaaddpath('poi_library/stax-api-1.0.1.jar');
        
        cd '/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/GKKOC_CODES/Sergio/Graphs/Pareto_Tail/'

        
%% Parameters 
    theta = 1.5 ;
    Threshold_Factor = 0.0;
    
    Progressive_Tax_Switch = 0 ;
    NSU_Switch = 1 ;

        
%% Fixed Parameters

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
    age_limit = [0, 5, 15, 25, 35, 45, 55, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;
    
% Percentiles
    prctl = [10, 25, 50, 75, 90, 99, 99.9];
    
% Other parameters
    if theta==1.0
        Params =[0.9415,  0.00,  0.50,  0.65,  0.34,  0.4494]; % tauL=0.224, tauC=0.075 calibration
    elseif theta==1.50
        if mu==0.85
        Params= [0.945, 0.00, 0.50, 0.7889, 0.34, 0.4494]; % mu=0.85 calibration, targetting 0.34, 0.69, vartheta1.5 
        else 
        Params= [0.9412, 0.0, 0.50, 0.640, 0.34, 0.4494]; % mu=0.9 calibration, targetting 0.34, 0.69, vartheta1.5
        end
    elseif theta==1.60
        Params= [0.9409, 0.0, 0.50, 0.640, 0.34, 0.4494]; % mu=0.9 calibration, targetting 0.34, 0.69, vartheta1.6
    elseif theta==2.00
        Params= [0.9405, 0.0, 0.50, 0.639, 0.34, 0.4494]; % mu=0.9 calibration, targetting 0.34, 0.69, vartheta2
    elseif theta==2.50
        Params= [0.9400, 0.0, 0.50, 0.639, 0.34, 0.4494]; % mu=0.9 calibration, targetting 0.34, 0.69, vartheta2.5
    else
        disp('No parameters for this theta, changing to default parameters (theta=1)')
        Params =[0.9436, 0.0, 0.50, 0.70444445, 0.34, 0.4494]; % tauL=0.224, tauC=0.075 calibration
    end

    beta   = Params(1) ;
    mu_z   = Params(2) ;
    rho_z  = Params(3) ;
    sigma_z_pdf      = Params(4) ;
    sigma_lambda_pdf = Params(5) ;
    gamma  = Params(6) ;

%% Excel file
    Tables_file = 'Tables_Pareto_Tail.xls' ;
    

%% Benchmark
    % Result_Folders
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        n_z=7;
    
    % Pareto Tail
        eval(['load ',Simul_Folder,'panela_bench']) ; 
        eval(['load ',Simul_Folder,'panela_exp'])   ;
        eval(['load ',Result_Folder,'agrid'])                      ;
        eval(['load ',Bench_Folder,'EBAR'])                        ;
        eval(['load ',Result_Folder,'Exp_Files/Exp_results_EBAR']) ;
        eval(['load ',Bench_Folder,'DBN'])            ;
        DBN_bench   = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ; clear DBN;
        DBN_a_bench = sum(sum(sum(sum(DBN_bench,5),4),3),1)  ; CCDF_bench = cumsum(DBN_a_bench,'reverse') ;
        n_w = 1000 ;
        N = numel(panela_bench) ;

        Wage_earnings = 47000;
        wealth_bench = sort(panela_bench) ; wealth_bench = Wage_earnings/EBAR * wealth_bench ;
        wealth_exp   = sort(panela_exp)   ; wealth_exp   = Wage_earnings/EBAR * wealth_exp   ;
        Exp_results_EBAR = Wage_earnings/EBAR * Exp_results_EBAR ; EBAR =  Wage_earnings ;
        agrid = Wage_earnings/EBAR * agrid ;
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
            ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
            set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
            file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_bench.fig'] ;
            file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_bench.pdf'] ;
            savefig(file_name); print('-dpdf',file_name_pdf) ;
                % Distribution graph
                [~,iii] = min(CCDF_bench>=ind(ii)/100) ; if ii==1; iii=1; end;
                DBN_aux = DBN_a_bench(iii:end)/sum(DBN_a_bench(iii:end)); CCDF_aux = cumsum(DBN_aux,'reverse') ; agrid_aux = agrid(iii:end) ;
                x_dbn = -log(agrid_aux/agrid_aux(1)) ; y_dbn = log(CCDF_aux) ;
                fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% bench - DBN'] ;
                figure; hold on; plot(-x,y,'k','linewidth',2); plot(-x_dbn,y_dbn,'linewidth',2); hold off;
                title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
                ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
                set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
                file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_DBN_bench.fig'] ;
                file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_DBN_bench.pdf'] ;
                savefig(file_name); print('-dpdf',file_name_pdf) ; 
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
                    EBAR,NaN,NaN,NaN,NaN,NaN];
        AA_exp   = [100 50 25 10 05 01;
                    alpha_exp;
                    SE_exp;
                    R2_exp;
                    alpha_exp_ml;
                    wealth_exp(nn);
                    wealth_exp(end),NaN,NaN,NaN,NaN,NaN;
                    Exp_results_EBAR,NaN,NaN,NaN,NaN,NaN];

        bench_name = {'Bench',' ',' ',' ',' ',' ',' '} ;
        Exp_name = {'Exp',' ',' ',' ',' ',' ',' '} ;
        blank_row = cell(1,7) ;
        row_name = {'top_x%';'Pareto_ind_reg';'SE';'R2';'Pareto_ind_ml';'Min_Wealth';'Max_Wealth';'EBAR'} ;
        Mat = [bench_name ; row_name num2cell(AA_bench) ; blank_row ; Exp_name; row_name num2cell(AA_exp) ]
        status = xlwrite(Tables_file,Mat,'Test_bench') ;

%% mu=0.95    
    % Result_Folders
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/mu_95/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/mu_95/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        n_z=7;

    % Pareto Tail
        eval(['load ',Simul_Folder,'panela_bench']) ; 
        eval(['load ',Simul_Folder,'panela_exp'])   ;
        eval(['load ',Result_Folder,'agrid'])                      ;
        eval(['load ',Bench_Folder,'EBAR'])                        ; EBAR_bench = EBAR;
        eval(['load ',Result_Folder,'Exp_Files/Exp_results_EBAR']) ;
        eval(['load ',Result_Folder,'Bench_Files/DBN'])            ;
        eval(['load ',Result_Folder,'Bench_Files/EBAR'])           ;
        DBN_bench   = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ; clear DBN;
        DBN_a_bench = sum(sum(sum(sum(DBN_bench,5),4),3),1)  ; CCDF_bench = cumsum(DBN_a_bench,'reverse') ;
        n_w = 1000 ;
        N = numel(panela_bench) ;

        Wage_earnings = 47000;
        wealth_bench = sort(panela_bench) ; wealth_bench = Wage_earnings/EBAR_bench * wealth_bench ;
        wealth_exp   = sort(panela_exp)   ; wealth_exp   = Wage_earnings/EBAR_bench * wealth_exp   ;
        Exp_results_EBAR = Wage_earnings/EBAR_bench * Exp_results_EBAR ; EBAR =  Wage_earnings/EBAR_bench * EBAR ;
        agrid = Wage_earnings/EBAR * agrid ;
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
            fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% mu=0.95'] ;
            figure; hold on; plot(-x,y,'linewidth',2); plot(xx,cons_bench(ii)-alpha_bench(ii)*xx,'-k'); hold off;
            title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
            ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
            set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
            file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_mu_95.fig'] ;
            file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_mu_95.pdf'] ;
            savefig(file_name); print('-dpdf',file_name_pdf) ;
                % Distribution graph
                [~,iii] = min(CCDF_bench>=ind(ii)/100) ; if ii==1; iii=1; end;
                DBN_aux = DBN_a_bench(iii:end)/sum(DBN_a_bench(iii:end)); CCDF_aux = cumsum(DBN_aux,'reverse') ; agrid_aux = agrid(iii:end) ;
                x_dbn = -log(agrid_aux/agrid_aux(1)) ; y_dbn = log(CCDF_aux) ;
                fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% mu=0.95 - DBN'] ;
                figure; hold on; plot(-x,y,'k','linewidth',2); plot(-x_dbn,y_dbn,'linewidth',2); hold off;
                title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
                ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
                set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
                file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_DBN_mu_95.fig'] ;
                file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_DBN_mu_95.pdf'] ;
                savefig(file_name); print('-dpdf',file_name_pdf) ; 
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
                    EBAR,NaN,NaN,NaN,NaN,NaN];
        AA_exp   = [100 50 25 10 05 01;
                    alpha_exp;
                    SE_exp;
                    R2_exp;
                    alpha_exp_ml;
                    wealth_exp(nn);
                    wealth_exp(end),NaN,NaN,NaN,NaN,NaN;
                    Exp_results_EBAR,NaN,NaN,NaN,NaN,NaN];

        bench_name = {'Bench',' ',' ',' ',' ',' ',' '} ;
        Exp_name = {'Exp',' ',' ',' ',' ',' ',' '} ;
        blank_row = cell(1,7) ;
        row_name = {'top_x%';'Pareto_ind_reg';'SE';'R2';'Pareto_ind_ml';'Min_Wealth';'Max_Wealth';'EBAR'} ;
        Mat = [bench_name ; row_name num2cell(AA_bench) ; blank_row ; Exp_name; row_name num2cell(AA_exp) ]
        status = xlwrite(Tables_file,Mat,'Test_mu_95') ;

%% rho=0.75    
    % Result_Folders
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/rho_75/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/rho_75/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;

    % Pareto Tail
        eval(['load ',Simul_Folder,'panela_bench']) ; 
        eval(['load ',Simul_Folder,'panela_exp'])   ;
        eval(['load ',Bench_Folder,'EBAR'])                        ; EBAR_bench = EBAR;
        eval(['load ',Result_Folder,'Exp_Files/Exp_results_EBAR']) ;
        n_w = 1000 ;
        N = numel(panela_bench) ;

        Wage_earnings = 47000 ;wealth_bench = sort(panela_bench) ; wealth_bench = Wage_earnings/EBAR_bench * wealth_bench ;
        wealth_exp   = sort(panela_exp)   ; wealth_exp   = Wage_earnings/EBAR_bench * wealth_exp   ;
        Exp_results_EBAR = Wage_earnings/EBAR_bench * Exp_results_EBAR ; EBAR =  Wage_earnings/EBAR_bench * EBAR ;
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
            fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% rho=0.75'] ;
            figure; hold on; plot(-x,y,'linewidth',2); plot(xx,cons_bench(ii)-alpha_bench(ii)*xx,'-k'); hold off;
            title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
            ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
            set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
            file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_rho_75.fig'] ;
            file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_rho_75.pdf'] ;
            savefig(file_name); print('-dpdf',file_name_pdf) ;
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
                    EBAR,NaN,NaN,NaN,NaN,NaN];
        AA_exp   = [100 50 25 10 05 01;
                    alpha_exp;
                    SE_exp;
                    R2_exp;
                    alpha_exp_ml;
                    wealth_exp(nn);
                    wealth_exp(end),NaN,NaN,NaN,NaN,NaN;
                    Exp_results_EBAR,NaN,NaN,NaN,NaN,NaN];

        bench_name = {'Bench',' ',' ',' ',' ',' ',' '} ;
        Exp_name = {'Exp',' ',' ',' ',' ',' ',' '} ;
        blank_row = cell(1,7) ;
        row_name = {'top_x%';'Pareto_ind_reg';'SE';'R2';'Pareto_ind_ml';'Min_Wealth';'Max_Wealth';'EBAR'} ;
        Mat = [bench_name ; row_name num2cell(AA_bench) ; blank_row ; Exp_name; row_name num2cell(AA_exp) ]
        status = xlwrite(Tables_file,Mat,'Test_rho_75') ;

        
%% Z_grid        
    % Result_Folders
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/z_grid/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/z_grid/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        n_z=9;

    % Pareto Tail
        eval(['load ',Simul_Folder,'panela_bench']) ; 
        eval(['load ',Simul_Folder,'panela_exp'])   ;
        eval(['load ',Result_Folder,'agrid'])                      ;
        eval(['load ',Bench_Folder,'EBAR'])                        ; EBAR_bench = EBAR;
        eval(['load ',Result_Folder,'Exp_Files/Exp_results_EBAR']) ;
        eval(['load ',Result_Folder,'Bench_Files/DBN'])            ;
        eval(['load ',Result_Folder,'Bench_Files/EBAR'])           ;
        DBN_bench   = reshape(DBN,[Max_Age,n_a,n_z,n_l,n_e]) ; clear DBN;
        DBN_a_bench = sum(sum(sum(sum(DBN_bench,5),4),3),1)  ; CCDF_bench = 1-cumsum(DBN_a_bench) ;
        n_w = 1000 ;
        N = numel(panela_bench) ;

        Wage_earnings = 47000 ;
        wealth_bench = sort(panela_bench) ; wealth_bench = Wage_earnings/EBAR_bench * wealth_bench ;
        wealth_exp   = sort(panela_exp)   ; wealth_exp   = Wage_earnings/EBAR_bench * wealth_exp   ;
        Exp_results_EBAR = Wage_earnings/EBAR_bench * Exp_results_EBAR ; EBAR =  Wage_earnings/EBAR_bench * EBAR ;
        agrid = Wage_earnings/EBAR * agrid ;
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
            fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% z_grid'] ;
            figure; hold on; plot(-x,y,'linewidth',2); plot(xx,cons_bench(ii)-alpha_bench(ii)*xx,'-k'); hold off;
            title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
            ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
            set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
            file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_z_grid.fig'] ;
            file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_z_grid.pdf'] ;
            savefig(file_name); print('-dpdf',file_name_pdf) ;
                % Distribution graph
                [~,iii] = min(CCDF_bench>=ind(ii)/100) ; if ii==1; iii=1; end;
                DBN_aux = DBN_a_bench(iii:end)/sum(DBN_a_bench(iii:end)); CCDF_aux = cumsum(DBN_aux,'reverse') ; agrid_aux = agrid(iii:end) ;
                x_dbn = -log(agrid_aux/agrid_aux(1)) ; y_dbn = log(CCDF_aux) ;
                fig_title = ['Pareto Tail Top ',num2str(ind(ii)),'% z_grid - DBN'] ;
                figure; hold on; plot(-x,y,'k','linewidth',2); plot(-x_dbn,y_dbn,'linewidth',2); hold off;
                title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([min(-x),max(-x)]);
                ww = linspace(log(min(wealth_bench(nn(ii):end)/w_min)),log(max(wealth_bench(nn(ii):end))/w_min),5) ;
                set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww)*w_min)
                file_name = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_DBN_z_grid.fig'] ;
                file_name_pdf = ['Pareto_Tail_Top_',num2str(ind(ii)),'%_DBN_z_grid.pdf'] ;
                savefig(file_name); print('-dpdf',file_name_pdf) ;     
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
                    EBAR,NaN,NaN,NaN,NaN,NaN];
        AA_exp   = [100 50 25 10 05 01;
                    alpha_exp;
                    SE_exp;
                    R2_exp;
                    alpha_exp_ml;
                    wealth_exp(nn);
                    wealth_exp(end),NaN,NaN,NaN,NaN,NaN;
                    Exp_results_EBAR,NaN,NaN,NaN,NaN,NaN];

        bench_name = {'Bench',' ',' ',' ',' ',' ',' '} ;
        Exp_name = {'Exp',' ',' ',' ',' ',' ',' '} ;
        blank_row = cell(1,7) ;
        row_name = {'top_x%';'Pareto_ind_reg';'SE';'R2';'Pareto_ind_ml';'Min_Wealth';'Max_Wealth';'EBAR'} ;
        Mat = [bench_name ; row_name num2cell(AA_bench) ; blank_row ; Exp_name; row_name num2cell(AA_exp) ]
        status = xlwrite(Tables_file,Mat,'Test_z_grid') ;

        
%% Vermuelen Data
    load /Users/s-ocampo/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Graphs/Vermeulen_Data/USA_IM1.dat
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
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        eval(['load ',Simul_Folder,'panela_bench']) ; panel_bench = panela_bench; 
        eval(['load ',Bench_Folder,'EBAR'])         ; EBAR_bench  = EBAR        ;
	% mu=0.95 Files
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/mu_95/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/mu_95/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        eval(['load ',Simul_Folder,'panela_bench'])      ; panel_mu = panela_bench; 
        eval(['load ',Result_Folder,'Bench_Files/EBAR']) ; EBAR_mu  = EBAR        ;
	% zgrid Files
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/z_grid/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/z_grid/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        eval(['load ',Simul_Folder,'panela_bench'])      ; panel_z = panela_bench; 
        eval(['load ',Result_Folder,'Bench_Files/EBAR']) ; EBAR_z  = EBAR        ;
          
    % Pareto Tail
        Wage_earnings = 47000 ;
        wealth_bench  = sort(panel_bench) ; wealth_bench = Wage_earnings/EBAR_bench * wealth_bench ;
        wealth_mu     = sort(panel_mu)    ; wealth_mu    = Wage_earnings/EBAR_bench * wealth_mu    ;
        wealth_z      = sort(panel_z)     ; wealth_z     = Wage_earnings/EBAR_bench * wealth_z     ;
        w_min = [1 500000 1000000 3000000 10000000];
        for ii=1:numel(w_min)
            % Bench
            %x_bench   = -log(wealth_bench(wealth_bench>w_min(ii)))     ;
            %y_bench   =  log(((numel(x_bench)):-1:1)/(numel(x_bench))) ;
            [y_bench,x_bench] = ecdf(wealth_bench(wealth_bench>w_min(ii))) ;
            x_bench = -log(x_bench(2:end)) ; y_bench = log(1-y_bench(1:end-1)) ;
            % mu=0.95
            %x_mu      = -log(wealth_mu(wealth_mu>w_min(ii)))           ;
            %y_mu      =  log(((numel(x_mu)):-1:1)/(numel(x_mu)))       ;
            [y_mu,x_mu] = ecdf(wealth_mu(wealth_mu>w_min(ii))) ;
            x_mu = -log(x_mu(2:end)) ; y_mu = log(1-y_mu(1:end-1)) ;
            % z grid
            %x_z       = -log(wealth_z(wealth_z>w_min(ii)))             ;
            %y_z       =  log(((numel(x_z)):-1:1)/(numel(x_z)))         ;
            [y_z,x_z] = ecdf(wealth_z(wealth_z>w_min(ii))) ;
            x_z = -log(x_z(2:end)) ; y_z = log(1-y_z(1:end-1)) ;
            % Vermuelen
            x_V_s     = -log(panel_V((panel_V>w_min(ii))&(type_V==1))) ;
            y_V_s     = log(cumsum(weights_V((panel_V>w_min(ii))&(type_V==1)),'reverse')/sum(weights_V(panel_V>w_min(ii)))) ;
            x_V_f     = -log(panel_V((panel_V>w_min(ii))&(type_V==0))) ;
            y_V_f     = log(cumsum(weights_V((panel_V>w_min(ii))&(type_V==0)),'reverse')/sum(weights_V(panel_V>w_min(ii)))) ;
            x_V       = [x_V_s x_V_f];      y_V = [y_V_s y_V_f];
            % Maximum Likelihood
            a_ml(ii)  = sum(weights_V(panel_V>w_min(ii))) / sum(weights_V(panel_V>w_min(ii)).*log(panel_V(panel_V>w_min(ii))/w_min(ii))) ;
            % Regression
            mdl_b     = fitlm(x_bench',y_bench')         ;
            a_b(ii)   = mdl_b.Coefficients.Estimate(2)   ; 
            mdl_mu    = fitlm(x_mu',y_mu')               ;
            a_mu(ii)  = mdl_mu.Coefficients.Estimate(2)  ;
            mdl_z     = fitlm(x_z',y_z')                 ;
            a_z(ii)   = mdl_z.Coefficients.Estimate(2)   ;
            mdl_V_s   = fitlm(x_V_s',y_V_s')             ;      
            C_V_s     = mdl_V_s.Coefficients.Estimate(1) ;
            a_V_s(ii) = mdl_V_s.Coefficients.Estimate(2) ;
            mdl_V     = fitlm(x_V',y_V')                 ;
            C_V       = mdl_V.Coefficients.Estimate(1)   ;
            a_V(ii)   = mdl_V.Coefficients.Estimate(2)   ;
            xx        = linspace(min(-x_V),max(-x_V),3)  ;
            % Figure
            fig_title = ['Pareto Tail Above $',num2str(w_min(ii)),' Vermeulen'] ;
            figure; hold on; 
            plot(-x_bench,y_bench,'linewidth',2); plot(-x_mu,y_mu,'linewidth',2); plot(-x_z,y_z,'linewidth',2); 
            scatter(-x_V_s,y_V_s); scatter(-x_V_f,y_V_f,'*'); 
            plot(xx,C_V_s-a_V_s(ii)*xx,'--k'); plot(xx,C_V-a_V(ii)*xx,':k'); plot(xx,a_ml(ii)*log(w_min(ii))-a_ml(ii)*xx,'-.k'); hold off;
            title(fig_title); xlabel('ln(Wealth/W(min))'); ylabel('ln(1-CDF(W))'); xlim([log(w_min(ii)),max(-x_V_f)]);
            ww = linspace(log(w_min(ii)),log(panel_V(end)),5) ;
            set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww));
            legend('Bench','mu=0.95','z-grid','Vermeulen Survey','Vermeulen Forbes','Reg Survey','Reg All','Pseudo Max Lik','location','SouthWest')
            file_name_pdf = ['Pareto_Tail_$',num2str(w_min(ii)),'_V.pdf'] ;
            print('-dpdf',file_name_pdf) ;
            % Table
            AA(:,ii) = [w_min(ii) ; a_b(ii) ; a_mu(ii) ; a_z(ii) ; a_V_s(ii) ; a_V(ii) ; a_ml(ii) ] ; 
        end
        
        bench_name = {'Bench',' ',' ',' ',' ',' '} ;
        mu_name    = {'mu=0.95',' ',' ',' ',' ',' '} ;
        z_name     = {'z_grid',' ',' ',' ',' ',' '} ;
        blank_row  = cell(1,6) ;
        
        row_name = {'w_min';'alpha_bench';'alpha_mu';'alpha_zgrid';'alpha_survey';'alpha_survey_forbes';'alpha_ml'} ;
        Mat = [blank_row ; row_name num2cell(AA)]
        status = xlwrite(Tables_file,Mat,'Test_Vermeulen') ;
        
    
    % Wealth Concentration - Vermeulen's Data
        CCDF_V         = cumsum(weights_V,'reverse')/sum(weights_V)       ;
        total_wealth_V = sum(panel_V.*weights_V)                          ; 
        CCDF_b         = (numel(wealth_bench):-1:1)/(numel(wealth_bench)) ;
        total_wealth_b = sum(wealth_bench)                                ;
        i=1;
        for prct   = [50 40 25 20 10 5 1 0.1]/100
            cum_wealth_V(i) = sum(panel_V(CCDF_V<=prct).*weights_V(CCDF_V<=prct))/total_wealth_V ;
            cum_wealth_b(i) = sum(wealth_bench(CCDF_b<=prct))/total_wealth_b   ;
            i=i+1;
        end 
        
        AA = [50 40 25 20 10 5 1 0.1; cum_wealth_b ; cum_wealth_V] ;
        
        rowname = {'top x%';'Bench';'Vermeulen'};
        Mat = [rowname num2cell(AA)]
        status = xlwrite(Tables_file,Mat,'Wealth_Concentration') ;
        
        % Fraction of observations above 1 million
        x_v = sum(weights_V(panel_V>1000000))/sum(weights_V) ;
        x_b = sum(wealth_bench>1000000)/numel(wealth_bench)  ;
        
        % Average wealth
        Ew_v = total_wealth_V/sum(weights_V)      ;
        Ew_b = total_wealth_b/numel(wealth_bench) ;
    
    
    

