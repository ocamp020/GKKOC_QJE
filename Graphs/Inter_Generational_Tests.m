% Intergenerational Mobility
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
        
        cd '/Users/s-ocampo/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Inter_Generation/'

        
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
    Tables_file = 'Tables_Inter_Generation.xls' ; 
    
%% Benchmark
    % Result_Folders
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        
    % Load Simulation
        eval(['load ',Simul_Folder,'panela_old_1']) ; 
        eval(['load ',Simul_Folder,'panela_old_2']) ; 
        eval(['load ',Simul_Folder,'panela_old_3']) ; eval(['load ',Simul_Folder,'panelage_old_3']) ;
        eval(['load ',Simul_Folder,'panela_new_1']) ; 
        eval(['load ',Simul_Folder,'panela_new_2']) ; 
        eval(['load ',Simul_Folder,'panela_new_3']) ; eval(['load ',Simul_Folder,'panelage_new_3']) ;
        eval(['load ',Bench_Folder,'EBAR'])         ;
        
    % Averages
        Wage_earnings = 47000;
        panela_old_1 = Wage_earnings/EBAR * panela_old_1 ;
        panela_old_2 = Wage_earnings/EBAR * panela_old_2 ;
        panela_old_3 = Wage_earnings/EBAR * panela_old_3 ;
        panela_new_1 = Wage_earnings/EBAR * panela_new_1 ;
        panela_new_2 = Wage_earnings/EBAR * panela_new_2 ;
        panela_new_3 = Wage_earnings/EBAR * panela_new_3 ;
        
        
        panela_old = log((panela_old_1+panela_old_2+panela_old_3)/3) ; panelage_old = panelage_old_3 ;
        panela_new = log((panela_new_1+panela_new_2+panela_new_3)/3) ; panelage_new = panelage_new_3 ;
        
    w_cut = log([1 1000000 10000000]);  
    for i=1:numel(w_cut)
	% Regression Log-Levels
        x           = panela_old(panela_old>=w_cut(i)) ;
        y           = panela_new(panela_old>=w_cut(i)) ;
        mdl         = fitlm(x',y')                 ;
        cons(i)     = mdl.Coefficients.Estimate(1) ;
        slope(i)    = mdl.Coefficients.Estimate(2) ;
        SE_slope(i) = mdl.Coefficients.SE(2)       ;
        R2(i)       = mdl.Rsquared.Ordinary        ;
        n(i)        = numel(x)                     ;
        xx          = linspace(min(x),max(x),3)    ;
        
    % Graph
        figure; hold on; 
        scatter(x,y); plot(xx,cons(i)+slope(i)*xx,'linewidth',2);
        hold off; 
        fig_title = ['Inter-Generational Wealth: W>',num2str(exp(w_cut(i))),' bench'] ;
        title(fig_title); xlabel('ln(Wealth Father)'); ylabel('ln(Wealth Son)'); xlim([min(xx),max(xx)]);
        ww = linspace(min(panela_old),max(panela_old),5) ;
        set(gca,'XTick',ww); set(gca,'XTickLabel',exp(ww))
        file_name = ['Inter_Generational_Wealth_',num2str(exp(w_cut(i))),'_bench.fig'] ;
        file_name_pdf = ['Inter_Generational_Wealth_',num2str(exp(w_cut(i))),'_bench.pdf'] ;
        savefig(file_name); print('-dpdf',file_name_pdf) ;
    
        
    % Regression Rank
        x           = panela_old(panela_old>=w_cut(i)) ; x_age = panelage_old(panela_old>=w_cut(i)) ; rank_x = NaN(numel(x)) ;
        y           = panela_new(panela_old>=w_cut(i)) ; y_age = panelage_new(panela_old>=w_cut(i)) ; rank_y = NaN(numel(x)) ;
        
        for age = unique(x_age)
           aux = x(x_age==age) ;
           [~, ~, rank_aux] = unique(aux) ; rank_aux = 100*(rank_aux-0.5)/numel(unique(aux)) ;
           rank_x(x_age==age) = rank_aux  ;
        end
        for age = unique(y_age)
           aux = y(y_age==age) ;
           [~, ~, rank_aux] = unique(aux) ; rank_aux = 100*(rank_aux-0.5)/numel(unique(aux)) ;
           rank_y(y_age==age) = rank_aux  ;
        end
        
        % Bins
        for j=1:100
            y_aux(j) = mean(rank_y( (rank_x>(j-1))&(rank_x<j) )) ;
            x_aux(j) = j                                         ;
        end
        % Regression
        mdl            = fitlm(x_aux',y_aux')       ;
        cons_R(i)      = mdl.Coefficients.Estimate(1) ;
        slope_R(i)     = mdl.Coefficients.Estimate(2) ;
        SE_slope_R(i)  = mdl.Coefficients.SE(2)       ;
        R2_R(i)        = mdl.Rsquared.Ordinary        ;
        xx             = linspace(0,100,3)             ;
	
    % Graph
        figure; hold on; 
        scatter(x_aux,y_aux); plot(xx,cons_R(i)+slope_R(i)*xx,'linewidth',2); plot(x_aux,x_aux,'--k')
        hold off; 
        fig_title = ['Inter-Generational Rank: W>',num2str(exp(w_cut(i))),' bench'] ;
        title(fig_title); xlabel('Rank Wealth Father'); ylabel('Average Rank Wealth Son'); xlim([min(xx),max(xx)]);
        file_name = ['Inter_Generational_Rank_',num2str(exp(w_cut(i))),'_bench.fig'] ;
        file_name_pdf = ['Inter_Generational_Rank_',num2str(exp(w_cut(i))),'_bench.pdf'] ;
        savefig(file_name); print('-dpdf',file_name_pdf) ;
        
        
	end 
        
	% Table
        col_name = {'w_cut','Cons','Slope','SE_slope','R2','obs'};
        AA       = [w_cut',cons',slope',SE_slope',R2',n'];
        Mat      = [col_name; num2cell(AA)]
        status   = xlwrite(Tables_file,Mat,'Log-Level-bench') ;
        col_name = {'w_cut','Cons','Slope','SE_slope','R2','obs'};
        AA       = [w_cut',cons_R',slope_R',SE_slope_R',R2_R',n'];
        Mat      = [col_name; num2cell(AA)]
        status   = xlwrite(Tables_file,Mat,'Rank-bench') ;
        
%% Benchmark

    % Result_Folders
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
        
    % Load Simulation
        eval(['load ',Simul_Folder,'panelz_bench'])       ;
        eval(['load ',Simul_Folder,'panela_bench'])       ;
        eval(['load ',Simul_Folder,'panelage_bench'])     ;
        eval(['load ',Simul_Folder,'panel_aprime_bench']) ;
        eval(['load ',Bench_Folder,'EBAR'])               ;
        
    % Averages
        Wage_earnings = 47000;
        assets   = Wage_earnings/EBAR * panela_bench       ;
        assets_p = Wage_earnings/EBAR * panel_aprime_bench ;
        
        for i=1:Max_Age
            assets_age(i) = mean(assets(panelage_bench==i))   ;
            saving_age(i) = mean(assets_p(panelage_bench==i)) ;
            for j=1:n_z
            assets_az(i,j) = mean(assets(   (panelage_bench==i)&(panelz_bench==j) )) ;
            saving_az(i,j) = mean(assets_p( (panelage_bench==i)&(panelz_bench==j) )) ;
            end
        end
        S_rate_age = 100*(saving_age./assets_age-1) ;
        S_rate_az  = 100*(saving_az./assets_az-1)   ;
        
	% Graph
        figure; 
        subplot(2,1,1); plot(20:19+Max_Age,assets_age); ylabel('Wealth'); title('Wealth Age Profile');
        subplot(2,1,2); plot(20:19+Max_Age,assets_az) ; xlabel('Age')   ; ylabel('Wealth'); legend('z1','z2','z3','z4','z5','z6','z7','location','northwest') 
        print('-dpdf','Wealth_Age_Profile.pdf') ;
        figure; 
        subplot(2,2,1); plot(20:19+Ret_Age,S_rate_age(1:Ret_Age)); ylabel('Saving Rate'); title('Wealth Saving Rate Profile Before Retirement'); xlim([20 19+Ret_Age])
        subplot(2,2,3); plot(20:19+Ret_Age,S_rate_az(1:Ret_Age,:)) ; xlabel('Age')   ; ylabel('Saving Rate'); legend('z1','z2','z3','z4','z5','z6','z7','location','best'); xlim([20 19+Ret_Age]) 
        subplot(2,2,2); plot(20+Ret_Age:19+Max_Age,S_rate_age(Ret_Age+1:end)); ylabel('Saving Rate'); title('Wealth Saving Rate Profile After Retirement'); xlim([20+Ret_Age Max_Age])
        subplot(2,2,4); plot(20+Ret_Age:19+Max_Age,S_rate_az(Ret_Age+1:end,:)) ; xlabel('Age')   ; ylabel('Saving Rate'); legend('z1','z2','z3','z4','z5','z6','z7','location','best'); xlim([20+Ret_Age Max_Age])
        print('-dpdf','Savings_Age_Profile.pdf') ;
        
        