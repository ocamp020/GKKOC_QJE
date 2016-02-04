% Simulation Graphs and Tables
% Sergio Ocampo Diaz

%% Clear and load required files
   % Clear and close
        close all;
        clear

   % Load files for excel in mac
        javaaddpath('poi_library/poi-3.8-20120326.jar');
        javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
        javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
        javaaddpath('poi_library/xmlbeans-2.3.0.jar');
        javaaddpath('poi_library/dom4j-1.6.1.jar');
        javaaddpath('poi_library/stax-api-1.0.1.jar');

        
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
    sigma_z_eps      = Params(4) ;
    sigma_lambda_eps = Params(5) ;
    gamma  = Params(6) ;
    
%% Result_Folders
cd '/Users/s-ocampo/Dropbox/ra_guvenen/wealth_tax/cggk_codes/Sergio/Graphs/Presentation/'

    if ((Progressive_Tax_Switch==0)&&(NSU_Switch==1))
        Result_Folder = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
    elseif ((Progressive_Tax_Switch==1)&&(NSU_Switch==1))
        Result_Folder = strcat('../../NSU_F_PT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../NSU_F_PT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../NSU_F_PT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
    elseif ((Progressive_Tax_Switch==0)&&(NSU_Switch==0))
        Result_Folder = strcat('../../SU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../SU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../SU_F_LT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
    elseif ((Progressive_Tax_Switch==1)&&(NSU_Switch==0))
        Result_Folder = strcat('../../SU_F_PT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
        Simul_Folder  = strcat('../../SU_F_PT_Results/Theta_',num2str(theta,'%.2f'),'/Factor_',num2str(Threshold_Factor,'%.2f'),'/Simul/') ;
        Bench_Folder  = strcat('../../SU_F_PT_Results/Theta_',num2str(theta,'%.2f'),'/Bench_Files/') ;
    end
        
%% Excel file
    Tables_file = 'Tables_Burhan.xls' ;
    
%% Composition of Z by welth percentile.
    % Each table shows the composition of the top x% into Z categories

    eval(['load ',Simul_Folder,'panela_bench']) ;    eval(['load ',Simul_Folder,'panelz_bench']);
    eval(['load ',Simul_Folder,'panela_exp'])   ;    eval(['load ',Simul_Folder,'panelz_exp'])  ;
    
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
        eval(['load ',Simul_Folder,'panelage_exp']) ; eval(['load ',Simul_Folder,'panela_exp']) ; 
        eval(['load ',Simul_Folder,'panel_return_exp']) ; eval(['load ',Simul_Folder,'panel_at_return_exp']) ;
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
    print('-dpdf','Wealth_Pareto.pdf') ;
    
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
        savefig(file_name); %print('-dpdf',file_name) ;
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
    print('-dpdf','Wealth_Prctl.pdf') ;
    

    
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
    axes1 = axes(...
      'FontName','Times New Roman',...
      'FontSize',18);
  
    G_K_Frac_k = Stats_by_tau_k(:,4) ;
    G_K_Frac_w = Stats_by_tau_w(:,4) ;
    

% Aggregate variabls as a function of capital/wealth tax revenue
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

    hgsave('1fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
    print -dpdf 1fig_KBAR_QBAR_by_CAP_TAX_REV.pdf
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


        hgsave('1.1fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
        print -dpdf 1.1fig_KBAR_QBAR_by_CAP_TAX_REV.pdf
        print -dps 1.1fig_KBAR_QBAR_by_CAP_TAX_REV.eps

        plot(G_K_Frac_k, 100*(Stats_by_tau_k(:,6)/Stats_by_tau_k(1,6)-1))
        plot(G_K_Frac_w, 100*(Stats_by_tau_w(:,6)/Stats_by_tau_w(1,6)-1))
        h2 = legend('$\bar k, \tau_k$' , '$\bar k, \tau_a$','$\bar Q, \tau_k$' , '$\bar Q, \tau_a$','Location','SouthWest');
        set(h2,'Interpreter','latex','FontSize',20)

        hgsave('1.2fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
        print -dpdf 1.2fig_KBAR_QBAR_by_CAP_TAX_REV.pdf
        print -dps 1.2fig_KBAR_QBAR_by_CAP_TAX_REV.eps


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


    hgsave('1fig_Wage_by_CAP_TAX_REV.fig')
    print -dpdf 1fig_Wage_by_CAP_TAX_REV.pdf
    print -dps  1fig_Wage_by_CAP_TAX_REV.eps
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
        print -dpdf 1fig_AT_Wage_by_CAP_TAX_REV.pdf
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
    print -dpdf 1fig_Average_Utility_by_CAP_TAX_REV.pdf
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
    print -dpdf 2fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
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
    print -dpdf 3fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
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

    annotation('textarrow',[0.76 0.73],[0.52 0.49], 'string','Cap. Income Tax Economy','linewidth',2,'FontSize',12)
    annotation('textarrow',[0.60 0.65],[0.5 0.55], 'string','Benchmark, \tau_k = 25%','linewidth',2,'FontSize',12)
    annotation('textarrow',[0.60 0.65],[0.17 0.12], 'string','0.28','linewidth',2,'FontSize',12)


    hgsave('1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -dpdf 1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.pdf
    print -dps 1.1.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps

    annotation('textarrow',[0.2 0.15],[0.71 0.76], 'string','Opt. \tau_k = 1.62%','linewidth',2,'FontSize',12)
    hgsave('1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -dpdf 1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.pdf
    print -dps 1.2.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.eps

    plot(G_K_Frac_w, (CE2_NB_w-CE2_NB_k(tauindx)),'b', 'linewidth',2)
    annotation('textarrow',[0.47 0.47],[0.81 0.86], 'string','Opt. \tau_a = 1.54%','linewidth',2,'FontSize',12)


    hgsave('1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.fig')
    print -dpdf 1.3.fig_Opt_Tax_Welfare_by_CAP_TAX_REV.pdf
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
    print -dpdf 1fig_diff_savings_rate_age1.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age1_z1.pdf
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

print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age1_z4.pdf
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

print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age1_z7.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age16_z1.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age16_z4.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age16_z7.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age31_z1.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age31_z4.pdf
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
print -dpdf 1fig_Welfare_Gain_by_Wealth_for_Age31_z7.pdf
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
print -dpdf 1fig_diff_savings_rate_age16.pdf
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
print -dpdf 1fig_diff_savings_rate_age31.pdf
print -dps  1fig_diff_savings_rate_age31.eps
print -dpng 1fig_diff_savings_rate_age31.png



