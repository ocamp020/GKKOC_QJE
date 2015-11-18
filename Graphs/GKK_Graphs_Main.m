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
    
    xls_file = 'AZ_Tables.xls' ;

%% Grids 
   
% A grid
    eval(['load ',Result_Folder,'agrid']);
    
    A_mat = repmat(agrid,[Max_Age,1,n_z,n_l,n_e]);
    
% Z grid
    eval(['load ',Result_Folder,'zgrid']);
    
    Z_mat = repmat(reshape(zgrid,[1,1,n_z,1,1]),[Max_Age,n_a,1,n_l,n_e]);
    
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
    
%% Consumption Equivalent matrices

    if ((tauPL==0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/CE_Files/') ;
        Bench_Folder  = '../NSU_LT_Results/Bench_Files/' ;
    elseif ((tauPL~=0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/CE_Files/') ;
        Bench_Folder  = '../NSU_PT_Results/Bench_Files/' ;
    elseif ((tauPL==0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/CE_Files/') ;
        Bench_Folder  = '../SU_LT_Results/Bench_Files/' ;
    elseif ((tauPL~=0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/CE_Files/') ;
        Bench_Folder  = '../SU_PT_Results/Bench_Files/' ;
    end

% Consumption equivalent
    CE = 100*((V_exp./V_bench).^(1/((1-sigma)*gamma)) - 1) ;
    
% Auxiliary consumption equivalents
    eval(['load ',Result_Folder,'CE_c'])
    eval(['load ',Result_Folder,'CE_cl'])
    eval(['load ',Result_Folder,'CE_cd'])
    eval(['load ',Result_Folder,'CE_h'])
    eval(['load ',Result_Folder,'CE_hl'])
    eval(['load ',Result_Folder,'CE_hd'])
    
    CE_c  = reshape(CE_c,[Max_Age,n_a,n_z,n_l,n_e]) ;
    CE_cl = reshape(CE_cl,[Max_Age,n_a,n_z,n_l,n_e]) ;
    CE_cd = reshape(CE_cd,[Max_Age,n_a,n_z,n_l,n_e]) ;
    CE_h  = reshape(CE_h,[Max_Age,n_a,n_z,n_l,n_e]) ;
    CE_hl = reshape(CE_hl,[Max_Age,n_a,n_z,n_l,n_e]) ;
    CE_hd = reshape(CE_hd,[Max_Age,n_a,n_z,n_l,n_e]) ;
    
    
    
%% Age - Z Tables

% Weights 
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
    
    Mat = [z_title;age_title,num2cell(Weights_AZ)] ;
   
    status = xlwrite(xls_file,Mat,'Weights') ;
    
% Wealth by AZ group
    A_AZ_bench  = zeros(n_age,n_z) ;
    A_AZ_exp    = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            A_AZ_bench(j,i) = A_AZ_bench(j,i) + sum(sum(sum(sum(A_mat(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            A_AZ_exp(j,i)   = A_AZ_exp(j,i)   + sum(sum(sum(sum(A_mat(l,:,i,:,:).*DBN_exp(l,:,i,:,:))))) ;
            end 
            if l==age(j); j=j+1; end 
        end
    end
    
    A_AZ_level = A_AZ_exp - A_AZ_bench ;
    A_AZ_prop  = 100*(A_AZ_exp - A_AZ_bench)./A_AZ_bench ;
    
    Weights_A_AZ_bench = 100*A_AZ_bench/sum(sum(A_AZ_bench)) ;
    Weights_A_AZ_exp   = 100*A_AZ_exp/sum(sum(A_AZ_exp))     ;
    
    Weights_A_AZ_diff  = Weights_A_AZ_exp-Weights_A_AZ_bench ;
    
    Mat_1 = [{'A_level',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_level)] ;
    Mat_2 = [{'A_prop',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_prop)] ;
    Mat_3 = [{'A_Weight_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Weights_A_AZ_bench)] ;
    Mat_4 = [{'A_Weight_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Weights_A_AZ_exp)] ;
    Mat_5 = [{'A_Weight_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Weights_A_AZ_diff)] ;
    Mat_6 = cell(n_age+2,n_z+1) ;
   
    Mat = [Mat_1,Mat_2;cell(2,2*n_z+2);Mat_3,Mat_4;cell(2,2*n_z+2);Mat_5,Mat_6] ;
    status = xlwrite(xls_file,Mat,'Assets') ;
    
    
    Average_A_AZ_bench = A_AZ_bench./Weights_AZ ;
    Mat_1 = [{'A_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_level)] ;
    Mat_2 = [{'A_prop',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(A_AZ_prop)] ;
    
    

% Capital by age group
    Q_AZ_bench  = zeros(n_age,n_z) ;
    Q_AZ_exp    = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            Q_AZ_bench(j,i) = Q_AZ_bench(j,i) + sum(sum(sum(sum(zgrid(i)*A_mat(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            Q_AZ_exp(j,i)   = Q_AZ_exp(j,i)   + sum(sum(sum(sum(zgrid(i)*A_mat(l,:,i,:,:).*DBN_exp(l,:,i,:,:))))) ;
            end 
            if l==age(j); j=j+1; end 
        end
    end
    
    Q_AZ_level = Q_AZ_exp - Q_AZ_bench ;
    Q_AZ_prop  = 100*(Q_AZ_exp - Q_AZ_bench)./Q_AZ_bench ;
    
    Weights_Q_AZ_bench = 100*Q_AZ_bench/sum(sum(Q_AZ_bench)) ;
    Weights_Q_AZ_exp   = 100*Q_AZ_exp/sum(sum(Q_AZ_exp))     ;
    
    Weights_Q_AZ_diff    = Weights_Q_AZ_exp-Weights_Q_AZ_bench ;
    
    
    Mat_1 = [{'Q_level',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Q_AZ_level)] ;
    Mat_2 = [{'Q_prop',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Q_AZ_prop)] ;
    Mat_3 = [{'Q_Weight_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Weights_Q_AZ_bench)] ;
    Mat_4 = [{'Q_Weight_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Weights_Q_AZ_exp)] ;
    Mat_5 = [{'Q_Weight_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(Weights_Q_AZ_diff)] ;
    Mat_6 = cell(n_age+2,n_z+1) ;
   
    Mat = [Mat_1,Mat_2;cell(2,2*n_z+2);Mat_3,Mat_4;cell(2,2*n_z+2);Mat_5,Mat_6] ;
    status = xlwrite(xls_file,Mat,'Capital') ;
    
    
% Savings 
    Ap_AZ_bench  = zeros(n_age,n_z) ;
    Ap_AZ_exp    = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            Ap_AZ_bench(j,i) = Ap_AZ_bench(j,i) + sum(sum(sum(sum(Ap_bench(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            Ap_AZ_exp(j,i)   = Ap_AZ_exp(j,i)   + sum(sum(sum(sum(Ap_exp(l,:,i,:,:).*DBN_exp(l,:,i,:,:)))))     ;
            end 
            if l==age(j); j=j+1; end 
        end
    end
    
    S_AZ_bench = 100*(Ap_AZ_bench-A_AZ_bench)./A_AZ_bench ;
    S_AZ_exp   = 100*(Ap_AZ_exp  -A_AZ_exp  )./A_AZ_exp   ;
    
    
    S_AZ_diff  = S_AZ_exp - S_AZ_bench  ;
    
    S_AZ_level = 100*((Ap_AZ_exp  -A_AZ_exp  ) - (Ap_AZ_bench-A_AZ_bench))./Weights_AZ  ;
    
    
    Mat_1 = [{'S_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_bench)] ;
    Mat_2 = [{'S_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_exp)] ;
    Mat_3 = [{'S_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_diff)] ;
    Mat_4 = [{'S_diff_level',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(S_AZ_level)] ;
   
    Mat = [Mat_1,Mat_2;cell(2,2*n_z+2);Mat_3,Mat_4] ;
    status = xlwrite(xls_file,Mat,'Savings') ;
    
% Hours
    H_AZ_bench  = zeros(n_age,n_z) ;
    H_AZ_exp    = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            H_AZ_bench(j,i) = H_AZ_bench(j,i) + sum(sum(sum(sum(H_bench(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            H_AZ_exp(j,i)   = H_AZ_exp(j,i)   + sum(sum(sum(sum(H_exp(l,:,i,:,:).*DBN_exp(l,:,i,:,:)))))     ;
            end 
            if l==age(j); j=j+1; end 
        end
    end

    H_AZ_bench  = 100*H_AZ_bench./Weights_AZ ;
    H_AZ_exp    = 100*H_AZ_exp./Weights_AZ   ;
    
    H_AZ_diff   = H_AZ_exp - H_AZ_bench      ;
    H_AZ_prop   = 100*(H_AZ_exp./H_AZ_bench-1) ;
    
    
    Mat_1 = [{'H_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_bench)] ;
    Mat_2 = [{'H_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_exp)] ;
    Mat_3 = [{'H_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_diff)] ;
    Mat_4 = [{'H_diff_prc',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(H_AZ_prop)] ;
   
    Mat = [Mat_1,Mat_2;cell(2,2*n_z+2);Mat_3,Mat_4] ;
    status = xlwrite(xls_file,Mat,'Hours') ;
    
% Consumption  
    C_AZ_bench  = zeros(n_age,n_z) ;
    C_AZ_exp    = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            C_AZ_bench(j,i) = C_AZ_bench(j,i) + sum(sum(sum(sum(C_bench(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            C_AZ_exp(j,i)   = C_AZ_exp(j,i)   + sum(sum(sum(sum(C_exp(l,:,i,:,:).*DBN_exp(l,:,i,:,:)))))     ;
            end 
            if l==age(j); j=j+1; end 
        end
    end

    C_AZ_bench  = 100*C_AZ_bench./Weights_AZ ;
    C_AZ_exp    = 100*C_AZ_exp./Weights_AZ   ;
    
    C_AZ_diff   = 100*(C_AZ_exp./C_AZ_bench-1) ;
    
    Mat_1 = [{'C_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_bench)] ;
    Mat_2 = [{'C_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_exp)] ;
    Mat_3 = [{'C_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(C_AZ_diff)] ;
    Mat_4 = cell(n_age+2,n_z+1) ;
   
    Mat = [Mat_1,Mat_2;cell(2,2*n_z+2);Mat_3,Mat_4] ;
    status = xlwrite(xls_file,Mat,'Consumption') ;


% Consumption Equivalent Welfare 
    CE_AZ  = zeros(n_age,n_z) ;
    CE_c_AZ  = zeros(n_age,n_z) ;
    CE_cl_AZ  = zeros(n_age,n_z) ;
    CE_cd_AZ  = zeros(n_age,n_z) ;
    CE_h_AZ  = zeros(n_age,n_z) ;
    CE_hl_AZ  = zeros(n_age,n_z) ;
    CE_hd_AZ  = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            CE_AZ(j,i)    = CE_AZ(j,i)    + sum(sum(sum(sum(CE(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            CE_c_AZ(j,i)  = CE_c_AZ(j,i)  + sum(sum(sum(sum(CE_c(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            CE_cl_AZ(j,i) = CE_cl_AZ(j,i) + sum(sum(sum(sum(CE_cl(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            CE_cd_AZ(j,i) = CE_cd_AZ(j,i) + sum(sum(sum(sum(CE_cd(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            CE_h_AZ(j,i)  = CE_h_AZ(j,i)  + sum(sum(sum(sum(CE_h(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            CE_hl_AZ(j,i) = CE_hl_AZ(j,i) + sum(sum(sum(sum(CE_hl(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            CE_hd_AZ(j,i) = CE_hd_AZ(j,i) + sum(sum(sum(sum(CE_hd(l,:,i,:,:).*DBN_bench(l,:,i,:,:))))) ;
            end 
            if l==age(j); j=j+1; end 
        end
    end

    CE_AZ     = 100*CE_AZ./Weights_AZ ;
    CE_c_AZ   = 100*CE_c_AZ./Weights_AZ ;
    CE_cl_AZ  = 100*CE_cl_AZ./Weights_AZ ;
    CE_cd_AZ  = 100*CE_cd_AZ./Weights_AZ ;
    CE_h_AZ   = 100*CE_h_AZ./Weights_AZ ;
    CE_hl_AZ  = 100*CE_hl_AZ./Weights_AZ ;
    CE_hd_AZ  = 100*CE_hd_AZ./Weights_AZ ;
    
    Mat_1 = [{'CE',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_AZ)] ;
    Mat_2 = [{'CE_cons',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_c_AZ)] ;
    Mat_3 = [{'CE_cons level',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_cl_AZ)] ;
    Mat_4 = [{'CE_cons dist',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_cd_AZ)] ;
    Mat_5 = [{'CE_hour',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_h_AZ)] ;
    Mat_6 = [{'CE_hour level',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_hl_AZ)] ;
    Mat_7 = [{'CE_hour dist',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(CE_hd_AZ)] ;
    Mat_8 = cell(n_age+2,n_z+1) ;
    Mat_9 = cell(n_age+2,1) ;
    Mat_10 = cell(1,3*(n_z+1)+3) ;
    
    
    Mat = [Mat_10 ; 
           Mat_9 Mat_1 Mat_9 Mat_2 Mat_9 Mat_5 ; 
           Mat_10 ; 
           Mat_9 Mat_8 Mat_9 Mat_3 Mat_9 Mat_6;
           Mat_10;
           Mat_9 Mat_8 Mat_9 Mat_4 Mat_9 Mat_7] ; 
    status = xlwrite(xls_file,Mat,'C-E') ;
    
    
%% Moments by age

% Distribution by age 
    DBN_age_bench = DBN_bench ;
    DBN_age_exp   = DBN_exp ;
    DBN_AZ_bench  = DBN_bench ;
    DBN_AZ_exp  = DBN_exp ;
    for i=1:Max_Age
        DBN_age_bench(i,:,:,:,:) = DBN_bench(i,:,:,:,:)/sum(sum(sum(sum(DBN_bench(i,:,:,:,:))))) ;
        DBN_age_exp(i,:,:,:,:)   = DBN_exp(i,:,:,:,:)/sum(sum(sum(sum(DBN_exp(i,:,:,:,:)))))     ;
        for j=1:n_z    
            DBN_AZ_bench(i,:,j,:,:) = DBN_bench(i,:,j,:,:)/sum(sum(sum(DBN_bench(i,:,j,:,:),5),4),2) ;
            DBN_AZ_exp(i,:,j,:,:)   = DBN_exp(i,:,j,:,:)/sum(sum(sum(DBN_exp(i,:,j,:,:),5),4),2)     ;
        end 
    end 

% Average consumption, labor and asset holdings by age
    mean_C_age_bench = NaN(Max_Age,1);
    mean_C_age_exp = NaN(Max_Age,1);
    mean_H_age_bench = NaN(Max_Age,1);
    mean_H_age_exp   = NaN(Max_Age,1);
    mean_A_age_bench = NaN(Max_Age,1);
    mean_A_age_exp   = NaN(Max_Age,1);
    mean_Ap_age_bench = NaN(Max_Age,1);
    mean_Ap_age_exp   = NaN(Max_Age,1);
    mean_S_age_bench = NaN(Max_Age,1);
    mean_S_age_exp   = NaN(Max_Age,1);
    
    mean_C_AZ_bench = NaN(Max_Age,n_z);
    mean_C_AZ_exp   = NaN(Max_Age,n_z);
    mean_H_AZ_bench = NaN(Max_Age,n_z);
    mean_H_AZ_exp   = NaN(Max_Age,n_z);
    mean_A_AZ_bench = NaN(Max_Age,n_z);
    mean_A_AZ_exp   = NaN(Max_Age,n_z);
    mean_Ap_AZ_bench = NaN(Max_Age,n_z);
    mean_Ap_AZ_exp   = NaN(Max_Age,n_z);
    mean_S_AZ_bench = NaN(Max_Age,n_z);
    mean_S_AZ_exp   = NaN(Max_Age,n_z);
    for i=1:Max_Age
        mean_C_age_bench(i) = sum(sum(sum(sum( C_bench(i,:,:,:,:).*DBN_age_bench(i,:,:,:,:) )))) ;
        mean_H_age_bench(i) = sum(sum(sum(sum( H_bench(i,:,:,:,:).*DBN_age_bench(i,:,:,:,:) )))) ;
        mean_A_age_bench(i) = sum(sum(sum(sum( A_mat(i,:,:,:,:).*DBN_age_bench(i,:,:,:,:) )))) ;
        mean_Ap_age_bench(i) = sum(sum(sum(sum(Ap_bench(i,:,:,:,:).*DBN_age_bench(i,:,:,:,:) )))) ;
        mean_S_age_bench(i) = 100*(mean_Ap_age_bench(i)/mean_A_age_bench(i)-1) ;
        
        mean_C_age_exp(i)   = sum(sum(sum(sum( C_exp(i,:,:,:,:).*DBN_age_exp(i,:,:,:,:) )))) ;
        mean_H_age_exp(i)  = sum(sum(sum(sum( H_exp(i,:,:,:,:).*DBN_age_exp(i,:,:,:,:) )))) ;
        mean_A_age_exp(i)  = sum(sum(sum(sum( A_mat(i,:,:,:,:).*DBN_age_exp(i,:,:,:,:) )))) ;
        mean_Ap_age_exp(i) = sum(sum(sum(sum(Ap_exp(i,:,:,:,:).*DBN_age_exp(i,:,:,:,:) )))) ;
        mean_S_age_exp(i)  = 100*(mean_Ap_age_exp(i)/mean_A_age_exp(i)-1) ;
        
        for j=1:n_z
            mean_C_AZ_bench(i,j)  = sum(sum(sum( C_bench(i,:,j,:,:).*DBN_AZ_bench(i,:,j,:,:) ))) ;
            mean_H_AZ_bench(i,j)  = sum(sum(sum( H_bench(i,:,j,:,:).*DBN_AZ_bench(i,:,j,:,:) ))) ;
            mean_A_AZ_bench(i,j)  = sum(sum(sum( A_mat(i,:,j,:,:).*DBN_AZ_bench(i,:,j,:,:) ))) ;
            mean_Ap_AZ_bench(i,j) = sum(sum(sum(Ap_bench(i,:,j,:,:).*DBN_AZ_bench(i,:,j,:,:) ))) ;
            mean_S_AZ_bench(i,j)  = 100*(mean_Ap_AZ_bench(i,j)/mean_A_AZ_bench(i,j)-1) ;

            mean_C_AZ_exp(i,j)   = sum(sum(sum( C_exp(i,:,j,:,:).*DBN_AZ_exp(i,:,j,:,:) ))) ;
            mean_H_AZ_exp(i,j)   = sum(sum(sum( H_exp(i,:,j,:,:).*DBN_AZ_exp(i,:,j,:,:) ))) ;
            mean_A_AZ_exp(i,j)   = sum(sum(sum( A_mat(i,:,j,:,:).*DBN_AZ_exp(i,:,j,:,:) ))) ;
            mean_Ap_AZ_exp(i,j)  = sum(sum(sum(Ap_exp(i,:,j,:,:).*DBN_AZ_exp(i,:,j,:,:) ))) ;
            mean_S_AZ_exp(i,j)   = 100*(mean_Ap_AZ_exp(i,j)/mean_A_AZ_exp(i,j)-1) ;
        end 
    end 
    
    figure; 
        subplot(1,4,1); plot(20:100,[mean_A_age_bench mean_A_age_exp]); xlim([19+1,19+Max_Age]); title('Mean Assets')
        subplot(1,4,2); plot(20:80,[mean_S_age_bench(1:end-20) mean_S_age_exp(1:end-20)]); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate')
        subplot(1,4,3); plot(20:100,[mean_C_age_bench mean_C_age_exp]); xlim([19+1,19+Max_Age]); title('Mean Consumption')
        subplot(1,4,4); plot(20:100,[mean_H_age_bench mean_H_age_exp]); xlim([19+1,19+Max_Age]); title('Mean Hours')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','Mean_Variables_Life_Cycle.pdf') ;
    
    z_vec = [2 4 6] ;
    figure; 
        subplot(2,4,1); plot(20:100,mean_A_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Assets - bench')
        subplot(2,4,2); plot(20:80,mean_S_AZ_bench(1:end-20,z_vec)); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate - bench')
        subplot(2,4,3); plot(20:100,mean_C_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Consumption - bench')
        subplot(2,4,4); plot(20:100,mean_H_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Hours - bench')
        legend('z2','z4','z6','location','southeast')
        subplot(2,4,5); plot(20:100,mean_A_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Assets - exp')
        subplot(2,4,6); plot(20:80,mean_S_AZ_exp(1:end-20,z_vec)); xlim([19+1,19+Max_Age-20]); title('Mean Saving Rate - exp')
        subplot(2,4,7); plot(20:100,mean_C_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Consumption - exp')
        subplot(2,4,8); plot(20:100,mean_H_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','Mean_Variables_Life_Cycle_by_Z.pdf') ;
    
% Std consumption, labor and asset holdings by age
    std_C_age_bench = NaN(Max_Age,1);
    std_C_age_exp   = NaN(Max_Age,1);
    std_H_age_bench = NaN(Max_Age,1);
    std_H_age_exp   = NaN(Max_Age,1);
    std_A_age_bench = NaN(Max_Age,1);
    std_A_age_exp   = NaN(Max_Age,1);
    
    std_C_AZ_bench = NaN(Max_Age,n_z);
    std_C_AZ_exp   = NaN(Max_Age,n_z);
    std_H_AZ_bench = NaN(Max_Age,n_z);
    std_H_AZ_exp   = NaN(Max_Age,n_z);
    std_A_AZ_bench = NaN(Max_Age,n_z);
    std_A_AZ_exp   = NaN(Max_Age,n_z);
    for i=1:Max_Age
        std_C_age_bench(i) = ( sum(sum(sum(sum( (C_bench(i,:,:,:,:)-repmat(mean_C_age_bench(i),[1,n_a,n_z,n_l,n_e])).^2 .*DBN_age_bench(i,:,:,:,:) )))) ).^0.5 ;
        std_H_age_bench(i) = ( sum(sum(sum(sum( (H_bench(i,:,:,:,:)-repmat(mean_H_age_bench(i),[1,n_a,n_z,n_l,n_e])).^2 .*DBN_age_bench(i,:,:,:,:) )))) ).^0.5 ;
        std_A_age_bench(i) = ( sum(sum(sum(sum( (A_mat(i,:,:,:,:)-repmat(mean_A_age_bench(i),[1,n_a,n_z,n_l,n_e])).^2 .*DBN_age_bench(i,:,:,:,:) )))) ).^0.5 ;
        
        std_C_age_exp(i) = ( sum(sum(sum(sum( (C_exp(i,:,:,:,:)-repmat(mean_C_age_exp(i),[1,n_a,n_z,n_l,n_e])).^2 .*DBN_age_exp(i,:,:,:,:) )))) ).^0.5 ;
        std_H_age_exp(i) = ( sum(sum(sum(sum( (H_exp(i,:,:,:,:)-repmat(mean_H_age_exp(i),[1,n_a,n_z,n_l,n_e])).^2 .*DBN_age_exp(i,:,:,:,:) )))) ).^0.5 ;
        std_A_age_exp(i) = ( sum(sum(sum(sum( (A_mat(i,:,:,:,:)-repmat(mean_A_age_exp(i),[1,n_a,n_z,n_l,n_e])).^2 .*DBN_age_exp(i,:,:,:,:) )))) ).^0.5 ;
        
        for j=1:n_z
            std_C_AZ_bench(i,j) = ( sum(sum(sum( (C_bench(i,:,j,:,:)-repmat(mean_C_AZ_bench(i,j),[1,n_a,1,n_l,n_e])).^2 .*DBN_AZ_bench(i,:,j,:,:) ))) ).^0.5 ;
            std_H_AZ_bench(i,j) = ( sum(sum(sum( (H_bench(i,:,j,:,:)-repmat(mean_H_AZ_bench(i,j),[1,n_a,1,n_l,n_e])).^2 .*DBN_AZ_bench(i,:,j,:,:) ))) ).^0.5 ;
            std_A_AZ_bench(i,j) = ( sum(sum(sum( (A_mat(i,:,j,:,:)-repmat(mean_A_AZ_bench(i,j),[1,n_a,1,n_l,n_e])).^2 .*DBN_AZ_bench(i,:,j,:,:) ))) ).^0.5 ;

            std_C_AZ_exp(i,j) = ( sum(sum(sum( (C_exp(i,:,j,:,:)-repmat(mean_C_AZ_exp(i,j),[1,n_a,1,n_l,n_e])).^2 .*DBN_AZ_exp(i,:,j,:,:) ))) ).^0.5 ;
            std_H_AZ_exp(i,j) = ( sum(sum(sum( (H_exp(i,:,j,:,:)-repmat(mean_H_AZ_exp(i,j),[1,n_a,1,n_l,n_e])).^2 .*DBN_AZ_exp(i,:,j,:,:) ))) ).^0.5 ;
            std_A_AZ_exp(i,j) = ( sum(sum(sum( (A_mat(i,:,j,:,:)-repmat(mean_A_AZ_exp(i,j),[1,n_a,1,n_l,n_e])).^2 .*DBN_AZ_exp(i,:,j,:,:) ))) ).^0.5 ;
        end 
    end 
    
    figure; 
        subplot(1,3,1); plot(20:100,[std_A_age_bench std_A_age_exp]); xlim([19+1,19+Max_Age]); title('Std Assets')
        subplot(1,3,2); plot(20:100,[std_C_age_bench std_C_age_exp]); xlim([19+1,19+Max_Age]); title('Std Consumption')
        subplot(1,3,3); plot(20:100,[std_H_age_bench std_H_age_exp]); xlim([19+1,19+Max_Age]); title('Std Hours')
        legend('Bench','Exp')
    print('-dpdf','Std_Variables_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,3,1); plot(20:100,std_A_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Assets - bench')
        subplot(2,3,2); plot(20:100,std_C_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Consumption - bench')
        subplot(2,3,3); plot(20:100,std_H_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Hours - bench')
        subplot(2,3,4); plot(20:100,std_A_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Assets - exp')
        subplot(2,3,5); plot(20:100,std_C_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Consumption - exp')
        subplot(2,3,6); plot(20:100,std_H_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Hours - exp')
        legend('z2','z4','z6','location','southeast')
    print('-dpdf','Std_Variables_Life_Cycle_by_Z.pdf') ;
    
%% Distribution Graphs
% Age brackets 
    age = [ 1,  6, 16, 36, 46, 56 ];
    %%%    20, 25, 35, 55, 65, 75
    n_age = numel(age) ;
    
    max_a = 30 ;
    
figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(agrid(1:max_a)',[sum(sum(sum(DBN_age_bench(age(i),1:max_a,:,:,:),5),4),3)' sum(sum(sum(DBN_age_exp(age(i),1:max_a,:,:,:),5),4),3)'])
       title(strcat('Age=',int2str(19+age(i)))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('bench','exp');
    print('-dpdf','Asset_DBN_by_age.pdf') ;
    
figure;
    for i=1:n_age
       subplot(2,n_age,i); 
       hold on;
       plot( agrid(1:max_a)',squeeze(sum(sum(DBN_age_bench(age(i),1:max_a,:,:,:),5),4)) )
       title(strcat('Age=',int2str(19+age(i)))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('z1','z2','z3','z4','z5','z6','z7');
    for i=1:n_age
        subplot(2,n_age,n_age+i); 
        hold on;
        plot(agrid(1:max_a)',squeeze(sum(sum(DBN_age_exp(age(i),1:max_a,:,:,:),5),4)) )
        title(strcat('Age=',int2str(19+age(i))));  xlim([0,floor(agrid(max_a))]);
        hold off;
    end
    legend('z1','z2','z3','z4','z5','z6','z7');
    print('-dpdf','Asset_DBN_by_z_age.pdf') ;
    
%% Policy Function Graphs
for z=[1,4,7]
figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(agrid(1:max_a)',[Ap_bench(age(i),1:max_a,z,3,3)' Ap_exp(age(i),1:max_a,z,3,3)'])
       plot(agrid(1:max_a),agrid(1:max_a))
       title(strcat('Age=',int2str(19+age(i)),' z',int2str(z))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('bench','exp','location','southeast'); 
    print('-dpdf',strcat('PolFun_Ap_by_age_z',int2str(z),'.pdf') ) ;
    
figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(agrid(1:max_a)',[C_bench(age(i),1:max_a,z,3,3)' C_exp(age(i),1:max_a,z,3,3)'])
       title(strcat('Age=',int2str(19+age(i)),' z',int2str(z))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('bench','exp','location','southeast'); 
    print('-dpdf',strcat('PolFun_C_by_age_z',int2str(z),'.pdf') ) ;
    
figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(agrid(1:max_a)',[H_bench(age(i),1:max_a,z,3,3)' H_exp(age(i),1:max_a,z,3,3)'])
       title(strcat('Age=',int2str(19+age(i)),' z',int2str(z))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('bench','exp','location','southeast'); 
    print('-dpdf',strcat('PolFun_H_by_age_z',int2str(z),'.pdf') ) ;
    
figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(agrid(1:max_a)',[V_bench(age(i),1:max_a,z,3,3)' V_exp(age(i),1:max_a,z,3,3)'])
       title(strcat('Age=',int2str(19+age(i)),' z',int2str(z))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('bench','exp','location','southeast'); 
    print('-dpdf',strcat('Value_by_age_z',int2str(z),'.pdf') ) ;
    
end 

max_a = 20 ;
for z=[2,3,4,5,6]
figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(agrid(1:max_a)',[Ap_bench(age(i),1:max_a,z,3,3)' Ap_exp(age(i),1:max_a,z,3,3)'])
       plot(agrid(1:max_a),agrid(1:max_a))
       title(strcat('Age=',int2str(19+age(i)),' z',int2str(z))); xlim([0,floor(agrid(max_a))]);
       hold off;
    end
    legend('bench','exp','45?','location','southeast'); 
    print('-dpdf',strcat('PolFun_Ap_by_age_z',int2str(z),'.pdf') ) ;
    
end 


%% Asset distribution
% Asset distribution by grid point
    DBN_agrid_bench = sum(sum(sum(sum(DBN_bench(:,:,:,:,:),5),4),3),1) ;
    result_agrid_1 = [114 agrid(114) sum(DBN_agrid_bench(114:end)) sum(agrid(114:end).*DBN_agrid_bench(114:end)) ...
              sum(agrid(114:end).*DBN_agrid_bench(114:end))/sum(agrid.*DBN_agrid_bench) ] ;
    result_agrid_2 = [170 agrid(170) sum(DBN_agrid_bench(170:end)) sum(agrid(170:end).*DBN_agrid_bench(170:end)) ...
              sum(agrid(170:end).*DBN_agrid_bench(170:end))/sum(agrid.*DBN_agrid_bench) ] ;

    DBN_agrid_exp = sum(sum(sum(sum(DBN_exp(:,:,:,:,:),5),4),3),1) ;
    result_agrid_3 = [101 agrid(101) sum(DBN_agrid_exp(101:end)) sum(agrid(101:end).*DBN_agrid_exp(101:end)) ...
              sum(agrid(101:end).*DBN_agrid_exp(101:end))/sum(agrid.*DBN_agrid_exp) ] ;
    result_agrid_4 = [151 agrid(151) sum(DBN_agrid_exp(151:end)) sum(agrid(151:end).*DBN_agrid_exp(151:end)) ...
              sum(agrid(151:end).*DBN_agrid_exp(151:end))/sum(agrid.*DBN_agrid_exp) ] ;

% Asset distribution NewBorns
    DBN_NB_bench = sum(sum(sum(DBN_bench(1,:,:,:,:),5),4),3)/sum(sum(sum(sum(DBN_bench(1,:,:,:,:))))) ;
    DBN_NB_az_bench = squeeze(sum(sum(DBN_bench(1,:,:,:,:),5),4)/sum(sum(sum(sum(DBN_bench(1,:,:,:,:)))))) ;
    DBN_NB_az_vec_bench = DBN_NB_az_bench(:) ;
    
% Assets by age-agrid-zgrid
    A_az_bench = NaN(Max_Age,n_a,n_z) ;
    Q_az_bench = NaN(Max_Age,n_a,n_z) ;
    A_az_exp   = NaN(Max_Age,n_a,n_z) ;
    Q_az_exp   = NaN(Max_Age,n_a,n_z) ;
    
    W_A_az_bench = NaN(Max_Age,n_a,n_z) ;
    W_Q_az_bench = NaN(Max_Age,n_a,n_z) ;
    W_A_az_exp   = NaN(Max_Age,n_a,n_z) ;
    W_Q_az_exp   = NaN(Max_Age,n_a,n_z) ;
    for l=1:Max_Age
        for i=1:n_a
            for j=1:n_z
                A_az_bench(l,i,j) =          agrid(i)*sum(sum(DBN_bench(l,i,j,:,:),5),4) ;
                Q_az_bench(l,i,j) = zgrid(j)*agrid(i)*sum(sum(DBN_bench(l,i,j,:,:),5),4) ;
                A_az_exp(l,i,j)   =          agrid(i)*sum(sum(DBN_exp(l,i,j,:,:),5),4) ;
                Q_az_exp(l,i,j)   = zgrid(j)*agrid(i)*sum(sum(DBN_exp(l,i,j,:,:),5),4) ;
            end 
        end 
        W_A_az_bench(l,:,:) = A_az_bench(l,:,:)./sum(sum(A_az_bench(l,:,:))) ;
        W_Q_az_bench(l,:,:) = Q_az_bench(l,:,:)./sum(sum(Q_az_bench(l,:,:))) ;
        W_A_az_exp(l,:,:)   = A_az_exp(l,:,:)./sum(sum(A_az_exp(l,:,:)))     ;
        W_Q_az_exp(l,:,:)   = Q_az_exp(l,:,:)./sum(sum(Q_az_exp(l,:,:)))     ;
    end 
    
    % Weight of each age
    W_age = sum(sum(sum(sum(DBN_bench,5),4),3),2) ;

%% Return on assets

    Return_bench_mat = 1 + ( mu * R_bench * Z_mat.^mu .* A_mat.^(mu-1) - delta )*(1-tauK) ;
    Return_exp_mat   = ( 1 + mu * R_bench * Z_mat.^mu .* A_mat.^(mu-1) - delta )*(1-tauW) ;
    Return_mat       =  1 + mu * R_bench * Z_mat.^mu .* A_mat.^(mu-1) - delta             ;
    
    Return_bench = squeeze(Return_bench_mat(1,:,:,1,1)) ;
    Return_vec_bench = Return_bench(:) ;
    [Return_vec_bench ,ind_bench] = sort(Return_vec_bench) ;
    
    Return_exp = squeeze(Return_exp_mat(1,:,:,1,1)) ;
    Return_vec_exp = Return_exp(:) ;
    [Return_vec_exp ,ind_exp] = sort(Return_vec_exp) ;
    
    Return_bt = squeeze(Return_mat(1,:,:,1,1)) ;
    Return_vec_bt = Return_bt(:) ;
    [Return_vec_bt ,ind_bt] = sort(Return_vec_bt) ;


%% Unweighted Return distribution

% Return distribution by age 
    DBN_return_bench     = NaN(n_a*n_z,Max_Age)       ;
    DBN_return_bench_bt  = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_bench     = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_bench_bt  = NaN(n_a*n_z,Max_Age)       ;
    mean_return_bench    = NaN(Max_Age,1)             ;
    mean_return_bench_bt = NaN(Max_Age,1)             ;
    std_return_bench     = NaN(Max_Age,1)             ;
    std_return_bench_bt  = NaN(Max_Age,1)             ;
    
    mean_return_AZ_bench    = NaN(Max_Age,n_z)        ;
    mean_return_AZ_bench_bt = NaN(Max_Age,n_z)        ;
    std_return_AZ_bench     = NaN(Max_Age,n_z)        ;
    std_return_AZ_bench_bt  = NaN(Max_Age,n_z)        ;
    
    DBN_return_exp     = NaN(n_a*n_z,Max_Age)       ;
    DBN_return_exp_bt  = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_exp     = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_exp_bt  = NaN(n_a*n_z,Max_Age)       ;
    mean_return_exp    = NaN(Max_Age,1)             ;
    mean_return_exp_bt = NaN(Max_Age,1)             ;
    std_return_exp     = NaN(Max_Age,1)             ;
    std_return_exp_bt  = NaN(Max_Age,1)             ;
    
    mean_return_AZ_exp    = NaN(Max_Age,n_z)        ;
    mean_return_AZ_exp_bt = NaN(Max_Age,n_z)        ;
    std_return_AZ_exp     = NaN(Max_Age,n_z)        ;
    std_return_AZ_exp_bt  = NaN(Max_Age,n_z)        ;
    
    for i=1:Max_Age
        DBN_az_age = squeeze(sum(sum(DBN_bench(i,:,:,:,:),5),4)/sum(sum(sum(sum(DBN_bench(i,:,:,:,:)))))) ;
        DBN_az_vec = DBN_az_age(:) ;
        DBN_return_bench(:,i)   = DBN_az_vec(ind_bench)   ;
        DBN_return_bench_bt(:,i)= DBN_az_vec(ind_bt)      ;
        CDF_return_bench(:,i)   = cumsum(DBN_return_bench(:,i)) ;
        CDF_return_bench_bt(:,i)= cumsum(DBN_return_bench_bt(:,i)) ;
        mean_return_bench(i)    = sum(Return_vec_bench.*DBN_return_bench(:,i)) ;
        mean_return_bench_bt(i) = sum(Return_vec_bt.*DBN_return_bench_bt(:,i)) ;
        std_return_bench(i)     = (sum((Return_vec_bench-mean_return_bench(i)).^2 .* DBN_return_bench(:,i)) )^.5 ;
        std_return_bench_bt(i)  = (sum((Return_vec_bt-mean_return_bench_bt(i)).^2 .* DBN_return_bench_bt(:,i)) )^.5 ;
        
        DBN_az_age = squeeze(sum(sum(DBN_exp(i,:,:,:,:),5),4)/sum(sum(sum(sum(DBN_exp(i,:,:,:,:)))))) ;
        DBN_az_vec = DBN_az_age(:) ;
        DBN_return_exp(:,i)   = DBN_az_vec(ind_exp)   ;
        DBN_return_exp_bt(:,i)= DBN_az_vec(ind_bt)      ;
        CDF_return_exp(:,i)   = cumsum(DBN_return_exp(:,i)) ;
        CDF_return_exp_bt(:,i)= cumsum(DBN_return_exp_bt(:,i)) ;
        mean_return_exp(i)    = sum(Return_vec_exp.*DBN_return_exp(:,i)) ;
        mean_return_exp_bt(i) = sum(Return_vec_bt.*DBN_return_exp_bt(:,i)) ;
        std_return_exp(i)     = (sum((Return_vec_exp-mean_return_exp(i)).^2 .* DBN_return_exp(:,i)) )^.5 ;
        std_return_exp_bt(i)  = (sum((Return_vec_bt-mean_return_exp_bt(i)).^2 .* DBN_return_exp_bt(:,i)) )^.5 ;
        
        for j=1:n_z
            Return_bench_aux = squeeze(Return_bench_mat(1,:,j,1,1)) ;
            Return_vec_AZ_bench = Return_bench_aux(:) ;
            [Return_vec_AZ_bench ,ind_AZ_bench] = sort(Return_vec_AZ_bench) ;
            
            Return_exp_aux = squeeze(Return_exp_mat(1,:,j,1,1)) ;
            Return_vec_AZ_exp = Return_exp_aux(:) ;
            [Return_vec_AZ_exp ,ind_AZ_exp] = sort(Return_vec_AZ_exp) ;
            
            Return_bt_aux = squeeze(Return_mat(1,:,j,1,1)) ;
            Return_vec_AZ_bt = Return_bt_aux(:) ;
            [Return_vec_AZ_bt ,ind_AZ_bt] = sort(Return_vec_AZ_bt) ;
            
            DBN_az_age = squeeze(sum(sum(DBN_bench(i,:,j,:,:),5),4)/sum(sum(sum(DBN_bench(i,:,j,:,:))))) ;
            DBN_az_vec = DBN_az_age(:) ;
            DBN_return_AZ_bench(:,i,j)   = DBN_az_vec(ind_AZ_bench)   ;
            DBN_return_AZ_bench_bt(:,i,j)= DBN_az_vec(ind_AZ_bt)      ;
            
            DBN_az_age = squeeze(sum(sum(DBN_exp(i,:,j,:,:),5),4)/sum(sum(sum(DBN_exp(i,:,j,:,:))))) ;
            DBN_az_vec = DBN_az_age(:) ;
            DBN_return_AZ_exp(:,i,j)   = DBN_az_vec(ind_AZ_exp)   ;
            DBN_return_AZ_exp_bt(:,i,j)= DBN_az_vec(ind_AZ_bt)      ;
            
            mean_return_AZ_bench(i,j)    = sum(Return_vec_AZ_bench.*DBN_return_AZ_bench(:,i,j)) ;
            mean_return_AZ_bench_bt(i,j) = sum(Return_vec_AZ_bt.*DBN_return_AZ_bench_bt(:,i,j)) ;
            std_return_AZ_bench(i,j)     = (sum((Return_vec_AZ_bench-mean_return_AZ_bench(i,j)).^2 .* DBN_return_AZ_bench(:,i,j)) )^.5 ;
            std_return_AZ_bench_bt(i,j)  = (sum((Return_vec_AZ_bt-mean_return_AZ_bench_bt(i,j)).^2 .* DBN_return_AZ_bench_bt(:,i,j)) )^.5 ;
            
            mean_return_AZ_exp(i,j)    = sum(Return_vec_AZ_exp.*DBN_return_AZ_exp(:,i,j)) ;
            mean_return_AZ_exp_bt(i,j) = sum(Return_vec_AZ_bt.*DBN_return_AZ_exp_bt(:,i,j)) ;
            std_return_AZ_exp(i,j)     = (sum((Return_vec_AZ_exp-mean_return_AZ_exp(i,j)).^2 .* DBN_return_AZ_exp(:,i,j)) )^.5 ;
            std_return_AZ_exp_bt(i,j)  = (sum((Return_vec_AZ_bt-mean_return_AZ_exp_bt(i,j)).^2 .* DBN_return_AZ_exp_bt(:,i,j)) )^.5 ;
        end 
    end 
    
    figure; 
        subplot(2,2,1); plot(20:100,[mean_return_bench_bt mean_return_exp_bt]); xlim([19+1,19+Max_Age]); title('Mean Return Before Tax')
        subplot(2,2,2); plot(20:100,[std_return_bench_bt  std_return_exp_bt ]); xlim([19+1,19+Max_Age]); title('Std Return Before Tax')
        subplot(2,2,3); plot(20:100,[mean_return_bench mean_return_exp]); xlim([19+1,19+Max_Age]); title('Mean Return After Tax')
        subplot(2,2,4); plot(20:100,[std_return_bench  std_return_exp ]); xlim([19+1,19+Max_Age]); title('Std Return After Tax')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','Return_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,mean_return_AZ_bench_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return Before Tax bench')
        subplot(2,2,2); plot(20:100,std_return_AZ_bench_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return Before Tax bench')
        subplot(2,2,3); plot(20:100,mean_return_AZ_exp_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return Before Tax exp')
        subplot(2,2,4); plot(20:100,std_return_AZ_exp_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return Before Tax exp')
        legend('z2','z4','z6','location','northeast')
    print('-dpdf','BT_Return_Life_Cycle_by_z.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,mean_return_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return After Tax bench')
        subplot(2,2,2); plot(20:100,std_return_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return After Tax bench')
        subplot(2,2,3); plot(20:100,mean_return_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return After Tax exp')
        subplot(2,2,4); plot(20:100,std_return_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return After Tax exp')
        legend('z2','z4','z6','location','northeast')
    print('-dpdf','AT_Return_Life_Cycle_by_z.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bench,DBN_return_bench(:,age(i)),Return_vec_exp,DBN_return_exp(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','Return_after_tax_DBN_by_age.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bench,CDF_return_bench(:,age(i)),Return_vec_exp,CDF_return_exp(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','Return_after_tax_CDF_by_age.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bt,DBN_return_bench_bt(:,age(i)),Return_vec_exp,DBN_return_exp_bt(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','Return_before_tax_DBN_by_age.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bt,CDF_return_bench_bt(:,age(i)),Return_vec_exp,CDF_return_exp_bt(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','Return_before_tax_CDF_by_age.pdf') ;
    
    prctl = [10:10:90, 95, 99, 99.9];
    Return_prctl_bench = NaN(Max_Age,numel(prctl)) ;
    Return_prctl_exp   = NaN(Max_Age,numel(prctl)) ;
    Return_prctl_bench_bt = NaN(Max_Age,numel(prctl)) ;
    Return_prctl_exp_bt   = NaN(Max_Age,numel(prctl)) ;
    for i=1:Max_Age
        for j=1:numel(prctl)
           [~,ind] = max(CDF_return_bench(:,i)>=prctl(j)/100);
           Return_prctl_bench(i,j) = 100*(Return_vec_bench(ind)-1) ;
           [~,ind] = max(CDF_return_exp(:,i)>=prctl(j)/100);
           Return_prctl_exp(i,j) = 100*(Return_vec_exp(ind)-1) ;
           
           [~,ind] = max(CDF_return_bench_bt(:,i)>=prctl(j)/100);
           Return_prctl_bench_bt(i,j) = 100*(Return_vec_bt(ind)-1) ;
           [~,ind] = max(CDF_return_exp_bt(:,i)>=prctl(j)/100);
           Return_prctl_exp_bt(i,j) = 100*(Return_vec_bt(ind)-1) ;
        end
    end 
    
    Return_prctl_diff = Return_prctl_exp - Return_prctl_bench ;
    Return_prctl_diff_bt = Return_prctl_exp_bt - Return_prctl_bench_bt ;
    
    col_title{1} = 'Age' ;
    for j=1:numel(prctl)
        col_title{j+1} = strcat('p',num2str(prctl(j))) ;
    end
    row_title = 20:100;
    
    Mat_bench = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_bench])] ;
            
    Mat_exp = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_exp])] ;
            
    Mat_diff = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_diff])] ;
            
    status = xlwrite('Return_Stat_after_tax.xls',Mat_bench,'Prctl_bench') ;
    status = xlwrite('Return_Stat_after_tax.xls',Mat_exp,'Prctl_exp') ;
    status = xlwrite('Return_Stat_after_tax.xls',Mat_diff,'Prctl_diff') ;
    
    
    Mat_bench = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_bench_bt])] ;
            
    Mat_exp = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_exp_bt])] ;
            
    Mat_diff = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_diff_bt])] ;
            
    status = xlwrite('Return_Stat_before_tax.xls',Mat_bench,'Prctl_bench') ;
    status = xlwrite('Return_Stat_before_tax.xls',Mat_exp,'Prctl_exp') ;
    status = xlwrite('Return_Stat_before_tax.xls',Mat_diff,'Prctl_diff') ;
    

    
%% Weighted Return distribution
% Return distribution by age 
    DBN_return_bench     = NaN(n_a*n_z,Max_Age)       ;
    DBN_return_bench_bt  = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_bench     = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_bench_bt  = NaN(n_a*n_z,Max_Age)       ;
    w_mean_return_bench    = NaN(Max_Age,1)             ;
    w_mean_return_bench_bt = NaN(Max_Age,1)             ;
    w_std_return_bench     = NaN(Max_Age,1)             ;
    w_std_return_bench_bt  = NaN(Max_Age,1)             ;
    
    w_mean_return_AZ_bench    = NaN(Max_Age,n_z)        ;
    w_mean_return_AZ_bench_bt = NaN(Max_Age,n_z)        ;
    w_std_return_AZ_bench     = NaN(Max_Age,n_z)        ;
    w_std_return_AZ_bench_bt  = NaN(Max_Age,n_z)        ;
    
    DBN_return_exp     = NaN(n_a*n_z,Max_Age)       ;
    DBN_return_exp_bt  = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_exp     = NaN(n_a*n_z,Max_Age)       ;
    CDF_return_exp_bt  = NaN(n_a*n_z,Max_Age)       ;
    w_mean_return_exp    = NaN(Max_Age,1)             ;
    w_mean_return_exp_bt = NaN(Max_Age,1)             ;
    w_std_return_exp     = NaN(Max_Age,1)             ;
    w_std_return_exp_bt  = NaN(Max_Age,1)             ;
    
    w_mean_return_AZ_exp    = NaN(Max_Age,n_z)        ;
    w_mean_return_AZ_exp_bt = NaN(Max_Age,n_z)        ;
    w_std_return_AZ_exp     = NaN(Max_Age,n_z)        ;
    w_std_return_AZ_exp_bt  = NaN(Max_Age,n_z)        ;
    
    for i=1:Max_Age
        DBN_az_age = squeeze(W_A_az_bench(i,:,:)) ;
        DBN_az_vec = DBN_az_age(:) ;
        DBN_return_bench(:,i)   = DBN_az_vec(ind_bench)   ;
        DBN_return_bench_bt(:,i)= DBN_az_vec(ind_bt)      ;
        CDF_return_bench(:,i)   = cumsum(DBN_return_bench(:,i)) ;
        CDF_return_bench_bt(:,i)= cumsum(DBN_return_bench_bt(:,i)) ;
        w_mean_return_bench(i)    = sum(Return_vec_bench.*DBN_return_bench(:,i)) ;
        w_mean_return_bench_bt(i) = sum(Return_vec_bt.*DBN_return_bench_bt(:,i)) ;
        w_std_return_bench(i)     = (sum((Return_vec_bench-w_mean_return_bench(i)).^2 .* DBN_return_bench(:,i)) )^.5 ;
        w_std_return_bench_bt(i)  = (sum((Return_vec_bt   -w_mean_return_bench_bt(i)).^2 .* DBN_return_bench_bt(:,i)) )^.5 ;
        
        DBN_az_age = squeeze(W_A_az_exp(i,:,:))  ;
        DBN_az_vec = DBN_az_age(:) ;
        DBN_return_exp(:,i)   = DBN_az_vec(ind_exp)   ;
        DBN_return_exp_bt(:,i)= DBN_az_vec(ind_bt)      ;
        CDF_return_exp(:,i)   = cumsum(DBN_return_exp(:,i)) ;
        CDF_return_exp_bt(:,i)= cumsum(DBN_return_exp_bt(:,i)) ;
        w_mean_return_exp(i)    = sum(Return_vec_exp.*DBN_return_exp(:,i)) ;
        w_mean_return_exp_bt(i) = sum(Return_vec_bt.*DBN_return_exp_bt(:,i)) ;
        w_std_return_exp(i)     = (sum((Return_vec_exp-w_mean_return_exp(i)).^2 .* DBN_return_exp(:,i)) )^.5 ;
        w_std_return_exp_bt(i)  = (sum((Return_vec_bt-w_mean_return_exp_bt(i)).^2 .* DBN_return_exp_bt(:,i)) )^.5 ;
        
        for j=1:n_z
            Return_bench_aux = squeeze(Return_bench_mat(1,:,j,1,1)) ;
            Return_vec_AZ_bench = Return_bench_aux(:) ;
            [Return_vec_AZ_bench ,ind_AZ_bench] = sort(Return_vec_AZ_bench) ;
            
            Return_exp_aux = squeeze(Return_exp_mat(1,:,j,1,1)) ;
            Return_vec_AZ_exp = Return_exp_aux(:) ;
            [Return_vec_AZ_exp ,ind_AZ_exp] = sort(Return_vec_AZ_exp) ;
            
            Return_bt_aux = squeeze(Return_mat(1,:,j,1,1)) ;
            Return_vec_AZ_bt = Return_bt_aux(:) ;
            [Return_vec_AZ_bt ,ind_AZ_bt] = sort(Return_vec_AZ_bt) ;
            
            DBN_az_age = squeeze(W_A_az_bench(i,:,j)/sum(W_A_az_bench(i,:,j))) ;
            DBN_az_vec = DBN_az_age(:) ;
            DBN_return_AZ_bench(:,i,j)   = DBN_az_vec(ind_AZ_bench)   ;
            DBN_return_AZ_bench_bt(:,i,j)= DBN_az_vec(ind_AZ_bt)      ;
            
            DBN_az_age = squeeze(W_A_az_exp(i,:,j)/sum(W_A_az_exp(i,:,j))) ;
            DBN_az_vec = DBN_az_age(:) ;
            DBN_return_AZ_exp(:,i,j)   = DBN_az_vec(ind_AZ_exp)   ;
            DBN_return_AZ_exp_bt(:,i,j)= DBN_az_vec(ind_AZ_bt)      ;
            
            w_mean_return_AZ_bench(i,j)    = sum(Return_vec_AZ_bench.*DBN_return_AZ_bench(:,i,j)) ;
            w_mean_return_AZ_bench_bt(i,j) = sum(Return_vec_AZ_bt.*DBN_return_AZ_bench_bt(:,i,j)) ;
            w_std_return_AZ_bench(i,j)     = (sum((Return_vec_AZ_bench-w_mean_return_AZ_bench(i,j)).^2 .* DBN_return_AZ_bench(:,i,j)) )^.5 ;
            w_std_return_AZ_bench_bt(i,j)  = (sum((Return_vec_AZ_bt-w_mean_return_AZ_bench_bt(i,j)).^2 .* DBN_return_AZ_bench_bt(:,i,j)) )^.5 ;
            
            w_mean_return_AZ_exp(i,j)    = sum(Return_vec_AZ_exp.*DBN_return_AZ_exp(:,i,j)) ;
            w_mean_return_AZ_exp_bt(i,j) = sum(Return_vec_AZ_bt.*DBN_return_AZ_exp_bt(:,i,j)) ;
            w_std_return_AZ_exp(i,j)     = (sum((Return_vec_AZ_exp-w_mean_return_AZ_exp(i,j)).^2 .* DBN_return_AZ_exp(:,i,j)) )^.5 ;
            w_std_return_AZ_exp_bt(i,j)  = (sum((Return_vec_AZ_bt-w_mean_return_AZ_exp_bt(i,j)).^2 .* DBN_return_AZ_exp_bt(:,i,j)) )^.5 ;
        end 
    end 
    
    figure; 
        subplot(2,2,1); plot(20:100,[w_mean_return_bench_bt w_mean_return_exp_bt]); xlim([19+1,19+Max_Age]); title('Mean Return Before Tax')
        subplot(2,2,2); plot(20:100,[w_std_return_bench_bt  w_std_return_exp_bt ]); xlim([19+1,19+Max_Age]); title('Std Return Before Tax')
        subplot(2,2,3); plot(20:100,[w_mean_return_bench w_mean_return_exp]); xlim([19+1,19+Max_Age]); title('Mean Return After Tax')
        subplot(2,2,4); plot(20:100,[w_std_return_bench  w_std_return_exp ]); xlim([19+1,19+Max_Age]); title('Std Return After Tax')
        legend('Bench','Exp','location','southeast')
    print('-dpdf','W_Return_Life_Cycle.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,w_mean_return_AZ_bench_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return Before Tax bench')
        subplot(2,2,2); plot(20:100,w_std_return_AZ_bench_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return Before Tax bench')
        subplot(2,2,3); plot(20:100,w_mean_return_AZ_exp_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return Before Tax exp')
        subplot(2,2,4); plot(20:100,w_std_return_AZ_exp_bt(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return Before Tax exp')
        legend('z2','z4','z6','location','northeast')
    print('-dpdf','W_BT_Return_Life_Cycle_by_z.pdf') ;
    
    figure; 
        subplot(2,2,1); plot(20:100,w_mean_return_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return After Tax bench')
        subplot(2,2,2); plot(20:100,w_std_return_AZ_bench(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return After Tax bench')
        subplot(2,2,3); plot(20:100,w_mean_return_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Mean Return After Tax exp')
        subplot(2,2,4); plot(20:100,w_std_return_AZ_exp(:,z_vec)); xlim([19+1,19+Max_Age]); title('Std Return After Tax exp')
        legend('z2','z4','z6','location','northeast')
    print('-dpdf','W_AT_Return_Life_Cycle_by_z.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bench,DBN_return_bench(:,age(i)),Return_vec_exp,DBN_return_exp(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','W_Return_after_tax_DBN_by_age.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bench,CDF_return_bench(:,age(i)),Return_vec_exp,CDF_return_exp(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','W_Return_after_tax_CDF_by_age.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bt,DBN_return_bench_bt(:,age(i)),Return_vec_exp,DBN_return_exp_bt(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','W_Return_before_tax_DBN_by_age.pdf') ;
    
    figure;
    for i=1:n_age
       subplot(2,3,i); 
       hold on;
       plot(Return_vec_bt,CDF_return_bench_bt(:,age(i)),Return_vec_exp,CDF_return_exp_bt(:,age(i)))
       title(strcat('Age=',int2str(19+age(i)))); xlim([1,1.5]);
       hold off;
    end
    legend('bench','exp','location','southeast');
    print('-dpdf','W_Return_before_tax_CDF_by_age.pdf') ;
    
    prctl = [10:10:90, 95, 99, 99.9];
    Return_prctl_bench = NaN(Max_Age,numel(prctl)) ;
    Return_prctl_exp   = NaN(Max_Age,numel(prctl)) ;
    Return_prctl_bench_bt = NaN(Max_Age,numel(prctl)) ;
    Return_prctl_exp_bt   = NaN(Max_Age,numel(prctl)) ;
    for i=1:Max_Age
        for j=1:numel(prctl)
           [~,ind] = max(CDF_return_bench(:,i)>=prctl(j)/100);
           Return_prctl_bench(i,j) = 100*(Return_vec_bench(ind)-1) ;
           [~,ind] = max(CDF_return_exp(:,i)>=prctl(j)/100);
           Return_prctl_exp(i,j) = 100*(Return_vec_exp(ind)-1) ;
           
           [~,ind] = max(CDF_return_bench_bt(:,i)>=prctl(j)/100);
           Return_prctl_bench_bt(i,j) = 100*(Return_vec_bt(ind)-1) ;
           [~,ind] = max(CDF_return_exp_bt(:,i)>=prctl(j)/100);
           Return_prctl_exp_bt(i,j) = 100*(Return_vec_bt(ind)-1) ;
        end
    end 
    
    Return_prctl_diff = Return_prctl_exp - Return_prctl_bench ;
    Return_prctl_diff_bt = Return_prctl_exp_bt - Return_prctl_bench_bt ;
    
    col_title{1} = 'Age' ;
    for j=1:numel(prctl)
        col_title{j+1} = strcat('p',num2str(prctl(j))) ;
    end
    row_title = 20:100;
    
    Mat_bench = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_bench])] ;
            
    Mat_exp = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_exp])] ;
            
    Mat_diff = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_diff])] ;
            
    status = xlwrite('W_Return_Stat_after_tax.xls',Mat_bench,'Prctl_bench') ;
    status = xlwrite('W_Return_Stat_after_tax.xls',Mat_exp,'Prctl_exp') ;
    status = xlwrite('W_Return_Stat_after_tax.xls',Mat_diff,'Prctl_diff') ;
    
    
    Mat_bench = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_bench_bt])] ;
            
    Mat_exp = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_exp_bt])] ;
            
    Mat_diff = [cell(1,numel(prctl)+1);
                col_title;
                num2cell([row_title' , Return_prctl_diff_bt])] ;
            
    status = xlwrite('W_Return_Stat_before_tax.xls',Mat_bench,'Prctl_bench') ;
    status = xlwrite('W_Return_Stat_before_tax.xls',Mat_exp,'Prctl_exp') ;
    status = xlwrite('W_Return_Stat_before_tax.xls',Mat_diff,'Prctl_diff') ;
    
%% 
% Age brackets 
    age = [5 , 15, 25, 35, 45, 55, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;
    
% Hours
    w_R_AZ_bench  = zeros(n_age,n_z) ;
    w_R_AZ_exp    = zeros(n_age,n_z) ;
    R_AZ_bench  = zeros(n_age,n_z) ;
    R_AZ_exp    = zeros(n_age,n_z) ;
    Weights_w_R_bench = zeros(n_age,n_z) ;
    Weights_w_R_exp = zeros(n_age,n_z) ;
    for i=1:n_z
        j=1;
        for l=1:Max_Age
            if l<=age(j)
            w_R_AZ_bench(j,i) = w_R_AZ_bench(j,i) + sum(Return_bt(:,i).*squeeze(A_az_bench(l,:,i))') ;
            w_R_AZ_exp(j,i)   = w_R_AZ_exp(j,i)   + sum(Return_bt(:,i).*squeeze(A_az_exp(l,:,i))') ;
            R_AZ_bench(j,i)   = R_AZ_bench(j,i)   + sum(Return_bt(:,i).*squeeze(sum(sum(DBN_bench(l,:,i,:,:),5),4))') ;
            R_AZ_exp(j,i)     = R_AZ_exp(j,i)     + sum(Return_bt(:,i).*squeeze(sum(sum(DBN_exp(l,:,i,:,:),5),4))') ;
            
            Weights_w_R_bench(j,i)  = Weights_w_R_bench(j,i) + sum(squeeze(A_az_bench(l,:,i))) ;
            Weights_w_R_exp(j,i)    = Weights_w_R_exp(j,i) + sum(squeeze(A_az_exp(l,:,i)))   ;
            end 
            if l==age(j); j=j+1; end 
        end
    end

    w_R_AZ_bench  = 100*(w_R_AZ_bench./Weights_w_R_bench-1) ;
    w_R_AZ_exp    = 100*(w_R_AZ_exp./Weights_w_R_exp-1)   ;
    R_AZ_bench    = 100*(100*R_AZ_bench./Weights_AZ-1)   ;
    R_AZ_exp      = 100*(100*R_AZ_exp./Weights_AZ-1)     ;
    
    w_R_AZ_diff   = w_R_AZ_exp - w_R_AZ_bench  ;
    R_AZ_diff     = R_AZ_exp   - R_AZ_bench    ;
    
    
    Mat_1 = [{'w_R_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(w_R_AZ_bench)] ;
    Mat_2 = [{'w_R_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(w_R_AZ_exp)] ;
    Mat_3 = [{'w_R_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(w_R_AZ_diff)] ;
    Mat_4 = [{'R_bench',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(R_AZ_bench)] ;
    Mat_5 = [{'R_exp',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(R_AZ_exp)] ;
    Mat_6 = [{'R_diff',' ',' ',' ',' ',' ',' ',' '};z_title;age_title,num2cell(R_AZ_diff)] ;
    Mat_7 = cell(1,8*3+3) ;
    Mat_8 = cell(9,1) ;
    
   
    Mat = [Mat_7;
           Mat_8 , Mat_1 , Mat_8 , Mat_2 , Mat_8 , Mat_3 ;
           Mat_7;
           Mat_8 , Mat_4 , Mat_8 , Mat_5 , Mat_8 , Mat_6 ] ;
    status = xlwrite(xls_file,Mat,'Return') ;
    
    
    
    
%% Tables by age-return
% For each age group I divide the population into segments according to the
% return on assets 

% Age brackets 
    age = [5 , 15, 25, 35, 45, 55, Max_Age ];
    %%%    24, 34, 44, 54, 64, 74, 80
    n_age = numel(age) ;

% Distribution of returns by age group 

    DBN_AG_bench = zeros(n_age,n_a,n_z) ;
    DBN_AG_exp   = zeros(n_age,n_a,n_z) ;
    
    j=1;
    for l=1:Max_Age
        if l<=age(j)
        DBN_AG_bench(j,:,:) = DBN_AG_bench(j,:,:) + sum(sum(DBN_bench(i,:,:,:,:),5),4) ;
        DBN_AG_exp(j,:,:)   = DBN_AG_exp(j,:,:)   + sum(sum(DBN_exp(i,:,:,:,:),5),4)   ;
        end 
        if l==age(j); j=j+1; end 
    end
    
    DBN_return_AG_bench  = NaN(n_a*n_z,n_age)       ;
    CDF_return_AG_bench  = NaN(n_a*n_z,n_age)       ;
    
    DBN_return_AG_exp    = NaN(n_a*n_z,n_age)       ;
    CDF_return_AG_exp    = NaN(n_a*n_z,n_age)       ;
    
    for i=1:n_age
        DBN_az_age = squeeze(DBN_AG_bench(i,:,:)/sum(sum(DBN_AG_bench(i,:,:)))) ;
        DBN_az_vec = DBN_az_age(:) ;
        DBN_return_AG_bench(:,i)   = DBN_az_vec(ind_bench)   ;
        CDF_return_AG_bench(:,i)   = cumsum(DBN_return_AG_bench(:,i)) ;

        DBN_az_age = squeeze(DBN_AG_exp(i,:,:)/sum(sum(DBN_AG_exp(i,:,:)))) ;
        DBN_az_vec = DBN_az_age(:) ;
        DBN_return_AG_exp(:,i)  = DBN_az_vec(ind_exp)   ;
        CDF_return_AG_exp(:,i)  = cumsum(DBN_return_AG_exp(:,i)) ;
    end 



