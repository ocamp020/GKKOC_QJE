
function [A_mat,C_mat,H_mat,Z_mat,R_mat,R_at_mat,e_mat,age_mat] = Model_Simulation(sw,rr,wage,C_pf,DBN,agrid,zgrid)

%% Parameters
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
	
% Set Tax Regime
    % tauPL = 0.185     ;
    % psi   = 0.77      ;
    tauPL = 0.0         ;
    psi   = 0.776       ;
    
if sw==1    
    tauK  = 0.25 ;
    tauW  = 0.00 ;
else
    tauK  = 0.00 ;
    tauW  = 0.017072675596579098 ;
end 
        
%% Aggregates 
    % Weights by z
        Weights_Z = NaN(1,n_z) ;
        for i=1:n_z
            Weights_Z(i) = sum(sum(sum(sum(DBN(:,:,i,:,:))))) ;
        end 
        Weights_Z = cumsum(Weights_Z) ;
        
    % Weights by lambda
        Weights_L = NaN(1,n_z) ;
        for i=1:n_l
            Weights_L(i) = sum(sum(sum(sum(DBN(:,:,:,i,:))))) ;
        end 
        Weights_L = cumsum(Weights_L) ;
        
    % Weights by initial A
        Weights_A = NaN(1,n_a) ;
        for i=1:n_a
            Weights_A(i) = sum(sum(sum(DBN(1,i,:,:,:)))) ;
        end 
        Weights_A = cumsum(Weights_A) ;
        
%% Result_Folders

    if ((tauPL==0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
    elseif ((tauPL~=0.0)&&(sigma~=1.0))  
        Result_Folder = strcat('../NSU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
    elseif ((tauPL==0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_LT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
    elseif ((tauPL~=0.0)&&(sigma==1.0))  
        Result_Folder = strcat('../SU_PT_Results/Factor_',num2str(Threshold_Factor,'%.2f'),'/') ;
    end
    
    eval(['load ',Result_Folder,'P_e']) ;
    P_e = reshape(P_e,[n_e,n_e]) ;
    
    eval(['load ',Result_Folder,'eff_un']) ;
    eff_un = reshape(eff_un,[Max_Age,n_l,n_e]) ;
    
    eval(['load ',Result_Folder,'survP']) ;
    survP = survP' ;
    
    eval(['load ',Result_Folder,'Ret_Y']) ;
    Ret_Y = reshape(Ret_Y,[n_l,n_e]) ;
    
    for i=1:Max_Age
        eff_un_age(i) = sum(sum((squeeze(eff_un(i,:,:)).*reshape(sum(sum(DBN_bench(i,:,:,:,:))),[5,5]))))/sum(sum(sum(sum(DBN_bench(i,:,:,:,:)))));
    end 
        

%% Variables for simulation
    n_agents = 10000 ; 
    T        = 1000  ;
    
    A_mat = NaN(T,n_agents)     ;
    C_mat = NaN(T,n_agents)     ;
    H_mat = NaN(T,n_agents)     ;
    R_mat = NaN(T,n_agents)     ;
    R_at_mat = NaN(T,n_agents)  ;
    Z_mat = NaN(T,n_agents)     ;
    e_mat = NaN(T,n_agents)     ;
    age_mat = NaN(T,n_agents)   ;
    
%% Random Numbers    
rng(1000)

	% Assignment of investment productivity 
    ran_Z = rand(n_agents,1) ;
    % Assignment of labor productivity
    ran_L = rand(n_agents,1) ;
    % Assignment of initial assets
    ran_A = rand(n_agents,1) ;
    % Assignment of transitory labor efficiency
    ran_e = rand(T,n_agents) ;
    % Death shock
    ran_d = rand(T,n_agents) ;


%% Simulation

for i=1:n_agents
    % Assgin initial conditions 
        % Agents are born with age 1
        
        % All agents are born with e=e(3)
            e_mat(1,i) = 3 ;
        % Investment productivity 
            [~,z_ind]  = max(ran_Z(i)<=Weights_Z) ;
            Z_mat(i) = zgrid(z_ind) ;
        % Labor productivity 
            [~,l_ind]  = max(ran_L(i)<=Weights_L) ;
        % Initial Assets
            [~,ind]    = max(ran_A(i)<=Weights_A) ;
            A_mat(1,i) = agrid(ind) ;
            
	% Consumption
        C_mat(1,i) = C_pf(1,ind,z_ind,l_ind,e_mat(1,i))   ;
    % Labor
        H_mat(1,i) = H_FOC(C_mat(1,i),1,l_ind,e_mat(1,i)) ;
        
    % Returns
        R_mat(1,i)    = MB(A_mat(1,i),Z_mat(i)) ;
        R_at_mat(1,i) = MB_at(A_mat(1,i),Z_mat(i)) ;
    
    % Future Assets     
        A_mat(2,i) = Y_a(A_mat(1,i),Z_mat(i)) + Y_h(H_mat(1,i),1,l_ind,e_mat(1,i)) - C_mat(1,i) ;
    
        if A_mat(2,i)<agrid(1)
            A_mat(2,i) = agrid(1) ;                   
            H_mat(1,i) = H_FOC_aux(A_mat(1,i),Z_mat(i),1,l_ind,e_mat(1,i),A_mat(2,i)) ;
            C_mat(1,i) = Y_a(A_mat(1,i),Z_mat(i)) + Y_h(H_mat(1,i),1,l_ind,e_mat(1,i)) - A_mat(2,i) ;
        end 
        if C_mat(1,i)<=0; error('Negative consumption'); end 
        
    % Working period    
    for t=2:Ret_Age-1
        % Death shock 
        if ran_d>=survP(t)
            
        else 
        % Labor Efficiency
        [~,e_mat(t,i)] = max( ran_e(t,i)<=P_e(e_mat(t-1,i),:) ) ;
        % Consumption
        C_mat(t,i) = interp1(agrid,C_pf(t,:,z_ind,l_ind,e_mat(t,i)),A_mat(t,i)) ;
        % Labor
        H_mat(t,i) = H_FOC(C_mat(t,i),t,l_ind,e_mat(t,i)) ;
        % Returns
        R_mat(t,i)    = MB(A_mat(t,i),Z_mat(i)) ;
        R_at_mat(t,i) = MB_at(A_mat(t,i),Z_mat(i)) ; 
        % Future Assetes
        A_mat(t+1,i) = Y_a(A_mat(t,i),Z_mat(i)) + Y_h(H_mat(t,i),t,l_ind,e_mat(t-1,i)) - C_mat(t,i) ;
        
        if A_mat(t+1,i)<agrid(1)
            A_mat(t+1,i) = agrid(1) ;                   
            H_mat(t,i) = H_FOC_aux(A_mat(t,i),Z_mat(i),t,l_ind,e_mat(t,i),A_mat(t+1,i)) ;
            C_mat(t,i) = Y_a(A_mat(t,i),Z_mat(i)) + Y_h(H_mat(t,i),t,l_ind,e_mat(t,i)) - A_mat(t+1,i) ;
        end 
        if C_mat(t,i)<=0; error('Negative consumption'); end
        end
        
    end 
    
    % Index of last level of labor efficiency
    e_ind = e_mat(Ret_Age-1,i) ;
    
    % Retirement period
    for t=Ret_Age:Max_Age-1
        
        % Consumption
        C_mat(t,i) = interp1(agrid,C_pf(t,:,z_ind,l_ind,e_ind),A_mat(t,i)) ;
        % Labor
        H_mat(t,i) = 0 ;
        % Returns
        R_mat(t,i)    = MB(A_mat(t,i),Z_mat(i)) ;
        R_at_mat(t,i) = MB_at(A_mat(t,i),Z_mat(i)) ; 
        % Future Assetes
        A_mat(t+1,i) = Y_a(A_mat(t,i),Z_mat(i)) + Ret_Y(l_ind,e_ind) - C_mat(t,i) ;
        
        if A_mat(t+1,i)<agrid(1)
            A_mat(t+1,i) = agrid(1) ;                   
            C_mat(t,i) = Y_a(A_mat(t,i),Z_mat(i)) + Ret_Y(l_ind,e_ind) - A_mat(t+1,i) ;
        end 
        if C_mat(t,i)<=0; error('Negative consumption'); end 
    end
    
    % Final period
        % Consumption
        C_mat(Max_Age,i) = Y_a(A_mat(Max_Age,i),Z_mat(i)) + Ret_Y(l_ind,e_ind) ;
        % Labor
        H_mat(Max_Age,i) = 0 ;
        % Returns
        R_mat(Max_Age,i)    = MB(A_mat(Max_Age,i),Z_mat(i)) ;
        R_at_mat(Max_Age,i) = MB_at(A_mat(Max_Age,i),Z_mat(i)) ; 
    
end 


%% Auxiliary functions
function Y = Y_a(a_in,z_in)
    Y = ( a_in + ( rr * (z_in * a_in )^mu - delta*a_in ) *(1-tauK) ) * (1-tauW) ;
end 

function Y = Y_h(h_in,age_in,lambda_in,e_in)
    Y = psi*( wage*eff_un(age_in,lambda_in,e_in)*h_in)^(1-tauPL) ;
end 

function MB = MB(a_in,z_in)
   MB =   1 + ( rr * (z_in)^mu * a_in^(mu-1)  - delta) ;
end

function MB = MB_at(a_in,z_in)
   MB =   (1 + ( rr * (z_in)^mu * a_in^(mu-1)  - delta)*(1-tauK))*(1-tauW) ;
end

function H = H_FOC(c_in,age_in,lambda_in,e_in)
   H = max( 0 , 1 - (1-gamma)*c_in/(gamma*psi*wage*eff_un(age_in,lambda_in,e_in)) ) ;   
end

function H = H_FOC_aux(a_in,z_in,age_in,lambda_in,e_in,ap)
   H = max( 0 , gamma - (1-gamma)*(Y_a(a_in,z_in)-ap)/(psi*wage*eff_un(age_in,lambda_in,e_in)) ) ;   
end

end 




% % 
% % % Bench
% %     [A_sim_ben,C_sim_ben,H_sim_ben,Z_sim_ben,R_sim_ben,R_at_sim_ben,e_sim_ben] = Model_Simulation(1,R_bench,W_bench,C_bench,DBN_bench,agrid,zgrid);
% % 
% % % Exp
% %     [A_sim_exp,C_sim_exp,H_sim_exp,Z_sim_exp,R_sim_exp,R_at_sim_exp,e_sim_exp] = Model_Simulation(1,R_exp,W_exp,C_exp,DBN_exp,agrid,zgrid);
% % 
% % % Life Cycle 
% %     mean_A_sim_ben = mean(A_sim_ben,2) ;
% %     std_A_sim_ben  = std(A_sim_ben,0,2);
% %     mean_C_sim_ben = mean(C_sim_ben,2) ;
% %     std_C_sim_ben  = std(C_sim_ben,0,2);
% %     mean_H_sim_ben = mean(H_sim_ben,2) ;
% %     std_H_sim_ben  = std(H_sim_ben,0,2);
% %     mean_R_sim_ben = mean(R_sim_ben,2) ;
% %     std_R_sim_ben  = std(R_sim_ben,0,2);
% %     mean_R_at_sim_ben = mean(R_at_sim_ben,2) ;
% %     std_R_at_sim_ben  = std(R_at_sim_ben,0,2);
% %     
% %     
% %     for i=1:n_z
% %         ind = (Z_sim_ben==zgrid(i)) ;
% %         mean_A_sim_ben_z(:,i) = mean(A_sim_ben(:,ind),2) ;
% %         std_A_sim_ben_z(:,i)  = std(A_sim_ben(:,ind),0,2);
% %         mean_C_sim_ben_z(:,i) = mean(C_sim_ben(:,ind),2) ;
% %         std_C_sim_ben_z(:,i)  = std(C_sim_ben(:,ind),0,2);
% %         mean_H_sim_ben_z(:,i) = mean(H_sim_ben(:,ind),2) ;
% %         std_H_sim_ben_z(:,i)  = std(H_sim_ben(:,ind),0,2);
% %         mean_R_sim_ben_z(:,i) = mean(R_sim_ben(:,ind),2) ;
% %         std_R_sim_ben_z(:,i)  = std(R_sim_ben(:,ind),0,2);
% %         mean_R_at_sim_ben_z(:,i) = mean(R_at_sim_ben(:,ind),2) ;
% %         std_R_at_sim_ben_z(:,i)  = std(R_at_sim_ben(:,ind),0,2);
% %     end 
% %     
% %     mean_A_sim_exp = mean(A_sim_exp,2) ;
% %     std_A_sim_exp  = std(A_sim_exp,0,2);
% %     mean_C_sim_exp = mean(C_sim_exp,2) ;
% %     std_C_sim_exp  = std(C_sim_exp,0,2);
% %     mean_H_sim_exp = mean(H_sim_exp,2) ;
% %     std_H_sim_exp  = std(H_sim_exp,0,2);
% %     mean_R_sim_exp = mean(R_sim_exp,2) ;
% %     std_R_sim_exp  = std(R_sim_exp,0,2);
% %     mean_R_at_sim_exp = mean(R_at_sim_exp,2) ;
% %     std_R_at_sim_exp  = std(R_at_sim_exp,0,2);
% %     
% %     
% %     for i=1:n_z
% %         ind = (Z_sim_exp==zgrid(i)) ;
% %         mean_A_sim_exp_z(:,i) = mean(A_sim_exp(:,ind),2) ;
% %         std_A_sim_exp_z(:,i)  = std(A_sim_exp(:,ind),0,2);
% %         mean_C_sim_exp_z(:,i) = mean(C_sim_exp(:,ind),2) ;
% %         std_C_sim_exp_z(:,i)  = std(C_sim_exp(:,ind),0,2);
% %         mean_H_sim_exp_z(:,i) = mean(H_sim_exp(:,ind),2) ;
% %         std_H_sim_exp_z(:,i)  = std(H_sim_exp(:,ind),0,2);
% %         mean_R_sim_exp_z(:,i) = mean(R_sim_exp(:,ind),2) ;
% %         std_R_sim_exp_z(:,i)  = std(R_sim_exp(:,ind),0,2);
% %         mean_R_at_sim_exp_z(:,i) = mean(R_at_sim_exp(:,ind),2) ;
% %         std_R_at_sim_exp_z(:,i)  = std(R_at_sim_exp(:,ind),0,2);
% %     end 
% %     
% %     
    
    
    