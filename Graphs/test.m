for i=1:Max_Age
        ind = (Sim_age_ben==i) ;
        n   = sum(ind)         ;
%         sim_mean_H_ben(i)  = sum(Sim_H_ben(ind))/n  ;
%         sim_std_H_ben(i)  = ( sum( (Sim_H_ben(ind) - sim_mean_H_ben(i)).^2 )/(n-1) )^0.5 ;
%         
%         
%         ind = (Sim_age_exp==i) ;
%         n   = sum(ind)         ;
%         sim_mean_H_exp(i)  = sum(Sim_H_exp(ind))/n  ;
%         sim_std_H_exp(i)  = ( sum( (Sim_H_exp(ind) - sim_mean_H_exp(i)).^2 )/(n-1) )^0.5 ;
%         
%         for j=1:n_z
%             ind = ((Sim_age_ben==i).*(Sim_Z_ben==j))==1 ;
%             n   = sum(ind)                         ;
%             sim_mean_H_ben_AZ(i,j)  = sum(Sim_H_ben(ind))/n  ;
%             sim_std_H_ben_AZ(i,j)  = ( sum( (Sim_H_ben(ind) - sim_mean_H_ben_AZ(i,j)).^2 )/(n-1) )^0.5 ;
%             
%             ind = ((Sim_age_exp==i).*(Sim_Z_exp==j))==1 ;
%             n   = sum(ind)                         ;
%             sim_mean_H_exp_AZ(i,j)  = sum(Sim_H_exp(ind))/n  ;
%             im_std_H_exp_AZ(i,j)  = ( sum( (Sim_H_exp(ind) - sim_mean_H_exp_AZ(i,j)).^2 )/(n-1) )^0.5 ;
%         end
%         
%         
%         eff_un_age(i) = sum(sum((squeeze(eff_un(i,:,:)).*reshape(sum(sum(DBN_bench(i,:,:,:,:))),[5,5]))))/sum(sum(sum(sum(DBN_bench(i,:,:,:,:)))));
        
        sim_var_lnC_ben(i) = var(log(Sim_C_ben(ind))) ;
        %sim_var_lnC_ben(i) = ( sum( ((Sim_C_ben(ind)) - (sim_mean_C_ben(i))).^2 )/(n-1) )^0.5 ;
        sim_var_lnC_ben(i) = mean((Sim_C_ben(ind))) ;
end

figure; plot(20:100,sim_var_lnC_ben); title('Variance log(C)'); xlim([19+1,19+Max_Age]);

% plot(eff_un_age)
% 
% mean(Sim_Yh_ben(Sim_age_ben==(25-19)))
% mean(Sim_Yh_ben(Sim_age_ben==(55-19)))
% mean(Sim_Yh_ben(Sim_age_ben==(64-19)))
% mean(Sim_Yh_ben(Sim_age_ben==(65-19)))