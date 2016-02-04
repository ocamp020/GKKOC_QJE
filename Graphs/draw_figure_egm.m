clear
clf
pause

cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CALIBRATION/ELASTIC_LABOR/SSC_Adjusted/NRHO04
cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CALIBRATION/ELASTIC_LABOR/Keep_SSC_at_Bench/NRHO04

%============================= OPTIMAL TAUK EXERCISE ======================
cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/
%cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/PosDep_rho05
%cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/Equal_bequest_rho05
%cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/same_z
%cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/Equal_Bequest_Same_z

cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/more_stats

load panela_bench
load panelz_bench

x=1;
for i=100-[1, 5, 10, 25,50,75, 90,95,99, 100]
   tempa=panela_bench(panela_bench>=prctile(panela_bench,i));
   tempz=panelz_bench(panela_bench>=prctile(panela_bench,i));
   disp([100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
          sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
          sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
          sum(tempz==7)/sum(tempz>-1) ] )
  table_bench(x,:)=[100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
          sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
          sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
          sum(tempz==7)/sum(tempz>-1) ] ;   
      x=x+1;
end

load panela_exp
load panelz_exp

x=1;
for i=100-[1, 5, 10, 25,50,75, 90,95,99, 100]
   tempa=panela_exp(panela_exp>=prctile(panela_exp,i));
   tempz=panelz_exp(panela_exp>=prctile(panela_exp,i));
   disp([100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
          sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
          sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
          sum(tempz==7)/sum(tempz>-1) ] )
       
      table_exp(x,:)=[100-i, sum(tempz==1)/sum(tempz>-1)  sum(tempz==2)/sum(tempz>-1) ...
          sum(tempz==3)/sum(tempz>-1)  sum(tempz==4)/sum(tempz>-1) ...
          sum(tempz==5)/sum(tempz>-1)  sum(tempz==6)/sum(tempz>-1) ...
          sum(tempz==7)/sum(tempz>-1) ] ;   
      x=x+1;

end

table_exp
table_bench
table_exp-table_bench


cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/
cd opt_tauK/brent_tauk/more_stats
% compute consumption increase for different percentiles
load panel_cons_bench  
load panel_cons_exp 
load panel_hours_bench  
load panel_hours_exp 
panel_cons_exp_K=panel_cons_exp;
panel_hours_exp_K=panel_hours_exp;


cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/
cd opt_tauW/brent_tauW/more_stats
load panel_cons_exp 
load panel_hours_exp 
panel_cons_exp_W=panel_cons_exp;
panel_hours_exp_W=panel_hours_exp;

display('consumption')
for i=[1, 10, 25,50,75, 90,95,99,100]
  disp([i, 100*(mean(panel_cons_exp_K(panel_cons_exp_K<prctile(panel_cons_exp_K,i)))/ mean(panel_cons_bench(panel_cons_bench<prctile(panel_cons_bench,i)))-1) ... 
  100*(mean(panel_cons_exp_W(panel_cons_exp_W<prctile(panel_cons_exp_W,i)))/ mean(panel_cons_bench(panel_cons_bench<prctile(panel_cons_bench,i)))-1)])

end

display('hours')
for i=[1, 10, 25,50,75, 90,95,99,100]
  disp([i, 100*(mean(panel_hours_exp_K(panel_cons_exp_K<prctile(panel_cons_exp_K,i)))/ mean(panel_hours_bench(panel_cons_bench<prctile(panel_cons_bench,i)))-1) ... 
           100*(mean(panel_hours_exp_W(panel_cons_exp_W<prctile(panel_cons_exp_W,i)))/ mean(panel_hours_bench(panel_cons_bench<prctile(panel_cons_bench,i)))-1)])

end

display('hours')
for i=[1, 10, 25,50,75, 90,95,99,100]
  disp([i, 100*(mean(panel_hours_exp_K(panel_hours_exp_K<prctile(panel_hours_exp_K,i)))/ mean(panel_hours_bench(panel_hours_bench<prctile(panel_hours_bench,i)))-1) ... 
           100*(mean(panel_hours_exp_W(panel_hours_exp_W<prctile(panel_hours_exp_W,i)))/ mean(panel_hours_bench(panel_hours_bench<prctile(panel_hours_bench,i)))-1)])

end


cd new_stats

tauK=0.25; tauW=0.0; % benchmark
load panel_return_bench
load panelage_bench
load panela_bench
panel_at_return_bench=(1+panel_return_bench*(1-tauK))-1;

tauK=0.0; tauW=  1.682506665068614E-002; % Tax reform
tauK=  0.014 ; tauW=0; % optimal capital income tax
tauW=  0.0166; tauK=0.0; % Optimal wealth tax


load panel_return_exp
load panelage_exp
load panela_exp
panel_at_return_exp=(1+panel_return_exp*(1-tauK))*(1-tauW)-1;

display('unweighted average')
[mean(panel_return_bench) std(panel_return_bench);
    mean(panel_return_exp) std(panel_return_exp)]

[mean(panel_at_return_bench) std(panel_at_return_bench);
    mean(panel_at_return_exp) std(panel_at_return_exp)]

display('weighted average')

[mean(panel_return_bench.*panela_bench)/mean(panela_bench) ;
    mean(panel_return_exp.*panela_exp)/mean(panela_exp) ]

[mean(panel_at_return_bench.*panela_bench)/mean(panela_bench) ;
    mean(panel_at_return_exp.*panela_exp)/mean(panela_exp) ]


disp([prctile(panel_return_bench,10) prctile(panel_return_bench,25) prctile(panel_return_bench,50) prctile(panel_return_bench,75) prctile(panel_return_bench,90) prctile(panel_return_bench,95) prctile(panel_return_bench,99)])

MaxAge=81;
age_limit = [0, 5, 15, 25, 35, 45, 55, MaxAge ];
for age_group_counter=1:7
    indx1 = ((panelage_bench>age_limit(age_group_counter)).*(panelage_bench<=age_limit(age_group_counter+1)));
    indx1=indx1>0;
    temp_return=panel_return_bench(indx1);
    disp([prctile(temp_return,10) prctile(temp_return,25) prctile(temp_return,50) prctile(temp_return,75) prctile(temp_return,90) prctile(temp_return,95) prctile(temp_return,99)])
    clear indx1 temp_return
end


disp([prctile(panel_at_return_bench,10) prctile(panel_at_return_bench,25) prctile(panel_at_return_bench,50) prctile(panel_at_return_bench,75) prctile(panel_at_return_bench,90) prctile(panel_at_return_bench,95) prctile(panel_at_return_bench,99)])
for age_group_counter=1:7
    indx1 = ((panelage_bench>age_limit(age_group_counter)).*(panelage_bench<=age_limit(age_group_counter+1)));
    indx1=indx1>0;
    temp_return=panel_at_return_bench(indx1);
    disp([prctile(temp_return,10) prctile(temp_return,25) prctile(temp_return,50) prctile(temp_return,75) prctile(temp_return,90) prctile(temp_return,95) prctile(temp_return,99)])
    clear indx1 temp_return
end



disp([prctile(panel_return_exp,10) prctile(panel_return_exp,25) prctile(panel_return_exp,50) prctile(panel_return_exp,75) prctile(panel_return_exp,90) prctile(panel_return_exp,95) prctile(panel_return_exp,99)])
age_limit = [0, 5, 15, 25, 35, 45, 55, MaxAge ];
for age_group_counter=1:7
    indx1 = ((panelage_exp>age_limit(age_group_counter)).*(panelage_exp<=age_limit(age_group_counter+1)));
    indx1=indx1>0;
    temp_return=panel_return_exp(indx1);
    disp([prctile(temp_return,10) prctile(temp_return,25) prctile(temp_return,50) prctile(temp_return,75) prctile(temp_return,90) prctile(temp_return,95) prctile(temp_return,99)])
    clear indx1 temp_return
end

disp([prctile(panel_at_return_exp,10) prctile(panel_at_return_exp,25) prctile(panel_at_return_exp,50) prctile(panel_at_return_exp,75) prctile(panel_at_return_exp,90) prctile(panel_at_return_exp,95) prctile(panel_at_return_exp,99)])
for age_group_counter=1:7
    indx1 = ((panelage_exp>age_limit(age_group_counter)).*(panelage_exp<=age_limit(age_group_counter+1)));
    indx1=indx1>0;
    temp_return=panel_at_return_exp(indx1);
    disp([prctile(temp_return,10) prctile(temp_return,25) prctile(temp_return,50) prctile(temp_return,75) prctile(temp_return,90) prctile(temp_return,95) prctile(temp_return,99)])
    clear indx1 temp_return
end

pause

cd /Users/bkuruscu/Dropbox/CGGK_CODES/kubu_fortran/EGM_NEW/DBN/CRRA_CALIBRATION/FLAT_TAX/calibration_tauL0224_tauC0075/

cd opt_tauK
%cd posDep_rho05_opt_tauK
%cd EQ_B_opt_tauK
%cd same_z_opt_tauK
 
load stat_all_tau
stats_all(1,:,:)=stat_all_tau;

cd ..
cd opt_tauW
%cd posDep_rho05_opt_tauW
%cd EQ_B_opt_tauW
%cd same_z_opt_tauW

load stat_all_tau
stats_all(2,:,:)=stat_all_tau;

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), 100*(stats_all(1,:,5)/stats_all(1,1,5)-1),'r')
plot(stats_all(2,:,4), 100*(stats_all(2,:,5)/stats_all(2,1,5)-1),'b')
plot(stats_all(1,:,4), 100*(stats_all(1,:,6)/stats_all(1,1,6)-1),'r-d')
plot(stats_all(2,:,4), 100*(stats_all(2,:,6)/stats_all(2,1,6)-1),'b-d')
plot(stats_all(1,:,4), 100*(stats_all(1,:,7)/stats_all(1,1,7)-1),'ro')
plot(stats_all(2,:,4), 100*(stats_all(2,:,7)/stats_all(2,1,7)-1),'bo')
plot(stats_all(1,:,4), 100*(stats_all(1,:,8)/stats_all(1,1,8)-1),'rd')
plot(stats_all(2,:,4), 100*(stats_all(2,:,8)/stats_all(2,1,8)-1),'bd')
grid on
legend(' KBAR, \tau_K', 'KBAR, \tau_W', 'QBAR, \tau_K', 'QBAR, \tau_W', ...
    'NBAR, \tau_K', 'NBAR, \tau_W', 'YBAR, \tau_K', 'YBAR, \tau_W','Location','SouthWest')
xlabel('GBAR_K')

hgsave('1fig_KBAR_QBAR_by_CAP_TAX_REV.fig')
print -dpdf 1fig_KBAR_QBAR_by_CAP_TAX_REV.pdf
print -dps 1fig_KBAR_QBAR_by_CAP_TAX_REV.eps

pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), (stats_all(1,:,10)-stats_all(1,1,10)),'r')
plot(stats_all(2,:,4), (stats_all(2,:,10)-stats_all(2,1,10)),'b')
grid on
legend('\tau_K', '\tau_W')
xlabel('GBAR_K')
title('CE NEWBORN')
[maxvalK indxK] = max((stats_all(1,:,10)-stats_all(1,1,10)));
text(stats_all(1,indxK,4)+0.0, maxvalK -0.15,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,10)-stats_all(2,1,10)));
text(stats_all(2,indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])


hgsave('1fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
print -dpdf 1fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
print -dps 1fig_CE_NEWBORN_by_CAP_TAX_REV.eps

pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), (stats_all(1,:,11)-stats_all(1,1,11)),'r')
plot(stats_all(2,:,4), (stats_all(2,:,11)-stats_all(2,1,11)),'b')
grid on
xlabel('GBAR_K')
title('Average Utility')

[maxvalK indxK] = max((stats_all(1,:,11)-stats_all(1,1,11)));
text(stats_all(1,indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,11)-stats_all(2,1,11)));
text(stats_all(2,indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])


hgsave('1fig_Average_Utility_by_CAP_TAX_REV.fig')
print -dpdf 1fig_Average_Utility_by_CAP_TAX_REV.pdf
print -dps 1fig_Average_Utility_by_CAP_TAX_REV.eps


clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), (stats_all(1,:,12)-stats_all(1,1,12)),'r')
plot(stats_all(2,:,4), (stats_all(2,:,12)-stats_all(2,1,12)),'b')
grid on
xlabel('GBAR_K')
title('Average Utility')

[maxvalK indxK] = max((stats_all(1,:,12)-stats_all(1,1,12)));
text(stats_all(1,indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,12)-stats_all(2,1,12)));
text(stats_all(2,indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), (stats_all(1,:,12) ),'r')
plot(stats_all(2,:,4), (stats_all(2,:,12) ),'b')
grid on
xlabel('GBAR_K')
title('CE Newborn computed at 2nd dbn')

[maxvalK indxK] = max((stats_all(1,:,12) ));
text(stats_all(1,indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,12) ));
text(stats_all(2,indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])


gamma=0.4494;
sigma=4;
stats_all(1,:,13) =100*( (stats_all(1,:,11)./stats_all(1,26,11)).^(1/(gamma*(1-sigma)))-1);
stats_all(2,:,13) =100*( (stats_all(2,:,11)./stats_all(1,26,11)).^(1/(gamma*(1-sigma)))-1);

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), (stats_all(1,:,13)-stats_all(1,1,13)),'r')
plot(stats_all(2,:,4), (stats_all(2,:,13)-stats_all(2,1,13)),'b')
grid on
xlabel('GBAR_K')
title('CE Newborn computed using average utilities')

[maxvalK indxK] = max((stats_all(1,:,13)-stats_all(1,1,13)));
text(stats_all(1,indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,13)-stats_all(2,1,13)));
text(stats_all(2,indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])

hgsave('2fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
print -dpdf 2fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
print -dps 2fig_CE_NEWBORN_by_CAP_TAX_REV.eps

pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4), (stats_all(1,:,13) ),'r')
plot(stats_all(2,:,4), (stats_all(2,:,13) ),'b')
grid on
xlabel('GBAR_K')
title('CE Newborn computed using average utilities')

[maxvalK indxK] = max((stats_all(1,:,13) ));
text(stats_all(1,indxK,4)+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,13) ));
text(stats_all(2,indxW,4)+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])


hgsave('3fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
print -dpdf 3fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
print -dps 3fig_CE_NEWBORN_by_CAP_TAX_REV.eps

pause

%--------------------------------------------------------------------------
%

GBAR_bench=  0.281709777634299;
SSC_Payments=  0.201224662623518;
TOT_budget= GBAR_bench+SSC_Payments;

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4)/TOT_budget, (stats_all(1,:,13)-stats_all(1,1,13)),'r')
plot(stats_all(2,:,4)/TOT_budget, (stats_all(2,:,13)-stats_all(2,1,13)),'b')
grid on
xlabel('GBAR_K / (GBAR+SSC)')
title('CE Newborn computed using average utilities')

[maxvalK indxK] = max((stats_all(1,:,13)-stats_all(1,1,13)));
text(stats_all(1,indxK,4)/TOT_budget+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,13)-stats_all(2,1,13)));
text(stats_all(2,indxW,4)/TOT_budget+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])

hgsave('2fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
print -dpdf 2fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
print -dps 2fig_CE_NEWBORN_by_CAP_TAX_REV.eps

pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
hold on
plot(stats_all(1,:,4)/TOT_budget, (stats_all(1,:,13) ),'r')
plot(stats_all(2,:,4)/TOT_budget, (stats_all(2,:,13) ),'b')
grid on
xlabel('GBAR_K / (GBAR+SSC)')
title('CE Newborn computed using average utilities')

[maxvalK indxK] = max((stats_all(1,:,13) ));
text(stats_all(1,indxK,4)/TOT_budget+0.01, maxvalK -0.10,['Opt. \tau_K=' ,num2str(stats_all(1,indxK,1))])

[maxvalW indxW] = max((stats_all(2,:,13) ));
text(stats_all(2,indxW,4)/TOT_budget+0.0, maxvalW +0.05,['Opt. \tau_W=' ,num2str(stats_all(2,indxW,2))])


hgsave('3fig_CE_NEWBORN_by_CAP_TAX_REV.fig')
print -dpdf 3fig_CE_NEWBORN_by_CAP_TAX_REV.pdf
print -dps 3fig_CE_NEWBORN_by_CAP_TAX_REV.eps

pause





%--------------------------------------------------------------------------
load mean_wealth_by_agegroup_z_bench
load mean_wealth_by_agegroup_z_exp

load agrid
load aprime_age1_bench
load aprime_age1_exp
load aprime_age16_bench
load aprime_age16_exp
load aprime_age31_bench
load aprime_age31_exp
load MeanAsset_by_z_med_E_Lambda_age1
load CE_by_Asset_z_med_E_Lambda_age1
load MeanAsset_by_z_med_E_Lambda_age16
load CE_by_Asset_z_med_E_Lambda_age16
load MeanAsset_by_z_med_E_Lambda_age31
load CE_by_Asset_z_med_E_Lambda_age31


alow=8;
ahigh=20;
clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_diff_savings_rate_age1.eps


pause


alow=9;
ahigh=20;
clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_diff_savings_rate_age16.eps


pause


alow=9;
ahigh=20;
clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
grid on


hgsave('1fig_diff_savings_rate_age31.fig')
print -dpdf 1fig_diff_savings_rate_age31.pdf
print -dps 1fig_diff_savings_rate_age31.eps


pause


%-------------------


naa=40;
nzz=1;

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age1_z1.eps
pause

naa=40;
nzz=4;

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age1_z4.eps
pause

naa=40;
nzz=6;

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age1_z7.eps
pause


clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age16_z1.eps
pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age16_z4.eps
pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age16_z7.eps
pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age31_z1.eps
pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age31_z4.eps
pause

clf
axes1 = axes(...
  'FontName','Times New Roman',...
  'FontSize',18);
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
print -dps 1fig_Welfare_Gain_by_Wealth_for_Age31_z7.eps
pause

%------------------------


