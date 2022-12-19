%% Use it or lose it: Efficiency and redistributional effects of wealth taxation
% Fatih Guvenen, Gueorgui Kambourov, Burhan Kuruscu, Sergio Ocampo, & Daphne Chen

% Graphs and Tables 


%% Set Up 
   % Clear and close
    close all; clear all; clc
 
    % Line style
    set(groot,'defaultAxesLineStyleOrder','-|--|-.|:')

    % Parameters for size
    n_a = 201 ; 
    n_z = 9   ;
    n_l = 5   ;
    n_e = 5   ;
    n_x = 3   ; 
    
    Max_Age = 81 ;
    Ret_Age = 45 ;
    
    % Parameters for scaling to dollars
    EBAR_data = 8.8891*10^12/(122.46*10^6) ; 
        % 2013 total compensation of employees' devided by number of HH's in 2013
    
    % Load wealth data (Vermuelen Data)
        load ../Graphs/Vermeulen_Data/USA_IM1.dat
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
        
        
    % Load agrid 
        load ../Model_2.1/agrid
        agrid_0 = agrid ;
        
    % Minimum value for pareto tail
        w_min = 1000000;
        
    % X Vec
        x_vec = [5 1 -inf] ;
        
    % Select min value for productivity
        prod_min = 2 ; 
      
    % Folder 
        mkdir('Model_2.1')
        mkdir('Model_2.1/png')
        mkdir('Model_2.1/fig')
        Tables_file = 'Model_2.1/Tables_Model_2.1.xls' ;
    
    % Values for pareto plot 
        a_plot        = 0.65                ;
        Data_Color    = [0 0.447 0.741]     ;
        Data_Color_a  = a_plot*Data_Color+(1-a_plot)*[1 1 1] ; 
        Model_Color   = [0.85 0.325 0.098]  ;
        Model_Color_a = a_plot*Model_Color+(1-a_plot)*[1 1 1] ; 



%% Figure 1.A (Pareto Tail, Benchmark Model) 


    % Simul folder 
    Simul_Folder = "../Model_2.1/Simul/" ;
    
    % Bench Files
        load ../Model_2.1/Simul/panel_Pareto
        wealth_bench = sort(panel_Pareto) ; clear panel_Pareto
        
    % Pareto Tail
        w_vec = [1000000 prctile(wealth_bench,[98 99])];
        for i=1:numel(w_vec)
        w_min = w_vec(i);%1000000;
        % Bench
            [y_bench,x_bench] = ecdf(wealth_bench(wealth_bench>=w_min)) ;
            x_bench = -log(x_bench(2:end)) ; y_bench = log(1-y_bench(1:end-1)) ;
        % Vermuelen
            x_V_s     = -log(panel_V((panel_V>w_min)&(type_V==1))) ;
            y_V_s     = log(cumsum(weights_V((panel_V>w_min)&(type_V==1)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V_f     = -log(panel_V((panel_V>w_min)&(type_V==0))) ;
            y_V_f     = log(cumsum(weights_V((panel_V>w_min)&(type_V==0)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V       = -log(panel_V(panel_V>w_min)) ;
            y_V       = log(cumsum(weights_V((panel_V>w_min)),'reverse')/sum(weights_V(panel_V>w_min))) ;
        % Maximum Likelihood
        a_ml  = sum(weights_V(panel_V>w_min)) / sum(weights_V(panel_V>w_min).*log(panel_V(panel_V>w_min)/w_min)) ;
        % Regression
        mdl_b     = fitlm(x_bench',y_bench')         ;
        C_b       = mdl_b.Coefficients.Estimate(1)   ;
        a_b   = mdl_b.Coefficients.Estimate(2)   ; 
        mdl_V_s   = fitlm(x_V_s',y_V_s')             ;      
        C_V_s     = mdl_V_s.Coefficients.Estimate(1) ;
        a_V_s = mdl_V_s.Coefficients.Estimate(2) ;
        mdl_V     = fitlm(x_V',y_V')                 ;
        C_V       = mdl_V.Coefficients.Estimate(1)   ;
        a_V   = mdl_V.Coefficients.Estimate(2)   ;
        xx        = linspace(min(-x_V),max(-x_V),3)  ;
        x_bench_2=x_bench(1:5:end); y_bench_2=y_bench(1:5:end);
        % Figure
        fig_title = ['Pareto Tail Above $',num2str(w_min)] ;
        figure; hold on; 
        plot(-x_V,y_V,'DisplayName','US Data','Marker','o','LineWidth',1,'LineStyle','none','Color',Data_Color_a,'MarkerSize',4);
        plot(xx,C_V-a_V*xx,'DisplayName','Regression Line','LineWidth',2.5,'LineStyle','-.','Color',Data_Color); 
        plot(-x_bench_2,y_bench_2,'DisplayName','Benchmark Model','Marker','diamond','LineWidth',1,'LineStyle','none','Color',Model_Color_a,'MarkerSize',4); 
        plot(xx,C_b-a_b*xx,'DisplayName','Regression Line','LineWidth',2.5,'LineStyle','--','Color',Model_Color); 
        hold off;
        xlabel('Wealth (log scale)','FontName','Times New Roman');
        ylabel('Log Counter-CDF','HorizontalAlignment','center','FontName','Times New Roman');
        grid(gca,'on');
        set(gca,'FontName','Times New Roman','FontSize',16);
        legend1 = legend(gca,'show');
        set(legend1,'Interpreter','latex','FontSize',13);
        xlim([log(w_min),max(-x_V_f)]);
        ylim([-16 0]); 
        ww = log([1e6 10e6 100e6 1e9 10e9 50e9]);
        set(gca,'XTick',ww,'YTick',[-16:2:0]); set(gca,'XTickLabel',{'$1M','$10M','$100M','$1B','$10B','$50B'}); 
        file_name_eps = ['Model_2.1/Pareto_Tail_$',num2str(w_min),'_Pooled_Cross_Sections.eps'] ;
        file_name_png = ['Model_2.1/png/Pareto_Tail_$',num2str(w_min),'_Pooled_Cross_Sections.png'] ;
        file_name_fig = ['Model_2.1/fig/Pareto_Tail_$',num2str(w_min),'_Pooled_Cross_Sections.fig'] ;
        print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);
        end


%% Figure 1.B (Pareto Tail, Alternative Models)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Low-Inequality Calibraton 

        load ../Model_2.1_Match_Return_Lambda/Simul/panel_Pareto
        wealth_bench = sort(panel_Pareto) ; clear panel_Pareto
        

        w_min = 1000000;
        % Bench
            [y_bench,x_bench] = ecdf(wealth_bench(wealth_bench>=w_min)) ;
            x_bench = -log(x_bench(2:end)) ; y_bench = log(1-y_bench(1:end-1)) ;
        % Vermuelen
            x_V_s     = -log(panel_V((panel_V>w_min)&(type_V==1))) ;
            y_V_s     = log(cumsum(weights_V((panel_V>w_min)&(type_V==1)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V_f     = -log(panel_V((panel_V>w_min)&(type_V==0))) ;
            y_V_f     = log(cumsum(weights_V((panel_V>w_min)&(type_V==0)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V       = -log(panel_V(panel_V>w_min)) ;
            y_V       = log(cumsum(weights_V((panel_V>w_min)),'reverse')/sum(weights_V(panel_V>w_min))) ;
        % Maximum Likelihood
        a_ml  = sum(weights_V(panel_V>w_min)) / sum(weights_V(panel_V>w_min).*log(panel_V(panel_V>w_min)/w_min)) ;
        % Regression
        mdl_b     = fitlm(x_bench',y_bench')         ;
        C_b       = mdl_b.Coefficients.Estimate(1)   ;
        a_b   = mdl_b.Coefficients.Estimate(2)   ; 
        mdl_V_s   = fitlm(x_V_s',y_V_s')             ;      
        C_V_s     = mdl_V_s.Coefficients.Estimate(1) ;
        a_V_s = mdl_V_s.Coefficients.Estimate(2) ;
        mdl_V     = fitlm(x_V',y_V')                 ;
        C_V       = mdl_V.Coefficients.Estimate(1)   ;
        a_V   = mdl_V.Coefficients.Estimate(2)   ;
        xx        = linspace(min(-x_V),max(-x_V),3)  ;
        % Select x_bench and y_bench
        x_bench_2=x_bench(1:5:end); y_bench_2=y_bench(1:5:end);
        % Figure
        fig_title = ['Pareto Tail Above $',num2str(w_min)] ;
        figure; hold on; 
        plot(-x_V,y_V,'DisplayName','US Data','Marker','o','LineWidth',1,'LineStyle','none','Color',Data_Color_a,'MarkerSize',4);
        plot(xx,C_V-a_V*xx,'DisplayName','Regression Line','LineWidth',2.5,'LineStyle','-.','Color',Data_Color); 
        plot(-x_bench_2,y_bench_2,'DisplayName','Benchmark Model','Marker','diamond','LineWidth',1,'LineStyle','none','Color',Model_Color_a,'MarkerSize',4); 
        plot(xx,C_b-a_b*xx,'DisplayName','Regression Line','LineWidth',2.5,'LineStyle','--','Color',Model_Color); 
        hold off;
        xlabel('Wealth (log scale)','FontName','Times New Roman');
        ylabel('Log Counter-CDF','HorizontalAlignment','center','FontName','Times New Roman');
        grid(gca,'on');
        set(gca,'FontName','Times New Roman','FontSize',16);
        legend1 = legend(gca,'show');
        set(legend1,'Interpreter','latex','FontSize',13);
        xlim([log(w_min),max(-x_V_f)]);
        ylim([-16 0]); 
        ww = log([1e6 10e6 100e6 1e9 10e9 50e9]);
        set(gca,'XTick',ww,'YTick',[-16:2:0]); set(gca,'XTickLabel',{'$1M','$10M','$100M','$1B','$10B','$50B'}); 
        file_name_eps = ['Model_2.1/Pareto_Tail_$',num2str(w_min),'_Low_Ineq.eps'] ;
        file_name_png = ['Model_2.1/png/Pareto_Tail_$',num2str(w_min),'_Low_Ineq.png'] ;
        file_name_fig = ['Model_2.1/fig/Pareto_Tail_$',num2str(w_min),'_Low_Ineq.fig'] ;
        print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);  
    


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Awesome State Income Shocks

        load ../Model_2.1_Super_Aiyagari/Simul/panela_bench
        load ../Model_2.1_Super_Aiyagari/Bench_Files/EBAR
        
        EBAR_bench   = EBAR*0.727853584919652 ;
        wealth_bench = EBAR_data/EBAR_bench * sort(panela_bench) ;
        
        % Pareto Tail
        w_min = 1000000;
        % Bench
            [y_bench,x_bench] = ecdf(wealth_bench(wealth_bench>w_min)) ;
            x_bench = -log(x_bench(2:end)) ; y_bench = log(1-y_bench(1:end-1)) ;
        % Vermuelen
            x_V_s     = -log(panel_V((panel_V>w_min)&(type_V==1))) ;
            y_V_s     = log(cumsum(weights_V((panel_V>w_min)&(type_V==1)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V_f     = -log(panel_V((panel_V>w_min)&(type_V==0))) ;
            y_V_f     = log(cumsum(weights_V((panel_V>w_min)&(type_V==0)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V       = -log(panel_V(panel_V>w_min)) ;
            y_V       = log(cumsum(weights_V((panel_V>w_min)),'reverse')/sum(weights_V(panel_V>w_min))) ;
        % Maximum Likelihood
        a_ml  = sum(weights_V(panel_V>w_min)) / sum(weights_V(panel_V>w_min).*log(panel_V(panel_V>w_min)/w_min)) ;
        % Regression
        mdl_b     = fitlm(x_bench',y_bench')         ;
        C_b       = mdl_b.Coefficients.Estimate(1)   ;
        a_b   = mdl_b.Coefficients.Estimate(2)   ; 
        mdl_V_s   = fitlm(x_V_s',y_V_s')             ;      
        C_V_s     = mdl_V_s.Coefficients.Estimate(1) ;
        a_V_s = mdl_V_s.Coefficients.Estimate(2) ;
        mdl_V     = fitlm(x_V',y_V')                 ;
        C_V       = mdl_V.Coefficients.Estimate(1)   ;
        a_V   = mdl_V.Coefficients.Estimate(2)   ;
        xx        = linspace(min(-x_V),max(-x_V),3)  ;
        % Select x_bench and y_bench
        x_bench_2 = [x_bench(1:200:end-600) ; x_bench(end-599:10:end-50) ; x_bench(end-50:end)];
        y_bench_2 = [y_bench(1:200:end-600) ; y_bench(end-599:10:end-50) ; y_bench(end-50:end)];
        % Figure
        fig_title = ['Pareto Tail Above $',num2str(w_min)] ;
        figure; hold on; 
        plot(-x_V,y_V,'DisplayName','US Data','Marker','o','LineWidth',1,'LineStyle','none','Color',Data_Color_a,'MarkerSize',4);
        plot(xx,C_V-a_V*xx,'DisplayName','Regression Line','LineWidth',2.5,'LineStyle','-.','Color',Data_Color); 
        plot(-x_bench_2,y_bench_2,'DisplayName','Benchmark Model','Marker','diamond','LineWidth',1,'LineStyle','none','Color',Model_Color_a,'MarkerSize',4); 
        plot(xx,C_b-a_b*xx,'DisplayName','Regression Line','LineWidth',2.5,'LineStyle','--','Color',Model_Color); 
        hold off;
        xlabel('Wealth (log scale)','FontName','Times New Roman');
        ylabel('Log Counter-CDF','HorizontalAlignment','center','FontName','Times New Roman');
        grid(gca,'on');
        set(gca,'FontName','Times New Roman','FontSize',16);
        legend1 = legend(gca,'show');
        set(legend1,'Interpreter','latex','FontSize',13);
        xlim([log(w_min),max(-x_V_f)]);
        ylim([-16 0]); 
        ww = log([1e5 1e6 10e6 100e6 1e9 10e9 50e9]);
        set(gca,'XTick',ww,'YTick',[-16:2:0]); set(gca,'XTickLabel',{'$100K','$1M','$10M','$100M','$1B','$10B','$50B'}); 
        file_name_eps = ['Model_2.1/Pareto_Tail_$',num2str(w_min),'_Simulation_Super_Aiyagari.eps'] ;
        file_name_png = ['Model_2.1/png/Pareto_Tail_$',num2str(w_min),'_Simulation_Super_Aiyagari.png'] ;
        file_name_fig = ['Model_2.1/fig/Pareto_Tail_$',num2str(w_min),'_Simulation_Super_Aiyagari.fig'] ;
        print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);




%% Figure 2 (Welfare Gains from Optimal Taxation)

% Read Tables 
    OTW_Table = readtable('/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Model_2.1/Opt_Tax_W/Wide_Grid/Stats_by_tau_w.txt') ;
    OTK_Table = readtable('/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Model_2.1/Opt_Tax_K/Wide_Grid/Stats_by_tau_k.txt') ;
    N         = numel(OTW_Table.Av_Util_NB) ;

% Capital Tax Revenue Over Total Tax Revenue
    K_Tax_Bench = 0.99*0.236992059442476 ;
    [~,OTK_grid]= max(OTK_Table.Av_Util_NB);
    K_Tax_OTK   = OTK_Table.GBAR_K_Tax_Rev_bench(OTK_grid); % -0.131283009707177 ;
    [~,OTW_grid]= max(OTW_Table.Av_Util_NB);
    K_Tax_OTW   = OTW_Table.GBAR_K_Tax_Rev_bench(OTW_grid); % 0.355360711880188 ;

% Figure: Welfare Gain
    graph_ind = 1:4:N ;
    x_lim = [-0.4  0.5];
    y_lim = [-6.0 10.0];
    figure; hold on; 
    % Graph welfare gain
    plot(OTK_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTK_Table.Av_Util_NB(graph_ind),'r','linewidth',2.5)
    plot(OTW_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTW_Table.Av_Util_NB(graph_ind),'linewidth',2.5,'color',[0.3,0.3,0.3])
    % Lines 
    plot(x_lim,[0,0],'k--','linewidth',2.0)
    plot([K_Tax_Bench,K_Tax_Bench],[y_lim(1),0],'k--','linewidth',2.0)
    % Arrows 
    % plot([K_Tax_OTK,K_Tax_OTK],y_lim,'k--','linewidth',2.0)
    % plot([K_Tax_OTW,K_Tax_OTW],y_lim,'k--','linewidth',2.0)
    annotation('textarrow',[0.30,0.33],[0.60,0.67],'String','Optimal $\tau_k$=-13.6\%','Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    annotation('textarrow',[0.76,0.78],[0.78,0.85],'String','Optimal $\tau_a$=3.03\%' ,'Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    annotation('textarrow',[0.62,0.67],[0.33,0.40],'String','U.S. Benchmark,  $\tau_k$=25\%','Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    % Labels 
    xlim(x_lim); ylim(y_lim); grid(gca,'on'); set(gca,'fontsize',12);
    xlabel('Tax Revenue from K $/$ Total Tax Revenue','Interpreter','latex','FontName','Times New Roman','FontSize',16);
    ylabel('$\overline{CE}_2$ Welfare Change from U.S. Benchmark','Interpreter','latex','FontName','Times New Roman','FontSize',16);
    legend({'Capital income tax economy','Wealth tax economy'},'Interpreter','latex','FontSize',13,'location','northwest','FontName','Times New Roman');  legend('boxoff');  
    % Save Figure 
    file_name_eps = 'Model_2.1/Opt_Tax_Welfare.eps' ;
    file_name_png = 'Model_2.1/png/Opt_Tax_Welfare.png' ;
    file_name_fig = 'Model_2.1/fig/Opt_Tax_Welfare.fig' ;
    print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);


%% Figure 3 (How K and Q Vary with Revenue Raised from Taxing Capital)

% Read Tables 
    OTW_Table = readtable('/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Model_2.1/Opt_Tax_W/Wide_Grid/Stats_by_tau_w.txt') ;
    OTK_Table = readtable('/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Model_2.1/Opt_Tax_K/Wide_Grid/Stats_by_tau_k.txt') ;
    N         = numel(OTW_Table.Av_Util_NB) ;

% Figure: K and Q
    % Get Change in K and Q 
        [~,ind_0] = min(abs(OTK_Table.tauK)) ; 
        OTK_K = 100*(OTK_Table.KBAR/OTK_Table.KBAR(ind_0)-1) ; 
        OTK_Q = 100*(OTK_Table.QBAR/OTK_Table.QBAR(ind_0)-1) ; 
        [~,ind_0] = min(abs(OTW_Table.tauW_at)) ; 
        OTW_K = 100*(OTW_Table.KBAR/OTW_Table.KBAR(ind_0)-1) ; 
        OTW_Q = 100*(OTW_Table.QBAR/OTW_Table.QBAR(ind_0)-1) ; 
    % Make Figure 
    graph_ind = 1:2:N ;
    x_lim = [-0.4  0.5];
    y_lim = [-40.0 40.0];
    figure; hold on; 
    % Graph Capital and Q 
    plot(OTK_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTK_K(graph_ind),'r'  ,'linewidth',2.5)
    plot(OTK_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTK_Q(graph_ind),'r--','linewidth',2.5)
    plot(OTW_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTW_K(graph_ind),'-'  ,'linewidth',2.5,'color',[0.3,0.3,0.3])
    plot(OTW_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTW_Q(graph_ind),'--' ,'linewidth',2.5,'color',[0.3,0.3,0.3])
    % Lines 
    plot(x_lim,[0,0],'--','linewidth',1.0,'color',[0.7,0.7,0.7])
    plot([0,0],y_lim,'--','linewidth',1.0,'color',[0.7,0.7,0.7])
    % Arrows 
    annotation('textarrow',[0.29,0.25],[0.77,0.77],'String','OCIT Economy','Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    annotation('textarrow',[0.69,0.69],[0.54,0.43],'String','OWT Economy' ,'Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    text(0.4,-15,'Q','FontSize',17,'FontWeight','bold','Color',[0.3,0.3,0.3])
    text(0.4,-25,'K','FontSize',17,'FontWeight','bold','Color',[0.3,0.3,0.3])
    text(0.30,-37,'Q','FontSize',17,'FontWeight','bold','Color',[1.0,0.0,0.0])
    text(0.35,-29,'K','FontSize',17,'FontWeight','bold','Color',[1.0,0.0,0.0])
    % Labels 
    xlim(x_lim); ylim(y_lim); grid(gca,'on'); set(gca,'fontsize',12);
    xlabel('Tax Revenue from K $/$ Total Tax Revenue','Interpreter','latex','FontName','Times New Roman','FontSize',16);
    ylabel('\% Deviation from Case with $\tau_k$=0 or $\tau_a$=0','Interpreter','latex','FontName','Times New Roman','FontSize',16);
    % legend({'Capital income tax economy','Wealth tax economy'},'Interpreter','latex','FontSize',13,'location','northwest','FontName','Times New Roman');  legend('boxoff');  
    % Save Figure 
    file_name_eps = 'Model_2.1/Opt_Tax_KQ.eps' ;
    file_name_png = 'Model_2.1/png/Opt_Tax_KQ.png' ;
    file_name_fig = 'Model_2.1/fig/Opt_Tax_KQ.fig' ;
    print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);



%% Figure C.1 (Pareto Tail, mu=0.8)

   load ../Model_2.1_mu80/Simul/panel_Pareto
    wealth_bench = sort(panel_Pareto) ; clear panel_Pareto
        
    % Pareto Tail
        w_min = 1000000;
        % Bench
            [y_bench,x_bench] = ecdf(wealth_bench(wealth_bench>w_min)) ;
            x_bench = -log(x_bench(2:end)) ; y_bench = log(1-y_bench(1:end-1)) ;
        % Vermuelen
            x_V_s     = -log(panel_V((panel_V>w_min)&(type_V==1))) ;
            y_V_s     = log(cumsum(weights_V((panel_V>w_min)&(type_V==1)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V_f     = -log(panel_V((panel_V>w_min)&(type_V==0))) ;
            y_V_f     = log(cumsum(weights_V((panel_V>w_min)&(type_V==0)),'reverse')/sum(weights_V(panel_V>w_min))) ;
            x_V       = -log(panel_V(panel_V>w_min)) ;
            y_V       = log(cumsum(weights_V((panel_V>w_min)),'reverse')/sum(weights_V(panel_V>w_min))) ;
        % Maximum Likelihood
        a_ml  = sum(weights_V(panel_V>w_min)) / sum(weights_V(panel_V>w_min).*log(panel_V(panel_V>w_min)/w_min)) ;
        % Regression
        mdl_b     = fitlm(x_bench',y_bench')         ;
        C_b       = mdl_b.Coefficients.Estimate(1)   ;
        a_b   = mdl_b.Coefficients.Estimate(2)   ; 
        mdl_V_s   = fitlm(x_V_s',y_V_s')             ;      
        C_V_s     = mdl_V_s.Coefficients.Estimate(1) ;
        a_V_s = mdl_V_s.Coefficients.Estimate(2) ;
        mdl_V     = fitlm(x_V',y_V')                 ;
        C_V       = mdl_V.Coefficients.Estimate(1)   ;
        a_V   = mdl_V.Coefficients.Estimate(2)   ;
        xx        = linspace(min(-x_V),max(-x_V),3)  ;
        % Select x_bench and y_bench
        x_bench_2=x_bench(1:5:end); y_bench_2=y_bench(1:5:end);
        % Figure
        fig_title = ['Pareto Tail Above $',num2str(w_min)] ;
        figure; hold on; 
        plot(-x_V,y_V,'Marker','o','LineWidth',1,'LineStyle','none','Color',Data_Color_a,'MarkerSize',4);
        plot(xx,C_V-a_V*xx,'LineWidth',2.5,'LineStyle','-.','Color',Data_Color); 
        plot(-x_bench_2,y_bench_2,'Marker','diamond','LineWidth',1,'LineStyle','none','Color',Model_Color_a,'MarkerSize',4); 
        plot(xx,C_b-a_b*xx,'LineWidth',2.5,'LineStyle','--','Color',Model_Color); 
        hold off;
        xlabel('Wealth (log scale)','FontName','Times New Roman');
        ylabel('Log Counter-CDF','HorizontalAlignment','center','FontName','Times New Roman');
        grid(gca,'on');
        set(gca,'FontName','Times New Roman','FontSize',16);
        legend('US Data','Regression Line','$\mu$=0.8 Model','Regression Line','Interpreter','latex','FontSize',13)
        xlim([log(w_min),max(-x_V_f)]);
        ylim([-16 0]); 
        if w_min~=1000000
            ww = linspace(log(w_min),log(panel_V(end)),5) ;
        else
            ww = log([1e6 10e6 100e6 1e9 10e9 50e9]);
        end 
        set(gca,'XTick',ww,'YTick',[-16:2:0]); set(gca,'XTickLabel',{'$1M','$10M','$100M','$1B','$10B','$50B'}); 
        file_name_eps = ['Model_2.1/Pareto_Tail_$',num2str(w_min),'_mu80.eps'] ;
        file_name_png = ['Model_2.1/png/Pareto_Tail_$',num2str(w_min),'_mu80.png'] ;
        file_name_fig = ['Model_2.1/fig/Pareto_Tail_$',num2str(w_min),'_mu80.fig'] ;
        print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);


%% Figure C.2 (Fraction of Entrepreneurs over the Life Cycle)

    % Load Benchmark results
    bench_profile = readmatrix('../Model_2.1/Entrepreneur_by_age.txt','NumHeaderLines',3) ;
    bench_profile = bench_profile(:,4) ; % Select entrepreneurs as business_inc/total_inc>50%

    % Load low inequality (lambda) results 
    low_ineq_profile = readmatrix('../Model_2.1_Match_Return_Lambda/Entrepreneur_by_age.txt','NumHeaderLines',3) ;
    low_ineq_profile = low_ineq_profile(:,4) ; % Select entrepreneurs as business_inc/total_inc>50%

    % Load low initial productivity results 
    low_z0_profile = readmatrix('../Model_2.1_X_Transition_Low_NB/Entrepreneur_by_age.txt','NumHeaderLines',3) ;
    low_z0_profile = low_z0_profile(:,4) ; % Select entrepreneurs as business_inc/total_inc>50%

    % Select valus every 5 years
    ind = 1:5:numel(bench_profile) ;

    % Separate graphs 

        % Benchmark model 
         figure; grid on; hold on;
            scatter(20:5:100,bench_profile(ind),'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            ylabel('Share of Entrepreneurs (%)','FontName','Times New Roman','FontSize',14);
            xlabel('Age','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([19,81]); xticks([20:5:100]); 
            ylim([0,17]); yticks([0:2:16]);
            % Save figure 
                file_name_eps = 'Model_2.1/Entrepreneur_Age_Profile_bench.eps'     ;
                file_name_png = 'Model_2.1/png/Entrepreneur_Age_Profile_bench.png' ;
                file_name_fig = 'Model_2.1/fig/Entrepreneur_Age_Profile_bench.fig'  ;
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)


        % Low inequality model 
         figure; grid on; hold on;
            scatter(20:5:100,low_ineq_profile(ind),'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            ylabel('Share of Entrepreneurs (%)','FontName','Times New Roman','FontSize',14);
            xlabel('Age','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([19,81]); xticks([20:5:100]); 
            ylim([0,17]); yticks([0:2:16]);
            % Save figure 
                file_name_eps = 'Model_2.1/Entrepreneur_Age_Profile_low_ineq.eps'     ;
                file_name_png = 'Model_2.1/png/Entrepreneur_Age_Profile_low_ineq.png' ;
                file_name_fig = 'Model_2.1/fig/Entrepreneur_Age_Profile_low_ineq.fig'  ;
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)


        % Low initial productivity model
         figure; grid on; hold on;
            scatter(20:5:100,low_z0_profile(ind),'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            ylabel('Share of Entrepreneurs (%)','FontName','Times New Roman','FontSize',14);
            xlabel('Age','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([19,81]); xticks([20:5:100]); 
            ylim([0,17]); yticks([0:2:16]);
            % Save figure 
                file_name_eps = 'Model_2.1/Entrepreneur_Age_Profile_low_initial_prod.eps'     ;
                file_name_png = 'Model_2.1/png/Entrepreneur_Age_Profile_low_initial_prod.png' ;
                file_name_fig = 'Model_2.1/fig/Entrepreneur_Age_Profile_low_initial_prod.fig'  ;
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)


    % Joint graph
        figure; grid on; hold on;
            scatter(20:5:100,low_z0_profile(ind),'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0.4 0.4 0.4],...
                'LineWidth',1.5,...
                'MarkerFaceColor',0.7*[1,1,1]+0.3*[0.4 0.4 0.4]);
            scatter(20:5:100,low_ineq_profile(ind),'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',Model_Color_a,...
                'LineWidth',1.5,...
                'MarkerFaceColor',0.7*[1,1,1]+0.3*Model_Color_a);
            scatter(20:5:100,bench_profile(ind),'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            text(61,3.9,'Baseline Model','color',[0 0.447058826684952 0.74117648601532],'FontName','Times New Roman','FontSize',14)
            text(56,13,'Low Inequality Calibration','color',Model_Color_a,'FontName','Times New Roman','FontSize',14)
            text(31,5.5,'Low Initial Productivity','color',[0.4 0.4 0.4],'FontName','Times New Roman','FontSize',14)
            ylabel('Share of Entrepreneurs (%)','FontName','Times New Roman','FontSize',14);
            xlabel('Age','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([19,81]); xticks([20:5:100]); 
            ylim([0,17]); yticks([0:2:16]);
            % Save figure 
                file_name_eps = 'Model_2.1/Entrepreneur_Age_Profile_all.eps'     ;
                file_name_png = 'Model_2.1/png/Entrepreneur_Age_Profile_all.png' ;
                file_name_fig = 'Model_2.1/fig/Entrepreneur_Age_Profile_all.fig'  ;
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)





%% Figure C.3.A (Intergenerational Rank-Rank Correlation, Model + Norwegian Data)

    % Load Simulation
        load ../Model_2.1/Simul/IGM_3050/panela_parents
        load ../Model_2.1/Simul/IGM_3050/panelr_parents
        load ../Model_2.1/Simul/IGM_3050/panelz_parents
        load ../Model_2.1/Simul/IGM_3050/panela_sons
        load ../Model_2.1/Simul/IGM_3050/panelr_sons
        load ../Model_2.1/Simul/IGM_3050/panelz_sons

        load ../Model_2.1/Bench_Files/EBAR ;
        EBAR_bench    = EBAR*0.727853584919652 ; 

    % Load Fagereng's data
        table = readtable('./Fagereng_Data/Fagereng_IGM_data.csv') ;
        
    % Dollar units
        panela_old  = log(panela_parents*EBAR_data/EBAR_bench)  ; panelr_old = panelr_parents ; panelz_old = panelz_parents ;
        panela_new  = log(panela_sons*EBAR_data/EBAR_bench)     ; panelr_new = panelr_sons    ; panelz_new = panelz_sons    ;
        
        clear panela_parents panela_sons panelr_parents panelr_sons panelz_parents panelz_sons
        
        panelage_old = ones(size(panela_old)) ; panelage_new = ones(size(panela_old)) ;
        
        x            = panela_old ; 
        y            = panela_new ;  
        N_y          = numel(y)   ;
       
	% Regression Log-Levels
        mdl          = fitlm(x',y')                 ;
        cons         = mdl.Coefficients.Estimate(1) ;
        slope        = mdl.Coefficients.Estimate(2) ;
        SE_slope     = mdl.Coefficients.SE(2)       ;
        R2           = mdl.Rsquared.Ordinary        ;
        n            = numel(x)                     ;
        xx           = linspace(min(x),max(x),3)    ;    
        
    % Regression Rank - Assets 
        % Bin cutoffs
        prct_x       = [min(x);prctile(x,[1:100]')];
        % Bins
        for j=1:100
            y_vec    = y( (x>=prct_x(j))&(x<=prct_x(j+1)) )  ; 
            y_aux(j) = 100*sum(y<=mean(y_vec))/N_y           ;
            x_aux(j) = j                                     ;
        end
        % Regression
        mdl             = fitlm(x_aux',y_aux')         ;
        cons_Ra         = mdl.Coefficients.Estimate(1) ;
        slope_Ra        = mdl.Coefficients.Estimate(2) ;
        SE_slope_Ra     = mdl.Coefficients.SE(2)       ;
        R2_Ra           = mdl.Rsquared.Ordinary        ;
        xx              = linspace(0,100,3)            ;

        mdl             = fitlm(x_aux',table.av_pct_child) ;
        cons_F         = mdl.Coefficients.Estimate(1) ;
        slope_F        = mdl.Coefficients.Estimate(2) ;
        
        % Graph Norwegian data
            figure; grid on; hold on; 
            daspect([1 1 1])
            scatter(x_aux,table.av_pct_child,'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            plot(xx, cons_F + slope_F*xx,'k--','LineWidth',1); % Regression line
            % plot(x_aux,x_aux,'--k') % 45º line
            ylabel('Offspring''s Average Wealth Percentile','FontName','Times New Roman','FontSize',14);
            xlabel('Father''s Wealth Percentile','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([0,101]); xticks([0:20:100]); 
            ylim([25,95]); yticks([30:10:90]);
            % Save figure 
                fig_title     = 'Inter-Generational Wealth Rank 30-50'                  ;
                file_name_eps = 'Model_2.1/Inter_Generational_Rank_Wealth_3050_Fagereng.eps'     ;
                file_name_png = 'Model_2.1/png/Inter_Generational_Rank_Wealth_3050_Fagereng.png' ;
                file_name_fig = 'Model_2.1/fig/Inter_Generational_Rank_Wealth_3050_Fagereng.fig'  ;
            %title(fig_title); 
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)


        % Graph Model 
            figure; grid on; hold on; 
            daspect([1 1 1])
            scatter(x_aux,y_aux,'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            plot(xx, cons_Ra + slope_Ra*xx,'k--','LineWidth',1); % Regression line
            % plot(x_aux,x_aux,'--k') % 45º line
            ylabel('Offspring''s Average Wealth Percentile','FontName','Times New Roman','FontSize',14);
            xlabel('Father''s Wealth Percentile','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([0,101]); xticks([0:20:100]); 
            ylim([25,95]); yticks([30:10:90]);
            % Save figure 
                fig_title     = 'Inter-Generational Wealth Rank 30-50'                  ;
                file_name_eps = 'Model_2.1/Inter_Generational_Rank_Wealth_3050.eps'     ;
                file_name_png = 'Model_2.1/png/Inter_Generational_Rank_Wealth_3050.png' ;
                file_name_fig = 'Model_2.1/fig/Inter_Generational_Rank_Wealth_305.fig'  ;
            %title(fig_title); 
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)

            % Add label 
            text(65,42,'Benchmark Model','color',[0 0.447058826684952 0.74117648601532],'FontName','Times New Roman','FontSize',14)
            text(70,39,'(slope=0.23)','color',[0 0.447058826684952 0.74117648601532],'FontName','Times New Roman','FontSize',12)

            % Add Fagereng data
            scatter(1:100,table.av_pct_child,'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',Model_Color_a,...
                'LineWidth',1.5,...
                'MarkerFaceColor',0.7*[1,1,1]+0.3*Model_Color_a);
            plot(xx, cons_F + slope_F*xx,'k--','LineWidth',1); % Regression line
            text(7,56,'Norway {\it (Fagereng, et al, 2020)}','color',Model_Color_a,'FontName','Times New Roman','FontSize',14)
            text(21,53,'(slope=0.16)','color',Model_Color_a,'FontName','Times New Roman','FontSize',12)
            
            % Save figure 
                fig_title     = 'Inter-Generational Wealth Rank 30-50'                  ;
                file_name_eps = 'Model_2.1/Inter_Generational_Rank_Wealth_3050_vs_Fagereng.eps'     ;
                file_name_png = 'Model_2.1/png/Inter_Generational_Rank_Wealth_3050_vs_Fagereng.png' ;
                file_name_fig = 'Model_2.1/fig/Inter_Generational_Rank_Wealth_3050_vs_Fagereng.fig'  ;
            %title(fig_title); 
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)

        % Graph Model vs Norwegian Data
            figure; grid on; hold on; 
            yyaxis left
            scatter(x_aux,y_aux,'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
                'LineWidth',1.5,...
                'MarkerFaceColor',[0.678431391716003 0.921568632125854 1]);
            plot(xx, cons_Ra + slope_Ra*xx,'k--','LineWidth',1); % Regression line
            % plot(x_aux,x_aux,'--k') % 45º line
            text(65,42,'Benchmark Model','color',[0 0.447058826684952 0.74117648601532],'FontName','Times New Roman','FontSize',14)
            text(69,39,'(slope=0.23)','color',[0 0.447058826684952 0.74117648601532],'FontName','Times New Roman','FontSize',12)
            % Labels 
            ylabel('Offspring''s Average Wealth Percentile','FontName','Times New Roman','FontSize',14);
            xlabel('Father''s Wealth Percentile','FontName','Times New Roman','FontSize',14);
            set(gca,'FontName','Times New Roman','FontSize',14);
            xlim([0,101]); xticks([0:20:100]); 
            ylim([25,95]); yticks([30:10:90]);
            
            % Add Fagereng data
            yyaxis right
            scatter(1:100,table.av_pct_child,'DisplayName','data1','SizeDataSource','6',...
                'MarkerEdgeColor',Model_Color_a,...
                'LineWidth',1.5,...
                'MarkerFaceColor',0.7*[1,1,1]+0.3*Model_Color_a);
            plot(xx, cons_F + slope_F*xx,'k--','LineWidth',1); % Regression line
            text(13,53,'Norway {\it (Fagereng, et al, 2020)}','color',Model_Color_a,'FontName','Times New Roman','FontSize',14)
            text(25,51,'(slope=0.16)','color',Model_Color_a,'FontName','Times New Roman','FontSize',12)
            ylim([40,75]); yticks([40:5:70]);
            ylabel('Offspring''s Average Wealth Percentile','FontName','Times New Roman','FontSize',14);
            
            hold off;
            % Save figure 
                fig_title     = 'Inter-Generational Wealth Rank 30-50'                  ;
                file_name_eps = 'Model_2.1/Inter_Generational_Rank_Wealth_3050_vs_Fagereng_yy.eps'     ;
                file_name_png = 'Model_2.1/png/Inter_Generational_Rank_Wealth_3050_vs_Fagereng_yy.png' ;
                file_name_fig = 'Model_2.1/fig/Inter_Generational_Rank_Wealth_305_vs_Fagereng_yy.fig'  ;
            %title(fig_title); 
            print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig)
        
 



%% Figure C.4 (Average After-Tax Labor and Capital Income vs Capital Tax Revenues)

% Read Tables 
    OTW_Table = readtable('/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Model_2.1/Opt_Tax_W/Wide_Grid/Stats_by_tau_w.txt') ;
    OTK_Table = readtable('/Users/ocamp020/Dropbox/RA_Guvenen/Wealth_Tax/CGGK_CODES/Sergio/Model_2.1/Opt_Tax_K/Wide_Grid/Stats_by_tau_k.txt') ;
    N         = numel(OTW_Table.Av_Util_NB) ;

% Capital Tax Revenue Over Total Tax Revenue
    K_Tax_Bench = 0.99*0.236992059442476 ;
    [~,OTK_grid]= max(OTK_Table.Av_Util_NB);
    K_Tax_OTK   = OTK_Table.GBAR_K_Tax_Rev_bench(OTK_grid); % -0.131283009707177 ;
    [~,OTW_grid]= max(OTW_Table.Av_Util_NB);
    K_Tax_OTW   = OTW_Table.GBAR_K_Tax_Rev_bench(OTW_grid); % 0.355360711880188 ;

    
% Figure: After tax labor and capital income 
    % Get after tax labor and capital income 
        OTK_W_Inc = OTK_Table.psi.*OTK_Table.wage.*OTK_Table.NBAR ; 
        OTK_K_Inc = (1-OTK_Table.tauK).*OTK_Table.Tot_Cap_Inc    ; 
        OTW_W_Inc = OTW_Table.psi.*OTW_Table.wage.*OTW_Table.NBAR ; 
        OTW_K_Inc = OTW_Table.Tot_Cap_Inc - OTW_Table.tauW_at.*OTW_Table.KBAR ; 
    % Get Change in after tax  and Q 
        [~,ind_0] = min(abs(OTK_Table.tauK)) ; 
        OTK_W_Inc = 100*(OTK_W_Inc/OTK_W_Inc(ind_0)-1) ; 
        OTK_K_Inc = 100*(OTK_K_Inc/OTK_K_Inc(ind_0)-1) ; 
        [~,ind_0] = min(abs(OTW_Table.tauW_at)) ; 
        OTW_W_Inc = 100*(OTW_W_Inc/OTW_W_Inc(ind_0)-1) ; 
        OTW_K_Inc = 100*(OTW_K_Inc/OTW_K_Inc(ind_0)-1) ; 
    % Make Figure 
    graph_ind = 1:2:N ;
    x_lim = [-0.4  0.5];
    y_lim_left = [-20.0 20.0]; y_lim_right = [-40.0 40.0];
    figure; hold on; 
    % Graph Income 
    plot(OTK_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTK_W_Inc(graph_ind),'r'  ,'linewidth',2.5)
    plot(OTW_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTW_W_Inc(graph_ind),'-'  ,'linewidth',2.5,'color',[0.3,0.3,0.3])
    yyaxis right
    plot(OTK_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTK_K_Inc(graph_ind),'r-.','linewidth',2.5)
    yyaxis right
    plot(OTW_Table.GBAR_K_Tax_Rev_bench(graph_ind),OTW_K_Inc(graph_ind),'--' ,'linewidth',2.5,'color',[0.3,0.3,0.3])
    % Lines 
    % plot([0,0]                ,y_lim,'--','linewidth',1.0,'color',[0.4,0.4,0.4])
    plot([K_Tax_OTK,K_Tax_OTK],y_lim,'--','linewidth',1.0,'color',[0.4,0.4,0.4])
    plot([K_Tax_OTW,K_Tax_OTW],y_lim,'--','linewidth',1.0,'color',[0.4,0.4,0.4])
    % Arrows 
    annotation('textarrow',[0.39,0.34],[0.87,0.87],'String','Opt. $\tau_k$','Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    annotation('textarrow',[0.72,0.77],[0.87,0.87],'String','Opt. $\tau_a$','Interpreter','latex','FontSize',15,'Linewidth',2,'FontWeight','bold','color',[0.2,0.2,0.2])
    % Labels 
    xlim(x_lim) ;                    xlabel('Tax Revenue from K $/$ Total Tax Revenue','Interpreter','latex','FontName','Times New Roman','FontSize',16);
    yyaxis left ; ylim(y_lim_left) ; ylabel('\% Deviation of After Tax Labor Income','Interpreter','latex','FontName','Times New Roman','FontSize',16); 
    yyaxis right; ylim(y_lim_right); ylabel('\% Deviation of After Tax Capital Income','Interpreter','latex','FontName','Times New Roman','FontSize',16); 
    grid(gca,'on'); set(gca,'fontsize',12,'YColor','k');
    legend({'Labor Income, $\tau_k$','Labor Income, $\tau_a$','Capital Income, $\tau_k$','Capital Income, $\tau_a$'},'Interpreter','latex','FontSize',13,'location','south','FontName','Times New Roman');  legend('boxoff');  
    % Save Figure 
    file_name_eps = 'Model_2.1/Opt_Tax_KvsW_Inc.eps' ;
    file_name_png = 'Model_2.1/png/Opt_Tax_KvsW_Inc.png' ;
    file_name_fig = 'Model_2.1/fig/Opt_Tax_KvsW_Inc.fig' ;
    print('-depsc',file_name_eps) ; print('-dpng',file_name_png) ; savefig(file_name_fig);
    
    


%% Table I (Summary of the Illustrative Example)

% Authors calculations 


%% Table II (Parameters for the Benchmark Model)

% Name and value of parameters reported in the paper


%% Table III (Targeted and Untargeted Moments: Model versus Data)

% Moments reported in "output.txt" from the Fortran code


%% Table IV (Distribution of Rates of Return (Untargeted) in the Model and the Data)

% Value of rates of return computed from the simulation of the model
% See output file "Simul/Asset_Return_Panel/Return_Stats_ben.txt"
% We report returns for 25-75 year olds (third row of the first table)


%% Table V (Tax Reform: Change in Macro Variables from Current US Benchmark)

% Value of aggregate variables reported in  in "output.txt" from the Fortran code
% See "Model_2.1/output.txt" and "Model_2.1/Tax_Reform/output.txt" 


%% Table VI (Average Welfare Gain from Tax Reform)

% Value reported at the bottom of "Model_2.1/Tax_Reform/output.txt" 
% Detailed results in "Model_2.1/Tax_Reform/CE.txt",
% "Model_2.1/Tax_Reform/CE_newborn.txt", "Model_2.1/Tax_Reform/CE_by_age.txt"  


%% Table VII (Welfare Gain/Loss by Age Group and Entrepreneurial Ability)

% Table reported in "Model_2.1/Tax_Reform/draft_group_CE.txt" 
% Last two columns correspond to partitions of the top 0.1% and are not
% reported in the paper



%% Table VIII (Optimal Taxation: Tax Rates and Average Welfare Effects)

% Value reported at the bottom of "Model_2.1/Opt_Tax_W/output.txt" 
% Detailed results in "Model_2.1/Opt_Tax_W/CE.txt",
% "Model_2.1/Opt_Tax_W/CE_newborn.txt", "Model_2.1/Opt_Tax_W/CE_by_age.txt"a

% Values for OKIT are in the "Opt_Tax_K" and for the wealth tax with
% exemption threshold in the "Opt_Tax_W_Threshold" folder.

% Values for the low inequality calibration are in the
% "Model_2.1_Match_Return_Lambda" folder and its sub-folders 



%% Table IX (Optimal Taxation: Changes in Macroeconomic Outcomes)

% Value of aggregate variables reported in  in "output.txt" from the Fortran code
% See "Model_2.1/output.txt" and "Model_2.1/Opt_Tax_W/output.txt" 
% See also the "Opt_Tax_K" and "Opt_Tax_W_Threshold" folders 


%% Table X (Decomposition of Welfare Gains)

% Table reported in "Model_2.1/Tax_Reform/CE_Decomposition",
% "Model_2.1/Opt_Tax_W/CE_Decomposition", and "Model_2.1/Opt_Tax_K/CE_Decomposition"
% The joint decomposition is reported in the paper (last three rows of each
% table)


%% Table XI (Policy Analysis Accounting for the Transition Path)

% Welfare gains are reported in "Model_2.1/Transition_OTK_tl/CE_Transition.txt" 
% and "Model_2.1/Transition_OTW_tl/CE_Transition.txt" 


%% Table XII (Robustness: Optimal Wealth Tax)

% Aggregates and welfare gains in the "ouput.txt" file from Fortran code
% for each of the extensions 


%% Table B.1 (Forbes Self-Made Index)

% Description provided in paper


%% Table B.2 (Concentration of Capital Income and Wealth)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load data 
    load ../Model_2.1_fine_grid/Simul/panela_bench ;
    load ../Model_2.1_fine_grid/Simul/panelage_bench ;
    load ../Model_2.1_fine_grid/Simul/panel_YL_bench ;
    load ../Model_2.1_fine_grid/Simul/panel_Ret_a_bench ;
    load ../Model_2.1_fine_grid/Simul/panelx_bench ;
    load ../Model_2.1_fine_grid/Simul/panelz_bench ;
    load ../Model_2.1_fine_grid/Bench_Files/EBAR ;
    panel_x     = panelx_bench          ; clear panelx_bench   ;
    panel_z     = panelz_bench          ; clear panelz_bench   ;
    panel_age   = panelage_bench        ; clear panelage_bench ;
    panel_L_inc = panel_YL_bench        ; clear panelYL_bench  ;
    panel_R     = panel_Ret_a_bench     ; clear panel_Ret_a_bench ;
    panel_K_inc = panel_R.*panela_bench';
    EBAR_bench  = EBAR*0.727853584919652 ;


    % Load auxiliary variables (taken form output.txt of Model_2.1)
    x_hi    = 1.50  ;
    tau_L   = 0.224 ; 
    DepRate = 0.05  ; 
    mu      = 0.90  ;  
    P       = 0.12476912139079619   ;
    R       = 0.026729010440549166 ; 
    theta   = [1.000 ; 1.225 ; 1.450 ; 1.675 ; 1.900 ; 2.125 ; 2.350 ; 2.575 ; 2.800 ] ;

    % Productivity 
    z_grid  = [0.42979057335122889       0.57304504679830393       0.75699738889794321        1.0000000000000000        1.3210085195350887        1.7450635086842872        2.3052437621017381        3.0452466493415149        4.0227967678658256]  ; 
    xz_grid = zeros(3,numel(z_grid)) ; 
        xz_grid(1,:) = exp(log(z_grid)*x_hi ) ;
		xz_grid(2,1:4) = xz_grid(1,1:4); xz_grid(2,5:end) = z_grid(5:end) ;
    for i=1:numel(panel_z)
    panel_xz(i) = xz_grid(panel_x(i),panel_z(i)) ; 
    end 

    % Pass labor income to before tax for comparisson with capital income
    panel_L_inc_aux = panel_L_inc ; 
    panel_L_inc_aux(panel_age<Ret_Age) = panel_L_inc(panel_age<Ret_Age)/(1-tau_L) ; 

    % Get capital demand 
    panel_K = min(theta(panel_z)'.*panela_bench,(mu*P*panel_xz.^mu/(R+DepRate)).^(1/(1-mu)) ) ; 

    % Modify K inc to only include profits 
    panel_Pr = panel_K_inc'-R*panela_bench+R*min(panela_bench,panel_K) ; 

    % Marginal Returns 
    panel_Mrg_R = R + (mu*P*((panel_xz'.*theta(panel_z)).^mu).*(panela_bench'.^(mu-1)) - (R+DepRate)*theta(panel_z)) ...
							.*(panel_K'>=(theta(panel_z).*panela_bench')) ;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Wealth Shares
    

    % Get percentiles 
    percentile = [0.999;0.995;0.99;0.90;0.50];
    prc        = prctile(panela_bench,percentile*100);
    
    % Get Shares
    total_assets = sum(panela_bench) ;
    total_K_inc  = sum(panel_K_inc)  ;
    share        = NaN(numel(prc),1) ;
    K_inc_share  = NaN(numel(prc),1) ;
    for i=1:numel(prc)
        share(i,1)       = 100*sum(panela_bench(panela_bench>=prc(i)))/total_assets ;
        K_inc_share(i,1) = 100*sum(panel_K_inc( panela_bench>=prc(i)))/total_K_inc  ;
    end 
    
    % Report Shares 
    Mat = [{'Top_x%','Wealth_Share','K_Inc_Share_of_Top_Wealth'};num2cell([100-100*percentile , share, K_inc_share])] ;
    disp(Mat)
    writecell(Mat,Tables_file,'Sheet','Wealth_Shares')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Top Capital Income Shares 

    % Get percentiles 
    percentile = [0.999;0.995;0.99;0.90;0.50];
    prc        = prctile(panel_K_inc,percentile*100);
    
    % Get Shares
    total_assets = sum(panela_bench) ;
    total_K_inc  = sum(panel_K_inc)  ;
    share   = NaN(numel(prc),1) ;
    A_share = NaN(numel(prc),1) ;
    for i=1:numel(prc)
        share(i,1)   = 100*sum(panel_K_inc(panel_K_inc>=prc(i)))/total_K_inc ;
        A_share(i,1) = 100*sum(panela_bench(panel_K_inc>=prc(i)))/total_assets ;
    end 
    
    % Report Shares 
    Mat = [{'Top x%','K_Inc Share','Wealth_share_of_top_K_inc'};num2cell([100-100*percentile , share, A_share])] ;
    disp(Mat)
    writecell(Mat,Tables_file,'Sheet','K_Inc_Shares')



%% Table B.3 (Optimal Tax Experiments: Distribution of Welfare Gains and Losses)

% Values reported in "draft_group_CE.txt" and "draft_group_fpos_welfare.txt"
% located in the folders "Model_2.1/Opt_Tax_W/" and "Model_2.1/Opt_Tax_K/"
% and "Model_2.1_Match_Return_Lambda/Opt_Tax_W_Threshold/"


%% Table B.4 (Welfare Change with Transition: Switch to Optimal Tax System with Transition)

% Values reported in "draft_group_CE.txt" and "draft_group_fpos_welfare.txt"
% located in the folders "Model_2.1/Transition_OTW_tl/" and "Model_2.1/Transition_OTK_tl/"


%% Table E.5 (Welfare Change: L-INEQ Calibration)

% Values reported in "draft_group_CE.txt" and "draft_group_fpos_welfare.txt"
% located in the folders "Model_2.1_Match_Return_Lambda/Tax_Reform/",
% "Model_2.1_Match_Return_Lambda/Opt_Tax_W/", and "Model_2.1/Opt_Tax_K/"


%% Table E.6 (Tax Reform: Change in Macro Variables from Low Inequality Calibration)

% Value of aggregate variables reported in  in "output.txt" from the Fortran code
% See "Model_2.1_Match_Returns_Lambda/output.txt" and 
% "Model_2.1_Match_Returns_Lambda/Tax_Reform/output.txt" 


%% Table E.7 (Corporate Sector: Optimal Wealth Tax)

% Value of aggregate variables reported in "output.txt" from the Fortran code
% See "Model_2.1_Corp_CD_1.5/output.txt" and "Model_2.1_Corp_CD_1.5/Opt_Tax_W/output.txt" 
% See also the "Opt_Tax_K" and "Tax_Reform" folders 

%% Table E.8 (Decomposition of the Change in Output)

% Values obtained from comparing aggregates from "output.txt" from the Fortran code
% from folders "Model_2.1" and "Model_2.1_Corp_CD_Large_Corp"


%% Table E.9 (Public Firms: Optimal Wealth Tax)

% See output in folder "Model_2.1_IPO_High_Entry"

%% Table E.10 (Additional Robustness and Extensions)

% Aggregates and welfare gains in the "ouput.txt" file from Fortran code
% for each of the extensions 

%% Table E.11 (Robustness Additional Results: Optimal Wealth Tax)

% Aggregates and welfare gains in the "ouput.txt" file from Fortran code
% for each of the extensions 


