
Module Simulation_Module
	use parameters
	use global
	use Toolbox
	    
	Contains

!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE  SIMULATION(bench_indx)
	use parameters
	use global
	use omp_lib
	IMPLICIT NONE	
	integer, intent(in) :: bench_indx
	! Break for comment
		integer  :: currentzi, currentlambdai, currentei, currentxi, thread_num
		REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempnox, tempno, currenta, currentY
		REAL(DP) :: start_timet, finish_timet, h_i
		INTEGER  :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi, paneli, simutime
		INTEGER , DIMENSION(MaxAge) :: requirednumberby_age, cdfrequirednumberby_age
		INTEGER , DIMENSION(totpop) :: panelage, panelz, panellambda, panele, panelx, panelz_old, panellambda_old
		REAL(SP), DIMENSION(totpop) :: panela, panelPV_a, panelK, panel_Y_L   

		! Intergenerational statistics
		INTEGER , DIMENSION(totpop) 			  :: eligible, death_count
		REAL(SP), DIMENSION(totpop) 			  :: panela_parents, panela_sons
		REAL(SP), DIMENSION(:)      , allocatable :: eligible_panela_parents, eligible_panela_sons
		INTEGER , DIMENSION(totpop) 			  :: panelage_parents, panelage_sons
		INTEGER , DIMENSION(:)      , allocatable :: eligible_panelage_parents, eligible_panelage_sons
		INTEGER                     			  :: n_eligible
		! Intergenerational statistics 30-50
		REAL(SP), DIMENSION(totpop) 	     :: assets_dad, assets_son, return_dad, return_son, PV_dad, PV_son
		INTEGER , DIMENSION(totpop) 	     :: age_dad, age_son, z_dad, z_son
		REAL(SP), DIMENSION(2,4000000)       :: IGM_a_matrix, IGM_r_matrix, IGM_pv_matrix
		INTEGER , DIMENSION(2,4000000) 		 :: IGM_z_matrix
		REAL(SP), DIMENSION(:) , allocatable :: panela_dad, panela_son, panelz_dad, panelz_son
		REAL(SP), DIMENSION(:) , allocatable :: panelr_dad, panelr_son, panelPV_dad, panelPV_son
		INTEGER 						     :: IGM_index
		! ! Intergenerational statistics 40-60
		! REAL(SP), DIMENSION(totpop) 	     :: assets_dad_2, assets_son_2, return_dad_2, return_son_2, PV_dad_2, PV_son_2
		! INTEGER , DIMENSION(totpop) 	     :: age_dad_2, age_son_2, z_dad_2, z_son_2
		! REAL(SP), DIMENSION(2,4000000)       :: IGM_a_matrix_2, IGM_r_matrix_2, IGM_pv_matrix_2
		! INTEGER , DIMENSION(2,4000000) 		 :: IGM_z_matrix_2
		! REAL(SP), DIMENSION(:) , allocatable :: panela_dad_2, panela_son_2, panelz_dad_2, panelz_son_2
		! REAL(SP), DIMENSION(:) , allocatable :: panelr_dad_2, panelr_son_2, panelPV_dad_2, panelPV_son_2
		! INTEGER 						     :: IGM_index_2
		! Average Return by age group
		INTEGER , PARAMETER  				 :: ret_size=100000
		REAL(DP), DIMENSION(ret_size) 	     :: ret_aux=0.0_sp, ret_20, ret_21_25 , ret_26_30 , ret_31_35 , ret_36_40 
		REAL(DP), DIMENSION(ret_size) 	     :: ret_41_45 , ret_46_50 , ret_51_55 , ret_56_60 , ret_61_65 , ret_66_70
		REAL(DP), DIMENSION(ret_size) 	     :: ret_w_aux=0.0_sp, ret_w_20, ret_w_21_25 , ret_w_26_30 , ret_w_31_35 , ret_w_36_40 
		REAL(DP), DIMENSION(ret_size) 	     :: ret_w_41_45 , ret_w_46_50 , ret_w_51_55 , ret_w_56_60 , ret_w_61_65 , ret_w_66_70
		REAL(DP), DIMENSION(ret_size) 	     :: ret_k_aux=0.0_sp, ret_k_20=0.0_sp, ret_k_21_25=0.0_sp, ret_k_26_30=0.0_sp
		REAL(DP), DIMENSION(ret_size) 	     :: ret_k_31_35=0.0_sp, ret_k_36_40=0.0_sp, ret_k_41_45=0.0_sp , ret_k_46_50=0.0_sp
		REAL(DP), DIMENSION(ret_size) 	     :: ret_k_51_55=0.0_sp , ret_k_56_60=0.0_sp , ret_k_61_65=0.0_sp , ret_k_66_70=0.0_sp
		REAL(DP), DIMENSION(ret_size) 	     :: ret_k_w_aux=0.0_sp, ret_k_w_20=0.0_sp, ret_k_w_21_25=0.0_sp, ret_k_w_26_30=0.0_sp
		REAL(DP), DIMENSION(ret_size) 	     :: ret_k_w_31_35=0.0_sp, ret_k_w_36_40=0.0_sp, ret_k_w_41_45=0.0_sp, ret_k_w_46_50=0.0_sp
		REAL(DP), DIMENSION(ret_size) 	     :: ret_k_w_51_55=0.0_sp, ret_k_w_56_60=0.0_sp, ret_k_w_61_65=0.0_sp, ret_k_w_66_70=0.0_sp
		INTEGER , DIMENSION(ret_size) 	     :: Ind_K=0, Ind_K_20, Ind_K_21_25 , Ind_K_26_30 , Ind_K_31_35 , Ind_K_36_40 
		INTEGER , DIMENSION(ret_size) 	     :: Ind_K_41_45 , Ind_K_46_50 , Ind_K_51_55 , Ind_K_56_60 , Ind_K_61_65 , Ind_K_66_70
		REAL(DP), DIMENSION(ret_size) 	     :: cum_assets, cum_K
		REAL(DP) 						 	 :: Std_Dev_Return_Age(11)    , Mean_Return_Age(11)    , prc_Return_Age(11,9)
		REAL(DP) 						 	 :: Std_Dev_Return_W_Age(11)  , Mean_Return_W_Age(11)  , prc_Return_W_Age(11,9)
		REAL(DP) 						 	 :: Std_Dev_Return_K_Age(11)  , Mean_Return_K_Age(11)  , prc_Return_K_Age(11,9)
		REAL(DP) 						 	 :: Std_Dev_Return_K_W_Age(11), Mean_Return_K_W_Age(11), prc_Return_K_W_Age(11,9)
		REAL(DP)  							 :: K_aux, prctile_ret(9)
		INTEGER 							 :: i_pct
		
		REAL :: k_igm

		! Top Agents 
		INTEGER       :: top_ind(80), panel_top_ind(totpop), top_ind_aux(80), n_top
		REAL(DP)      :: top_A(80), A_cut, A_hi, A_low
		character(10) :: top_folder

		call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/' )

		!$ call omp_set_num_threads(20)

			! Set Seeds for each thread
		!$OMP parallel
	    !$OMP critical
			thread_num = omp_get_thread_num()
			newiseed   = -5 - thread_num
			tempno     = omp_ran1(newiseed)
	       	! write(*,'("inside critical myt=",i4,f12.8)') thread_num,tempno,newiseed
	    !$OMP end critical
	    !$OMP end parallel

		print*,'SIMULATION STARTED'

	    age=1
	    requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
	    cdfrequirednumberby_age(age) = requirednumberby_age(age)
	    DO age=2,MaxAge
	        requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
	        cdfrequirednumberby_age(age) = requirednumberby_age(age) + cdfrequirednumberby_age(age-1)
	    ENDDO
	    ! If the total number of people are not equal to the total population, then I will add the remainder to the last age
	    requirednumberby_age(MaxAge) =  requirednumberby_age(MaxAge)-cdfrequirednumberby_age(MaxAge) + totpop
	    cdfrequirednumberby_age(MaxAge) = totpop

		!=====================================================================
		!                     GENERATE   INITIAL   PANEL
		!=====================================================================


		newiseed=-1

		!$omp parallel do private(tempnoage,age,tempnoz,zi,tempnolambda,lambdai,tempnoe,ei,xi)
		DO paneli=1,totpop

		! AGE
		   tempnoage = omp_ran1() ! ran1(newiseed)
		   age=1
		   DO WHILE (tempnoage*totpop .gt. cdfrequirednumberby_age(age))
		       age=age+1
		   ENDDO

		! Z   
		   tempnoz = omp_ran1() ! ran1(newiseed)
		   zi=1
		   DO WHILE (tempnoz .gt. cdf_Gz(zi))
		       zi=zi+1
		   ENDDO

		! X   
		   xi=1
		 
		! LAMBDA  
		   tempnolambda = omp_ran1() ! ran1(newiseed) 
		   lambdai=1
		   DO WHILE (tempnolambda .gt. cdf_Glambda(lambdai))
		       lambdai=lambdai+1
		   ENDDO

		! E   
		   tempnoe = omp_ran1() ! ran1(newiseed)   
		   ei=1
		   DO WHILE (tempnoe .gt. cdf_Ge_byage(age,ei))
		       ei=ei+1
		   ENDDO


		   panelage(paneli)		= age
		   panelz(paneli)		= zi
		   panellambda(paneli)	= lambdai
		   panele(paneli)		= ei
		   panelx(paneli)		= xi
		   
		ENDDO
		print*, sum(panelage)/real(totpop,8), sum(panelz)/real(totpop,8), sum(panele)/real(totpop,8)


		! SET INITIAL ASSET DISTRIBUTION
		panela            = 1.0_DP

		!=============================================================================
		!
		! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
		!
		!=============================================================================

		!call cpu_time(start_timet) 
		
		eligible    = 1 
		death_count = 0

		age_dad = 0 ; age_son = 0 ; assets_dad = 0.0_dp ; assets_son = 0.0_dp ; PV_dad = 0.0_dp ; PV_son = 0.0_dp;
		IGM_index = 1 ; IGM_a_matrix = 0.0_dp ; IGM_r_matrix = 0.0_dp ; IGM_z_matrix = 0 ; 
		IGM_pv_matrix = 0.0_dp ; IGM_pv_matrix = 0 ; 
		! age_dad_2 = 0 ; age_son_2 = 0 ; assets_dad_2 = 0.0_dp ; assets_son_2 = 0.0_dp ; PV_dad_2 = 0.0_dp ; PV_son_2 = 0.0_dp;
		! IGM_index_2 = 1 ; IGM_a_matrix_2 = 0.0_dp ; IGM_r_matrix_2 = 0.0_dp ; IGM_z_matrix_2 = 0 ; 
		! IGM_pv_matrix_2 = 0.0_dp ; IGM_pv_matrix_2 = 0 ; 
		
		print*, 'Starting Simutime loop'
		DO simutime=1, MaxSimuTime
			!$omp parallel do private(tempnoage,age,tempnoz,zi,tempnolambda,lambdai,tempnoe,ei,xi, &
			!$omp& currenta,currentzi,currentlambdai,currentei,currentxi,tklo,tkhi,tempno,k_igm)
		   	DO paneli=1,totpop
		    
		       	currenta  		= panela(paneli)
		       	age       		= panelage(paneli)
		       	currentzi      	= panelz(paneli)
		       	currentlambdai 	= panellambda(paneli) 
		       	currentei 		= panele(paneli)
		       	currentxi 		= panelx(paneli)
		        
				!  COMPUTE NEXT PERIOD'S ASSET
		       	if (age .lt. MaxAge) then
		            ! do linear interpolation here to find aprime. calling the function takes much more time
		            if (currenta .ge. amax) then
	                	tklo = na-1
	                elseif (currenta .lt. amin) then
	                    tklo = 1
	                else
	                    tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
		            endif 
		            tkhi = tklo + 1        

		            panela(paneli) = ((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai,currentei,currentxi) 		&
		                           	&  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei,currentxi)) &
		                            &  / ( agrid(tkhi) - agrid(tklo) )    
		            
		            if (panela(paneli)  .ge. amax) then
		                panela(paneli) = min(panela(paneli), amax) 
		            endif      
		            if (panela(paneli)  .lt. amin) then
		                panela(paneli) = max(panela(paneli), amin) 
		            endif      

		       	endif !age .lt. MaxAge
				!  NEXT PERIOD'S ASSET IS COMPUTED

		             
				! DRAW NEXT PERIOD'S AGE DBN
		      	tempnoage = omp_ran1() ! ran1(newiseed)  
		  
		      	IF (tempnoage .gt. survP(age)) THEN
					panelage(paneli)   = 1
				ELSE
					panelage(paneli)   = age+1
					panelz(paneli)     = currentzi
					panellambda(paneli)= currentlambdai   
		      	ENDIF
		    
				! DRAW Z and LAMBDA DISTRIBUTION FOR ONE-YEAR OLDS
		     	age = panelage(paneli)   
		     	IF (age .eq. 1) THEN   
					! Z      
			       	tempnoz = omp_ran1() ! ran1(newiseed) 
			       	zi=1
					DO WHILE (tempnoz .gt. cdf_pr_z(currentzi,zi))
						zi=zi+1
					ENDDO

					! X     
			       	xi = 1
		       
					! LAMBDA  
					tempnolambda = omp_ran1() ! ran1(newiseed) 
					lambdai=1
					DO WHILE (tempnolambda .gt. cdf_pr_lambda(currentlambdai,lambdai))
						lambdai=lambdai+1
					ENDDO

					! E       
					currentei = panele(paneli)
					ei = ne/2+1 ! ALL NEWBORN START FROM THE MEDIAN E  but if too many people died and started from median E, draw a new E for them

					panelz(paneli)		= zi    
					panellambda(paneli)	= lambdai
					panele(paneli) 		= ei 
					panelx(paneli)		= xi

	       		else  IF (age .gt. 1) THEN
	       			! !$omp critical
	       			! print*,'Test 1'
	       			currentxi = panelx(paneli)
	       			tempno 	  = omp_ran1() ! ran1(newiseed)   
		            xi 		  = 1
		            DO WHILE (tempno .gt. cdf_pr_x(currentxi,xi,currentzi,age-1))
		            	! print*, 'tempno',tempno,'xi',xi
		               xi = xi+1
		            ENDDO            
		            panelx(paneli)=xi          
		            ! print*,'Test 1.1'
		            IF (age.lt.RetAge) THEN
			            currentei = panele(paneli)   
			            tempno 	  = omp_ran1() ! ran1(newiseed)   
			            ei 		  = 1
			            DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
			               ei = ei+1
			            ENDDO            
			            panele(paneli)=ei 
		            ENDIF 
		            ! !$omp end critical          
		     	ENDIF ! new age==1

		     	! Inter-Generation Mobility 30-50
		     	if (IGM_index.le.4000000) then
		     		! Reset variables if son dies before 50
			     	if ((age.eq.1).and.(age_son(paneli).lt.31)) then 
			     		! !$omp critical
			     		! print*, ' Agent died', IGM_index, 'age_son',age_son(paneli), 'agent', paneli
			     		! !$omp end critical
			     		age_dad(paneli)    = 0 		; age_son(paneli)    = 0 
			     		assets_dad(paneli) = 0.0_dp ; assets_son(paneli) = 0.0_dp
			     		z_dad(paneli)      = 0      ; z_son(paneli)  	 = 0
			     		return_dad(paneli) = 0.0_dp ; return_son(paneli) = 0.0_dp
			     		PV_dad(paneli)     = 0.0_dp ; PV_son(paneli)     = 0.0_dp
			     	endif 
		     		! Update age of current "son"
		     			age_son(paneli)    = age
		     		! Update variables for agents between 30-50 
		     		if ((age.ge.11).and.(age.le.31)) then 
		     			k_igm = min(theta(panelz(paneli))*currenta,&
		     					&(mu*P*xz_grid(panelx(paneli),panelz(paneli))**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
			     		assets_son(paneli) = panela(paneli) + assets_son(paneli)
			     		return_son(paneli) = ( P*(xz_grid(panelx(paneli),panelz(paneli))*k_igm)**mu - (R+DepRate)*k_igm +&
		     								&   R*panela(paneli) )/panela(paneli) + return_son(paneli)
			     		if (panela(paneli) .ge. amax) then
					        tklo = na-1
					    elseif (panela(paneli) .lt. amin) then
				            tklo = 1
				        else
				            tklo = ((panela(paneli)- amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
					    endif    
					    tkhi = tklo + 1        

					    PV_son(paneli)    = (    (agrid(tkhi) - panela(paneli)) * & 
					    					&		V_Pr(age,tklo,panelz(paneli),panellambda(paneli),panele(paneli), panelx(paneli))    &
					                       	&  + (panela(paneli) - agrid(tklo)) * &
					                       	&		V_Pr(age,tkhi,panelz(paneli),panellambda(paneli),panele(paneli), panelx(paneli)) )  &
					                       	&  / ( agrid(tkhi) - agrid(tklo) )  + (1.0_dp+R)*panela(paneli)

			     		! !$omp critical
			     		! print*, ' Potential Agent', IGM_index, 'age_son',age_son(paneli), 'agent', paneli
			     		! !$omp end critical
			     	endif 
			     	! Generation change and Save results 
			     	if (age.eq.31) then 
			     		z_son(paneli) = panelz(paneli)
			     		!$omp critical
			     		!print*, ' Son is 50:', IGM_index, 'age_son',age_son(paneli), 'age_dad',age_dad(paneli)
			     		if ((age_dad(paneli).eq.31).and.(simutime.gt.1800)) then  
				     		IGM_a_matrix(1,IGM_index)  = assets_dad(paneli)
				     		IGM_a_matrix(2,IGM_index)  = assets_son(paneli)
				     		IGM_r_matrix(1,IGM_index)  = return_dad(paneli)
				     		IGM_r_matrix(2,IGM_index)  = return_son(paneli)
				     		IGM_pv_matrix(1,IGM_index) = PV_dad(paneli)
				     		IGM_pv_matrix(2,IGM_index) = PV_son(paneli)
				     		IGM_z_matrix(1,IGM_index)  = z_dad(paneli)
				     		IGM_z_matrix(2,IGM_index)  = z_son(paneli)
				     		IGM_index = IGM_index + 1
				     		! print*, ' Save result', IGM_index-1
			     		endif 
			     		!$omp end critical
			     		age_dad(paneli)    = 31
			     		assets_dad(paneli) = assets_son(paneli)
			     		return_dad(paneli) = return_son(paneli)
			     		PV_dad(paneli)     = PV_son(paneli)
			     		z_dad(paneli)      = panelz(paneli)
			     		assets_son(paneli) = 0.0_dp 
			     		return_son(paneli) = 0.0_dp  
			     		PV_son(paneli)     = 0.0_dp    		
			     		z_son(paneli)      = 0
			     	endif 
		     	endif

		     	! Inter-Generation Mobility 40-60
		     	! if (IGM_index_2.le.4000000) then
		     	! 	! Reset variables if son dies before 60
		     	! 	if ((age.eq.1).and.(age_son_2(paneli).lt.41)) then 
			     ! 		! !$omp critical
			     ! 		! print*, ' Agent died', IGM_index, 'age_son',age_son(paneli), 'agent', paneli
			     ! 		! !$omp end critical
			     ! 		age_dad_2(paneli)    = 0 	  ; age_son_2(paneli)    = 0 
			     ! 		assets_dad_2(paneli) = 0.0_dp ; assets_son_2(paneli) = 0.0_dp
			     ! 		z_dad_2(paneli)      = 0      ; z_son_2(paneli)  	 = 0
			     ! 		return_dad_2(paneli) = 0.0_dp ; return_son_2(paneli) = 0.0_dp
			     ! 		PV_dad_2(paneli)     = 0.0_dp ; PV_son_2(paneli)     = 0.0_dp
			     ! 	endif 
			     ! 	! Update age of current "son"
			     ! 		age_son_2(paneli)    = age 
		     	! 	! Update variables for agents between 40-60 
		     	! 	if ((age.ge.21).and.(age.le.41)) then 
		     	! 		k_igm = min(theta(panelz(paneli))*currenta,&
		     	! 				& (mu*P*xz_grid(panelx(paneli),panelz(paneli))**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
			     ! 		assets_son_2(paneli) = panela(paneli) + assets_son_2(paneli)
			     ! 		return_son_2(paneli) = ( P*(xz_grid(panelx(paneli),panelz(paneli))*k_igm)**mu - (R+DepRate)*k_igm +&
		     	! 							&   R*panela(paneli) )/panela(paneli) + return_son_2(paneli)
			     ! 		if (panela(paneli) .ge. amax) then
					   !      tklo = na-1
					   !  elseif (panela(paneli) .lt. amin) then
				    !         tklo = 1
				    !     else
				    !         tklo = ((panela(paneli)- amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
					   !  endif    
					   !  tkhi = tklo + 1        

					   !  PV_son_2(paneli)  = (   (agrid(tkhi) - panela(paneli)) * & 
					   !  					&		V_Pr(age,tklo,panelz(paneli),panellambda(paneli),panele(paneli), panelx(paneli))  &
					   !                     	&  + (panela(paneli) - agrid(tklo)) * &
					   !                     	&		V_Pr(age,tkhi,panelz(paneli),panellambda(paneli),panele(paneli), panelx(paneli)) ) &
					   !                     	&  / ( agrid(tkhi) - agrid(tklo) )  + (1.0_dp+R)*panela(paneli)
			     ! 		! !$omp critical
			     ! 		! print*, ' Potential Agent', IGM_index, 'age_son',age_son(paneli), 'agent', paneli
			     ! 		! !$omp end critical
			     ! 	endif 
			     ! 	! Generation change and Save results 
			     ! 	if (age.eq.41) then 
			     ! 		z_son_2(paneli) = panelz(paneli)
			     ! 		!$omp critical
			     ! 		!print*, ' Son is 50:', IGM_index, 'age_son',age_son(paneli), 'age_dad',age_dad(paneli)
			     ! 		if ((age_dad_2(paneli).eq.41).and.(simutime.gt.1800)) then  
				    !  		IGM_a_matrix_2(1,IGM_index_2) = assets_dad_2(paneli)
				    !  		IGM_a_matrix_2(2,IGM_index_2) = assets_son_2(paneli)
				    !  		IGM_r_matrix_2(1,IGM_index_2) = return_dad_2(paneli)
				    !  		IGM_r_matrix_2(2,IGM_index_2) = return_son_2(paneli)
				    !  		IGM_pv_matrix_2(1,IGM_index)  = PV_dad_2(paneli)
				    !  		IGM_pv_matrix_2(2,IGM_index)  = PV_son_2(paneli)
				    !  		IGM_z_matrix_2(1,IGM_index_2) = z_dad_2(paneli)
				    !  		IGM_z_matrix_2(2,IGM_index_2) = z_son_2(paneli)
				    !  		IGM_index_2 = IGM_index_2 + 1
				    !  		! print*, ' Save result', IGM_index-1
			     ! 		endif 
			     ! 		!$omp end critical
			     ! 		age_dad_2(paneli)    = 41
			     ! 		assets_dad_2(paneli) = assets_son_2(paneli)
			     ! 		return_dad_2(paneli) = return_son_2(paneli)
			     ! 		PV_dad_2(paneli)     = PV_son_2(paneli)
			     ! 		z_dad_2(paneli)      = panelz(paneli)
			     ! 		assets_son_2(paneli) = 0.0_dp    
			     ! 		return_son_2(paneli) = 0.0_dp
			     ! 		PV_son_2(paneli)     = 0.0_dp  
			     ! 		z_son_2(paneli)      = 0  		
			     ! 	endif 
		     	! endif

		     	! Average Return by age group
		     	if (paneli.le.ret_size) then 
		     	K_aux = min(theta(panelz(paneli))*panela(paneli),&
		     			&(mu*P*xz_grid(panelx(paneli),panelz(paneli))**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )
	     		ret_aux(paneli)     = ( P*(xz_grid(panelx(paneli),panelz(paneli))*K_aux)**mu - (R+DepRate)*K_aux +&
	     								&   R*panela(paneli) )/panela(paneli) + ret_aux(paneli)
	     		ret_w_aux(paneli)   = ( P*(xz_grid(panelx(paneli),panelz(paneli))*K_aux)**mu - (R+DepRate)*K_aux +&
	     								&   R*panela(paneli) ) + ret_w_aux(paneli)
	     		cum_assets(paneli)  = panela(paneli) + cum_assets(paneli)
	     		if (K_aux.gt.0.0_dp) then 
	     		ret_k_aux(paneli)   = ( P*(xz_grid(panelx(paneli),panelz(paneli))*K_aux)**mu - (R+DepRate)*K_aux +&
	     								&   R*panela(paneli) )/K_aux + ret_k_aux(paneli)
	     		ret_k_w_aux(paneli) = ( P*(xz_grid(panelx(paneli),panelz(paneli))*K_aux)**mu - (R+DepRate)*K_aux +&
	     								&   R*panela(paneli) ) + ret_k_w_aux(paneli)
	     		cum_K(paneli)       = K_aux + cum_K(paneli)
	     		ind_K(paneli) 		= 1 + ind_K(paneli)
	     		endif 

	     		if (age.eq.1) then
	     			ret_20(paneli)        = ret_aux(paneli) 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_20(paneli)      = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_20(paneli)      = ret_k_aux(paneli)/ind_K(paneli)
	     			ret_k_w_20(paneli)    = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_20(paneli)      = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
	     		elseif (age.eq.6) then
	     			ret_21_25(paneli)     = ret_aux(paneli)/6.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_21_25(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_21_25(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_21_25(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_21_25(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.11) then 
     				ret_26_30(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_26_30(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_26_30(paneli)   = ret_k_aux(paneli)/ind_K(paneli)
	     			ret_k_w_26_30(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_26_30(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.15) then 
     				ret_31_35(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_31_35(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_31_35(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_31_35(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_31_35(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.21) then 
     				ret_36_40(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_36_40(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_36_40(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_36_40(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_36_40(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.26) then 
     				ret_41_45(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_41_45(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_41_45(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_41_45(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_41_45(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.31) then 
     				ret_46_50(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_46_50(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_46_50(paneli)   = ret_k_aux(paneli)/ind_K(paneli)
	     			ret_k_w_46_50(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_46_50(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.36) then 
     				ret_51_55(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_51_55(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_51_55(paneli)   = ret_k_aux(paneli)/ind_K(paneli)
	     			ret_k_w_51_55(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_51_55(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.41) then 
     				ret_56_60(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_56_60(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_56_60(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_56_60(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_56_60(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.46) then 
     				ret_61_65(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_61_65(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp 
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_61_65(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_61_65(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_61_65(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			elseif (age.eq.51) then 
     				ret_66_70(paneli)     = ret_aux(paneli)/5.0_sp 
	     			ret_aux(paneli)       = 0.0_sp 
	     			ret_w_66_70(paneli)   = ret_w_aux(paneli)/cum_assets(paneli)
	     			ret_w_aux(paneli)     = 0.0_sp  
	     			cum_assets(paneli)    = 0.0_sp 
	     			if (cum_K(paneli).gt.0.0_sp) then 
	     			ret_k_66_70(paneli)   = ret_k_aux(paneli)/ind_K(paneli) 
	     			ret_k_w_66_70(paneli) = ret_k_w_aux(paneli)/cum_K(paneli)
	     			ind_K_66_70(paneli)   = 1
	     			endif 
	     			ret_k_aux(paneli)     = 0.0_sp 
	     			ret_k_w_aux(paneli)   = 0.0_sp 
	     			cum_K(paneli)         = 0.0_sp
	     			ind_K(paneli)         = 0
     			endif 
				endif 

			ENDDO ! paneli
			

			! Save state of fathers in last period
		     	if (simutime.eq.MaxSimuTime-1) then 
		     		panelz_old = panelz 
		     		panellambda_old = panellambda
	     		endif 

			! Save data on assets for the last periods
			! Agents are eligible if:
				! 1) They don't die during the first two recording periods
				! 2) They they die between the third recording period and the recording periods for the next generation
				! 3) They don't die again
				if (simutime.eq.(MaxSimuTime-54)) then 
			    	panela_parents   = panela
			    	! panelage_parents = panelage
		        endif 
		        if ((simutime.ge.(MaxSimuTime-(5+53))).and.(simutime.le.(MaxSimuTime-(5+51)))) then 
			    	panela_parents   = panela_parents  + panela
			    	! panelage_parents = panelage
			    	where(panelage==1) eligible = 0 
		        endif 
		        if (simutime.eq.(MaxSimuTime-(5+50))) then 
			    	panela_parents   = panela_parents   + panela
			    	panelage_parents = panelage
			    	where(panelage==1) eligible = 0 
		        endif 
		        if ((simutime.ge.(MaxSimuTime-(5+49))).and.(simutime.le.(MaxSimuTime-4))) then
		        	where(panelage==1) death_count = death_count + 1
		        endif
		        if (simutime.eq.(MaxSimuTime-4)) then 
			    	panela_sons   = panela
			    	! panelage_sons = panelage 
		        endif 
		        if ((simutime.ge.(MaxSimuTime-3)).and.(simutime.le.(MaxSimuTime-1))) then 
			    	panela_sons   = panela_sons   + panela
			    	! panelage_sons = panelage 
			    	where(panelage==1) eligible = 0 
		        endif 
		        if (simutime.eq.(MaxSimuTime)) then 
			    	! panela_sons   = panela_sons   + panela
			    	panelage_sons = panelage 
			    	where(panelage==1) eligible = 0 
		        endif 

		 		
		 		! print*, "Simulation period", simutime

		 		
			ENDDO ! simutime
			print*,' '
			print*,'Averages'
			print*, sum(panelage)/real(totpop,8), sum(panelz)/real(totpop,8), sum(panele)/real(totpop,8), sum(panela)/real(totpop,8)

			! IGM 30-50
				! Get mean of assets and return
				IGM_a_matrix  = IGM_a_matrix/real(21,8) 
				IGM_r_matrix  = IGM_r_matrix/real(21,8) 
				IGM_pv_matrix = IGM_pv_matrix/real(21,8) 
				! Get number of eligibles
				n_eligible = count(IGM_a_matrix(1,:).gt.0.0_dp)
				! Allocate variables
				allocate(panela_dad(n_eligible) , panela_son(n_eligible) )
				allocate(panelr_dad(n_eligible) , panelr_son(n_eligible) )
				allocate(panelpv_dad(n_eligible), panelpv_son(n_eligible))
				allocate(panelz_dad(n_eligible) , panelz_son(n_eligible) )
				panela_dad  = pack(IGM_a_matrix(1,:)  , (IGM_a_matrix(1,:).gt.0.0_dp) )
				panela_son  = pack(IGM_a_matrix(2,:)  , (IGM_a_matrix(2,:).gt.0.0_dp) )
				panelr_dad  = pack(IGM_r_matrix(1,:)  , (IGM_r_matrix(1,:).gt.0.0_dp) )
				panelr_son  = pack(IGM_r_matrix(2,:)  , (IGM_r_matrix(2,:).gt.0.0_dp) )
				panelpv_dad = pack(IGM_pv_matrix(1,:) , (IGM_pv_matrix(1,:).gt.0.0_dp))
				panelpv_son = pack(IGM_pv_matrix(2,:) , (IGM_pv_matrix(2,:).gt.0.0_dp))
				panelz_dad  = pack(IGM_z_matrix(1,:)  , (IGM_z_matrix(1,:).gt.0.0_dp) )
				panelz_son  = pack(IGM_z_matrix(2,:)  , (IGM_z_matrix(2,:).gt.0.0_dp) )
			! Print
				print*, 'IGM 30-50'
				print*, 'n_eligible', n_eligible, 'mean_panel_dad', sum(panela_dad)/n_eligible, 'mean_panel_son', sum(panela_son)/n_eligible
			! ! IGM 40-60
			! 	! Get mean of assets and return
			! 	IGM_a_matrix_2  = IGM_a_matrix_2/real(21,8) 
			! 	IGM_r_matrix_2  = IGM_r_matrix_2/real(21,8) 
			! 	IGM_pv_matrix_2 = IGM_pv_matrix_2/real(21,8) 
			! 	! Get number of eligibles
			! 	n_eligible = count(IGM_a_matrix_2(1,:).gt.0.0_dp)
			! 	! Allocate variables
			! 	allocate(panela_dad_2(n_eligible) , panela_son_2(n_eligible))
			! 	allocate(panelr_dad_2(n_eligible) , panelr_son_2(n_eligible))
			! 	allocate(panelpv_dad_2(n_eligible), panelpv_son_2(n_eligible))
			! 	allocate(panelz_dad_2(n_eligible) , panelz_son_2(n_eligible))
			! 	panela_dad_2  = pack(IGM_a_matrix_2(1,:)  , (IGM_a_matrix_2(1,:).gt.0.0_dp) )
			! 	panela_son_2  = pack(IGM_a_matrix_2(2,:)  , (IGM_a_matrix_2(2,:).gt.0.0_dp) )
			! 	panelr_dad_2  = pack(IGM_r_matrix_2(1,:)  , (IGM_r_matrix_2(1,:).gt.0.0_dp) )
			! 	panelr_son_2  = pack(IGM_r_matrix_2(2,:)  , (IGM_r_matrix_2(2,:).gt.0.0_dp) )
			! 	panelpv_dad_2 = pack(IGM_pv_matrix_2(1,:) , (IGM_pv_matrix_2(1,:).gt.0.0_dp))
			! 	panelpv_son_2 = pack(IGM_pv_matrix_2(2,:) , (IGM_pv_matrix_2(2,:).gt.0.0_dp))
			! 	panelz_dad_2  = pack(IGM_z_matrix_2(1,:)  , (IGM_z_matrix_2(1,:).gt.0.0_dp) )
			! 	panelz_son_2  = pack(IGM_z_matrix_2(2,:)  , (IGM_z_matrix_2(2,:).gt.0.0_dp) )
			! ! Print
			! 	print*, 'IGM 20-40'
			! 	print*, 'n_eligible', n_eligible, 'mean_panel_dad', sum(panela_dad_2)/n_eligible, 'mean_panel_son', sum(panela_son_2)/n_eligible


			! Mean of assets 
				panela_parents = panela_parents/5.0_dp 
				panela_sons    = panela_sons/5.0_dp 

			! Clean eligibles 
				where(death_count/=1) eligible = 0
				where(panelage_sons.lt.25) eligible = 0
				where(panelage_sons.gt.41) eligible = 0
				where(panelage_parents.lt.25) eligible = 0
				where(panelage_parents.gt.41) eligible = 0


			! Get data on intergenerational mobility
				n_eligible = sum(eligible)

				allocate( eligible_panela_parents(n_eligible), eligible_panela_sons(n_eligible) )

				eligible_panela_parents 	= pack(panela_parents   , (eligible.eq.1) )
				eligible_panela_sons    	= pack(panela_sons      , (eligible.eq.1) )

				allocate( eligible_panelage_parents(n_eligible), eligible_panelage_sons(n_eligible) )

				eligible_panelage_parents 	= pack(panelage_parents , (eligible.eq.1) )
				eligible_panelage_sons    	= pack(panelage_sons	, (eligible.eq.1) )
				

			print*, ' '
			print*, 'n_eligible', sum(eligible)
			print*, 'panela_parents', sum(eligible_panela_parents)/n_eligible, 'panela_sons', sum(eligible_panela_sons)/n_eligible
			print*, 'panelage_parents', sum(eligible_panelage_parents)/n_eligible, 'panelage_sons', sum(eligible_panelage_sons)/n_eligible
			print*, ' '

			print*, 'Start of std dev of return by age'
			! Std Dev of return by age
			Std_Dev_Return_Age(1)  = sqrt( sum( (ret_20   -sum(ret_20)   /ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(2)  = sqrt( sum( (ret_21_25-sum(ret_21_25)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(3)  = sqrt( sum( (ret_26_30-sum(ret_26_30)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(4)  = sqrt( sum( (ret_31_35-sum(ret_31_35)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(5)  = sqrt( sum( (ret_36_40-sum(ret_36_40)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(6)  = sqrt( sum( (ret_41_45-sum(ret_41_45)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(7)  = sqrt( sum( (ret_46_50-sum(ret_46_50)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(8)  = sqrt( sum( (ret_51_55-sum(ret_51_55)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(9)  = sqrt( sum( (ret_56_60-sum(ret_56_60)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(10) = sqrt( sum( (ret_61_65-sum(ret_61_65)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_Age(11) = sqrt( sum( (ret_66_70-sum(ret_66_70)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )

			Std_Dev_Return_W_Age(1)  = sqrt( sum( (ret_w_20   -sum(ret_w_20)   /ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(2)  = sqrt( sum( (ret_w_21_25-sum(ret_w_21_25)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(3)  = sqrt( sum( (ret_w_26_30-sum(ret_w_26_30)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(4)  = sqrt( sum( (ret_w_31_35-sum(ret_w_31_35)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(5)  = sqrt( sum( (ret_w_36_40-sum(ret_w_36_40)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(6)  = sqrt( sum( (ret_w_41_45-sum(ret_w_41_45)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(7)  = sqrt( sum( (ret_w_46_50-sum(ret_w_46_50)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(8)  = sqrt( sum( (ret_w_51_55-sum(ret_w_51_55)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(9)  = sqrt( sum( (ret_w_56_60-sum(ret_w_56_60)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(10) = sqrt( sum( (ret_w_61_65-sum(ret_w_61_65)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )
			Std_Dev_Return_W_Age(11) = sqrt( sum( (ret_w_66_70-sum(ret_w_66_70)/ret_size)**2.0_sp ) /real(ret_size-1,SP)  )

			Std_Dev_Return_K_Age(1)  = sqrt( sum( (ret_k_20   -sum(ret_k_20   , Ind_K_20   ==1)/sum(Ind_K_20   )**2.0_sp) , Ind_K_20   ==1) & 
										&	/real(sum(Ind_K_20)   -1,SP)  )
			Std_Dev_Return_K_Age(2)  = sqrt( sum( (ret_k_21_25-sum(ret_k_21_25, Ind_K_21_25==1)/sum(Ind_K_21_25)**2.0_sp) , Ind_K_21_25==1) & 
										&	/real(sum(Ind_K_21_25)-1,SP)  )
			Std_Dev_Return_K_Age(3)  = sqrt( sum( (ret_k_26_30-sum(ret_k_26_30, Ind_K_26_30==1)/sum(Ind_K_26_30)**2.0_sp) , Ind_K_26_30==1) & 
										&	/real(sum(Ind_K_26_30)-1,SP)  )
			Std_Dev_Return_K_Age(4)  = sqrt( sum( (ret_k_31_35-sum(ret_k_31_35, Ind_K_31_35==1)/sum(Ind_K_31_35)**2.0_sp) , Ind_K_31_35==1) & 
										&	/real(sum(Ind_K_31_35)-1,SP)  )
			Std_Dev_Return_K_Age(5)  = sqrt( sum( (ret_k_36_40-sum(ret_k_36_40, Ind_K_36_40==1)/sum(Ind_K_36_40)**2.0_sp) , Ind_K_36_40==1) & 
										&	/real(sum(Ind_K_36_40)-1,SP)  )
			Std_Dev_Return_K_Age(6)  = sqrt( sum( (ret_k_41_45-sum(ret_k_41_45, Ind_K_41_45==1)/sum(Ind_K_41_45)**2.0_sp) , Ind_K_41_45==1) & 
										&	/real(sum(Ind_K_41_45)-1,SP)  )
			Std_Dev_Return_K_Age(7)  = sqrt( sum( (ret_k_46_50-sum(ret_k_46_50, Ind_K_46_50==1)/sum(Ind_K_46_50)**2.0_sp) , Ind_K_46_50==1) & 
										&	/real(sum(Ind_K_46_50)-1,SP)  )
			Std_Dev_Return_K_Age(8)  = sqrt( sum( (ret_k_51_55-sum(ret_k_51_55, Ind_K_51_55==1)/sum(Ind_K_51_55)**2.0_sp) , Ind_K_51_55==1) & 
										&	/real(sum(Ind_K_51_55)-1,SP)  )
			Std_Dev_Return_K_Age(9)  = sqrt( sum( (ret_k_56_60-sum(ret_k_56_60, Ind_K_56_60==1)/sum(Ind_K_56_60)**2.0_sp) , Ind_K_56_60==1) & 
										&	/real(sum(Ind_K_56_60)-1,SP)  )
			Std_Dev_Return_K_Age(10) = sqrt( sum( (ret_k_61_65-sum(ret_k_61_65, Ind_K_61_65==1)/sum(Ind_K_61_65)**2.0_sp) , Ind_K_61_65==1) & 
										&	/real(sum(Ind_K_61_65)-1,SP)  )
			Std_Dev_Return_K_Age(11) = sqrt( sum( (ret_k_66_70-sum(ret_k_66_70, Ind_K_66_70==1)/sum(Ind_K_66_70)**2.0_sp) , Ind_K_66_70==1) & 
										&	/real(sum(Ind_K_66_70)-1,SP)  )

			Std_Dev_Return_K_W_Age(1)  = sqrt( sum( (ret_k_w_20   -sum(ret_k_w_20   , Ind_K_20   ==1)/sum(Ind_K_20   )**2.0_sp) , &
										& Ind_K_20   ==1) / real(sum(Ind_K_20)   -1,SP)  )
			Std_Dev_Return_K_W_Age(2)  = sqrt( sum( (ret_k_w_21_25-sum(ret_k_w_21_25, Ind_K_21_25==1)/sum(Ind_K_21_25)**2.0_sp) , &
										& Ind_K_21_25==1) / real(sum(Ind_K_21_25)-1,SP)  )
			Std_Dev_Return_K_W_Age(3)  = sqrt( sum( (ret_k_w_26_30-sum(ret_k_w_26_30, Ind_K_26_30==1)/sum(Ind_K_26_30)**2.0_sp) , &
										& Ind_K_26_30==1) / real(sum(Ind_K_26_30)-1,SP)  )
			Std_Dev_Return_K_W_Age(4)  = sqrt( sum( (ret_k_w_31_35-sum(ret_k_w_31_35, Ind_K_31_35==1)/sum(Ind_K_31_35)**2.0_sp) , &
										& Ind_K_31_35==1) / real(sum(Ind_K_31_35)-1,SP)  )
			Std_Dev_Return_K_W_Age(5)  = sqrt( sum( (ret_k_w_36_40-sum(ret_k_w_36_40, Ind_K_36_40==1)/sum(Ind_K_36_40)**2.0_sp) , &
										& Ind_K_36_40==1) / real(sum(Ind_K_36_40)-1,SP)  )
			Std_Dev_Return_K_W_Age(6)  = sqrt( sum( (ret_k_w_41_45-sum(ret_k_w_41_45, Ind_K_41_45==1)/sum(Ind_K_41_45)**2.0_sp) , &
										& Ind_K_41_45==1) / real(sum(Ind_K_41_45)-1,SP)  )
			Std_Dev_Return_K_W_Age(7)  = sqrt( sum( (ret_k_w_46_50-sum(ret_k_w_46_50, Ind_K_46_50==1)/sum(Ind_K_46_50)**2.0_sp) , &
										& Ind_K_46_50==1) / real(sum(Ind_K_46_50)-1,SP)  )
			Std_Dev_Return_K_W_Age(8)  = sqrt( sum( (ret_k_w_51_55-sum(ret_k_w_51_55, Ind_K_51_55==1)/sum(Ind_K_51_55)**2.0_sp) , &
										& Ind_K_51_55==1) / real(sum(Ind_K_51_55)-1,SP)  )
			Std_Dev_Return_K_W_Age(9)  = sqrt( sum( (ret_k_w_56_60-sum(ret_k_w_56_60, Ind_K_56_60==1)/sum(Ind_K_56_60)**2.0_sp) , &
										& Ind_K_56_60==1) / real(sum(Ind_K_56_60)-1,SP)  )
			Std_Dev_Return_K_W_Age(10) = sqrt( sum( (ret_k_w_61_65-sum(ret_k_w_61_65, Ind_K_61_65==1)/sum(Ind_K_61_65)**2.0_sp) , &
										& Ind_K_61_65==1) / real(sum(Ind_K_61_65)-1,SP)  )
			Std_Dev_Return_K_W_Age(11) = sqrt( sum( (ret_k_w_66_70-sum(ret_k_w_66_70, Ind_K_66_70==1)/sum(Ind_K_66_70)**2.0_sp) , &
										& Ind_K_66_70==1) / real(sum(Ind_K_66_70)-1,SP)  )
			print*, 'End of std dev of return by age'

			print*, ' '
			print*, 'Start of mean of return by age'
			! Mean of return by age
			Mean_Return_Age(1)  = sum(ret_20)   /ret_size
			Mean_Return_Age(2)  = sum(ret_21_25)/ret_size
			Mean_Return_Age(3)  = sum(ret_26_30)/ret_size
			Mean_Return_Age(4)  = sum(ret_31_35)/ret_size
			Mean_Return_Age(5)  = sum(ret_36_40)/ret_size
			Mean_Return_Age(6)  = sum(ret_41_45)/ret_size
			Mean_Return_Age(7)  = sum(ret_46_50)/ret_size
			Mean_Return_Age(8)  = sum(ret_51_55)/ret_size
			Mean_Return_Age(9)  = sum(ret_56_60)/ret_size
			Mean_Return_Age(10) = sum(ret_61_65)/ret_size
			Mean_Return_Age(11) = sum(ret_66_70)/ret_size

			Mean_Return_W_Age(1)  = sum(ret_w_20)   /ret_size
			Mean_Return_W_Age(2)  = sum(ret_w_21_25)/ret_size
			Mean_Return_W_Age(3)  = sum(ret_w_26_30)/ret_size
			Mean_Return_W_Age(4)  = sum(ret_w_31_35)/ret_size
			Mean_Return_W_Age(5)  = sum(ret_w_36_40)/ret_size
			Mean_Return_W_Age(6)  = sum(ret_w_41_45)/ret_size
			Mean_Return_W_Age(7)  = sum(ret_w_46_50)/ret_size
			Mean_Return_W_Age(8)  = sum(ret_w_51_55)/ret_size
			Mean_Return_W_Age(9)  = sum(ret_w_56_60)/ret_size
			Mean_Return_W_Age(10) = sum(ret_w_61_65)/ret_size
			Mean_Return_W_Age(11) = sum(ret_w_66_70)/ret_size

			Mean_Return_K_Age(1)  = sum(ret_k_20   , Ind_K_20   ==1)/sum(Ind_K_20   )
			Mean_Return_K_Age(2)  = sum(ret_k_21_25, Ind_K_21_25==1)/sum(Ind_K_21_25)
			Mean_Return_K_Age(3)  = sum(ret_k_26_30, Ind_K_26_30==1)/sum(Ind_K_26_30)
			Mean_Return_K_Age(4)  = sum(ret_k_31_35, Ind_K_31_35==1)/sum(Ind_K_31_35)
			Mean_Return_K_Age(5)  = sum(ret_k_36_40, Ind_K_36_40==1)/sum(Ind_K_36_40)
			Mean_Return_K_Age(6)  = sum(ret_k_41_45, Ind_K_41_45==1)/sum(Ind_K_41_45)
			Mean_Return_K_Age(7)  = sum(ret_k_46_50, Ind_K_46_50==1)/sum(Ind_K_46_50)
			Mean_Return_K_Age(8)  = sum(ret_k_51_55, Ind_K_51_55==1)/sum(Ind_K_51_55)
			Mean_Return_K_Age(9)  = sum(ret_k_56_60, Ind_K_56_60==1)/sum(Ind_K_56_60)
			Mean_Return_K_Age(10) = sum(ret_k_61_65, Ind_K_61_65==1)/sum(Ind_K_61_65)
			Mean_Return_K_Age(11) = sum(ret_k_66_70, Ind_K_66_70==1)/sum(Ind_K_66_70)

			Mean_Return_K_W_Age(1)  = sum(ret_k_w_20   , Ind_K_20   ==1)/sum(Ind_K_20   )
			Mean_Return_K_W_Age(2)  = sum(ret_k_w_21_25, Ind_K_21_25==1)/sum(Ind_K_21_25)
			Mean_Return_K_W_Age(3)  = sum(ret_k_w_26_30, Ind_K_26_30==1)/sum(Ind_K_26_30)
			Mean_Return_K_W_Age(4)  = sum(ret_k_w_31_35, Ind_K_31_35==1)/sum(Ind_K_31_35)
			Mean_Return_K_W_Age(5)  = sum(ret_k_w_36_40, Ind_K_36_40==1)/sum(Ind_K_36_40)
			Mean_Return_K_W_Age(6)  = sum(ret_k_w_41_45, Ind_K_41_45==1)/sum(Ind_K_41_45)
			Mean_Return_K_W_Age(7)  = sum(ret_k_w_46_50, Ind_K_46_50==1)/sum(Ind_K_46_50)
			Mean_Return_K_W_Age(8)  = sum(ret_k_w_51_55, Ind_K_51_55==1)/sum(Ind_K_51_55)
			Mean_Return_K_W_Age(9)  = sum(ret_k_w_56_60, Ind_K_56_60==1)/sum(Ind_K_56_60)
			Mean_Return_K_W_Age(10) = sum(ret_k_w_61_65, Ind_K_61_65==1)/sum(Ind_K_61_65)
			Mean_Return_K_W_Age(11) = sum(ret_k_w_66_70, Ind_K_66_70==1)/sum(Ind_K_66_70)
			print*, 'End of mean of return by age'

			print*, ' '
			print*, 'Start of pcrt of return by age'

			! Percentiles of return by age
			prctile_ret = (/0.999_sp, 0.99_sp, 0.95_sp, 0.90_sp, 0.75_sp, 0.50_sp, 0.25_sp, 0.10_sp, 0.01_sp/)

			do i_pct=1,9 
				print*, 'Ret prc=', prctile_ret(i_pct)
				prc_Return_Age(1 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_20)
				prc_Return_Age(2 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_21_25)
				prc_Return_Age(3 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_26_30)
				prc_Return_Age(4 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_31_35)
				prc_Return_Age(5 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_36_40)
				prc_Return_Age(6 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_41_45)
				prc_Return_Age(7 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_46_50)
				prc_Return_Age(8 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_51_55)
				prc_Return_Age(9 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_56_60)
				prc_Return_Age(10,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_61_65)
				prc_Return_Age(11,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_66_70)
				print*, 'Ret W prc=', prctile_ret(i_pct)
				prc_Return_W_Age(1 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_20)
				prc_Return_W_Age(2 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_21_25)
				prc_Return_W_Age(3 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_26_30)
				prc_Return_W_Age(4 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_31_35)
				prc_Return_W_Age(5 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_36_40)
				prc_Return_W_Age(6 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_41_45)
				prc_Return_W_Age(7 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_46_50)
				prc_Return_W_Age(8 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_51_55)
				prc_Return_W_Age(9 ,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_56_60)
				prc_Return_W_Age(10,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_61_65)
				prc_Return_W_Age(11,i_pct) = Percentile(prctile_ret(i_pct),ret_size,ret_w_66_70)
				print*, 'Ret K prc=', prctile_ret(i_pct)
				prc_Return_K_Age(1 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_20   ),pack(ret_20   ,Ind_K_20   ==1))
				prc_Return_K_Age(2 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_21_25),pack(ret_21_25,Ind_K_21_25==1))
				prc_Return_K_Age(3 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_26_30),pack(ret_26_30,Ind_K_26_30==1))
				prc_Return_K_Age(4 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_31_35),pack(ret_31_35,Ind_K_31_35==1))
				prc_Return_K_Age(5 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_36_40),pack(ret_36_40,Ind_K_36_40==1))
				prc_Return_K_Age(6 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_41_45),pack(ret_41_45,Ind_K_41_45==1))
				prc_Return_K_Age(7 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_46_50),pack(ret_46_50,Ind_K_46_50==1))
				prc_Return_K_Age(8 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_51_55),pack(ret_51_55,Ind_K_51_55==1))
				prc_Return_K_Age(9 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_56_60),pack(ret_56_60,Ind_K_56_60==1))
				prc_Return_K_Age(10,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_61_65),pack(ret_61_65,Ind_K_61_65==1))
				prc_Return_K_Age(11,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_66_70),pack(ret_66_70,Ind_K_66_70==1))
				print*, 'Ret K W prc=', prctile_ret(i_pct)
				prc_Return_K_W_Age(1 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_20   ),pack(ret_w_20   ,Ind_K_20   ==1))
				prc_Return_K_W_Age(2 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_21_25),pack(ret_w_21_25,Ind_K_21_25==1))
				prc_Return_K_W_Age(3 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_26_30),pack(ret_w_26_30,Ind_K_26_30==1))
				prc_Return_K_W_Age(4 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_31_35),pack(ret_w_31_35,Ind_K_31_35==1))
				prc_Return_K_W_Age(5 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_36_40),pack(ret_w_36_40,Ind_K_36_40==1))
				prc_Return_K_W_Age(6 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_41_45),pack(ret_w_41_45,Ind_K_41_45==1))
				prc_Return_K_W_Age(7 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_46_50),pack(ret_w_46_50,Ind_K_46_50==1))
				prc_Return_K_W_Age(8 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_51_55),pack(ret_w_51_55,Ind_K_51_55==1))
				prc_Return_K_W_Age(9 ,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_56_60),pack(ret_w_56_60,Ind_K_56_60==1))
				prc_Return_K_W_Age(10,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_61_65),pack(ret_w_61_65,Ind_K_61_65==1))
				prc_Return_K_W_Age(11,i_pct) = Percentile(prctile_ret(i_pct),sum(Ind_K_66_70),pack(ret_w_66_70,Ind_K_66_70==1))
			enddo 
			print*, 'End of prc of return by age'


			if (bench_indx.eq.1) then
			OPEN(UNIT=10, FILE=trim(Result_Folder)//'Simul/Return_Stats_by_Age_bench.txt', STATUS='replace')
			else
			OPEN(UNIT=10, FILE=trim(Result_Folder)//'Simul/Return_Stats_by_Age_exp.txt'  , STATUS='replace')
			endif 

			WRITE(UNIT=10, FMT=*) ' '
			WRITE(UNIT=10, FMT=*) 'Std Dev of Return by Age'
			WRITE(UNIT=10, FMT=*) '20 ','21-25 ','26-30 ','31-35 ','36-40 ','41-45 ','46-50 ','51-55 ','56-60 ','61-65 ','66-70 '
			WRITE(UNIT=10, FMT=*) Std_Dev_Return_Age
			WRITE(UNIT=10, FMT=*) Mean_Return_Age
			do i_pct=1,9 
				WRITE(UNIT=10, FMT=*) prc_Return_Age(:,i_pct)
			enddo 
			WRITE(UNIT=10, FMT=*) ' '
			WRITE(UNIT=10, FMT=*) 'Std Dev of Return W by Age'
			WRITE(UNIT=10, FMT=*) ' 20 ','21-25 ','26-30 ','31-35 ','36-40 ','41-45 ','46-50 ','51-55 ','56-60 ','61-65 ','66-70 '
			WRITE(UNIT=10, FMT=*) Std_Dev_Return_W_Age
			WRITE(UNIT=10, FMT=*) Mean_Return_W_Age
			do i_pct=1,9 
				WRITE(UNIT=10, FMT=*) prc_Return_W_Age(:,i_pct)
			enddo 
			WRITE(UNIT=10, FMT=*) ' '
			WRITE(UNIT=10, FMT=*) 'Std Dev of Return K by Age'
			WRITE(UNIT=10, FMT=*) ' 20 ','21-25 ','26-30 ','31-35 ','36-40 ','41-45 ','46-50 ','51-55 ','56-60 ','61-65 ','66-70 '
			WRITE(UNIT=10, FMT=*) Std_Dev_Return_K_Age
			WRITE(UNIT=10, FMT=*) Mean_Return_K_Age
			do i_pct=1,9 
				WRITE(UNIT=10, FMT=*) prc_Return_K_Age(:,i_pct)
			enddo 
			WRITE(UNIT=10, FMT=*) ' '
			WRITE(UNIT=10, FMT=*) 'Std Dev of Return K W by Age'
			WRITE(UNIT=10, FMT=*) ' 20 ','21-25 ','26-30 ','31-35 ','36-40 ','41-45 ','46-50 ','51-55 ','56-60 ','61-65 ','66-70 '
			WRITE(UNIT=10, FMT=*) Std_Dev_Return_K_W_Age
			WRITE(UNIT=10, FMT=*) Mean_Return_K_W_Age
			do i_pct=1,9 
				WRITE(UNIT=10, FMT=*) prc_Return_K_W_Age(:,i_pct)
			enddo 
			WRITE(UNIT=10, FMT=*) ' '

			CLOSE(UNIT=10)

		!$omp parallel do private(currenta,age,currentzi,currentlambdai,currentei,tklo,tkhi,h_i)
		DO paneli=1,totpop
		    currenta 		= panela(paneli)
		    age 			= panelage(paneli)
		    currentzi 		= panelz(paneli)
		    currentlambdai 	= panellambda(paneli) 
		    currentei 		= panele(paneli)
		    currentxi 		= panelx(paneli)

		    if (currenta .ge. amax) then
		        tklo = na-1
		    elseif (currenta .lt. amin) then
	            tklo = 1
	        else
	            tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
		    endif    
		    tkhi = tklo + 1   

		    panelK(paneli)    = min(theta(currentzi)*currenta,(mu*P*xz_grid(currentxi,currentzi)**mu/(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )

			if (age.lt.RetAge) then 
		        h_i  = ((agrid(tkhi) - currenta)*hours(age,tklo,currentzi,currentlambdai,currentei,currentxi) &
		           &  + (currenta - agrid(tklo))*hours(age,tkhi,currentzi,currentlambdai,currentei,currentxi) ) &
		                                &  / ( agrid(tkhi) - agrid(tklo) )  

				panel_Y_L(paneli) = psi*( Wage*eff_un(age,currentlambdai,currentei)*h_i)**(1.0_dp-tauPL)
			else 
				panel_Y_L(paneli) = RetY_lambda_e(currentlambdai,currentei)
			endif 

			if (age.eq.1) then 
				currentzi = panelz_old(paneli)
				currentlambdai = panellambda_old(paneli)
				panelPV_a(paneli) = (   (agrid(tkhi) - currenta) * V_Pr_nb(tklo,currentzi,currentlambdai)  &
		                       &  + (currenta - agrid(tklo)) * V_Pr_nb(tkhi,currentzi,currentlambdai)) &
		                       &  / ( agrid(tkhi) - agrid(tklo) )
		    else 
		    	
		    	panelPV_a(paneli) = (   (agrid(tkhi) - currenta) * V_Pr(age,tklo,currentzi,currentlambdai, currentei, currentxi)  &
		                       &  + (currenta - agrid(tklo)) * V_Pr(age,tkhi,currentzi,currentlambdai, currentei, currentxi)) &
		                       &  / ( agrid(tkhi) - agrid(tklo) )  + (1.0_dp+R)*currenta 
		    endif 
		           
		ENDDO ! paneli


		print*, ' '
		print*, 'Writing simulation results'

		if (bench_indx.eq.1) then
			OPEN(UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_bench'		, STATUS='replace')
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_bench'	  	, STATUS='replace')
			OPEN(UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_bench'		, STATUS='replace')
			OPEN(UNIT=13, FILE=trim(Result_Folder)//'Simul/panellambda_bench'   , STATUS='replace')
			OPEN(UNIT=14, FILE=trim(Result_Folder)//'Simul/panele_bench'        , STATUS='replace')
			OPEN(UNIT=26, FILE=trim(Result_Folder)//'Simul/panelPV_a_bench'     , STATUS='replace')
			OPEN(UNIT=27, FILE=trim(Result_Folder)//'Simul/panelK_bench'        , STATUS='replace')
			OPEN(UNIT=28, FILE=trim(Result_Folder)//'Simul/panelx_bench'        , STATUS='replace')
			OPEN(UNIT=24, FILE=trim(Result_Folder)//'Simul/panel_YL_bench'    	, STATUS='replace')
		else 
			OPEN(UNIT=10, FILE=trim(Result_Folder)//'Simul/panela_exp'		 	, STATUS='replace')
			OPEN(UNIT=11, FILE=trim(Result_Folder)//'Simul/panelage_exp'		, STATUS='replace')
			OPEN(UNIT=12, FILE=trim(Result_Folder)//'Simul/panelz_exp'		 	, STATUS='replace')
			OPEN(UNIT=13, FILE=trim(Result_Folder)//'Simul/panellambda_exp'   	, STATUS='replace')
			OPEN(UNIT=14, FILE=trim(Result_Folder)//'Simul/panele_exp'        	, STATUS='replace')
			OPEN(UNIT=26, FILE=trim(Result_Folder)//'Simul/panelPV_a_exp'       , STATUS='replace')
			OPEN(UNIT=27, FILE=trim(Result_Folder)//'Simul/panelK_exp' 	        , STATUS='replace')
			OPEN(UNIT=28, FILE=trim(Result_Folder)//'Simul/panelx_exp'	        , STATUS='replace')
			OPEN(UNIT=24, FILE=trim(Result_Folder)//'Simul/panel_YL_exp'    	, STATUS='replace')
		endif 


		WRITE  (UNIT=10, FMT='(F12.4)') panela
		WRITE  (UNIT=11, FMT=*) panelage 
		WRITE  (UNIT=12, FMT=*) panelz 
		WRITE  (UNIT=13, FMT=*) panellambda 
		WRITE  (UNIT=14, FMT=*) panele 
		WRITE  (UNIT=26, FMT='(F12.4)') panelPV_a
		WRITE  (UNIT=27, FMT='(F12.4)') panelK
		WRITE  (UNIT=28, FMT=*) panelx
		WRITE  (UNIT=24, FMT='(F12.4)') panel_Y_L

		close (unit=10); close (unit=11); close (unit=12); close (unit=13); close (unit=14)
		close (unit=26); close (unit=27); close (unit=28); close (unit=24) 

		if (bench_indx==1) then
			OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/panela_parents' 	, STATUS='replace')
			OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/panela_sons'    	, STATUS='replace')
			OPEN(UNIT=22, FILE=trim(Result_Folder)//'Simul/panelage_parents' 	, STATUS='replace')
			OPEN(UNIT=23, FILE=trim(Result_Folder)//'Simul/panelage_sons'    	, STATUS='replace')
			WRITE (UNIT=20, FMT='(F12.4)') eligible_panela_parents
			WRITE (UNIT=21, FMT='(F12.4)') eligible_panela_sons
			WRITE (UNIT=22, FMT=*) eligible_panelage_parents
			WRITE (UNIT=23, FMT=*) eligible_panelage_sons
			close (unit=20); close (unit=21); close (unit=22); close (unit=23)

			call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/IGM_3050' )
			OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_3050/panela_parents' , STATUS='replace')
			OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_3050/panela_sons'     , STATUS='replace')
			WRITE (UNIT=20, FMT='(F12.4)') panela_dad
			WRITE (UNIT=21, FMT='(F12.4)') panela_son
			close (unit=20); close (unit=21); 
			OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_3050/panelr_parents' , STATUS='replace')
			OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_3050/panelr_sons'     , STATUS='replace')
			WRITE (UNIT=20, FMT='(F12.4)') panelr_dad
			WRITE (UNIT=21, FMT='(F12.4)') panelr_son
			close (unit=20); close (unit=21); 
			OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_3050/panelpv_parents' , STATUS='replace')
			OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_3050/panelpv_sons'     , STATUS='replace')
			WRITE (UNIT=20, FMT='(F12.4)') panelpv_dad
			WRITE (UNIT=21, FMT='(F12.4)') panelpv_son
			close (unit=20); close (unit=21); 
			OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_3050/panelz_parents' , STATUS='replace')
			OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_3050/panelz_sons'     , STATUS='replace')
			WRITE (UNIT=20, FMT=*) panelz_dad
			WRITE (UNIT=21, FMT=*) panelz_son
			close (unit=20); close (unit=21); 

			! call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/IGM_4060' )
			! OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_4060/panela_parents' , STATUS='replace')
			! OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_4060/panela_sons'     , STATUS='replace')
			! WRITE (UNIT=20, FMT='(F12.4)') panela_dad_2
			! WRITE (UNIT=21, FMT='(F12.4)') panela_son_2
			! close (unit=20); close (unit=21); 
			! OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_4060/panelr_parents' , STATUS='replace')
			! OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_4060/panelr_sons'     , STATUS='replace')
			! WRITE (UNIT=20, FMT='(F12.4)') panelr_dad
			! WRITE (UNIT=21, FMT='(F12.4)') panelr_son
			! close (unit=20); close (unit=21); 
			! OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_4060/panelpv_parents' , STATUS='replace')
			! OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_4060/panelpv_sons'     , STATUS='replace')
			! WRITE (UNIT=20, FMT='(F12.4)') panelpv_dad
			! WRITE (UNIT=21, FMT='(F12.4)') panelpv_son
			! close (unit=20); close (unit=21); 
			! OPEN(UNIT=20, FILE=trim(Result_Folder)//'Simul/IGM_4060/panelz_parents' , STATUS='replace')
			! OPEN(UNIT=21, FILE=trim(Result_Folder)//'Simul/IGM_4060/panelz_sons'     , STATUS='replace')
			! WRITE (UNIT=20, FMT=*) panelz_dad
			! WRITE (UNIT=21, FMT=*) panelz_son
			! close (unit=20); close (unit=21); 
			
		endif


		if (bench_indx==1) then
		print*, 'Identifying top agents'
		! Top Agents - Assets
			n_top = totpop 
			A_hi  = maxval(panela)
			A_low = minval(panela)
			A_cut = (A_hi+A_low)/2.0_dp 
			do while (n_top.ne.80)
				n_top = count(panela.ge.A_cut)
				if (n_top.gt.80) then
					A_low = A_cut 
				endif 
				if (n_top.lt.80) then
					A_hi  = A_cut 
				endif 
				A_cut = (A_hi+A_low)/2.0_dp 
				print*, A_cut, n_top, count(panela.ge.A_cut) 
			enddo 
			print*,' '
			print*,'A_cut final:'
			print*, A_cut, n_top, count(panela.ge.A_cut) 

			panel_top_ind = (/(paneli, paneli=1,totpop, 1)/)
			top_ind = pack(panel_top_ind  , (panela.ge.A_cut) )
			top_A   = pack(panela         , (panela.gt.A_cut) )
			! Sort by assets
			call Sort(80,top_A,top_A,top_ind_aux)
			top_ind = top_ind(top_ind_aux)
			! Test print
			print*,' '
			print*,'Top Agents and Top Assets'
			do paneli=1,80
				print*, top_ind(paneli),top_A(paneli)
			enddo
			print*,'Max(A)', maxval(panela)

			top_folder = 'Top_A/'
			call SIMULATION_TOP(bench_indx,top_ind,top_folder)

		! Top Agents - PV
			n_top = totpop 
			A_hi  = maxval(panelPV_a)
			A_low = minval(panelPV_a)
			A_cut = (A_hi+A_low)/2.0_dp 
			do while (n_top.ne.80)
				n_top = count(panelPV_a.ge.A_cut)
				if (n_top.gt.80) then
					A_low = A_cut 
				endif 
				if (n_top.lt.80) then
					A_hi  = A_cut 
				endif 
				A_cut = (A_hi+A_low)/2.0_dp 
				print*, A_cut, n_top, count(panelPV_a.ge.A_cut) 
			enddo 
			print*,' '
			print*,'A_cut final:'
			print*, A_cut, n_top, count(panelPV_a.ge.A_cut) 

			panel_top_ind = (/(paneli, paneli=1,totpop, 1)/)
			top_ind = pack(panel_top_ind  , (panelPV_a.ge.A_cut) )
			top_A   = pack(panelPV_a      , (panelPV_a.gt.A_cut) )
			! Sort by assets
			call Sort(80,top_A,top_A,top_ind_aux)
			top_ind = top_ind(top_ind_aux)
			! Test print
			print*,' '
			print*,'Top Agents and Top PV'
			do paneli=1,80
				print*, top_ind(paneli),top_A(paneli)
			enddo
			print*,'Max(PV)', maxval(panelPV_a)

			top_folder = 'Top_PV/'
			call SIMULATION_TOP(bench_indx,top_ind,top_folder)

		endif 



END SUBROUTINE SIMULATION


!========================================================================================
!========================================================================================
!========================================================================================

SUBROUTINE  SIMULATION_TOP(bench_indx,top_ind,folder)
	use parameters
	use global
	use omp_lib
	IMPLICIT NONE
	integer      , intent(in) :: bench_indx, top_ind(80)
	character(10), intent(in) :: folder
	integer  :: currentzi, currentlambdai, currentei, currentxi
	REAL(DP) :: tempnoage, tempnoz, tempnolambda, tempnoe, tempnox, tempno, currenta, currentY
	REAL(DP) :: start_timet, finish_timet, h_i
	INTEGER  :: agecounter, agesign, tage, tzi, tlambdai, tei, tklo, tkhi, paneli, simutime
	INTEGER , DIMENSION(MaxAge) :: requirednumberby_age, cdfrequirednumberby_age
	INTEGER , DIMENSION(totpop) :: panelage, panelz, panellambda, panele, panelx, panelz_old, panellambda_old 
	REAL(SP), DIMENSION(totpop) :: panela, panelPV_a, panelK 
	INTEGER , DIMENSION(150,80) :: panelage_top, panelz_top, panelx_top, panel_lambda_top, panele_top
	REAL(SP), DIMENSION(150,80) :: panela_top, panelK_top, panel_YL_top, panel_PV_top
	REAL(SP), DIMENSION(150,80) :: prc_all_top, prc_cohort_top, prc_PV_all_top, prc_PV_cohort_top
	INTEGER 				    :: ii, age_top, thread_num

	!$ call omp_set_num_threads(20)
	!$OMP parallel
    !$OMP critical
		thread_num = omp_get_thread_num()
		newiseed   = -5 - thread_num
		tempno     = omp_ran1(newiseed)
       	! write(*,'("inside critical myt=",i4,f12.8)') thread_num,tempno,newiseed
    !$OMP end critical
    !$OMP end parallel


	print*,'SIMULATION STARTED (for top agents)'

    age=1
    requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
    cdfrequirednumberby_age(age) = requirednumberby_age(age)
    DO age=2,MaxAge
        requirednumberby_age(age)    = NINT(totpop*pop(age)/sum(pop))
        cdfrequirednumberby_age(age) = requirednumberby_age(age) + cdfrequirednumberby_age(age-1)
    ENDDO
    ! If the total number of people are not equal to the total population, then I will add the remainder to the last age
    requirednumberby_age(MaxAge) =  requirednumberby_age(MaxAge)-cdfrequirednumberby_age(MaxAge) + totpop
    cdfrequirednumberby_age(MaxAge) = totpop

	!=====================================================================
	!                     GENERATE   INITIAL   PANEL
	!=====================================================================


	newiseed=-1

	!$omp parallel do private(tempnoage,age,tempnoz,zi,xi,tempnolambda,lambdai,tempnoe,ei)
	DO paneli=1,totpop

	! AGE
	   tempnoage = omp_ran1() ! ran1(newiseed)
	   age=1
	   DO WHILE (tempnoage*totpop .gt. cdfrequirednumberby_age(age))
	       age=age+1
	   ENDDO

	! Z   
	   tempnoz = omp_ran1() ! ran1(newiseed)
	   zi=1
	   DO WHILE (tempnoz .gt. cdf_Gz(zi))
	       zi=zi+1
	   ENDDO


	! X   
	   xi=1
	 
	! LAMBDA  
	   tempnolambda = omp_ran1() ! ran1(newiseed) 
	   lambdai=1
	   DO WHILE (tempnolambda .gt. cdf_Glambda(lambdai))
	       lambdai=lambdai+1
	   ENDDO

	! E   
	   tempnoe = omp_ran1() ! ran1(newiseed)   
	   ei=1
	   DO WHILE (tempnoe .gt. cdf_Ge_byage(age,ei))
	       ei=ei+1
	   ENDDO

	   panelage(paneli)		= age
	   panelz(paneli)		= zi
	   panellambda(paneli)	= lambdai
	   panele(paneli)		= ei
	   panelx(paneli)		= xi
	   
	ENDDO
	print*, sum(panelage)/real(totpop,8), sum(panelz)/real(totpop,8), sum(panele)/real(totpop,8)

	print*, ' Initial states ready'

	! SET INITIAL ASSET DISTRIBUTION
	panela            = 1.0_DP

	!=============================================================================
	!
	! SIMULATE FROM THE SECOND PERIOD SHOCKS AND UPDATE NEW DISTRIBUTIONS
	!
	!=============================================================================

	!call cpu_time(start_timet) 

	DO simutime=1, MaxSimuTime

		! Record satates of father
		panelz_old = panelz 
		panellambda_old = panellambda 

		!$omp parallel do private(tempnoage,age,tempnoz,zi,xi,tempnolambda,lambdai,tempnoe,ei, &
		!$omp& currenta,currentzi,currentlambdai,currentei,currentxi,tklo,tkhi,tempno)
	    DO paneli=1,totpop
	    
	       	currenta  		= panela(paneli)
	       	age       		= panelage(paneli)
	       	currentzi      	= panelz(paneli)
	       	currentlambdai 	= panellambda(paneli) 
	       	currentei 		= panele(paneli)
	       	currentxi 		= panelx(paneli)
	        
			!  COMPUTE NEXT PERIOD'S ASSET
	       	if (age .lt. MaxAge) then
	            ! do linear interpolation here to find aprime. calling the function takes much more time
	            if (currenta .ge. amax) then
                	tklo = na-1
                elseif (currenta .lt. amin) then
                    tklo = 1
                else
                    tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
	            endif 
	            tkhi = tklo + 1        

	            panela(paneli) = ((agrid(tkhi) - currenta)*Aprime(age,tklo,currentzi,currentlambdai,currentei,currentxi) 		&
	                           	&  + (currenta - agrid(tklo))*Aprime(age,tkhi,currentzi,currentlambdai, currentei,currentxi)) &
	                            &  / ( agrid(tkhi) - agrid(tklo) )    
	            
	            if (panela(paneli)  .ge. amax) then
	                panela(paneli) = min(panela(paneli), amax) 
	            endif      
	            if (panela(paneli)  .lt. amin) then
	                panela(paneli) = max(panela(paneli), amin) 
	            endif      

	       	endif !age .lt. MaxAge
			!  NEXT PERIOD'S ASSET IS COMPUTED

	             
			! DRAW NEXT PERIOD'S AGE DBN
	      	tempnoage = omp_ran1() ! ran1(newiseed)  
	  
	      	IF (tempnoage .gt. survP(age)) THEN
				panelage(paneli)   = 1
			ELSE
				panelage(paneli)   = age+1
				panelz(paneli)     = currentzi
				panellambda(paneli)= currentlambdai   
	      	ENDIF
	    
			! DRAW Z and LAMBDA DISTRIBUTION FOR ONE-YEAR OLDS
	     	age = panelage(paneli)   
	     	IF (age .eq. 1) THEN    
				! Z      
		       	tempnoz = omp_ran1() ! ran1(newiseed) 
		       	zi=1
				DO WHILE (tempnoz .gt. cdf_pr_z(currentzi,zi))
					zi=zi+1
				ENDDO

				! X     
		       	xi = 1
	       
				! LAMBDA  
				tempnolambda = omp_ran1() ! ran1(newiseed) 
				lambdai=1
				DO WHILE (tempnolambda .gt. cdf_pr_lambda(currentlambdai,lambdai))
					lambdai=lambdai+1
				ENDDO

				! E       
				currentei = panele(paneli)
				ei = ne/2+1 ! ALL NEWBORN START FROM THE MEDIAN E  but if too many people died and started from median E, draw a new E for them

				panelz(paneli)		= zi    
				panellambda(paneli)	= lambdai
				panele(paneli) 		= ei 
				panelx(paneli)		= xi

       		else  IF (age .gt. 1) THEN
       			currentxi = panelx(paneli)
       			tempno 	  = omp_ran1() ! ran1(newiseed)   
	            xi 		  = 1
	            DO WHILE (tempno .gt. cdf_pr_x(currentxi,xi,currentzi,age-1))
	               xi = xi+1
	            ENDDO            
	            panelx(paneli)=xi          

	            IF (age.lt.RetAge) THEN
		            currentei = panele(paneli)   
		            tempno 	  = omp_ran1() ! ran1(newiseed)   
		            ei 		  = 1
		            DO WHILE (tempno .gt. cdf_pr_e(currentei,ei))
		               ei = ei+1
		            ENDDO            
		            panele(paneli)=ei 
	            ENDIF           
	     	ENDIF ! new age==1
		ENDDO ! paneli



	    ! Record data of top agents
     	if (simutime.ge.(MaxSimuTime-149)) then 

     		age_top = simutime-MaxSimuTime+150

			!$omp parallel do private(currenta,age,currentzi,currentlambdai,currentei,tklo,tkhi)
			DO ii=1,totpop
			    currenta 		= panela(ii)

			    if (currenta .ge. amax) then
			        tklo = na-1
			    elseif (currenta .lt. amin) then
		            tklo = 1
		        else
		            tklo = ((currenta - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
			    endif    
			    tkhi = tklo + 1        

			    if (panelage(ii).eq.1) then 
			    panelPV_a(ii) = ( (agrid(tkhi) - currenta) * V_Pr_nb(tklo,panelz_old(ii),panellambda_old(ii)) &
			                &  + (currenta - agrid(tklo)) * V_Pr_nb(tkhi,panelz_old(ii),panellambda_old(ii))) &
			                &  / ( agrid(tkhi) - agrid(tklo) )  + (1.0_dp+R)*currenta 
			    else 
			    panelPV_a(ii) = ( (agrid(tkhi) - currenta) * V_Pr(panelage(ii),tklo,panelz(ii),panellambda(ii),panele(ii),panelx(ii)) &
			                &  + (currenta - agrid(tklo)) * V_Pr(panelage(ii),tkhi,panelz(ii),panellambda(ii),panele(ii),panelx(ii))) &
			                &  / ( agrid(tkhi) - agrid(tklo) )  + (1.0_dp+R)*currenta 
			    endif 
		    enddo 


     		! print*, "Selecting top agents - Period", age_top
     		panelage_top(age_top,:) = panelage(top_ind)
     		panelz_top(age_top,:)   = panelz(top_ind)
     		panela_top(age_top,:)   = panela(top_ind)
     		panelx_top(age_top,:)   = panelx(top_ind)
     		panele_top(age_top,:)   = panele(top_ind)
     		panel_lambda_top(age_top,:)   = panellambda(top_ind)
     		panel_PV_top(age_top,:) = panelPV_a(top_ind)


     		!$omp parallel do private(tklo,tkhi,h_i)
     		do ii=1,80   
     		panelk_top(age_top,ii)  = min( theta(panelz_top(age_top,ii))*panela_top(age_top,ii) ,&
     			& (mu*P*xz_grid(panelx_top(age_top,ii),panelz_top(age_top,ii))**mu & 
     				& /(R+DepRate))**(1.0_dp/(1.0_dp-mu)) )

			prc_all_top(age_top,ii)    = 100.0_dp*real(count(panela.le.panela_top(age_top,ii)),8)/real(totpop,8)
     		prc_cohort_top(age_top,ii) = 100.0_dp* &
     							& real(count((panela.le.panela_top(age_top,ii)).and.(panelage.eq.panelage_top(age_top,ii))),8) &
     							& /real(count(panelage.eq.panelage_top(age_top,ii)),8)

			prc_PV_all_top(age_top,ii)    = 100.0_dp*real(count(panelPV_a.le.panel_PV_top(age_top,ii)),8)/real(totpop,8)
     		prc_PV_cohort_top(age_top,ii) = 100.0_dp* &
     							& real(count((panelPV_a.le.panel_PV_top(age_top,ii)).and.(panelage.eq.panelage_top(age_top,ii))),8) &
     							& /real(count(panelage.eq.panelage_top(age_top,ii)),8)

			if (panelage_top(age_top,ii).lt.RetAge) then 

     			if (panela_top(age_top,ii) .ge. amax) then
		            tklo = na-1
		        else if (panela_top(age_top,ii) .lt. amin) then
		            tklo = 1
		        else
		            tklo = ((panela_top(age_top,ii) - amin)/(amax-amin))**(1.0_DP/a_theta)*(na-1)+1          
		        endif 
		        
		        tkhi = tklo + 1     

		        h_i  = ((agrid(tkhi) - panela_top(age_top,ii))*&
		        	& hours(panelage_top(age_top,ii),tklo,&
		        	&	    panelz_top(age_top,ii),panel_lambda_top(age_top,ii),&
		        	&	    panele_top(age_top,ii),panelx_top(age_top,ii)) &
		            &  + (panela_top(age_top,ii) - agrid(tklo))*&
		            & hours(panelage_top(age_top,ii),tkhi,&
		        	&	    panelz_top(age_top,ii),panel_lambda_top(age_top,ii),&
		        	&	    panele_top(age_top,ii),panelx_top(age_top,ii))  ) &
		                                &  / ( agrid(tkhi) - agrid(tklo) )  

				panel_YL_top(age_top,ii) = psi*( Wage*&
					& eff_un(panelage_top(age_top,ii),panel_lambda_top(age_top,ii),&
		        	&	    panele_top(age_top,ii))&
					& *h_i)**(1.0_dp-tauPL)
			else 
				panel_YL_top(age_top,ii) = RetY_lambda_e(panel_lambda_top(age_top,ii),&
		        	&	    panele_top(age_top,ii))
			endif 

     		enddo
     	endif
	    ! print*, "Simulation period", simutime
	ENDDO ! simutime
	print*,' '
	print*,'Averages'
	print*, sum(panelage)/real(totpop,8), sum(panelz)/real(totpop,8), sum(panele)/real(totpop,8), sum(panela)/real(totpop,8)


	print*, ' '
	print*, 'Writing simulation results for top agents'
	call system( 'mkdir -p ' // trim(Result_Folder) // 'Simul/' // trim(folder) )


	if (bench_indx.eq.1) then
		OPEN(UNIT=10, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panela_top_bench'			, STATUS='replace')
		OPEN(UNIT=11, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelage_top_bench'		, STATUS='replace')
		OPEN(UNIT=12, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelz_top_bench'			, STATUS='replace')
		OPEN(UNIT=27, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelK_top_bench'    		, STATUS='replace')
		OPEN(UNIT=28, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelx_top_bench'   	 	, STATUS='replace')
		OPEN(UNIT=29, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panele_top_bench'    		, STATUS='replace')
		OPEN(UNIT=30, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panel_lambda_top_bench'	, STATUS='replace')
		OPEN(UNIT=31, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panel_YL_top_bench'	    , STATUS='replace')
		OPEN(UNIT=32, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_all_top_bench'			, STATUS='replace')
		OPEN(UNIT=33, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_cohort_top_bench'	    , STATUS='replace')
		OPEN(UNIT=34, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panel_PV_top_bench'    	, STATUS='replace')
		OPEN(UNIT=35, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_PV_all_top_bench'		, STATUS='replace')
		OPEN(UNIT=36, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_PV_cohort_top_bench'	, STATUS='replace')
	else 
		OPEN(UNIT=10, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panela_top_exp'			, STATUS='replace')
		OPEN(UNIT=11, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelage_top_exp'		    , STATUS='replace')
		OPEN(UNIT=12, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelz_top_exp'			, STATUS='replace')
		OPEN(UNIT=27, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelK_top_exp' 	    	, STATUS='replace')
		OPEN(UNIT=28, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panelx_top_exp'	   	 	, STATUS='replace')
		OPEN(UNIT=29, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panele_top_exp'    		, STATUS='replace')
		OPEN(UNIT=30, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panel_lambda_top_exp'	    , STATUS='replace')
		OPEN(UNIT=31, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panel_YL_top_exp'	        , STATUS='replace')
		OPEN(UNIT=32, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_all_top_exp'			, STATUS='replace')
		OPEN(UNIT=33, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_cohort_top_exp'		, STATUS='replace')
		OPEN(UNIT=34, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'panel_PV_top_exp'    	    , STATUS='replace')
		OPEN(UNIT=35, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_PV_all_top_exp'		, STATUS='replace')
		OPEN(UNIT=36, FILE=trim(Result_Folder)//'Simul/'//trim(folder)//'prc_PV_cohort_top_exp'		, STATUS='replace')
	endif 


	WRITE  (UNIT=10, FMT='(F12.4)') panela_top
	WRITE  (UNIT=11, FMT=*) panelage_top  
	WRITE  (UNIT=12, FMT=*) panelz_top 
	WRITE  (UNIT=27, FMT='(F12.4)') panelK_top
	WRITE  (UNIT=28, FMT=*) panelx_top
	WRITE  (UNIT=29, FMT=*) panele_top
	WRITE  (UNIT=30, FMT=*) panel_lambda_top
	WRITE  (UNIT=31, FMT='(F12.4)') panel_YL_top
	WRITE  (UNIT=32, FMT='(F12.4)') prc_all_top
	WRITE  (UNIT=33, FMT='(F12.4)') prc_cohort_top
	WRITE  (UNIT=34, FMT='(F12.4)') panel_PV_top
	WRITE  (UNIT=35, FMT='(F12.4)') prc_PV_all_top
	WRITE  (UNIT=36, FMT='(F12.4)') prc_PV_cohort_top

	close (unit=10); close (unit=11); close (unit=12); close (unit=27)
	close (unit=28); close (unit=29); close (unit=30); close (unit=31)
	close (unit=32); close (unit=33); close (unit=34); close (unit=35); close (unit=36)


END SUBROUTINE SIMULATION_TOP




!========================================================================================
!========================================================================================
!========================================================================================




    function omp_ran1(idum)
        ! returns a uniform random number between 0 and 1
        ! see numerical recipes
        ! press,flannery,teukolsky & vetterling
        ! cambridge university press 1986 pp 191-2-3
        ! http://geco.mines.edu/prototype/Show_me_some_local_HPC_tutorials/examples/openmp/ranmod.f90
        ! http://inside.mines.edu/~tkaiser/fortran//
        ! https://people.sc.fsu.edu/~jburkardt/f_src/random_openmp/random_openmp.html
        ! http://jblevins.org/log/openmp
        implicit none
        integer, intent(inout), optional :: idum
        real(DP) :: omp_ran1
        integer  :: m1,m2,m3,ia1,ia2,ia3,ic1,ic2,ic3
        integer  :: iff
        integer  :: ix1,ix2,ix3,j
        real(DP) :: r(97),rm1,rm2
        
        parameter (m1=259200,ia1=7141,ic1=54773)
        parameter (m2=134456,ia2=8121,ic2=28411)
        parameter (m3=243000,ia3=4561,ic3=51349)
        data iff /0/
        
        !$OMP THREADPRIVATE(iff,ix1,ix2,ix3,j,r,rm1,rm2) 
        save iff,ix1,ix2,ix3,j,r,rm1,rm2
        if(present(idum))then
            if (idum<0.or.iff.eq.0)then
                rm1=1.0_dp/m1
                rm2=1.0_dp/m2
                iff=1
                ix1=mod(ic1-idum,m1)
                ix1=mod(ia1*ix1+ic1,m1)
                ix2=mod(ix1,m2)
                ix1=mod(ia1*ix1+ic1,m1)
                ix3=mod(ix1,m3)
                do  j=1,97
                    ix1=mod(ia1*ix1+ic1,m1)
                    ix2=mod(ia2*ix2+ic2,m2)
                    r(j)=(float(ix1)+float(ix2)*rm2)*rm1
                enddo 
                idum=1
            endif
        endif
        
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ia2*ix2+ic2,m2)
        ix3=mod(ia3*ix3+ic3,m3)
        j=1+(97*ix3)/m3

        if(j>97.or.j<1)then
            write(*,*)' error in omp_ran1 j=',j
            stop
        endif

        omp_ran1=r(j)
        r(j)=(float(ix1)+float(ix2)*rm2)*rm1
        return
    end function omp_ran1



!========================================================================================
!========================================================================================
!========================================================================================



end Module Simulation_Module