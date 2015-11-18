function V = Value_Function(agrid,Cons,Hours,par)



% Compute value function at grid nodes
%% Final age
	age=MaxAge ;
	for ai=1:na    
	    for zi=1:nz
	        for lambdai=1:nlambda          
	              for ei=1:ne
	              	if (sigma==1.0) 
	                  	V(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) ; 
	                else 
	                	V(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)^gamma) ...
                   		   * (1.0-Hours(age,ai,zi,lambdai,ei))^(1.0-gamma))^(1.0-sigma)/(1.0-sigma) 
	                end 
	                  % print*,Cons(age, ai, zi, lambdai, ei),  V(age, ai, zi, lambdai, ei) 
	                  % pause
	              end % ei          
	        end % lambdai
	    end % zi
	end % ai

%% Retirement Period
	for age=MaxAge-1,RetAge,-1
	    for zi=1:nz
	        for lambdai=1:nlambda          
	              for ei=1:ne            
	          		if (sigma==1.0) 
	                    CALL spline( agrid, V(age+1, :, zi, lambdai, ei) , na , ...
	                     1.0/Cons(age+1, 1, zi, lambdai,ei) , 1.0/Cons(age+1, na, zi, lambdai,ei) , ValueP2)  
	                else 
	                	CALL spline( agrid, V(age+1, :, zi, lambdai, ei) , na , ...
                  			 gamma*MBGRID(1,zi) *Cons(age+1, 1, zi, lambdai, ei) ^((1.0-sigma)*gamma-1.0)/(1+tauC), ...
                    		 gamma*MBGRID(na,zi)*Cons(age+1, na, zi, lambdai, ei)^((1.0-sigma)*gamma-1.0)/(1+tauC), ValueP2)  
                  	end 
	                  
	                    for ai=1:na    
	                        call splint( agrid, V(age+1, :, zi, lambdai, ei), ...
	                                 ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))  ; 
	    					if sigma==1.0 
                                V(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) ...
                                     + beta*survP(age)* ValueP(ai);
	                        else 
                                V(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)^gamma) ...
                                     * (1.0-Hours(age,ai,zi,lambdai,ei))^(1.0-gamma))^(1.0-sigma)/(1.0-sigma) ...
                                     + beta*survP(age)* ValueP(ai);
                           	end
	                    end % ai
	              
	            end % ei          
	        end % lambdai
	    end % zi
	end % age
   

%% Working Period
	for age=RetAge-1,1,-1
	    for zi=1:nz
	        for lambdai=1:nlambda          
	              for ei=1:ne
	                    for ai=1:na    
	                          ExpValueP(ai) = sum(V(age+1, ai, zi, lambdai, :) * pr_e(ei,:))
	                    end

	                    if (sigma==1.0)  
	                    CALL spline( agrid, ExpValueP , na , ...
	                     sum(pr_e(ei,:)/Cons(age+1, 1, zi, lambdai,:)) , sum(pr_e(ei,:)/Cons(age+1, na, zi, lambdai,:)) , ValueP2)  
	                    else 
	                    CALL spline( agrid, ExpValueP , na , ...
		                     (gamma*MBGRID(1,zi)/(1.0+tauC)) * sum(pr_e(ei,:)* ...
		                     Cons(age+1, 1, zi, lambdai, :)^((1.0-sigma)*gamma-1.0) * ...
		                     (1.0-Hours(age+1,1,zi,lambdai,:))^((1.0-gamma)*(1.0-sigma))),...                    
		                     (gamma*MBGRID(na,zi)/(1.0+tauC)) * sum(pr_e(ei,:)* ...
		                     Cons(age+1, na, zi, lambdai, :)^((1.0-sigma)*gamma-1.0) * ...
		                     (1.0-Hours(age+1,na,zi,lambdai,:))^((1.0-gamma)*(1.0-sigma))),...
		                     ValueP2)  
	                    end 

	                    for ai=1:na 
	                    	if (sigma==1.0) 
	                         call splint( agrid, ExpValueP, ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
	                         V(age, ai, zi, lambdai, ei) = log(Cons(age, ai, zi, lambdai, ei)) ...
	                             + ((1.0-gamma)/gamma) * log(1.0-Hours(age, ai, zi, lambdai, ei)) + beta*survP(age)*ValueP(ai)
	                        else 
	                          call splint( agrid, ExpValueP, ValueP2, na, Aprime(age,ai,zi,lambdai, ei), ValueP(ai))   
		                         V(age, ai, zi, lambdai, ei) = ((Cons(age,ai,zi,lambdai,ei)^gamma) ...
		                            * (1.0-Hours(age,ai,zi,lambdai,ei))^(1.0-gamma))^(1.0-sigma)/(1.0-sigma) ...
		                            + beta*survP(age)* ValueP(ai)
	                        end 
	                    end % ai
	               end % ei          
	        end % lambdai
	    end % zi
	end % age

end