MODULE Opt_Tax_Parameters
    USE parameters
    USE global

    IMPLICIT NONE 

	! Variables for optimal taxe 
	INTEGER  :: opt_tax_switch, tauindx
	REAL(DP) :: Opt_TauK, Opt_TauW, maxbrentvaluet, brentvaluet, GBAR_K

END MODULE Opt_Tax_Parameters