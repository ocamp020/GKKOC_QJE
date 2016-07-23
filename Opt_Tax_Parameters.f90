MODULE Opt_Tax_Parameters
    USE parameters
    USE global

    IMPLICIT NONE 

	! Variables for optimal taxe 
	INTEGER  :: opt_tax_switch, tauindx
	REAL(DP) :: Opt_TauK, Opt_TauW, Opt_psi, maxbrentvaluet, brentvaluet

END MODULE Opt_Tax_Parameters