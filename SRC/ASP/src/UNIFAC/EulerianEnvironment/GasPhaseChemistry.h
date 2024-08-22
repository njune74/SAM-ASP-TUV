!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! GasPhaseChemistry.h
!! This file contains routines for inputing and outputing gas phase
!! chemical concnetrations to the grid points, allocating the 
!! gas phase chemical grid, conveting units, and calculating 
!! gas phase reaction rates.
!!
!! WARNING!: MELAM considers concentrations different from burdens,
!! with Burden = Concentration*GridPointVolume. ASP ignores the 
!! distinction, so burden and concentration are equivalent, and 
!! in units of molecules/cm^{3}. 
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History	                     !!
!! 07/30  2007   Matt Alvarado     Removed spurious eps factor from RH calc  !!
!!                                 in SetWaterVaporField		     !!
!! 02/15  2012   Matt Alvarado     Removing Eulerian grids, making ASP       !!
!!                                 a one-box model or subroutine.            !!
!! 01/03  2013   Matt Alvarado     Fixed second order "+ M" reactions from   !!
!!                                   MCM v3.2 to act as pseudo-1st order     !!
!!                                   and added IUPAC style 3-body rxn rates  !!
!!                                   to ResetGasPhaseReactionRates           !!
!! 01/18  2013   Matt Alvarado     Removed hardwired CO + OH reaction and    !!
!!                                  replaced with a K1 + K2 reaction         !!
!! 01/18  2013   Matt Alvarado    Added  thermal 3rd order rate constant     !!
!!                                 and fixed N2O5 and O1D hardwired rxns     !!
!!                                 and made one call to GetM per chem step   !!
!! 01/25  2013   Matt Alvarado    Added GetTotalPeroxy to top of             !!
!!                                 ResetGasPhaseReactionRates and converted  !!
!!                                 + RO2_T reactions to special first order  !!
!! 01/25  2013   Matt Alvarado    Added reactions like HO2 + HO2 + H2O       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. SUBROUTINE ResetGasPhaseReactionRates (ReactionRates,ArraySize)
!!    NOTE: Any changes made to rate constant functions here must be reflected
!!    in subroutine MakeODESet in Chemistry/ChemistryParametersInitRoutines.h
!! 2. SUBROUTINE AllocateGasChemistryGrid (NumbChems)
!! 3. SUBROUTINE SetHomogeneousChemFieldPPB (Index, Value)
!! 4. SUBROUTINE SetHomogeneousChemFieldMGperM3 (Index, Value, MolecMass)
!! 5. FUNCTION	 GetGasChemVector () RESULT (Y)
!! 6. SUBROUTINE UpdateChemicalBurdens (Y)
!! 7. SUBROUTINE UpdateChemicalConcentrations (, Y)
!! 8. FUNCTION   GetGasChems () RESULT (Y)
!! 9. FUNCTION   GetGasChem (SingleChemIndex) RESULT (Y)
!!10. FUNCTION   GetGasBurden (SingleChemIndex) RESULT (Y)
!!11. FUNCTION   GetGasChemFromBurden (Burden)
!!12. FUNCTION GetGridCellChemBurden (SingleChemIndex)
!!13. SUBROUTINE AddToGridCellChemBurden ( SingleChemIndex, DeltaChemBurden)
!!14. SUBROUTINE ReplaceGridCellChemBurden ( SingleChemIndex, ChemBurden)
!!15. FUNCTION   GetPartialPressure (SingleChemIndex)
!!16. SUBROUTINE SetWaterVaporField (RH)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ************************************************************************
! Modified by CMB for use by gfortran: dlog10() needs a double precision
! value, so I forced an explicit cast each time
! ************************************************************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Each timestep the reaction rates are recalculated !!
!! This subroutine performs that service.	     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ResetGasPhaseReactionRates (ReactionRates,ArraySize, TotalPeroxy, TotalAcylPeroxy, H2Oconc)

  Use ModelParameters
  USE InfrastructuralCode, ONLY : Transcript, ERROR, Int2Str, Real2Str
  IMPLICIT NONE

  REAL*8,  INTENT(in) :: ReactionRates(:,:)
  INTEGER, INTENT(in) :: ArraySize(2)
  REAL*8              :: TotalPeroxy, TotalAcylPeroxy, H2OConc

  INTEGER :: L
  REAL*8	:: Temp,Mgas,X(3),k0,ki, Y_o, Y_o_term, Y_inf_term, Z, F, &
  	   alpha, beta, m_o, m_inf, Ratio, Scale, Pres, k1, k2, k3, kf, Keq, &
	   KINF, NCD, FD, RO2_T, RCO3_T, C_H2O, Knaught, N2gas, O2gas

  Temp = GetTemp()
  Mgas = GetM   ()
!  WRITE(*,*) "Mgas", Mgas
  Pres = GetPress() !mbar
  RO2_T = TotalPeroxy
  RCO3_T = TotalAcylPeroxy  
  C_H2O = H2Oconc

  !! Loop over all reactions
  DO L = 1, ArraySize(1)
      !! Use Select to Differentiate Reaction Types
      SELECT CASE (INT(ReactionRates(L,1)))	
      ! Number of Bodies in the reaction

      !! One body reactions (units are per-second)
      CASE(1)
        
	!Photolysis Case
	IF (INT(ReactionRates(L,3)) .EQ. 0) THEN
	  GasPhaseChemicalRates(L) = ReactionRates(L,3)

        !PAN Degradation case (NOT USED! Use type 4!)
        ELSE IF (INT(ReactionRates(L,3)) .EQ. 1) THEN
	  k0   = ReactionRates(L,4)*EXP(ReactionRates(L,5)/Temp)* Mgas
	  ki   = ReactionRates(L,6)*EXP(ReactionRates(L,7)/Temp)
	  GasPhaseChemicalRates(L) = (k0*ki)/(k0+ki) * ReactionRates(L,8)** &
						(1./(1.+(DLOG10(dble(k0/ki))**2.)))
								
        !Isomerization Case
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 2) THEN
	  GasPhaseChemicalRates(L) = ReactionRates(L,4) *		&
					 Temp**ReactionRates(L,5)    *	&
					 EXP(ReactionRates(L,6)/Temp)**	&
					 ReactionRates(L,7)					
	!Heterogeneous Case
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 3) THEN
	  !Temporary - call ResetHeteroRates to set rate!
	  GasPhaseChemicalRates(L) = 0.0 

	!Equilibrium Degradation case
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 4) THEN			
	  k0  = ReactionRates(L,4)*(Temp/300.)**(-1.*ReactionRates(L,5)) * Mgas
	  ki  = ReactionRates(L,6)*(Temp/300.)**(-1.*ReactionRates(L,7))
	  kf  =	(k0*ki)/(k0+ki) * ReactionRates(L,8)**			&
				(1./(1.+(DLOG10(dble(k0/ki))**2.)))
	  Keq = ReactionRates(L,9)*EXP(ReactionRates(L,10)/TEMP)
								
	  GasPhaseChemicalRates(L) = kf/Keq
						
	!Fixed Rate Case
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 5) THEN
	  GasPhaseChemicalRates(L) = ReactionRates(L,4)
							
						
	!HNO4 Thermal degradation (obsolete!)
        ELSE IF (INT(ReactionRates(L,3)) .EQ. 6) THEN
          Knaught = 4.1E-5*EXP(-10649.2/TEMP)*Mgas
	  KINF = 5.7E+15*EXP(-11172.6/TEMP)
	  GasPhaseChemicalRates(L) = Knaught/(1.0+(Knaught/KINF)) &
				 *0.5**(1.0+(LOG10(Knaught/KINF))**2)**(-1)
							
	!N2O5 Thermal degradation - IUPAC 2012
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 7) THEN
	  Knaught = Mgas*1.3E-3*((300.0/TEMP)**3.5)*EXP(-11000.0/TEMP)
	  !Knaught = Mgas*1.3E-3*((TEMP/300.0)**3.5)*EXP(-11000.0/TEMP)
	  KINF = 9.7E+14*((TEMP/300.0)**0.1)*EXP(-11080.0/TEMP)
          NCD = 0.75-1.27*(DLOG10(dble(0.35)))
!          FD = 10.**(DLOG10(0.35)/(1.+(DLOG10(Knaught/KINF)/NCD)**2))
          FD = 10.**(DLOG10(dble(0.35))/(1.+(DLOG10(dble(Knaught/KINF))/NCD)**2))
          GasPhaseChemicalRates(L) = FD*(Knaught*KINF)/(Knaught+KINF)	
!          FD = 10.**(DLOG10(0.35)/(1.+(DLOG10(Knaught/KINF)/NCD)**2))
          FD = 10.**(DLOG10(dble(0.35))/(1.+(DLOG10(dble(Knaught/KINF))/NCD)**2))
          GasPhaseChemicalRates(L) = FD*(Knaught*KINF)/(Knaught+KINF)	

!          WRITE(*,*) GasPhaseChemicalRates(L), (1.+(DLOG10(Knaught/KINF)/NCD)**2)
!          CALL TRANSCRIPT("N2O5 deg: "//REAL2str(GasPhaseChemicalRates(L)))
!          WRITE(*,*) "Mgas  ", Mgas, EXP(-11000.0/TEMP)
!          CALL TRANSCRIPT("k_o: "//REAL2str(Knaught))
!          CALL TRANSCRIPT("k_inf: "//REAL2str(KINF))
!          CALL TRANSCRIPT("NCD: "//REAL2str(NCD))
!          CALL TRANSCRIPT("FD: "//REAL2str(FD)) 							

        !O1D Thermal deactivation - IUPAC 2012, 
        !Assuming N2 = 0.79*Mgas and O2 = 0.21*Mgas
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 8) THEN

	   GasPhaseChemicalRates(L) = Mgas *(0.21*3.2e-11*exp(67./TEMP) &
                                            +0.79*2.0e-11*exp(130./TEMP))
				
	!O + O2 + M => O3 + M 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 9) THEN
	   N2gas = 0.79*Mgas
           O2gas = 0.21*Mgas
           GasPhaseChemicalRates(L) = (5.6e-34*(TEMP/300.)**(-2.6)*N2gas  &
                                     +6.0e-34*(TEMP/300.)**(-2.6)*O2gas)*O2gas
							
	!RAD* + O2  
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 10) THEN
	   GasPhaseChemicalRates(L) = 0.9*4.62e7/Temp/60
  
        !reaction form k0 = a*exp(B/T) (cm3/molec/s) and ki = c*exp(D/T) (1/s)
        !From MCM v3.2 for PAN degradation
        !See KBPAN on http://mcm.leeds.ac.uk/MCM/parameters/complex.htt
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 11) THEN
	  k0 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP) * Mgas
	  ki = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
	  NCD = 0.75-1.27*(DLOG10(dble(ReactionRates(L,8))))
          FD = 10.**(DLOG10(dble(ReactionRates(L,8)))/(1.+(DLOG10(dble(k0/ki))/NCD)**2))
          GasPhaseChemicalRates(L) = FD*(k0*ki)/(k0+ki)

        !+ RO2_T pseudo-1st order Case (MJA, 01-25-2013)
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 12) THEN
	  GasPhaseChemicalRates(L) = RO2_T*ReactionRates(L,4) *		&
					 Temp**ReactionRates(L,5)    *	&
					 EXP(ReactionRates(L,6)/Temp)**	&
					 ReactionRates(L,7)	
	ELSE
	   CALL ERROR("Error in ResetGasPhaseReactionRates - your 1st order rate constant doesn't make sense")								
	END IF
						
					
      !! Two body reactions (units are 1 / ppb / s)
      CASE(2)
						
	!Normal 2-Body Case
	IF (INT(ReactionRates(L,3)) .EQ. 0) THEN
           GasPhaseChemicalRates(L) = ReactionRates(L,3)          *	&
					Temp**ReactionRates(L,4)    *	&
					EXP(ReactionRates(L,5)/Temp)**	&
					ReactionRates(L,6)
	!Organic Nitrate Case
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 1 .OR. &
         INT(ReactionRates(L,3)) .EQ. 2) THEN
								
           !!Formula taken from Carter and Atkinson (1989)
	   alpha = 1.94E-22
	   beta = 0.97
	   Y_o = alpha*EXP(beta*ReactionRates(L,8))
	   m_o = 0.
	   m_inf = 8.1

	   Y_o_term = Y_o*Mgas*(Temp/300.0)**(-m_o)
	   Y_inf_term = 0.826*(Temp/300)**(-m_inf)

	   Z = (1+(log10(Y_o_term/Y_inf_term))**(2.0))**(-1.0)
           F = 0.411
	   Ratio=ReactionRates(L,9)*(Y_o_term/(1+(Y_o_term/Y_inf_term)))*F**(Z)

	   !Special scaling for a few reactions
	   Ratio = ReactionRates(L,10)*Ratio + ReactionRates(L,11)
								
	   IF(INT(ReactionRates(L,3)) .EQ. 1) THEN
		Scale = Ratio/(1.0 + Ratio)
	   ELSE
		Scale = 1.0 - Ratio/(1.0 + Ratio)
	   END IF
	   GasPhaseChemicalRates(L) = Scale*ReactionRates(L,4) *  &
				   Temp**ReactionRates(L,5)    *  &
				   EXP(ReactionRates(L,6)/Temp)**	&
				   ReactionRates(L,7)
							
	!This is a k = k1 + k2 case
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 3) THEN	
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
	  k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
          GasPhaseChemicalRates(L) = k1 + k2
							
	!This is a k = k1 + k2*M case 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 4) THEN
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
	  k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
	  GasPhaseChemicalRates(L) = k1 + k2*Mgas
							
	!This is a k = k1 + k3*M*(1+k3*M/k2) case 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 5) THEN
	   k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
	   k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
	   k3 = ReactionRates(L,8)*exp(ReactionRates(L,9)/TEMP)
	  GasPhaseChemicalRates(L) = k1 + k3*Mgas*(1 + k3*Mgas/k2)
							
	!This is a heterogeneous reaction 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 6) THEN
	  !Temporary - call ResetHeteroRates to set rate!
	  GasPhaseChemicalRates(L) = 0.0

	!k = A*exp(B/T)*(1-(1/(1+C*exp(D/T)))) 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 7) THEN
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
          k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
          GasPhaseChemicalRates(L) = k1*(1 - (1/(1 + k2))) 

	!k = A*exp(B/T)*(1/(1+C*exp(D/T))) 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 8) THEN
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
          k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
          GasPhaseChemicalRates(L) = k1*(1/(1 + k2)) 

        !k = A*exp(B/T)*(1-C*exp(D/T)) 
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 9) THEN
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
          k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
          GasPhaseChemicalRates(L) = k1*(1 - k2) 

	!k = C_H2O*(A1*exp(C1/T)+Mgas*A2*exp(C2/T)
       ELSE IF (INT(ReactionRates(L,3)) .EQ. 10) THEN
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
          k2 = ReactionRates(L,6)*exp(ReactionRates(L,7)/TEMP)
          GasPhaseChemicalRates(L) = C_H2O*(k1+Mgas*k2)

	!!     k = TMP3*TMP4/100
        !!     TMP3 = A*exp(B/T)
        !!     TMP4 = [C/T + D*PRES+F] 
        !!     PRES = ATMOSPHERIC PRESSURE IN PASCAL
       ELSE IF (INT(ReactionRates(L,3)) .EQ. 11) THEN
	  k1 = ReactionRates(L,4)*exp(ReactionRates(L,5)/TEMP)
          k2 = ReactionRates(L,6)/TEMP + ReactionRates(L,7)*PRES*100.0 + ReactionRates(L,8)
          GasPhaseChemicalRates(L) = k1*k2/100.0
        END IF

     !! Three body reactions, using the troe pressure-broadening 
     !! From from p.1303 of Seinfeld & Pandis
     !! This yields an effective second-order rate constant.
     !! Rates are set in 1 / ppb / ppb / sec 
      CASE(3)
     
        !Case for third-order rate constant
	!(like for O + O2 + M => O3 + M) (NOT USED - see first order class 9)
	IF (INT(ReactionRates(L,3)) .EQ. 1) THEN
	   GasPhaseChemicalRates(L) = ReactionRates(L,4)*(Temp/300.)**(-1.*ReactionRates(L,5)) * Mgas
						
        !Thermal 3rd order rate constant
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 2) THEN
           GasPhaseChemicalRates(L) = Mgas*ReactionRates(L,3)*	&
					Temp**ReactionRates(L,4)    *	&
					EXP(ReactionRates(L,5)/Temp)**	&
					ReactionRates(L,6)

        !From MCM v3.2 for PAN formation (Also IUPAC format?)
        !See KFPAN on http://mcm.leeds.ac.uk/MCM/parameters/complex.htt
	ELSE IF (INT(ReactionRates(L,3)) .EQ. 3) THEN
	  k0 = ReactionRates(L,4)*(Temp/300.)**(-1.*ReactionRates(L,5)) * Mgas
	  ki = ReactionRates(L,6)*(Temp/300.)**(-1.*ReactionRates(L,7))
	  NCD = 0.75-1.27*(DLOG10(dble(ReactionRates(L,8))))
          FD = 10.**(DLOG10(dble(ReactionRates(L,8)))/(1.+(DLOG10(dble(k0/ki))/NCD)**2))
          GasPhaseChemicalRates(L) = FD*(k0*ki)/(k0+ki)													
	!Calculate pseudo-2nd order rate constant
	ELSE !JPL Style							
	  k0 = ReactionRates(L,3)*(Temp/300.)**(-1.*ReactionRates(L,4)) * Mgas
	  ki = ReactionRates(L,5)*(Temp/300.)**(-1.*ReactionRates(L,6))
	  GasPhaseChemicalRates(L) = (k0*ki)/(k0+ki) * ReactionRates(L,7)** &
					(1./(1.+(DLOG10(dble(k0/ki))**2.)))
	END IF
      END SELECT

	
      END DO

END SUBROUTINE ResetGasPhaseReactionRates

!!! --- Allocate the gas phase chemical grids --- !!!
!!! --- Called from  Chemistry Parameters.f90 --- !!!
SUBROUTINE AllocateGasChemistryGrid (NumbChems)

	USE InfrastructuralCode, ONLY : ERROR
	IMPLICIT NONE

	INTEGER :: NumbChems, Allocation_Error 
		
	!! Bring the monster to life
	ALLOCATE (GridGasChem(NumbChems),stat = allocation_error)

	! Allocation Error !
	IF (allocation_error > 0) &
	CALL ERROR ("Allocation of GridGasChem could not proceed in AllocateGasChemistryGrid()")

        !! Water Mixing Ratios have been temporarily stored and are now put
	!! in their final holding place

	GridGasChem(1) = TemporaryWaterArray

	RETURN
END SUBROUTINE AllocateGasChemistryGrid

SUBROUTINE AllocateGasReactionGrid (NumbRates)

	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	INTEGER :: NumbRates, Allocation_Error

	!! Bring the monster to life
	ALLOCATE (GasPhaseChemicalRates(NumbRates),stat = allocation_error)

	! Allocation Error !
	IF (allocation_error > 0) CALL ERROR ("Allocation of GridGasChem could not proceed in AllocateGasReactionGrid()")

	RETURN
END SUBROUTINE AllocateGasReactionGrid
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialize all of the grid points with the same initial value !!
!! input value of the chemical concentration is in ppbv			 !!
!! Stored value is in total abundance.                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetHomogeneousChemFieldPPB (Index, Value)

	USE InfrastructuralCode, ONLY : ERROR, Transcript
	USE ModelParameters,	 ONLY : ppbscale,chemscale
	implicit none

	integer ::  Index
	real*8  :: Value, Mgas

        !! Set each chemical to the specified value.
        Mgas = GetM ()
	GridGasChem(Index) = Value*ppbscale*Mgas/ChemScale

	RETURN
END SUBROUTINE SetHomogeneousChemFieldPPB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialize all of the grid points with the same initial value !!
!! input value of the chemical concentration is in ppbv	         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetHomogeneousChemFieldMGperM3 (Index, Value, MolecMass)

	USE InfrastructuralCode, ONLY : ERROR, Transcript
	USE ModelParameters,	 ONLY : chemscale,grams,avogadro
	implicit none

	integer ::  Index
	real*8  :: Value, Mgas
	real*8  :: MolecMass			
        !! Molecular mass of gas chem # Index -- 
        !! circular ref to access chemistry module here.

	GridGasChem(Index) = Value/1.e12/	&
				MolecMass*grams/chemscale*Avogadro

	RETURN
END SUBROUTINE SetHomogeneousChemFieldMGperM3


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Give a vector of the current chemical concentrations to the outside world.!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GetGasChemVector () RESULT (Y)

	IMPLICIT NONE

	!! External Variables
	REAL*8  :: Y(ylen)

	!! Internal Variables
	INTEGER   :: I

	DO I = 1, HowManyGasChems
		Y(I) = GridGasChem(I)
	END DO

	IF (ylen > HowManyGasChems) THEN
		DO I = HowManyGasChems+1, ylen
			Y(I) = 0.
		END DO
	END IF
END FUNCTION GetGasChemVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Once a routine calculates new values for the chemicals at a grid point !!
!! (which it does locally), it places those values back in the grid.      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Call it like this (assuming the vector you're using is longer than the !!
!! evolving gas chemical length:					  !!
!!    CALL UpdateChemicalConcentrations(Y(1:HowManyEvolveGasChems))       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE UpdateChemicalBurdens (Y)

	IMPLICIT NONE

	!! External Variables
        REAL*8  :: Y(HowManyEvolveGasChems)

	!! Internal Variables
	INTEGER :: I

	!! We only replace those chemicals that are allowed to evolve 
	!! (which are listed first, by construction)
	DO I = 1, HowManyEvolveGasChems
		IF (Y(I) .EQ. -1) CYCLE
		CALL ReplaceGridCellChemBurden (I, Y(I))
	END DO

END SUBROUTINE

SUBROUTINE UpdateChemicalConcentrations (Y)

	IMPLICIT NONE

	!! External Variables
	REAL*8  :: Y(HowManyEvolveGasChems)

	!! Internal Variables
	INTEGER :: I

	!! We only replace those chemicals that are allowed to evolve 
	!! (which are listed first, by construction)
	DO I = 1, HowManyEvolveGasChems
		GridGasChem(I) = Y(I)
	END DO

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retrieve the Chemistry Vector Y() Interpolated to Any Point !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GetGasChems () RESULT (Y)

	IMPLICIT NONE

	REAL*8				:: Y(HowManyGasChems)	! Output
	INTEGER				:: i

	DO i = 1,HowManyGasChems
		Y(i) = GridGasChem(i)
	END DO

	RETURN
END FUNCTION GetGasChems

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retrieve a single Gas Phase Chemical Interpolated to Any Point. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetGasChem (SingleChemIndex) RESULT (Y)

	implicit none

	INTEGER, intent(in) :: SingleChemIndex

	Y = GridGasChem(SingleChemIndex)

	RETURN
END FUNCTION GetGasChem

REAL*8 FUNCTION GetGasBurden (SingleChemIndex) RESULT (Y)

	implicit none

	INTEGER, intent(in) :: SingleChemIndex
		
	Y = GridGasChem(SingleChemIndex)

	RETURN
END FUNCTION GetGasBurden

REAL*8 FUNCTION GetGasChemFromBurden (Burden)

	IMPLICIT NONE

	REAL*8, intent (in) :: Burden

	GetGasChemFromBurden = Burden

	RETURN
END FUNCTION GetGasChemFromBurden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retrieve a single Gas Phase Chemical Total Burden within a Grid Cell !!
!! (this is just the concentration - I've removed the scaling by volume !!
!! (Matt Alvarado, 4/12/06).						!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetGridCellChemBurden (SingleChemIndex)

	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	INTEGER, intent(in) :: SingleChemIndex
		
	GetGridCellChemBurden = GridGasChem(SingleChemIndex)

	RETURN
END FUNCTION GetGridCellChemBurden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This one replaces that same burden number when condensation is over !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE AddToGridCellChemBurden (SingleChemIndex, DeltaChemBurden)

	implicit none

	REAL*8,  intent(in) :: DeltaChemBurden
	INTEGER, intent(in) :: SingleChemIndex

		GridGasChem(SingleChemIndex) = GridGasChem(SingleChemIndex) + &
						DeltaChemBurden

	RETURN
END SUBROUTINE AddToGridCellChemBurden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This one replaces that same burden number when condensation is over !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReplaceGridCellChemBurden (SingleChemIndex, ChemBurden)

	implicit none

	REAL*8,  intent(in) :: ChemBurden
	INTEGER, intent(in) :: SingleChemIndex

	GridGasChem(SingleChemIndex) = ChemBurden

	RETURN
END SUBROUTINE ReplaceGridCellChemBurden


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retrieve a single Gas Phase Chemical Partial Pressure.                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetPartialPressure (SingleChemIndex)

	USE ModelParameters, ONLY : kB

	implicit none

	INTEGER, intent(in) :: SingleChemIndex

	!! Calculate partial pressure using Dalton's Law
	GetPartialPressure = kB * GetTemp() * GetM() * &
				GridGasChem(SingleChemIndex)

	RETURN
END FUNCTION GetPartialPressure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Populate the Gas Phase Water Vapor Fields !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetWaterVaporField (RH)

	USE InfrastructuralCode, ONLY : ERROR
	USE ModelParameters,     ONLY : ppbscale, chemscale

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: RH

	!! Internal Variables
	INTEGER :: allocation_error
	REAL*8  :: SPress !, MM
				
	SPress = GetSatVapPress()
	TemporaryWaterArray  = ChemScale * SPress * RH * GetM() / GetPress()  

END SUBROUTINE SetWaterVaporField

