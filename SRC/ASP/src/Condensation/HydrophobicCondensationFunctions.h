!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! HydrophobicCondensationFunctions.h
!! Contains many functions needed by the numerical integration routines
!! for organic dissolution in the aqueous and hydrophobic phases.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						             !!
!!								   	     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History		             !!
!! 08/28  2006   Matt Alvarado     1. Changed zero flag for	             !!
!!				   EquilibrateOrganicReaction		     !!
!!				   2. Set a minimum value for		     !!
!!				      TotalMolesOM (10^-30)		     !!
!! 08/29  2006   Matt Alvarado	   Reset minimum value for		     !!
!!				      TotalMolesOM (10^-25)		
!! 09/08  2006   Matt Alvarado     Created BinaryWaterActivityOrganic
!! 07/10  2007   Matt Alvarado     Rewrote CalculateSurrogateGasConc
!!				   to set all negative values to 0.
!! 07/17  2007   Matt Alvarado     Edited TotalMolesOM to apply
!!				  minimum flag correctly.
!! 09/10  2007   Matt Alvarado     Edited EquilibrateOrganicParticle
!!				  to just exit after 500 iterations
!! 05/26  2010   Matt Alvarado    Just use unity activity coefficients
!! 10/15  2010   Matt Alvarado    Remove Optional Arguments
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 01/02  2013   Matt Alvarado    Replaced UR21 with CH3COCO2H               !!
!! 01/23  2013   Matt Alvarado    Added flag to SaturationVaporPressure to   !!
!!                                  decide between two vaopr pressure methods!!
!! 01/24  2013   Matt Alvarado     Added flag to UpdateHydrophobicUNIFAC to  !!
!!                                  decide between two org act. coeff.methods!!
!! 01/24  2013   Matt Alvarado    Removed calls to CalculateSurrogateGasConc !!
!!                                 and UpdateGasPhaseOrganics from           !!
!!                                 EquilibrateHydrophobicPhaseAtGridPoint    !!
!!                                 and EquilibrateAllOrganicsAtGridPoint     !!
!! 01/25  2013   Matt Alvarado    Created wrapper subroutines                !!
!!                                 UpdateHydrophobicUNIFAC_all and           !!
!!                                 UpdateHydrophilicUNIFAC_all               !!
!! 01/30  2013   Matt Alvarado    Added flag to UpdateHydrophilicUNIFAC to   !!
!!                                  decide between two org act. coeff.methods!!
!!                                Removed CalculateSurrogateGasConc,         !!
!!                                UpdateGasPhaseOrganics, TotalAqMoles and   !!
!!                                CalcInfiniteDilutionHydrophilicUNIFAC      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		     !!
!! 1. FUNCTION SaturationVaporPressure (Temperature, OrganicDissolutionRxn)  !!
!! 2. FUNCTION TotalMolesOM(InParticle)                                      !!
!! 3. SUBROUTINE UpdateHydrophobicUNIFAC(InParticle, Temperature)
!! 4. SUBROUTINE EquilibrateHydrophobicPhaseAtGridPoint()
!! 5. SUBROUTINE EquilibrateHydrophobicReaction (GasPhaseBurden, WhichRxn, &
!!               InParticle, Temperature, ReturnType, OptErrorTolerance)
!! 6. FUNCTION HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, & 
!!             WhichRxn, Temperature, InParticle)
!! 7. SUBROUTINE EquilibrateHydrophilicReaction
!! 8. FUNCTION HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, & 
!!              WhichRxn, Temperature, InParticle)		
!! 9. FUNCTION EffectiveOrganicHenrysLaw (Temperature, ProtonMolality, & 
!!              AqOrganicDissolutionRxn)
!! 10. SUBROUTINE UpdateHydrophilicUNIFAC(InParticle, Temperature)
!! 12. SUBROUTINE EquilibrateAllOrganicsAtGridPoint
!! 13. SUBROUTINE EquilibrateOrganicReaction
!! 14. FUNCTION TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, InParticle, & 
!!                 Temperature)
!! 15. SUBROUTINE EquilibrateOrganicParticle
!! 16. FUNCTION BinaryWaterActivityOrganic
!! 17. SUBROUTINE UpdateHydrophobicUNIFAC_all
!! 18. SUBROUTINE UpdateHydrophilicUNIFAC_all
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate SaturationVaporPressure !!
!! This is indexed by the organic    !!
!! dissolution reaction. 	     !!
!! Puts out in MBar, so make         !!
!! sure to convert the weird         !!
!! pressure units.		     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SaturationVaporPressure  (Temperature, OrganicDissolutionRxn)

	USE ModelParameters, ONLY : EqCoeffTref, RstarMB

	IMPLICIT NONE

	REAL*8  :: Temperature, TempRatio
	INTEGER :: OrganicDissolutionRxn


        IF(OrganicDissolutionData(OrganicDissolutionRxn,2) .LT. 0.005) THEN
	    !Boiling point divided by Temperature
	    TempRatio = OrganicDissolutionData(OrganicDissolutionRxn,5)/Temperature

	    !Units of mbar, see Myrdal and Yalkowski and Schwarzenbach et al.,2003
	    !"Environmental Organic Chemistry"
            SaturationVaporPressure  = 1000*EXP ( &
                 -1.0*(21.2 + 0.3*OrganicDissolutionData(OrganicDissolutionRxn,6) &
                 + 177.0*OrganicDissolutionData(OrganicDissolutionRxn,7) ) &
                 *(TempRatio - 1.0) + (10.8 + 0.25*OrganicDissolutionData(OrganicDissolutionRxn,6)) &
                 *LOG(TempRatio) ) / OrganicDissolutionData(OrganicDissolutionRxn,8)
        ELSE
        
!            SaturationVaporPressure = OrganicDissolutionData(OrganicDissolutionRxn,14)&
!                 *(300.0/Temperature) & 
!                 *exp(OrganicDissolutionData(OrganicDissolutionRxn,15)*(1000.0/8.314) &
!                 *(1.0/300.0-1.0/Temperature))
            SaturationVaporPressure = OrganicDissolutionData(OrganicDissolutionRxn,14)&
                 *exp(OrganicDissolutionData(OrganicDissolutionRxn,15)*(1000.0/8.3145) &
                 *(1.0/300.0-1.0/Temperature))
!            WRITE(*,*) RstarMB, SaturationVaporPressure, OrganicDissolutionRxn
!            IF(OrganicDissolutionRxn .EQ. 21) STOP
        END IF   

	RETURN
END FUNCTION SaturationVaporPressure 

!Moles per particle of hydrophobic organics
!Note that this value is forced to be above 1.0e-25 to
!keep hydrophobic equilibrium constants from dividing by ~0, and
!giving an error in LSODES
REAL*8 FUNCTION TotalMolesOM(InParticle) 

	USE Chemistry,				ONLY :	HowManyOrgChems
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	Total = 0.0
	!Skip I = 1 since that is hard-wired for BC
	DO I = 2, HowManyOrgChems
		Total = Total + InParticle%OrgChems(I)
	END DO

	!Force total moles of OM to be above 1.0e-25 moles/particle
	!(prevents zeros in equilibrium coefficients)
	IF(Total .GE. 1.0E-25) THEN
		TotalMolesOM = Total
	ELSE
		TotalMolesOM = 1.0E-25
	END IF
	

	RETURN

END FUNCTION TotalMolesOM

!! This subroutine provides the link between the main model and the 
!! UNIFAC hydrophobic driver file soa_unidriv.f
SUBROUTINE UpdateHydrophobicUNIFAC(InParticle, Temperature)

        USE InfrastructuralCode, ONLY: ERROR

	USE Chemistry,		ONLY :	HowManyOrgChems, OrgActivityFlag, OrgActivityIndex

	IMPLICIT NONE

	INTEGER :: I, N, Q, status
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: TotalOM, Temperature, TotalOM_unifac
	


	N = 23 !Number of compounds in UNIDRIV and unifacparams.h
	ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1), STAT = status)
	IF (status > 0) THEN
		  CALL ERROR("Allocation error for XPASS and GAMMA in " &
		   //" UpdateHydrophobicUNIFAC of HydrophobicCondensationFunctions.h")
	ENDIF

        DO Q = 0, N-1
           XPASS(Q) = 0.0
	ENDDO
	!Calculate mole fractions, putting them in an array ordered from 0:N-1
	TotalOM = TotalMolesOM(InParticle)
		
	!If minimum organic is present, set activity coefficients to 1.0
	IF(TotalOM .LE. 1.0E-25) THEN
		DO I = 2, HowManyOrgChems !Skip BC
	 	    InParticle%HydrophobicActivityCoeffs(I) = 1.0
		END DO
	        RETURN
        ELSE !Calculate Activity Coefficients
                !WRITE(*,*) "Start UNIFAC"	
                DO I = 2, HowManyOrgChems 
                   IF (OrgActivityFlag(I) .LT. 0.005) THEN !UNIFAC
                       Q = NINT(OrgActivityIndex(I)-1.0) !Index in XPASS and GAMMA, starting with 0
                       !WRITE(*,*) Q
                       XPASS(Q) = InParticle%OrgChems(I)
                       !WRITE(*,*) XPASS(q)
                   ELSE !Constant Activity Coefficients
                          InParticle%HydrophobicActivityCoeffs(I) = OrgActivityIndex(I) 
                   ENDIF
                ENDDO  
                 
                !Normalize XPASS based only on UNIFAC compounds
                TotalOM_unifac = 0.0
                !WRITE(*,*) "TotalOM_unifac: ", TotalOM_unifac
                DO Q = 0, N-1
                    !WRITE(*,*) XPASS(q)
                    TotalOM_unifac = TotalOM_unifac+XPASS(Q)
                ENDDO
                !WRITE(*,*) "TotalOM_unifac: ", TotalOM_unifac

                !!MJA, 20170922, protect against zero OM from UNIFAC compounds 
                !!by setting activity coefficient to 1.0 in this case
                !!Note TotalOM_unifac is in mol/particle at this point
                IF (TotalOM_unifac .lt. 1e-40) THEN
                    DO Q = 0, N-1
                        GAMMA(Q) = 1.0
                    ENDDO
                ELSE
                    DO Q = 0, N-1
                        XPASS(Q) = XPASS(Q)/TotalOM_unifac
                    ENDDO

		    !Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
                    !WRITE(*,*) "Before UNIDRIV", XPASS, TotalOM_unifac
		    CALL UNIDRIV(XPASS,GAMMA,N,Temperature)
                    !WRITE(*,*) "After UNIDRIV"  
                ENDIF
	
		!Place UNIFAC activity coefficients into particle array
		DO I = 2, HowManyOrgChems !Skip BC
		    IF (OrgActivityFlag(I) .LT. 0.005) THEN !UNIFAC	
                        Q = NINT(OrgActivityIndex(I)-1.0)
                        InParticle%HydrophobicActivityCoeffs(I) = GAMMA(Q)
                    ENDIF
		END DO
	ENDIF
!WRITE(*,*) "Done!"
        DEALLOCATE(XPASS, GAMMA, STAT = status)
	RETURN
END SUBROUTINE UpdateHydrophobicUNIFAC


SUBROUTINE EquilibrateHydrophobicPhaseAtGridPoint()

	USE Aerosols,    ONLY : SortAerosolAtGridPointForCondensation, &
				Particles,			&
				Particle, &
				DumpParticleContentsAtError
	
	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				GasMolecularMass
	
	USE GridPointFields, ONLY : GetTemp,	&
			            GetGridCellChemBurden,	&
				    ReplaceGridCellChemBurden

	USE ModelParameters,     ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					Avogadro

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR
	
	IMPLICIT NONE

	INTEGER ::  HowManyLinks
	REAL*8  :: S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
			S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
			S7MassFracs(8), S8MassFracs(4), S9MassFracs(4), &
			S10MassFracs(4), InitialBurden(10) 
	TYPE(Particle),POINTER :: Current

	LOGICAL :: Equilibrated, RxnEquilibrated
	REAL*8  :: GasPhaseBurden, Temperature,  II
	INTEGER :: I, J, K, L, M, Q, ReturnType

	
	!Calculate surrogate concentrations and 
	!mass fractions of organics contributing to surrogates, and
	!Store initial surrogate gas phase burdens in a dummy vector
	!CALL CalculateSurrogateGasConc(S1MassFracs, S2MassFracs, &
        !        S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	!S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)

	!! Sort the gridcell, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()

	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	DO WHILE (ASSOCIATED(Current))
           IF (Current%NumberOfParticles .GT. 0 ) &
                HowManyLinks = HowManyLinks + 1
           Current => Current%Next
	END DO

	Temperature = GetTemp()

	Equilibrated = .FALSE.
	K = 0

	!! The main loop calculates whether it's equilibrated
	DO WHILE (.NOT.Equilibrated)

		K = K + 1
		Equilibrated = .TRUE.

		!! Loop over the reactions
		DO I = 1, HowManyOrganicDissolutionReactions						 
			
			!WRITE (*,*) "Reaction: ", I
			L = 0
			RxnEquilibrated = .FALSE.
				
			GasPhaseBurden = GetGridCellChemBurden ( INT(OrganicDissolutionData(I,1)))
		!WRITE(*,*) GasPhaseBurden
			!! Loop over all of the particles until 
			!! this particular reaction is happy
			DO WHILE (.NOT.RxnEquilibrated)

				RxnEquilibrated = .TRUE.

				Current => Particles%First
				DO J = 1, HowManyLinks
									
					!Equilibrate one reaction
					!WRITE(*,*) "Rxn: ", I
					CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)							
					
					!! Make a more refined guess if still equilibrating late in the game
					!! (That is, force it to equilibrate much closer to 1 each time through)
					IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
						CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
					END IF

					!WRITE(*,*) "Past Equil"
					IF (ReturnType .NE. 2 ) THEN
						RxnEquilibrated = .FALSE.
						Equilibrated    = .FALSE.
					END IF

					!WRITE(*,*) "Call UNIFAC"
					!Recalculate activity coefficients
					CALL UpdateHydrophobicUNIFAC(Current, Temperature)
					!WRITE(*,*) "Past UNIFAC"

					L = L+1
					Current => Current%Next
				END DO

				IF (L .GT. AqThermoNumbAllowedIterations) THEN
					CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
					CALL ERROR ("EquilibrateHydrophobicPhaseAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to organic dissolution reaction "//TRIM(INT2STR(I))//" (As indexed "// &
								"as input in OrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophobic.txt")
				END IF
			END DO
			
			!! Replace the surrogate gas phase chemical concentrations
			IF (INT(OrganicDissolutionData(I,1)) .LE. HowManyEvolveGasChems)										&
				CALL ReplaceGridCellChemBurden ( INT(OrganicDissolutionData(I,1)), GasPhaseBurden)			
		END DO



		IF (K .GT. AqThermoNumbAllowedIterations) THEN
			CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
			CALL ERROR ("EquilibrateHydrophobicPhaseAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all of the hydrophobic dissolution "// &
					    "reactions and failed.  Usually this implies an "// &
						"unrealistic starting condition.  You may try to fix this by increasing AqThermoNumbAllowedIterations,"// &
						" the maximum  allowed iterations of this type of routine.")
		END IF

	END DO

	!WRITE(*,*) "Before org update"
	!Update Concentrations of individual organics
	!CALL UpdateGasPhaseOrganics( S1MassFracs, S2MassFracs, &
        !        S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	!S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	!WRITE(*,*) "aFTER org update"

	RETURN

END SUBROUTINE EquilibrateHydrophobicPhaseAtGridPoint




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- EQUILIBRATE ONE HYDROPHOBIC DISSOLUTION REACTION ---    !!
!!							       !!
!! Iterate a particular reaction until it reaches equilibrium. !!
!!							       !!
!! ReturnType: 1 is normal, 2 is "already was in eq"	       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateHydrophobicReaction (GasPhaseBurden, WhichRxn, &
            InParticle, Temperature, ReturnType, OptErrorTolerance)

	USE Chemistry,		ONLY :	AqEquilibriaList,		&
					HowManyOrgChems,		&
					OrgPhaseChemicalNames,	&
					HowManyEvolveGasChems

	USE ModelParameters,	ONLY :	AqThermoEquilibriumError,	&
					WaterContentPrecision,		&
					AqThermoNumbAllowedIterations,	&
					AerosolWaterEquilibriumRHThreshold,&
					Avogadro,&
					moles, grams, WaterMolecMass, micron 

	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,ONLY :	INT2STR, REAL2STR,		&
					Transcript,			&
					ERROR, WARN

	USE Aerosols,		ONLY :	DumpParticleContentsAtError

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, ActivityCoefficients,&
                   ErrorTolerance, DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	! CMB: floating-point equality
	if (inparticle%numberofparticles .le. 1.0e-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN


	ErrorTolerance = OptErrorTolerance
	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine generally follows the Mass Flux Iteration method reviewed!! 
!! in Jacobson, 1999 (the textbook) and reviewed earlier in              !!
!! Jacobson et al. 1996, and Villars 1959                                !!
!!                                                                       !!
!! Jacobson, M.Z., A. Tabazadeh, and R.P. Turco, Simulating equilibrium  !!
!!       within aerosols and nonequilibrium between gases and aerosols,  !!
!!       Journal of Geophysical Research, 101 (D4), 9079-9091, 1996.     !!
!!									 !!
!! Villars, D.S., A method of successive approximations for computing    !!
!!    combustion equilibria on a high speed digital computer,            !!
!!    Journal of Physical Chemistry, 63, 521-5, 1959.			 !!
!!									 !!
!! STEP 1: Find the most aberrant ratio of concentration to              !!
!!         stoicheometric coefficient  (cf., Jacobson 1999, p.497)	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check for already equilibrated (either by none of 
!! the species being present or explicit equilibrium)
!!
!! And then pre-calculate the Qnumer and Qdenomer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is only one reaction type:					     !!
!!									     !!
!! 8. Direct dissolution of a single species not inclorporating dissociation.!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Reaction Check 1", GasPhaseBurden, WhichRxn, Temperature
	EqRatio = HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) "Reaction Check 2"
	IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has condensed, but org-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The gas-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		ReturnType = 2
		RETURN
	END IF

	Qnumer = InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3)))
	Qdenom = GasPhaseBurden / Avogadro  / InParticle%NumberOfParticles
	!WRITE(*,*) "Reaction Check 3", Qnumer, Qdenom

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	
	!! Loop over Corrections until convergence if there is anything to do
	! CMB: floating-point inequality
	if (dabs(qnumer) .gt. 1.0-40 .or. dabs(qdenom) .gt. 1.0e-40) then
	!IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!

	InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) = InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) + dX
	IF (OrganicDissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems) &
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles

	!WRITE(*,*) InParticle%OrgChems(OrganicDissolutionData(WhichRxn,3)), GasPhaseBurden
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio

	!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
	EqRatio = HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) WhichRxn, EqRatio
	!WRITE(*,*) "Check 9"

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	!IF (WhichRxn .EQ. 9) WRITE(*,*) EqRatio, InParticle%OrgChems(OrganicDissolutionData(WhichRxn,3)), GasPhaseBurden

	!! Count the number of attempts before convergence
	I = I+1

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		RETURN
	END IF

	!Check if all organic has condensed, but org-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The gas-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		RETURN
	END IF
	
	!Error if no convergence
	IF (I .GT. 500) THEN
		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle-HydrophobicDissolutionError.txt", &
											 InDoSurfaceTension=.FALSE.,InDoWaterActivity=.FALSE., InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., &
                                                     InDensity=.TRUE.,  InIonicStrength=.TRUE.)
		WRITE(*,*) "Org Reaction #: ", WhichRxn, "Eq. Ratio: ", EqRatio
		!WRITE(*,*) GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles
		CALL ERROR("Apparent Non Convergence in EquilibrateHydrophobicReaction().  Trying to converge "// &
				   "towards equilibrium in hydrophobic dissolution reaction number "//TRIM(INT2STR(WhichRxn))// ".")

	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	RETURN
END SUBROUTINE EquilibrateHydrophobicReaction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check the Ratio between the Existing and Ideal Equilibrium   !!
!! Constants for hydrophobic reactions.                         !!
!! This is used by several other routines.		        !!
!!								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This should only be used in the context of EquilibrateHydrophobicReaction!  
!! Otherwise there is no check for lack of the species.
REAL*8 FUNCTION HydrophobicEquilibriumConstantsRatio (GasPhaseBurden, & 
                 WhichRxn, Temperature, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyOrgChems

	USE ModelParameters, ONLY : RStarMB, moles, grams, &
                                    WaterMolecMass, Avogadro

	USE InfrastructuralCode, ONLY : Error

	USE Aerosols, ONLY : OrgCurvatureCorrection

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReactionType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature, ActCoeff, OrgMolFrac

	!! Internal Variables
	REAL*8	:: Pvap, P1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is only one reaction type:					     !!
!!									     !!
!! 8. Direct dissolution of a single species not incorporating dissociation. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! See if the reaction is empty !!
	! CMB (AER, Inc): floating-point equality check
	if (inparticle%orgchems(int(organicdissolutionData(WhichRxn, 3))) .le. 1.0e-40 &
			.and. GasPhaseBurden .le. 1.0e-40) then
	!IF (InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3))) .EQ. 0. & 
    !        .AND. GasPhaseBurden .EQ. 0) THEN
		
		HydrophobicEquilibriumConstantsRatio = 1.
		RETURN
	END IF
	!WRITE(*,*) "Flag 2"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	!! Get the mole fraction of the condensed organic,
	!! the activity coefficient of the condensed organic,
	!! the saturation vapor pressure,
	!! and the partial pressure and the equilibrium at this temperature
	Pvap = SaturationVaporPressure (Temperature, WhichRxn)
	
	!WRITE(*,*) "Flag 3a"
	ActCoeff = InParticle%HydrophobicActivityCoeffs(INT(OrganicDissolutionData(WhichRxn,3)))
	!WRITE(*,*) WhichRxn, ActCoeff
	!WRITE(*,*) "Flag 3b"
	OrgMolFrac = InParticle%OrgChems(INT(OrganicDissolutionData(WhichRxn,3)))/TotalMolesOM(InParticle)
	!WRITE(*,*) "Flag 3c"
	P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro
	!WRITE(*,*) "Flag 3d"	
	!! If the solution is set to blow up, then establish a large artificial value.
	! CMB: float equality check
	if (gasphaseburden .le. 1.0e-40) then
	!IF (GasPhaseBurden .LE. 0.) THEN
		HydrophobicEquilibriumConstantsRatio = 1.e12
		RETURN
	END IF
	!WRITE(*,*) "Flag 4"

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!	
		
	HydrophobicEquilibriumConstantsRatio = OrgCurvatureCorrection(InParticle)*OrgMolFrac*ActCoeff*Pvap/P1
        !WRITE(*,*) WhichRxn, OrgMolFrac, Pvap, P1
	!IF(WhichRxn .EQ. 21) STOP			
	!IF(WhichRxn .EQ. 23) WRITE(*,*) "Phobic Ratio", HydrophobicEquilibriumConstantsRatio, GasPhaseBurden
	!WRITE(*,*) "Ratio", HydrophobicEquilibriumConstantsRatio
	!WRITE(*,*) "Flag 8"
	RETURN
END FUNCTION HydrophobicEquilibriumConstantsRatio


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! --- EQUILIBRATE ONE HYDROPHIlIC DISSOLUTION REACTION ---    !!
!!							       !!
!! Iterate a particular reaction until it reaches equilibrium. !!
!!							       !!
!! ReturnType: 1 is normal, 2 is "already was in eq"	       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateHydrophilicReaction (GasPhaseBurden, WhichRxn, & 
            InParticle, Temperature, ReturnType, OptErrorTolerance)

	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				HowManyOrgChems,		&
				OrgPhaseChemicalNames,	&
				HowManyEvolveGasChems

	USE ModelParameters,	ONLY :	AqThermoEquilibriumError,	&
					WaterContentPrecision,		&
					AqThermoNumbAllowedIterations,	&
					AerosolWaterEquilibriumRHThreshold,&
					Avogadro,			&
					moles, grams, WaterMolecMass, micron 

	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,	ONLY :	INT2STR, REAL2STR,	&
						Transcript,		&
						ERROR, WARN

	USE Aerosols,		ONLY :	DumpParticleContentsAtError

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, ActivityCoefficients, ErrorTolerance, DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	! CMB (AER, Inc): floating-point equality check
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN
	IF (InParticle%NumberOfParticles .le. 1.0e-40 .OR. InParticle%Dry) RETURN
	
	!! Set the error tolerance
	ErrorTolerance = AqThermoEquilibriumError / 2.

        EqRatio = HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .LT. 1e-28) THEN !There is no organic left              
                ReturnType = 2
		RETURN
	END IF

	!Check if all organic has condensed, but aq-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The aq-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-28) THEN !There is no organic left in gas
                ReturnType = 2
		RETURN
	END IF
	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine generally follows the Mass Flux Iteration method reviewed!! 
!! in Jacobson, 1999 (the textbook) and reviewed earlier in              !!
!! Jacobson et al. 1996, and Villars 1959                                !!
!!                                                                       !!
!! Jacobson, M.Z., A. Tabazadeh, and R.P. Turco, Simulating equilibrium  !!
!!       within aerosols and nonequilibrium between gases and aerosols,  !!
!!       Journal of Geophysical Research, 101 (D4), 9079-9091, 1996.     !!
!!									 !!
!! Villars, D.S., A method of successive approximations for computing    !!
!!    combustion equilibria on a high speed digital computer,            !!
!!    Journal of Physical Chemistry, 63, 521-5, 1959.			 !!
!!									 !!
!! STEP 1: Find the most aberrant ratio of concentration to              !!
!!         stoicheometric coefficient  (cf., Jacobson 1999, p.497)	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check for already equilibrated (either by none of 
!! the species being present or explicit equilibrium)
!!
!! And then pre-calculate the Qnumer and Qdenomer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is only one reaction type:					     !!
!!									     !!
!! 9. Direct dissolution of a single species not inclorporating dissociation.!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Reaction Check 1", GasPhaseBurden, WhichRxn, Temperature
	EqRatio = HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) "Reaction Check 2"
	IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .LT. 1e-30) THEN !There is no organic left
		ReturnType = 2
		RETURN
	END IF

	!Check if all organic has condensed, but aq-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The aq-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-30) THEN !There is no organic left in gas
		ReturnType = 2
		RETURN
	END IF


	Qnumer = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3)))
	Qdenom = GasPhaseBurden / Avogadro / InParticle%NumberOfParticles
	!WRITE(*,*) "Reaction Check 3"

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	!WRITE(*,*) "dx: ", dX, " Z: ", Z
	!! Loop over Corrections until convergence if there is anything to do
	IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!

	InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) + dX
	IF (AqOrganicDissolutionData(WhichRxn,1) .LE. HowManyEvolveGasChems) &
			GasPhaseBurden = GasPhaseBurden - dX * Avogadro * InParticle%NumberOfParticles

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio

	!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
	EqRatio = HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, WhichRxn, Temperature, InParticle)
	!WRITE(*,*) "Check 9"

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	!WRITE(*,*) EqRatio, Z, dX

	!! Count the number of attempts before convergence
	I = I+1

	!Check if all organic has evaporated, but gas-phase is still unsaturated
	IF (EqRatio .GT. 1.0 & !The gas-phase is sub-saturated
		.AND. InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .LT. 1e-28) THEN !There is no organic left
		RETURN
	END IF

	!Check if all organic has condensed, but aq-phase is still unsaturated
	IF (EqRatio .LT. 1.0 & !The aq-phase is sub-saturated
		.AND. GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles .LT. 1e-28) THEN !There is no organic left in gas
		RETURN
	END IF
	
	!Exit if no convergence after 200 iterations
	IF (I .GT. 200) THEN
		!WRITE(*,*) "Aq Org Reaction #: ", WhichRxn, "Eq. Ratio: ", EqRatio, InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3)))
		!WRITE(*,*) GasPhaseBurden/ Avogadro / InParticle%NumberOfParticles
		!CALL WARN("Apparent Non Convergence in EquilibrateHydrophilicReaction().  Trying to converge "// &
		!		   "towards equilibrium in hydrophilic dissolution reaction number "//TRIM(INT2STR(WhichRxn))// ".")
                ReturnType = 2
                RETURN
	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	!WRITE(*,*) "Check 12"
	RETURN
END SUBROUTINE EquilibrateHydrophilicReaction


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check the Ratio between the Existing and Ideal Equilibrium   !!
!! Constants for hydrophilic reactions.                         !!
!! This is used by several other routines.		        !!
!!								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This should only be used in the context of EquilibrateHydrophilicReaction!  
!! Otherwise there is no check for lack of the species.
REAL*8 FUNCTION HydrophilicEquilibriumConstantsRatio (GasPhaseBurden, & 
                  WhichRxn, Temperature, InParticle)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyAqOrgChems

	USE ModelParameters, ONLY : RStarMB, moles, grams, & 
                                    WaterMolecMass, Avogadro, ProtonIndex

	USE InfrastructuralCode, ONLY : Error
	USE Aerosols, ONLY : Molality, CurvatureCorrection

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichRxn, ReactionType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: GasPhaseBurden, Temperature, ActCoeff, OrgMolality, & 
                  ProtonMolality

	!! Internal Variables
	REAL*8	:: Henry,  P1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There is only one reaction type:					     !!
!!									     !!
!! 9. Direct dissolution of a single species not incorporating dissociation. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! See if the reaction is empty !!

	!WRITE(*,*) "Flag 1"
	!! If there is none of the chemical to deal with
	!WRITE(*,*) WhichRxn
	!WRITE(*,*) InParticle%AqOrgChems(AqOrganicDissolutionData(WhichRxn,3))
	!WRITE(*,*) GasPhaseBurden
	! CMB: floating-point equality checks
	IF (InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .le. 1.0e-40 .AND.	&
			GasPhaseBurden .le. 1.0e-40) THEN
	!IF (InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) .EQ. 0. .AND.	&
	!	GasPhaseBurden .EQ. 0) THEN
		
		HydrophilicEquilibriumConstantsRatio = 1.
		RETURN
	END IF
	!WRITE(*,*) "Flag 2"

	!! If the solution is set to blow up, then establish a large artificial value.
	! CMB: floating-point equality check
	if (GasPhaseBurden .le. 1.0e-40) then
	!IF (GasPhaseBurden .LE. 0.) THEN
		HydrophilicEquilibriumConstantsRatio = 1.e12
		RETURN
	END IF
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	
	!Calculate the proton molality
	!WRITE(*,*) ProtonIndex
	ProtonMolality = Molality (ProtonIndex, InParticle)
	
	!Calculate the effective Henry's Law Constant (mol/kg H2O/mbar)
	Henry = EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, WhichRxn)
	
	!Retrieve the hydrophilic activity coefficient (corrected for infinite dilution convention)
	ActCoeff = InParticle%HydrophilicActivityCoeffs(INT(AqOrganicDissolutionData(WhichRxn,3)))

	!Partial pressure in mbar
	P1 = GasPhaseBurden*RstarMB*Temperature/Avogadro

	!Calculate molality of species in aqueous phase
	OrgMolality = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichRxn,3))) / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!
		
	HydrophilicEquilibriumConstantsRatio = CurvatureCorrection(InParticle)*ActCoeff*OrgMolality/(P1*Henry)
				
	!IF(WhichRxn .EQ. 23) WRITE(*,*) "Ratio", HydrophilicEquilibriumConstantsRatio, GasPhaseBurden
	!WRITE(*,*) "Flag 8"
	RETURN
END FUNCTION HydrophilicEquilibriumConstantsRatio


!Calculates effective henry's law for hydrophilic organics 
!in units of mol/kg/mbar 
REAL*8 FUNCTION EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, AqOrganicDissolutionRxn)

	USE ModelParameters, ONLY : EqCoeffTref
	USE InfrastructuralCode, ONLY : Error, warn	! cmb add warn 

	IMPLICIT NONE

	REAL*8  :: Temperature, ProtonMolality, KA1, KA2, Henry, AcidFactor
	INTEGER :: AqOrganicDissolutionRxn, I
	
	I = AqOrganicDissolutionRxn
	
	!in mol/kg
	KA1 = AqOrganicDissolutionData(i,7)*EXP((AqOrganicDissolutionData(i,8)/0.001097) &
											*((1.0/Temperature)-(1.0/298.15)))
	!in mol/kg
	KA2 = AqOrganicDissolutionData(i,9)*EXP((AqOrganicDissolutionData(i,12)/0.001097) &
											*((1.0/Temperature)-(1.0/298.15)))
	
	!Correction to Henry's law due to acid dissociation
	AcidFactor = 1 + KA1/ProtonMolality + KA1*KA2/(ProtonMolality**2)
	
	! CMB (AER, Inc): Temporarily, let's comment this out to allow an aerosol run to finish
	!				  and replace it with a less-strict version
	IF(AcidFactor .LT. 1.0) then
		write(*,*) 'AcidFactor < 1; ProtonMolality = ', ProtonMolality
		! CMB: seg fault?
		write(*,*) "Error in EffectiveOrganicHenrysLaw - corrected value is less than uncorrected."
		!CALL ERROR("Error in EffectiveOrganicHenrysLaw - corrected value is less than uncorrected.")
	endif
	!if (acidFactor .lt. 1.0) acidFactor = 1.0

	!in (ug i/ ug H2O) / (ug i / m3 air)
	Henry = AqOrganicDissolutionData(i,5)
	
	!Conversion to units of mol i/ kg H20 / mbar
	Henry = (Henry*1e6)/8.314e-5/298.15
	
	!Temperature Dependence
	Henry = Henry*EXP((AqOrganicDissolutionData(i,6)/0.001097) &
						*((1.0/Temperature)-(1.0/298.15)))
		
	EffectiveOrganicHenrysLaw  = Henry*AcidFactor
	
	RETURN
END FUNCTION EffectiveOrganicHenrysLaw

!! This subroutine provides the link between the main model and the 
!! UNIFAC hydrophilic driver file soa_unidriva.f
!! Note that the order of the hydrophilic organic compounds in 
!! unifacparama.h and in HydrophilicOrgChems.in must be the same! 
SUBROUTINE UpdateHydrophilicUNIFAC(InParticle, Temperature)

        USE InfrastructuralCode, ONLY: ERROR
	
        USE Chemistry,		ONLY :	HowManyAqOrgChems, GammaInf, &
                                      AqOrgActivityFlag, AqOrgActivityIndex  

	IMPLICIT NONE

	INTEGER :: I, N, Q, status
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: TotalMoles, Temperature
		
	N = 24 !Number of compounds
	ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1), STAT =status)
	IF (status > 0) THEN
		  CALL ERROR("Allocation error for XPASS and GAMMA in " &
		   //" UpdateHydrophobicUNIFAC of HydrophobicCondensationFunctions.h")
	ENDIF

        DO Q = 0, N-1
           XPASS(Q) = 0.0
	ENDDO

        !Calculate Activity Coefficients
        !WRITE(*,*) "Start UNIFAC"	
        DO I = 1, HowManyAqOrgChems 
              IF (AqOrgActivityFlag(I) .LT. 0.005) THEN !UNIFAC
                  Q = NINT(AqOrgActivityIndex(I)-1.0) !Index in XPASS and GAMMA, starting with 0
                  !WRITE(*,*) Q
                  XPASS(Q) = InParticle%AqOrgChems(I)
                  !WRITE(*,*) XPASS(q)
              ELSE !Constant Activity Coefficients
                  InParticle%HydrophilicActivityCoeffs(I) = AqOrgActivityIndex(I) 
              ENDIF
        ENDDO  
                 
        !Normalize XPASS based only on UNIFAC compounds and water
        TotalMoles = InParticle%AqChems(1)

        DO Q = 0, N-2
             !WRITE(*,*) XPASS(q)
             TotalMoles = TotalMoles+XPASS(Q)
        ENDDO

        DO Q = 0, N-2
             XPASS(Q) = XPASS(Q)/TotalMoles
        ENDDO

	!Water is the last compound (index 23)
	XPASS(N-1) = InParticle%AqChems(1)/TotalMoles

	!Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
				!WRITE(*,*) "Before UNIDRIVA", XPASS 
                                !WRITE(*,*)  TotalMoles, InParticle%AqChems(1), InParticle%Dry
	CALL UNIDRIVA(XPASS,GAMMA,N,Temperature)
				!WRITE(*,*) "After UNIDRIVA"
	!Place UNIFAC activity coefficients into particle array
	DO I = 1, HowManyAqOrgChems
	   IF (AqOrgActivityFlag(I) .LT. 0.005) THEN !UNIFAC	
               Q = NINT(AqOrgActivityIndex(I)-1.0)
               InParticle%HydrophilicActivityCoeffs(I) = GAMMA(Q)/GammaInf(I)
           ENDIF
           !WRITE(*,*) "Henry's law Activity ",I, InParticle%HydrophilicActivityCoeffs(I) 
        END DO
        
        DEALLOCATE(XPASS, GAMMA, STAT = status)
	RETURN

END SUBROUTINE UpdateHydrophilicUNIFAC

SUBROUTINE EquilibrateAllOrganicsAtGridPoint()

	USE Aerosols,    ONLY : SortAerosolAtGridPointForCondensation, &
				Particles,				&
				Particle, &
				DumpParticleContentsAtError
	
	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				GasMolecularMass
	
	USE GridPointFields, ONLY : GetTemp,	&
				    GetGridCellChemBurden,		&
				    ReplaceGridCellChemBurden

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
					AqThermoEquilibriumError, &
					ThermoBulkMode, &
					Avogadro

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR
	
	IMPLICIT NONE

	INTEGER ::  HowManyLinks
	REAL*8  :: S1MassFracs(2), S2MassFracs(6), S3MassFracs(2), &
			S4MassFracs(3), S5MassFracs(2), S6MassFracs(3), &
			S7MassFracs(8), S8MassFracs(4), S9MassFracs(4), &
			S10MassFracs(4), InitialBurden(10) 
	TYPE(Particle),POINTER :: Current

	LOGICAL :: Equilibrated, RxnEquilibrated
	REAL*8  :: GasPhaseBurden, Temperature, II
	INTEGER :: I, J, K, L, M, Q, ReturnType, AqRxnIndex

	
	!Calculate surrogate concentrations and 
	!mass fractions of organics contributing to surrogates, and
	!Store initial surrogate gas phase burdens in a dummy vector
	!CALL CalculateSurrogateGasConc( S1MassFracs, S2MassFracs, &
        !        S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	!S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)

	!! Sort the gridcell, pushing all empty sections to the end
	CALL SortAerosolAtGridPointForCondensation ()

	!! Count the number of links
	HowManyLinks = 0
	Current => Particles%First
	DO WHILE (ASSOCIATED(Current))
           IF (Current%NumberOfParticles .GT. 0 ) &
                HowManyLinks = HowManyLinks + 1
           Current => Current%Next
	END DO

	Temperature = GetTemp()

	Equilibrated = .FALSE.
	K = 0

	!! The main loop calculates whether it's equilibrated
	DO WHILE (.NOT.Equilibrated)

		K = K + 1
		Equilibrated = .TRUE.

		!! Loop over the reactions
		DO I = 1, HowManyOrganicDissolutionReactions+HowManyAqOrganicDissolutionReactions
			
			!WRITE (*,*) "Reaction: ", I
			L = 0
			RxnEquilibrated = .FALSE.
				
			!Hydrophobic Reactions
			IF (I .LE. HowManyOrganicDissolutionReactions) THEN
			
				GasPhaseBurden = GetGridCellChemBurden ( INT(OrganicDissolutionData(I,1)))
				!WRITE(*,*) GasPhaseBurden
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
									
						!Equilibrate one reaction
						!WRITE(*,*) "Rxn: ", I
						CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)							
					
						!! Make a more refined guess if still equilibrating late in the game
						!! (That is, force it to equilibrate much closer to 1 each time through)
						IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
							CALL EquilibrateHydrophobicReaction (GasPhaseBurden, I, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
						END IF

						!WRITE(*,*) "Past Equil"
						IF (ReturnType .NE. 2 ) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!WRITE(*,*) "Call UNIFAC"
						!Recalculate activity coefficients
						CALL UpdateHydrophobicUNIFAC(Current, Temperature)
						!WRITE(*,*) "Past UNIFAC"

							L = L+1
					Current => Current%Next
					END DO

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
						CALL ERROR ("EquilibrateAllOrganicsAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to hydrophobic organic dissolution reaction "//TRIM(INT2STR(I))//" (As indexed "// &
								"as input in OrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophobic.txt")
					END IF
				END DO
			
				!! Replace the surrogate gas phase chemical concentrations
				IF (INT(OrganicDissolutionData(I,1)) .LE. HowManyEvolveGasChems)										&
					CALL ReplaceGridCellChemBurden ( INT(OrganicDissolutionData(I,1)), GasPhaseBurden)			
		
			!Hydrophilic Reactions
			ELSE
		
				AqRxnIndex = I - HowManyOrganicDissolutionReactions
				GasPhaseBurden = GetGridCellChemBurden ( INT(AqOrganicDissolutionData(AqRxnIndex,1)))
				!WRITE(*,*) AqRxnIndex
				!IF (AqRxnIndex .EQ. 1) WRITE(*,*) GasPhaseBurden
				!! Loop over all of the particles until 
				!! this particular reaction is happy
				DO WHILE (.NOT.RxnEquilibrated)

					RxnEquilibrated = .TRUE.

					Current => Particles%First
					DO J = 1, HowManyLinks
									
						!Equilibrate one reaction
						!WRITE(*,*) "Rxn: ", I, AqRxnIndex
						CALL EquilibrateHydrophilicReaction (GasPhaseBurden, AqRxnIndex, Current, Temperature, ReturnType, OptErrorTolerance=AqThermoEquilibriumError/2.)							
						!WRITE(*,*) "Equilibrate Okay"

						!! Make a more refined guess if still equilibrating late in the game
						!! (That is, force it to equilibrate much closer to 1 each time through)
						IF ((K .GT. 25 .OR. L .GT. 10) .AND. ReturnType .NE. 2) THEN
							CALL EquilibrateHydrophilicReaction (GasPhaseBurden, AqRxnIndex, Current, Temperature, ReturnType, &
																	 OptErrorTolerance=AqThermoEquilibriumError/10.)
						END IF

						!WRITE(*,*) "Past Equil"
						IF (ReturnType .NE. 2 ) THEN
							RxnEquilibrated = .FALSE.
							Equilibrated    = .FALSE.
						END IF

						!WRITE(*,*) "Call UNIFAC"
						!Recalculate activity coefficients
						CALL UpdateHydrophilicUNIFAC(Current, Temperature)
						!WRITE(*,*) "Past UNIFAC"

							L = L+1
					Current => Current%Next
					END DO

					IF (L .GT. AqThermoNumbAllowedIterations) THEN
						CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
						CALL ERROR ("EquilibrateAllOrganicsAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all "     // &
					            "of the particles (stored as "//TRIM(INT2STR(HowManyLinks))//" links) at a grid point "  // &
								"to hydrophilic organic dissolution reaction "//TRIM(INT2STR(AqRxnIndex))//" (As indexed "// &
								"as input in HydrophilicOrganicDissolution.in) and failed.  The contents of the first particle in that "   // &
								"section are printed as UnEquilableParticle-Hydrophilic.txt")
					END IF
				END DO
			
				!! Replace the surrogate gas phase chemical concentrations
				IF (INT(AqOrganicDissolutionData(AqRxnIndex,1)) .LE. HowManyEvolveGasChems)										&
					CALL ReplaceGridCellChemBurden ( INT(AqOrganicDissolutionData(AqRxnIndex,1)), GasPhaseBurden)			
				
			END IF
		END DO



		IF (K .GT. AqThermoNumbAllowedIterations) THEN
			CALL DumpParticleContentsAtError(Particles%First, "UnEquilableParticle-Hydrophobic.txt", InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE., InDoWaterActivity=.TRUE., &
                                                     InDensity=.TRUE., InDoSurfaceTension=.TRUE., InIonicStrength=.TRUE.) 
			CALL ERROR ("EquilibrateHydrophobicPhaseAtGridPoint() in HydrophobicCondensationFunctions.h tried to equilibrate all of the hydrophobic dissolution "// &
					    "reactions and failed.  Usually this implies an "// &
						"unrealistic starting condition.  You may try to fix this by increasing AqThermoNumbAllowedIterations,"// &
						" the maximum  allowed iterations of this type of routine.")
		END IF

	END DO

	!WRITE(*,*) "Before org update"
	!Update Gas Phase Concentrations of individual organics
	!CALL UpdateGasPhaseOrganics(S1MassFracs, S2MassFracs, &
        !        S3MassFracs, S4MassFracs, S5MassFracs, S6MassFracs, &
	!S7MassFracs, S8MassFracs, S9MassFracs, S10MassFracs, InitialBurden)
	!WRITE(*,*) "After org update"

	RETURN

END SUBROUTINE EquilibrateAllOrganicsAtGridPoint

REAL*8 function dev(xin, Q, RH, Temp)
   
   IMPLICIT NONE
   
   !External Variables
   real*8, INTENT(IN) :: xin, RH, Temp
   integer, INTENT(IN) :: Q

   !Internal Variables
   REAL*8 :: GAMMA(0:23), XPASS(0:23)
	 
	 XPASS(Q) = xin
	 XPASS(23) = 1.0 - XPASS(Q)
	 
	 CALL UNIDRIVA(XPASS,GAMMA,24,Temp)
	 
	 dev=GAMMA(23)*XPASS(23)-RH/100.
	 return

end function dev

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate until a single organic in two phases (org and aqueous)!! 
!! is in equilibrium.                                            !!
!!								 !!
!! ReturnType: 1 is normal, 2 is "already was in eq"		 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EquilibrateOrganicReaction (WhichAqRxn, WhichOrgRxn, &
               InParticle, Temperature, ReturnType, OptErrorTolerance)

        ! CB: Ifort complains when HowMany is included
	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				!HowManyOrgChems,		&
				OrgPhaseChemicalNames,	&
				HowManyEvolveGasChems

	USE ModelParameters
			
	USE GridPointFields,		ONLY : 	GetRelativeHumidity

	USE InfrastructuralCode,	ONLY :	INT2STR, REAL2STR,	&
						Transcript,	&
						ERROR, WARN

	USE Aerosols,		ONLY :	DumpParticleContentsAtError

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichAqRxn, WhichOrgRxn, ReturnType
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I
	REAL*8  :: II, Qdenom, Qnumer, Z, dX, EqRatio, &
                    ActivityCoefficients, ErrorTolerance, &
                    DissociatedScaling, DissociatedExp
	LOGICAL :: WaterRxn

	!! SCSCSCSCSC
	LOGICAL :: Scaffolding = .TRUE. ! .TRUE.

	!If there are no particles, dissolution isn't allowed into one phase
	!or org phase is empty, skip this step
	IF (InParticle%NumberOfParticles .EQ. 0.) RETURN
	IF (.NOT.(DoHydrophobicOrgDissolutionFlag) .OR. .NOT.(DoHydrophilicOrgDissolutionFlag)) RETURN
	IF (TotalMolesOM(InParticle) .LE. 0.0) RETURN 

	!! Set the error tolerance according to the 
	!! subroutine's input.
	ErrorTolerance = OptErrorTolerance

	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine generally follows the Mass Flux Iteration method reviewed!! 
!! in Jacobson, 1999 (the textbook) and reviewed earlier in              !!
!! Jacobson et al. 1996, and Villars 1959                                !!
!!                                                                       !!
!! Jacobson, M.Z., A. Tabazadeh, and R.P. Turco, Simulating equilibrium  !!
!!       within aerosols and nonequilibrium between gases and aerosols,  !!
!!       Journal of Geophysical Research, 101 (D4), 9079-9091, 1996.     !!
!!									 !!
!! Villars, D.S., A method of successive approximations for computing    !!
!!    combustion equilibria on a high speed digital computer,            !!
!!    Journal of Physical Chemistry, 63, 521-5, 1959.			 !!
!!									 !!
!! STEP 1: Find the most aberrant ratio of concentration to              !!
!!         stoicheometric coefficient  (cf., Jacobson 1999, p.497)	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check for already equilibrated (either by none of 
!! the species being present or explicit equilibrium)
!!
!! And then pre-calculate the Qnumer and Qdenomer

	!If none of the organic is present, return
	IF(InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn, 3))) .EQ. 0.0 &
	   .AND. InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn, 3))) .EQ. 0.0) THEN
		ReturnType = 2
		RETURN
	END IF

	!Check if reaction is already in equilibrium
	EqRatio = TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, InParticle, Temperature)

	IF  (ABS(EqRatio-1.) .LT. AqThermoEquilibriumError) THEN
		ReturnType = 2
		RETURN
	END IF

	!If the concentration in either phase is going to 0, and the
	!other phase is subsaturated, return.
	IF((InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .GT. 1.0) .OR. &
		(InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .LT. 1.0) ) THEN
		ReturnType = 2
		RETURN
	END IF

	Qnumer = InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3)))
	Qdenom = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3)))
	!WRITE(*,*) "Reaction Check 3", Qnumer, Qdenom

	!WRITE(*,*) "Check 4"
	!! I counts the number of iterations
	I = 1

	!! STEP 2: Calculate the Mass Step Sizes by which to Correct Concentrations
	Z  = (Qdenom + Qnumer) / 2.
	dX = Qdenom - Z				!! The Mass Flux Factor
	
	!! Loop over Corrections until convergence if there is anything to do
	IF (.NOT. (Qnumer .EQ. 0 .AND. Qdenom .EQ. 0)) THEN
	DO WHILE ((ABS(EqRatio-1) .GT. ErrorTolerance))

	!WRITE(*,*) "Check 5"
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 3: Adjust Molalities: !!
	!!							  !!
	!! -- Ions First			  !!

	InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) = InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) + dX
	InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) - dX
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! STEP 4: Recalculate Z and dX for a new iteration !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!WRITE(*,*) "Check 7"
	Z = 0.5*Z

	!! The recalculation is based on the equilibrium ratio

	!WRITE(*,*) "Check 8", GasPhaseBurden, WhichRxn
	EqRatio = TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, InParticle, Temperature)
	!WRITE(*,*) WhichAqRxn, EqRatio
	!WRITE(*,*) "Check 9"

	!WRITE(*,*) "Check 10"
	IF (EqRatio .GE. 1.+ErrorTolerance) dX = -1.*Z
	IF (EqRatio .LE. 1.-ErrorTolerance) dX = Z
	
	!If the concentration in either phase is going to 0, and the
	!other phase is subsaturated, return.
	IF((InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .GT. 1.0) .OR. &
		(InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) .LT. 1.0e-30 &
		.AND. EqRatio .LT. 1.0) ) THEN
		RETURN
	END IF
	
	!! Count the number of attempts before convergence
	I = I+1

	!Error if no convergence
	IF (I .GT. 500) THEN
		WRITE(*,*) EqRatio, InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3))), &
		            InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3)))
		CALL DumpParticleContentsAtError(InParticle, "UnEquilableParticle-HydrophobicDissolutionError.txt", &
											 InDoSurfaceTension=.FALSE.,InDoWaterActivity=.FALSE., InRelativeHumidity=.TRUE., &
                                                     InRadius=.TRUE., InpH=.TRUE.,  &
                                                     InDensity=.TRUE., InIonicStrength=.TRUE.)
		CALL ERROR("Apparent Non Convergence in EquilibrateOrganicReaction.  Trying to converge "// &
				   "towards equilibrium in hydrophilic internal reaction number "//TRIM(INT2STR(WhichAqRxn))// ".")

	END IF

	!WRITE(*,*) "Check 11"
	END DO ; END IF

	RETURN
END SUBROUTINE EquilibrateOrganicReaction

REAL*8 FUNCTION TwoPhaseOrgEqRatio(WhichAqRxn, WhichOrgRxn, & 
                                   InParticle, Temperature)

	USE Chemistry,       ONLY : AqEquilibriaList, HowManyOrgChems

	USE ModelParameters, ONLY : RStarMB, moles, grams, WaterMolecMass, &
                                    Avogadro, &
				    ProtonIndex

	USE InfrastructuralCode, ONLY : Error

	USE Aerosols, ONLY : Molality, CurvatureCorrection, &
                             OrgCurvatureCorrection

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichAqRxn, WhichOrgRxn
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: Temperature

	!! Internal Variables
	REAL*8	:: Pvap, OrgActCoeff, OrgMolFrac, ProtonMolality, &
				Henry, AqOrgMolality, AqActCoeff
	

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Prepare the needed quantities !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!WRITE(*,*) "Flag 3"
	!! Get the mole fraction of the condensed organic,
	!! the activity coefficient of the condensed organic,
	!! the saturation vapor pressure,
	!! and the partial pressure and the equilibrium at this temperature
	Pvap = SaturationVaporPressure (Temperature, WhichOrgRxn)
	
	OrgActCoeff = InParticle%HydrophobicActivityCoeffs(INT(OrganicDissolutionData(WhichOrgRxn,3)))

	OrgMolFrac = InParticle%OrgChems(INT(OrganicDissolutionData(WhichOrgRxn,3)))/TotalMolesOM(InParticle)

	!Calculate the proton molality
	!WRITE(*,*) "ProtonIndex", ProtonIndex

	ProtonMolality = Molality (ProtonIndex, InParticle)
	
	!Calculate the effective Henry's Law Constant (mol/kg H2O/mbar)
	Henry = EffectiveOrganicHenrysLaw  (Temperature, ProtonMolality, WhichAqRxn)
	
	!Retrieve the hydrophilic activity coefficient (corrected for infinite dilution convention)
	AqActCoeff = InParticle%HydrophilicActivityCoeffs(INT(AqOrganicDissolutionData(WhichAqRxn,3)))

	!Calculate molality of species in aqueous phase
	AqOrgMolality = InParticle%AqOrgChems(INT(AqOrganicDissolutionData(WhichAqRxn,3))) / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	!Prevents a return of infinity if there is no organic in aqueous phase
	IF(AqOrgMolality .EQ. 0.0) THEN
		TwoPhaseOrgEqRatio = 1.0e10
		RETURN
	END IF
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make the primary calculation !!	
		
	TwoPhaseOrgEqRatio = OrgCurvatureCorrection(InParticle)*OrgMolFrac*OrgActCoeff*Pvap*Henry/(CurvatureCorrection(InParticle)*AqActCoeff*AqOrgMolality)
				
	!WRITE(*,*) "Ratio", TwoPhaseOrgEqRatio
	!WRITE(*,*) "Flag 8"
	RETURN


END FUNCTION TwoPhaseOrgEqRatio

SUBROUTINE EquilibrateOrganicParticle (InParticle, Temperature, &
                                      ReturnType, OptErrorTolerance)

        ! CB: Ifort complains when HowMany is included
	USE Chemistry,	ONLY :	AqEquilibriaList,		&
				!HowManyOrgChems,		&
				OrgPhaseChemicalNames,	&
				HowManyEvolveGasChems

	USE ModelParameters
			
	USE GridPointFields,	ONLY :  GetRelativeHumidity

	USE InfrastructuralCode,ONLY :	INT2STR, REAL2STR,	&
					Transcript,	&
					ERROR, WARN

	USE Aerosols,		ONLY :	DumpParticleContentsAtError

	IMPLICIT NONE

	!! External Variables
	INTEGER :: ReturnType, allocation_error
	TYPE(PARTICLE),POINTER :: InParticle
	REAL*8 :: Temperature
	REAL*8 :: OptErrorTolerance

	!! Local Variables
	INTEGER :: I, K, C, InternalReturnType
	REAL*8  :: II, ErrorTolerance
	INTEGER, ALLOCATABLE :: OrgIndex(:)

	!! SCSCSCSCSC
	LOGICAL :: Eq, Scaffolding = .TRUE. ! .TRUE.

	ReturnType = 2
	!If there are no particles, dissolution isn't allowed into one phase
	!or org phase is empty, skip this step
	IF (InParticle%NumberOfParticles .EQ. 0.) RETURN
    ! CMB (AER): Make conditional play nice with gfortran
	IF (.NOT.(DoHydrophobicOrgDissolutionFlag) .OR. .NOT.(DoHydrophilicOrgDissolutionFlag)) RETURN
	!WRITE(*,*) "Total OM in EquilibrateOrganicParticle ", TotalMolesOM(InParticle)
	IF (TotalMolesOM(InParticle) .LE. 0.0) RETURN 

	!! Set the error tolerance according to the input
	ErrorTolerance = OptErrorTolerance
	
	!! Presume innocence in the form of a successful return.
	ReturnType = 1

	!! Allocate the LSODES working arrays
	ALLOCATE (OrgIndex(HowManyAqOrganicDissolutionReactions), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of OrgIndex could not proceed in EquilibrateOrganicParticle")

	
	!Find matching hydrophobic reaction for hydrophilic reaction
	DO I = 1, HowManyAqOrganicDissolutionReactions
		OrgIndex(I) = 0
		DO C = 1, HowManyOrganicDissolutionReactions
			IF(OrganicDissolutionData(C,1) .EQ. AqOrganicDissolutionData(I,1)) THEN
				OrgIndex(I) = C
			END IF
		END DO
	END DO

	
	K = 0
	Eq = .FALSE.

	DO WHILE(.NOT.(Eq))

		ReturnType = 2
		Eq = .TRUE.
		DO I = 1, HowManyAqOrganicDissolutionReactions
			CALL EquilibrateOrganicReaction (I, OrgIndex(I), InParticle, Temperature, InternalReturnType, ErrorTolerance)
			IF(InternalReturnType .NE. 2) Eq = .FALSE.
			IF(K .EQ. 0 .AND. InternalReturnType .NE. 2) ReturnType = 1
		END DO
	
	
	!! Count the number of attempts before convergence
	K = K+1

	!Error if no convergence
	!Actually, just exit (MJA, 091107)
	IF (K .GT. 500) THEN
		Eq = .TRUE.
	END IF

	END DO
        DEALLOCATE(OrgIndex, STAT = allocation_error)
	RETURN
END SUBROUTINE EquilibrateOrganicParticle

REAL*8 FUNCTION BinaryWaterActivityOrganic(InParticle, WhichRxn, Temp)
!! This function calculates the water activity for a binary solution 
!! of aqueous organic and water (using UNIFAC) when the organic 
!! is at the same 
!! molality as in the complex solution. This is later used to 
!! calculate the total water activity of the 
!! organic-inorganic solution.
!! Created by Matt Alvarado, 09/08/2006
        USE InfrastructuralCode, ONLY: ERROR

	IMPLICIT NONE

	!Input Variables
	TYPE(PARTICLE),POINTER :: InParticle
	INTEGER :: WhichRxn
	REAL*8 :: Temp !Temperature

	!Internal Variables
	INTEGER :: I, J, N, status
	REAL*8, ALLOCATABLE :: XPASS(:), GAMMA(:)
	REAL*8 :: xorg, Ratio

	!I is the index of the organic in AqOrgChems
	I = AqOrganicDissolutionData(WhichRxn,3)
	Ratio = InParticle%AqOrgChems(I) &
	           / InParticle%AqChems(1)

	!Mole fraction of organic in binary solution at same molality
	!as complex solution
	xorg = Ratio/(1 + Ratio)

	N = 24 !Number of compounds
	ALLOCATE(XPASS(0:N-1), GAMMA(0:N-1), STAT = status)
	IF (status > 0) THEN
		  CALL ERROR("Allocation error for XPASS and GAMMA in " &
		   //" BinaryWaterActivityOrganic of HydrophobicCondensationFunctions.h")
	ENDIF
	!Assign mole fractions
	DO J = 1,N
		IF (J == I) THEN !Organic Compound
			XPASS(J-1) = xorg
		ELSE IF (J == N) THEN !Water
			XPASS(J-1) = 1.0 - xorg
		ELSE
			XPASS(J-1) = 0.0
		END IF
	END DO

	!Call UNIFAC, placing activity coefficients in GAMMA(0:N-1)
	CALL UNIDRIVA(XPASS,GAMMA,N,Temp)
		
	BinaryWaterActivityOrganic = GAMMA(N-1)*XPASS(N-1)

        DEALLOCATE(XPASS, GAMMA, STAT = status)
	RETURN
	
END FUNCTION BinaryWaterActivityOrganic

!!This just updates the UNIFAC activity coefficients
!!for all aprticles before a condensation step
SUBROUTINE UpdateHydrophobicUNIFAC_all

    USE GridPointFields, ONLY: GetTemp
    USE Aerosols,		 ONLY : Particle, Particles

    IMPLICIT NONE 
    
    TYPE(Particle), POINTER :: Current
    REAL*8 :: Temp

    Temp = GetTemp()

    current => particles%first
    DO WHILE(associated(current)) 
         CALL UpdateHydrophobicUNIFAC(current, Temp)
         current => current%next
    END DO
END SUBROUTINE UpdateHydrophobicUNIFAC_all

!!This just updates the UNIFAC activity coefficients
!!for all aprticles before a condensation step
SUBROUTINE UpdateHydrophilicUNIFAC_all

    USE GridPointFields, ONLY: GetTemp
    USE Aerosols,        ONLY : Particle, Particles

    IMPLICIT NONE 
    
    TYPE(Particle), POINTER :: Current
    REAL*8 :: Temp
    INTEGER :: Count

    Temp = GetTemp()

    current => particles%first
    Count = 1
    DO WHILE(associated(current)) 
         IF(current%AqChems(1) .GT. 0.0) CALL UpdateHydrophilicUNIFAC(current, Temp)
!         WRITE(*,*) Count
         current => current%next
         Count = Count+1
    END DO
END SUBROUTINE UpdateHydrophilicUNIFAC_all
