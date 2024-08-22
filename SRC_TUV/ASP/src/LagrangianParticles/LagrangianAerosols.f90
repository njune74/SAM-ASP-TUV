!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! LagrangianAerosols.f90
!! This is the main source file for particle structure, properties,
!! and aqueous-phase inorganic thermodynamics. It also includes the 
!! particle sorting routines needed for coagulation, condensation, 
!! and the moving-center distribution.
!! It mainly links all the particle header files together.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 09/07  2006   Matt Alvarado     Added array BoundaryParticles             !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 03/01  2012   Matt Alvarado     Added SetInitialDens to                   !!
!!                                  SUBROUTINE InitializeParticlesSectional  !!
!! 07/19  2013   Matt Alvarado     Moved IsNaN to InfrastructuralCode        !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES				!!
!! 1. ModelParameters			!!
!! 2. InfrastructuralCode		!!
!! 3. GridPointFields			!!
!! 4. Time				!!
!! 5. Chemistry				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES					            !!
!! 1. ParticleStructure.h					    !!
!! 2. InitializationAndSamplingRoutines.h			    !!
!! 3. ParticleAttributes.h					    !!
!! 4. Thermodynamics.h						    !!
!! 5. SortAerosol.h						    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE InitializeParticles()
!! 2. SUBROUTINE ClearParticles ()
!! 3. SUBROUTINE SCInitializeWaterParticles()
!! 4. SUBROUTINE ReInitializeParticles()
!! 5. FUNCTION Molality (ChemIndex, InParticle)
!! 6. FUNCTION CationMolality (ChemIndex, InParticle)
!! 7. FUNCTION AnionMolality (ChemIndex, InParticle)
!! 8. SUBROUTINE ChangeMolality (MolalityChange, ChemIndex, InParticle)
!! 9. SUBROUTINE ChangeCationMolality (MolalityChange, ChemIndex, InParticle)
!!10. SUBROUTINE ChangeAnionMolality (MolalityChange, ChemIndex, InParticle)
!!11. FUNCTION CopyLagrangianParticle (InParticle)
!!12. SUBROUTINE InitializeParticlesSectional()
!!13. SUBROUTINE DissociateAllElectrolytes()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE Aerosols

IMPLICIT NONE

PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Now specify which subroutines may be accessed by the outside world	 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PUBLIC	:: AerosolModeNames,			&
           AerosolModes,                        &
           AerosolPH,				&
           CurvatureCorrection,			&
           OrgCurvatureCorrection,		&
           DumpParticleContentsAtError,         &
           FindElectrolyteEquilibrium,	        &
           FindAqElectrolyteEquilibriumForGridPoint, &
           HypotheticalElectrolyteConcentrations, &
           InitializeParticlesSectional,        &
           KusikMeissner,			&
           Molality,				&
           ParticleKnudsen,			&
           Particle,				&
           Particles, ParticleArray,	        &
           ParticleKnudsenForEnergy,	        &   
           ParticleMass,			&
           ReynoldsNumber,			&
           ReformElectrolyte,			&
           RecalculateRadius,			&
           SortAerosolAtGridPointForCondensation, &
           SortAerosolAtGridPointForCoagulation, &
           TerminalVelocity,			&
           UpdateWaterActivity,                 &
           RegridAerosol,                       & 
           EquilibriumWaterContentAmount,       &
           ReadBoundaryDistributionsOrganic,    &
           BoundaryParticles,                   &
           SortBoundaryPart,                    &
           AerosolOpticalDepth,                 &
           ShellRefIndAndRad,                   &
           AerosolOptProp,                      &
           AerosolOptPropFASTTUV,                      &
           PopulateParticlesSectionsRightAway,  &
           ReadDistributionsOrganic,            &
           InputFlagRatioOrMass

INCLUDE "ParticleStructure.h"

!! The flag from the input file AerosolModes.in
INTEGER :: InputFlagRatioOrMass

!! The windspeed for the O'Dowd 1997 Salt Parameterization 
!! (0 means don't use it)
!REAL*8  :: ODowdSaltWindSpeed
!INTEGER :: ODowdSaltNumbBins

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Aerosol Modes array holds, in order (first index is mode identifier): !!
!!   Nt        (the number concentration of the mode per cc)		     !!
!!   Dbar_g_N  (the mean of the mode)					     !!
!!   Rho_g	   (the standard deviation)				     !!
!!   C(1...N)  (fractional mass concentration of all aq. chemicals,          !!
!!             using aqueous indexing)	                                     !!
!!   Corg(1..L) (fractional mass concentration of all org. chemicals,        !!
!!              using org. indexing)                                         !!
!!   IC(1..M)  (Number of each type of insoluble core aerosol per input)     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8, ALLOCATABLE :: AerosolModes(:,:)
INTEGER	            :: HowManyAerosolModes
LOGICAL             :: InsolCoresYesOrNo
CHARACTER (len=64), ALLOCATABLE  :: AerosolModeNames(:)

!!These are for the environmental (boundary) aerosol distribution
REAL*8, ALLOCATABLE :: EnvAerosolModes(:,:)
INTEGER		    :: HowManyEnvAerosolModes
CHARACTER (len=64), ALLOCATABLE  :: EnvAerosolModeNames(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Particles array is the linked list                   !!
!! of particles observable to the outside world.	    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TYPE (ParticleArray) :: Particles

!!BoundaryParticles is the environmental particle list
TYPE (ParticleArray) :: BoundaryParticles

!

CONTAINS

INCLUDE "InitializationAndSamplingRoutines.h"
INCLUDE "ParticleAttributes.h"
INCLUDE "Thermodynamics.h"
INCLUDE "SortAerosol.h"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Erase all particles so we can start again !!
SUBROUTINE ClearParticles ()
	
        ! CB: what happens if I comment this out?
	!USE Chemistry,       ONLY : HowManyAqChems, &
        !                            HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : HowManyBins
	USE InfrastructuralCode, ONLY : TRANSCRIPT, INT2STR

	IMPLICIT NONE

	INTEGER :: L
	TYPE(Particle), POINTER :: Leader, Current

	Current => Particles%First
	Leader  => Current%Next
	Particles%First => Null()

10	CALL TRANSCRIPT("Clearing particle "//TRIM(INT2STR(Current%ParticleID)))

	DEALLOCATE(Current%GammaMixed, stat=L)
	DEALLOCATE(Current%AqChems, stat=L)
	DEALLOCATE(Current%OrgChems, stat=L)
	DEALLOCATE(Current, stat=L)


	IF (ASSOCIATED(Leader)) THEN
		Current => Leader
		Leader  => Leader%Next
		GOTO 10
	END IF

	RETURN
END SUBROUTINE ClearParticles


!!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC!!
!! Make particles that fit a pre-specified set of sizes, all dry. !!
!!SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC!!
SUBROUTINE SCInitializeWaterParticles()

	USE InfrastructuralCode, ONLY : ERROR
	USE ModelParameters,     ONLY : micron

	IMPLICIT NONE

	!! Local Variables
	INTEGER :: Allocation_Error,I,J,K,HowManyRadii
	REAL*8  :: Radius(1000)
	TYPE(Particle), POINTER :: CurrentParticle



	Particles%First => Null()
	
	HowManyRadii = 2000

        DO I = 1, HowManyRadii
	    Radius(I) = 2.*micron/2.
        END DO


	CurrentParticle        => SCTailorMakeWaterParticle(Radius(1))
	Particles%First => CurrentParticle

	DO I = 2,HowManyRadii

		CurrentParticle%Next => SCTailorMakeWaterParticle(Radius(I))
		CurrentParticle      => CurrentParticle%Next

	END DO

	RETURN
END SUBROUTINE SCInitializeWaterParticles

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Molality is the unit of choice for calculating	  !!
!! anything that goes on in the aqueous phase.  But	  !!
!! the amount of each chemical in each particle is stored !!
!! in moles, to make conservation more exact, and so this !!
!! function converts moles to Molality.			  !!
!!							  !!
!! Molality is MOLES OF SOLUTE / KG of SOLVENT		  !!
!!							  !!
!! BUT THIS ISN'T SCALED BY MOLES, SO BE CAREFUL!         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! CMB (AER, Inc): None of these molality functions initialize the value of the 
!                 returned function, nor do they handle floating-point equality.
!                 Revisions made...
REAL*8 FUNCTION Molality (ChemIndex, InParticle)

	USE ModelParameters, ONLY : moles, grams,WaterMolecMass

	IMPLICIT NONE

	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    Molality = 0.0
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	Molality = InParticle%AqChems(ChemIndex) / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

END FUNCTION Molality

!! -- Ion molalities are calculated in a pass-along format
REAL*8 FUNCTION CationMolality (ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems

	IMPLICIT NONE

	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    CationMolality = 0.0
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	CationMolality = Molality (ChemIndex + HowManyAqChems, InParticle)

END FUNCTION CationMolality

!! -- Ion molalities are calculated in a pass-along format
REAL*8 FUNCTION AnionMolality (ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems, HowManyAqCations

	IMPLICIT NONE

	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    AnionMolality = 0.0
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	AnionMolality = Molality (ChemIndex + HowManyAqChems + HowManyAqCations, InParticle)

END FUNCTION AnionMolality

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This function adds amount MolalityChange to the molality	   !!
!! for a certain chemical in a particular particle.  Chemicals !!
!! are stored in terms of moles, so conversion is necessary,   !!
!! hence this subroutine.									   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ChangeMolality (MolalityChange, ChemIndex, InParticle)

	USE ModelParameters, ONLY : moles, WaterMolecMass

	IMPLICIT NONE

	REAL*8		   :: MolalityChange
	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	!! Change the molality by the requested amount.
	InParticle%AqChems(ChemIndex) = InParticle%AqChems(ChemIndex) + MolalityChange / 1000. *	&
									(InParticle%AqChems(1) / moles * WaterMolecMass)

END SUBROUTINE  ChangeMolality

!! -- Changes in Ion molalities are calculated in a pass-along format
SUBROUTINE ChangeCationMolality (MolalityChange, ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems

	IMPLICIT NONE

	REAL*8		   :: MolalityChange
	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

!print *, "Before Cation Molality:", MOLALITY(ChemIndex+ HowManyAqChems,InParticle)," + ",MolalityChange

	CALL ChangeMolality (MolalityChange, ChemIndex + HowManyAqChems, InParticle)

!print *, "After Cation Molality:", MOLALITY(ChemIndex+ HowManyAqChems,InParticle)


END SUBROUTINE ChangeCationMolality

!! -- Changes in Ion molalities are calculated in a pass-along format
SUBROUTINE ChangeAnionMolality (MolalityChange, ChemIndex, InParticle)

	USE ModelParameters, ONLY : HowManyAqChems, HowManyAqCations

	IMPLICIT NONE

	REAL*8		   :: MolalityChange
	INTEGER		   :: ChemIndex
	TYPE(Particle) :: InParticle

    ! CMB: mods for floating-point equality and function initialization
    if (inparticle%numberofparticles .le. 1.0E-40 .or. inparticle%dry) return
	!IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) RETURN

	CALL ChangeMolality (MolalityChange, ChemIndex + HowManyAqChems + HowManyAqCations, InParticle)

END SUBROUTINE ChangeAnionMolality 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Copy a Lagrangian Particle Exactly !!
FUNCTION CopyLagrangianParticle (InParticle)

	USE Chemistry,			 ONLY : HowManyAqChems,		&
						HowManyAqCations,	&
						HowManyAqAnions,	&
						AqMolecularMass,	&
						HowManyAqEqReactions

	USE InfrastructuralCode, ONLY : ERROR
	IMPLICIT NONE

	TYPE(Particle), POINTER :: CopyLagrangianParticle
	TYPE(Particle) :: InParticle

	INTEGER :: I, allocation_error

	ALLOCATE (CopyLagrangianParticle, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle Failed in CopyLagrangianParticle()")
	ALLOCATE (CopyLagrangianParticle%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%AqChems Failed in CopyLagrangianParticle()")
	ALLOCATE (CopyLagrangianParticle%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%GammaMixed Failed in CopyLagrangianParticle()")

	CopyLagrangianParticle%ParticleID = InParticle%ParticleID
	CopyLagrangianParticle%ParticleDistribution = InParticle%ParticleDistribution

	CopyLagrangianParticle%Sectional = InParticle%Sectional
	CopyLagrangianParticle%NumberOfParticles = InParticle%NumberOfParticles
	CopyLagrangianParticle%Edges(1) = InParticle%Edges(1)
	CopyLagrangianParticle%Edges(2) = InParticle%Edges(2)

	CopyLagrangianParticle%velocity(1) = InParticle%velocity(1)
	CopyLagrangianParticle%velocity(2) = InParticle%velocity(2)
	CopyLagrangianParticle%velocity(3) = InParticle%velocity(3)

	CopyLagrangianParticle%Temperature = InParticle%Temperature
	CopyLagrangianParticle%EmbryoRadius = InParticle%EmbryoRadius
	CopyLagrangianParticle%InsolubleRadius = InParticle%InsolubleRadius
	CopyLagrangianParticle%EffectiveRadius = InParticle%EffectiveRadius

	CopyLagrangianParticle%IonicStr = InParticle%IonicStr
	CopyLagrangianParticle%WaterActivity = InParticle%WaterActivity
	CopyLagrangianParticle%SolutionDensity = InParticle%SolutionDensity
	CopyLagrangianParticle%ParticleDensity = InParticle%ParticleDensity
	CopyLagrangianParticle%SurfaceTension = InParticle%SurfaceTension

	DO I = 1, HowManyAqEqReactions
	CopyLagrangianParticle%GammaMixed(I) = InParticle%GammaMixed(I)
	END DO

	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		CopyLagrangianParticle%AqChems(I) = InParticle%AqChems(I)
	END DO
	CopyLagrangianParticle%Next => Null()

	RETURN
END FUNCTION CopyLagrangianParticle

SUBROUTINE InitializeParticlesSectional()

	USE InfrastructuralCode, ONLY : ERROR
        USE ModelParameters,     ONLY : SetInitialDens

	IMPLICIT NONE

	!! Local Variables
	INTEGER :: Allocation_Error

	Particles%First => Null()
        CALL SetInitialDens(1.5000) !MJA 03-01-2012
        !Moved here so that the initial density can be set and reset at will

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Open the Particle Modes Input Deck and Read It In !!
	!! And then use it to populate the AerosolModes Array!!
	CALL ReadDistributionsOrganic() !!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Fill the Array with a Selection of Particles, !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!This puts all particles into sections
	Call PopulateParticlesSectionsRightAway(.FALSE.)
	RETURN
END SUBROUTINE InitializeParticlesSectional

SUBROUTINE DissociateAllElectrolytes()
	
	USE GridPointFields, ONLY : GetTemp,		        &
				    GetGridCellChemBurden,	&
				    ReplaceGridCellChemBurden,	&
				    GetRelativeHumidity

	USE Chemistry,	 ONLY : HowManyEvolveGasChems, &
				HowManyAqEqReactions, &
				AqEquilibriaList, &
				HowManyAqChems, &
				HowManyAqCations, &
				HowManyAqAnions

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterations,	&
				    AqThermoEquilibriumError, &
				    ThermoBulkMode, &
				    WaterContentPrecision

	USE Time, ONLY : CurrentTime, BeginningTime

	USE InfrastructuralCode, ONLY : ERROR, WARN, INT2STR, REAL2STR

	IMPLICIT NONE

	!! Internal Variables
	REAL*8 :: II
	INTEGER :: I, J, K, L, M, Q, HowManyLinks, ReturnType, OrgRxnIndex, AqRxnIndex, StoreReturn
	TYPE(Particle),POINTER :: Current	
	
!!Matt Alvarado added this to dissociate all solid electrolytes at 
!!initialization of the particles, creating a metastable aerosol initially 
!!that will then dry to equilibrium. 

	Current => Particles%First
	IF (ASSOCIATED(Current)) THEN
20		DO I = 1, HowManyAqEqReactions
			IF (AqEquilibriaList(I,6) .EQ. 1.0 ) THEN 
                        !!Is a solid Eq. Rxn. 
		
			!! The change amount, two cases: 
                        !! equations with and without water
				IF (AqEquilibriaList(I,8) .EQ. 0) THEN
				  II = Current%AqChems(INT(AqEquilibriaList(I,1)))
				ELSE
				  II = MIN(Current%AqChems(INT(AqEquilibriaList(I,1))),	& ! water limited, or
					Current%AqChems(1)/AqEquilibriaList(I,8))	  ! electrolyte limited?
				END IF		
				
				!! Add to ions:
				!! CATION:
				Current%AqChems(HowManyAqChems+INT(AqEquilibriaList(I,2)))	=	&
				Current%AqChems(HowManyAqChems+INT(AqEquilibriaList(I,2))) +	&				
				II * AqEquilibriaList(I,4)

				!! ANION:
				Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,3))) =	&
				Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,3))) +	&
				II * AqEquilibriaList(I,5)

				IF(AqEquilibriaList(I,22) .NE. 0) THEN
					!!Levitocite case, extra anion
					Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,22))) =	&
					Current%AqChems(HowManyAqChems+HowManyAqCations+INT(AqEquilibriaList(I,22))) +	&
					II * AqEquilibriaList(I,23)
				END IF
				

				! And Adjust WATER CONTENT appropriately
				Current%AqChems(1) = Current%AqChems(1) - II * AqEquilibriaList(I,8)

				!! -- Then the Electrolyte
				!IF (AqEquilibriaList(I,1) .GT. 0.) THEN					
				Current%AqChems(INT(AqEquilibriaList(I,1)))	= &
				Current%AqChems(INT(AqEquilibriaList(I,1))) - II

			END IF
		END DO !!Solid Dissociation loop
		
		!Set the particle to WET initially, 
                !so dissolution reactions will happen
		Current%Dry = .FALSE. 
	
		IF (ASSOCIATED(Current%Next)) THEN
			Current => Current%Next
			GOTO 20
		END IF
	END IF

END SUBROUTINE DissociateAllElectrolytes

END MODULE Aerosols
