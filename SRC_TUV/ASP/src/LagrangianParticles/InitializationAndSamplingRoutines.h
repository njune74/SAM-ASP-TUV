!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! InitializationAndSamplingRoutines.h
!! These routines are to do with the initialization	
!! of the aerosol particle distribution
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History		             !!
!! 07/24  2006   Matt Alvarado     Change "Pointer => NULL()" to	     !!
!!				   NULLIFY(POINTER) to fit pgf90	     !!
!! 09/07  2006   Matt Alvarado     Created PopBoundaryPart and		     !!
!!				       MakeBoundaryParticleExponential	     !!
!! 09/21  2006   Matt Alvarado     Created MakeBoundaryMonodisperse and      !!
!!				       MakeBoundaryBulk                      !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/29  2012   Matt Alvarado     Revised MakeParticle to allow particules 
!!                                  to have a specified amount of initial 
!!                                  water and to allow particles with only BC.
!! 05/29  2012   Matt Alvarado     Added option to 
!!                                    PopulateParticlesSectionsRightAway
!!                                    to follow Jacobson formulas for 
!!                                    initialization
!! 06/04  2012   Matt Alvarado     Modified the Jacobson init in 
!!                                    PopulateParticlesSectionsRightAway
!!                                    to ignore smallest bin
!! 01/24  2013   Matt Alvarado     Modified the Jacobson init in 
!!                                    PopulateParticlesSectionsRightAway
!!                                    to set smallest bin number and conc to 0
!! 01/25  2013   Matt Alvarado     Added Particle%Dry = .FALSE. to 
!!                                    MakeMonodisperse
!! 01/28  2013   Matt Alvarado     Removed calls to EquilibriumWaterContent
!! 01/29  2013   Matt Alvarado     Added Boundary variable to MakeBulk and 
!!                                    removed MakeBoundaryBulk
!!                           Added Boundary variable to MakeMonodisperse and 
!!                                    removed MakeBoundaryMonodisperse       !!
!!                           Added Boundary variable to MakeParticle and 
!!                                    removed MakeBoundaryParticle
!!                    Added Boundary variable to MakeParticleExponential and 
!!                                    removed MakeBoundaryParticleExponential
!!          Added Boundary variable to PopulateParticlesSectionsRightAway and 
!!                                    removed PopBoundaryPart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. FUNCTION MakeSection(WhichPopulation, WhichBin)
!! 2. SUBROUTINE MakeParticle(WhichPopulation, ResultInt, &
!!      LagrangianOnly, StructuredDistribution, Boundary)
!! 3. FUNCTION MakeBulk(Boundary)
!! 4. SUBROUTINE PopulateParticlesSectionsRightAway()
!! 5. SUBROUTINE ReadDistributionsOrganic ()
!! 6. FUNCTION MakeMonodisperse(Boundary)
!! 7. SUBROUTINE MakeParticleExponential(WhichPopulation, Boundary)
!! 8. SUBROUTINE ReadBoundaryDistributionsOrganic ()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! Modified by CMB (AER): gfortran doesn't like string formatting...
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initializes and returns a sectional pseudo-particle following !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION MakeSection(WhichPopulation, WhichBin)

	! 3/30/2016 CMB (AER, Inc): In order to figure out why I am getting
	!							a massive conversion into Wall_NOx, let's
	!							change "0" for the initialization value
	!						    of floating-point quantities to 1e-30
	!							or some small value.
	
	USE GridPointFields,	 ONLY : GetTemp, GetPress

	USE Chemistry,		ONLY : HowManyAqChems,			&
				HowManyAqCations,			&
				HowManyAqAnions,			&
				HowManyAqEqReactions, &
				HowManyOrgChems, &
				HowManyAqOrgChems

	USE InfrastructuralCode, ONLY : ERROR

	USE ModelParameters,     ONLY : micron, BinEdges

	IMPLICIT NONE

	!! External Variables
	INTEGER :: WhichPopulation, WhichBin
	TYPE(Particle), POINTER :: MakeSection


	!! Local Variables
	INTEGER :: allocation_error, I,J
    
	! CMB add
	real*8 :: filler_value
	filler_value = 1.0e-30
	! end CMB add
	
	!! If particle is large enough to be tracked, 
	!! create the New Particle and allocate its arrays
	ALLOCATE (MakeSection, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection Failed in MakeSection()")
	ALLOCATE (MakeSection%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection%AqChems Failed in MakeSection()")
	ALLOCATE (MakeSection%OrgChems(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection%OrgChems Failed in MakeSection()")
	ALLOCATE (MakeSection%AqOrgChems(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection%AqOrgChems Failed in MakeSection()")
	
		
	!! Allocate the Appropriate Electrolytic Arrays
	ALLOCATE (MakeSection%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection%GammaMixed Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyAqEqReactions
		MakeSection%GammaMixed(I) = 0
	END DO

	!Allocate the hydrophobic acitivty coefficients arrays
	ALLOCATE (MakeSection%HydrophobicActivityCoeffs(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection%HydrophobicActivityCoeffs Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyOrgChems
		MakeSection%HydrophobicActivityCoeffs(I) = 1.0
	END DO

	!Allocate the hydrophilic acitivty coefficients arrays
	ALLOCATE (MakeSection%HydrophilicActivityCoeffs(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeSection%HydrophilicActivityCoeffs Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyAqOrgChems
		MakeSection%HydrophilicActivityCoeffs(I) = 1.0
	END DO
	
	MakeSection%IonicStr        = 0.
	MakeSection%EmbryoRadius	= 0. 
	MakeSection%InsolubleRadius = 0.
	MakeSection%EffectiveRadius = 0.
	!MakeSection%IonicStr        = filler_value
	!MakeSection%EmbryoRadius	= filler_value
	!MakeSection%InsolubleRadius = filler_value
	!MakeSection%EffectiveRadius = filler_value
	
	MakeSection%Dry  = .FALSE.
	NULLIFY(MakeSection%Next)
	!MakeSection%Next => Null()

	MakeSection%ParticleID           = 0.
	!MakeSection%ParticleID = filler_value
	MakeSection%ParticleDistribution = WhichPopulation
	MakeSection%NumberOfParticles    = 0.
	!MakeSection%NumberOfParticles = filler_value
	MakeSection%Sectional			 = .TRUE.

	MakeSection%Edges(1) = BinEdges(WhichBin,1)
	MakeSection%Edges(2) = BinEdges(WhichBin,2)

	MakeSection%SurfaceTension  = 0.
	MakeSection%ParticleDensity = 0.
	MakeSection%SolutionDensity = 0.
	MakeSection%InsolubleDensity = 0.
	!MakeSection%SurfaceTension  = filler_value
	!MakeSection%ParticleDensity = filler_value
	!MakeSection%SolutionDensity = filler_value
	!MakeSection%InsolubleDensity = filler_value
		
	!! For Now, the Initial Temperature of the Particle is that of the 
	!! Local Environment (1/2002)
	MakeSection%Temperature = GetTemp()
	MakeSection%OriginPressureLevel = GetPress()

	!! Use a fake initial water activity:
	MakeSection%WaterActivity = 1.

	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		MakeSection%AqChems(I) = 0.
		!MakeSection%AqChems(i) = filler_value
	END DO
	
	DO I = 1, HowManyOrgChems
		MakeSection%OrgChems(I) = 0.
		!MakeSection%OrgChems(i) = filler_value
	END DO
		
	DO I = 1, HowManyAqOrgChems
		MakeSection%AqOrgChems(I) = 0.
		!MakeSection%AqOrgChems(i) = filler_value
	END DO

	RETURN
END FUNCTION MakeSection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Create a New Particle From A Given !!
!! log-normal mode and Hang it in the !!
!! Particle Array Structure.	      !!
!!                                    !!
!! MakeLooseParticle() is a replica of!!
!! this which does not place the      !!
!! aerosol in the domain, rather      !!
!! returns it.  Any changes here need !!
!! to be replicated there.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MakeParticle(WhichPopulation, ResultInt, LagrangianOnly, StructuredDistribution, Boundary)
	
	! 3/30/2016 CMB (AER, Inc): In order to figure out why I am getting
	!							a massive conversion into Wall_NOx, let's
	!							change "0" for the initialization value
	!						    of floating-point quantities to 1e-30
	!							or some small value.
	
	USE GridPointFields,	 ONLY : GetTemp, GetPress

	USE ModelParameters,	 ONLY : LagrSectSizeCutoff,		&
					lengthscale,			&
					Pi, EulerMass,			&
					micron,	moles,			&
					AverageAerosolDensity,		&
					DomainX, DomainY, DomainZ,	&
					Monodisperse,			&
					BinEdges, HowManyBins,		&
					NumbInsolTypes,			&
					InsolRadii,			&
					InsolDensity,			&
					InsolContactAngle

	USE Chemistry,		ONLY : HowManyAqChems,			&
				       HowManyAqCations,		&
				       HowManyAqAnions,			&
				       AqMolecularMass,			&
				       HowManyAqEqReactions,            &
				       HowManyOrgChems,                 &
				       OrgMolecularMass,                &
				       HowManyAqOrgChems,               &
				       AqOrgMolecularMass

	USE InfrastructuralCode, ONLY : REAL2STR, INT2STR, Transcript,  &
                                        Warn, Error, GetParticleID

	IMPLICIT NONE

	!! Input Variables
	INTEGER :: WhichPopulation, ResultInt
	LOGICAL :: LagrangianOnly 
!! If this is true, then return always as lagrangian

	!! If we want to draw from a distribution in a structured way, use 
        !! this input
	!! and it will go this far up the cumulative distribution function.
	REAL*8 :: StructuredDistribution
        LOGICAL :: Boundary

	!! Local Variables
	REAL*8 :: II, JJ, KK, RR(1)
	TYPE(Particle), POINTER :: NewParticle, CurrentParticle
	INTEGER			:: allocation_error, I,J,X,Y,Z
	LOGICAL, PARAMETER	:: Scaffolding = .FALSE. !.TRUE.

    !!This holds the AerosolModes or EnvAerosolModes arrays, depending on the value of Boundary
    REAL*8, ALLOCATABLE :: ModesArray(:,:)
	
	! CMB add
	real*8 :: filler_value
	filler_value = 1.0e-30
	! end CMB add
	
	IF (SCAFFOLDING) CALL Transcript(">>>Entering MakeParticle()<<<")
	IF (SCAFFOLDING) CALL Transcript("")


        IF (.NOT.(BOUNDARY)) THEN
           ALLOCATE (ModesArray(HowManyAerosolModes,3+HowManyAqChems+HowManyOrgChems) , &
                stat = allocation_error)
        ELSE
           ALLOCATE (ModesArray(HowManyEnvAerosolModes,3+HowManyAqChems+HowManyOrgChems) , &
                stat = allocation_error)
        ENDIF

	if (allocation_error > 0) CALL ERROR("Allocation of ModesArray Failed in MakeParticle()")

        IF (.NOT.(BOUNDARY)) THEN
           ModesArray = AerosolModes
        ELSE
           ModesArray = EnvAerosolModes
        ENDIF          

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For sampling, we want to use a random number to find a number of standard !!
!! deviations away from the mean value.  To do this, we must invert the error!!
!! function.  To do this, we use the following numerical trick, which inverts!!
!!an error function for a series with a mean at the Euler-Masscheroni Constant!
!! and a standard deviation of pi/sqrt(6):			             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    ! CMB (AER): Rewrite this conditional to play nice with gfortran
    if (.not. monodisperse) then
!	IF (MONODISPERSE .EQ. .FALSE.) THEN

           CALL PPND (StructuredDistribution, I, II)

           II = EXP(LOG(ModesArray(WhichPopulation,3)) * II + LOG(ModesArray(WhichPopulation,2))) / 2.*micron

	ELSE
           II = ModesArray(WhichPopulation,2)/2.*micron
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! PUSH PARTICLE INTO THE SECTIONAL ARRAYS !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (II .LE. BinEdges(HowManyBins,2) .AND. .NOT.LagrangianOnly) THEN

           IF (.NOT.(Boundary)) THEN
              CurrentParticle => Particles%First
           ELSE
              CurrentParticle => BoundaryParticles%First
           ENDIF

		!! Find the appropriate section to put it in, or err
           IF (ASSOCIATED(CurrentParticle)) THEN

               ! CMB (AER): Rewrite this conditional for gfortran
 10            if (currentparticle%sectional .and. &
!10            IF ((CurrentParticle%Sectional .EQ. .TRUE.) .AND.	&
                   (II .GT. CurrentParticle%Edges(1)) .AND.		&
                    (II .LE. CurrentParticle%Edges(2)))	THEN
			
                 !! We've found the proper section
                 !! So add the chemicals

                 CurrentParticle%NumberOfParticles = CurrentParticle%NumberOfParticles + 1
                 !DO I = 2, HowManyAqChems
                 DO I = 1, HowManyAqChems !MJA, 02-28-2012
                    !WRITE(*,*) "MakeParticle: ", AverageAerosolDensity
                    CurrentParticle%AqChems(I) = (ModesArray(WhichPopulation,3+I) *	&
                         4./3.*Pi*II*II*II*AverageAerosolDensity/		&
                         AqMolecularMass(I) +				&
                         (CurrentParticle%NumberOfParticles-1.)*CurrentParticle%AqChems(I))/ &
                         CurrentParticle%NumberOfParticles
                 END DO
			
                 DO I = 1, HowManyOrgChems
                    CurrentParticle%OrgChems(I) = (ModesArray(WhichPopulation,3+HowManyAqChems+I) *	&
                         4./3.*Pi*II*II*II*AverageAerosolDensity/		&
                         OrgMolecularMass(I) +					&
                         (CurrentParticle%NumberOfParticles-1.)*CurrentParticle%OrgChems(I))/ &
                         CurrentParticle%NumberOfParticles 
                 END DO
			
                 !NOTE: We initialize with all of the organic compounds in 
                 !the organic phase, and none in the hydrophilic phase
                 DO I = 1, HowManyAqOrgChems
                    CurrentParticle%AqOrgChems(I) = 0.0
					!CurrentParticle%AqOrgChems(i) = filler_value
                 END DO
                        

              ELSE
                 IF(ASSOCIATED(CurrentParticle%Next)) THEN
                    CurrentParticle => CurrentParticle%Next
                    GOTO 10
                 ELSE
                    CALL ERROR ("Couldn't Find the appropriate section to put a too-small particle of mode "  &
                         //TRIM(INT2STR(WhichPopulation))//" in subroutine MakeParticle.")
                 END IF
              END IF
           END IF

           IF (SCAFFOLDING) CALL WARN("Particle too small")
           ResultInt = 0 
           RETURN
	END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! If particle is large enough to be tracked        !!
!!(Lagrangian, above section cutoff),               !!
!! create the New Particle and allocate its arrays  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	ALLOCATE (NewParticle, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle Failed in MakeParticle()")
	ALLOCATE (NewParticle%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%AqChems Failed in MakeParticle()")
	ALLOCATE (NewParticle%OrgChems(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%OrgChems Failed in MakeParticle()")
	ALLOCATE (NewParticle%AqOrgChems(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%AqOrgChems Failed in MakeParticle()")

	!! Allocate the Appropriate Electrolytic Arrays
	ALLOCATE (NewParticle%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%GammaMixed Failed in MakeParticle()")

	!! Initialize
	DO I = 1, HowManyAqEqReactions
           NewParticle%GammaMixed(I) = 0
	END DO

	!Allocate the hydrophobic acitivty coefficients arrays
	ALLOCATE (NewParticle%HydrophobicActivityCoeffs(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%HydrophobicActivityCoeffs Failed in MakeParticle()")

	!! Initialize
	DO I = 1, HowManyOrgChems
		NewParticle%HydrophobicActivityCoeffs(I) = 1.0
	END DO

	!Allocate the hydrophilic acitivty coefficients arrays
	ALLOCATE (NewParticle%HydrophilicActivityCoeffs(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of NewParticle%HydrophilicActivityCoeffs Failed in MakeParticle()")

	!! Initialize
	DO I = 1, HowManyAqOrgChems
           NewParticle%HydrophilicActivityCoeffs(I) = 1.0
	END DO

	NewParticle%IonicStr		= 1.
	NewParticle%InsolubleRadius = 0.
	!NewParticle%InsolubleRadius = filler_value
	NewParticle%EmbryoRadius	= II
	NewParticle%EffectiveRadius = II

	NewParticle%ParticleID           = GetParticleID()
	NewParticle%ParticleDistribution = WhichPopulation
	NewParticle%NumberOfParticles    = 1.
	NewParticle%Sectional			 = .FALSE.
	NewParticle%Edges(1)			 = 0.
	NewParticle%Edges(2)			 = 0.
	!NewParticle%Edges(1)			 = filler_value
	!NewParticle%Edges(2)			 = filler_value
		
	!! Budge the new particle in at the head of the list
        IF (.NOT.(BOUNDARY)) THEN	
            IF (ASSOCIATED(Particles%First)) THEN
		CurrentParticle		   => Particles%First

		!! Sectional Bins are at the head of the lists
        ! CMB (AER): Rewrite this conditional to play well with gfortran
 20     if (currentparticle%sectional .and. associated(currentparticle%next)) then
!20		IF(CurrentParticle%Sectional .EQ. .TRUE. .AND. ASSOCIATED(CurrentParticle%Next)) THEN
			CurrentParticle => CurrentParticle%Next
			GOTO 20	
		END IF

		IF (ASSOCIATED(CurrentParticle%Next)) THEN
			NewParticle%Next	 => CurrentParticle%Next
			CurrentParticle%Next => NewParticle
		ELSE
			NULLIFY(NewParticle%Next)
			!NewParticle%Next	 => Null()
			CurrentParticle%Next => NewParticle
		END IF
	    ELSE
		Particles%First => NewParticle
		NULLIFY(NewParticle%Next)
		!NewParticle%Next	   => Null()
	    END IF
        ELSE
            IF (ASSOCIATED(BoundaryParticles%First)) THEN
		CurrentParticle		   => BoundaryParticles%First

		!! Sectional Bins are at the head of the lists
        ! CMB (AER): Rewrite this conditional to play well with gfortran
 88     if (currentparticle%sectional .and. &
            associated(currentparticle%next)) then
! 88		IF(CurrentParticle%Sectional .EQ. .TRUE. .AND. ASSOCIATED(CurrentParticle%Next)) THEN
			CurrentParticle => CurrentParticle%Next
			GOTO 88	
		END IF

		IF (ASSOCIATED(CurrentParticle%Next)) THEN
			NewParticle%Next	 => CurrentParticle%Next
			CurrentParticle%Next => NewParticle
		ELSE
			NULLIFY(NewParticle%Next)
			!NewParticle%Next	 => Null()
			CurrentParticle%Next => NewParticle
		END IF
	    ELSE
		BoundaryParticles%First => NewParticle
		NULLIFY(NewParticle%Next)
		!NewParticle%Next	   => Null()
	    END IF
        ENDIF
	!! For Now, the Initial Temperature of the Particle is that of the 
	!! Local Environment
	NewParticle%Temperature = GetTemp()
	NewParticle%OriginPressureLevel  = GetPress()

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Create a Velocity for the Particle.  It is some random	!!
	!! percentage of the terminal velocity, randomly apportioned	!!
	!! between the three directions.				!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	CALL Random_Number(II)
	CALL Random_Number(JJ)
	CALL Random_Number(KK)

	!! Set initial velocity at half or more of terminal velocity
	KK = (1.+KK)/2.  ! * TerminalVelocity(NewParticle%Radius)

	!! Set initial angles
	II = 2.*Pi*II	 ! Flat Angle is out of 360 degrees
	JJ = Pi*JJ-Pi/2. ! Azimuthal angle is out of 180 degrees

	!! Change this into two random angles 
        !! and then move to cosines and sines from there:
	NewParticle%velocity(1) = KK*DCOS(II)*DCOS(JJ)
	NewParticle%velocity(2) = KK*DSIN(II)*DCOS(JJ)
	NewParticle%velocity(3) = KK*DSIN(JJ)

	IF (SCAFFOLDING) CALL Transcript("Effective Radius : ("//TRIM(REAL2STR(NewParticle%EffectiveRadius))//")")
	IF (SCAFFOLDING) CALL Transcript("Velocity (X,Y,Z): ("//TRIM(REAL2STR(NewParticle%velocity(1)))//","//			&
				TRIM(REAL2STR(NewParticle%velocity(2)))//","//TRIM(REAL2STR(NewParticle%velocity(3)))//")")

	II = 0.
	DO I = 2, HowManyAqChems
 
	
			!! Load in the chemistry assignments, stored as Moles.
			!! Calculate a Mass of each chemical type (using ModesArray, which is simply a mass fraction value)
			!! Assume, for this calculation, that the mean desity is AverageAerosolDensity (=1).
			NewParticle%AqChems(I) = ModesArray(WhichPopulation,3+I) *	&
						4./3.*Pi*(NewParticle%EffectiveRadius)**3.		&
						*AverageAerosolDensity/				&
						AqMolecularMass(I) !*moles
			II = II + NewParticle%AqChems(I)

	END DO

	DO I = 1, HowManyOrgChems
 
	
			!! Load in the chemistry assignments, stored as Moles.
			!! Calculate a Mass of each chemical type (using ModesArray, which is simply a mass fraction value)
			!! Assume, for this calculation, that the mean desity is AverageAerosolDensity (=1).
			NewParticle%OrgChems(I) = ModesArray(WhichPopulation,3+HowManyAqChems+I) *	&
						4./3.*Pi*(NewParticle%EffectiveRadius)**3.		&
						*AverageAerosolDensity/				&
						OrgMolecularMass(I) !*moles
	END DO

	DO I = 1, HowManyAqOrgChems
 
			!! Initialize all the organic compounds in the hydrophobic phase
			NewParticle%AqOrgChems(I) = 0.0
			!NewParticle%AqOrgChems(i) = filler_value
	END DO
	
	
	!! Use a fake initial water activity:
	NewParticle%WaterActivity = 1.

	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		NewParticle%AqChems(I) = 0.
		!NewParticle%AqChems(i) = filler_value
	END DO

	!! If didn't just return, then it can deal with water
	NewParticle%Dry = .FALSE.

	!! As a token amount of water to avoid initial numerical problems, set the 
	!! water mixing ratio at 0.5
	NewParticle%AqChems(1) = II*1.

	IF (SCAFFOLDING) CALL Transcript(">>>Exiting MakeParticle()<<<")
	
	ResultInt = 1
	RETURN

END SUBROUTINE MakeParticle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initializes and returns a sectional pseudo-particle following !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION MakeBulk(Boundary)

	USE GridPointFields,	 ONLY : GetTemp, GetPress

	USE Chemistry,		 ONLY : HowManyAqChems,			&
					HowManyAqCations,		&
					HowManyAqAnions,		&
					HowManyAqEqReactions,           &
					HowManyOrgChems,                &
					HowManyAqOrgChems

	USE InfrastructuralCode, ONLY : ERROR

	USE ModelParameters,     ONLY : micron,NumbInsolTypes

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle), POINTER :: MakeBulk
        LOGICAL :: Boundary !TRUE if boundary distribution, FALSE if main distribution


	!! Local Variables
	INTEGER :: allocation_error, I,J
	REAL*8  :: II

	!! If particle is large enough to be tracked, 
	!! create the New Particle and allocate its arrays
	ALLOCATE (MakeBulk, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk Failed in MakeBulk()")
	ALLOCATE (MakeBulk%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk%AqChems Failed in MakeBulk()")
	ALLOCATE (MakeBulk%OrgChems(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk%OrgChems Failed in MakeBulk()")
	ALLOCATE (MakeBulk%AqOrgChems(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk%AqOrgChems Failed in MakeBulk()")



	!! Allocate the Appropriate Electrolytic Arrays
	ALLOCATE (MakeBulk%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk%GammaMixed Failed in MakeBulk()")

	!! Initialize
	DO I = 1, HowManyAqEqReactions
		MakeBulk%GammaMixed(I) = 0
	END DO

	!Allocate the hydrophobic acitivty coefficients arrays
	ALLOCATE (MakeBulk%HydrophobicActivityCoeffs(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk%HydrophobicActivityCoeffs Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyOrgChems
		MakeBulk%HydrophobicActivityCoeffs(I) = 1.0
	END DO

	!Allocate the hydrophobic acitivty coefficients arrays
	ALLOCATE (MakeBulk%HydrophilicActivityCoeffs(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeBulk%HydrophilicActivityCoeffs Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyAqOrgChems
		MakeBulk%HydrophilicActivityCoeffs(I) = 1.0
	END DO


	MakeBulk%IonicStr = 1.
	MakeBulk%EmbryoRadius	 = 0.1*micron
	MakeBulk%InsolubleRadius = 0.
	MakeBulk%EffectiveRadius = 0.

	MakeBulk%ParticleID           = 0.
	MakeBulk%ParticleDistribution = 1
	MakeBulk%NumberOfParticles    = 1
	MakeBulk%Sectional	      = .TRUE.
	MakeBulk%Edges(1)	      = 0.
	MakeBulk%Edges(2)	      = 1.0e9 

	!! For Now, the Initial Temperature of the Particle is that of the 
	!! Local Environment (1/2002)
	MakeBulk%Temperature = GetTemp()
	MakeBulk%OriginPressureLevel  = GetPress()


	!! Use a fake initial water activity:
	MakeBulk%WaterActivity = 1.

	DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		MakeBulk%AqChems(I) = 0.
	END DO

	II = 0.
	DO I = 2, HowManyAqChems
		IF(.NOT.(Boundary)) THEN
                   MakeBulk%AqChems(I) = AerosolModes(1,3+I)
                ELSE
                   MakeBulk%AqChems(I) = EnvAerosolModes(1,3+I)
                ENDIF
		II = II + MakeBulk%AqChems(I)
	END DO

	DO I = 1, HowManyOrgChems
		IF(.NOT.(Boundary)) THEN
                    MakeBulk%OrgChems(I) = AerosolModes(1,3+HowManyAqChems+I)
                ELSE
                    MakeBulk%OrgChems(I) = EnvAerosolModes(1,3+HowManyAqChems+I)
                ENDIF
	END DO
	
	!Initialize all organic species in hydrophobic phase
	DO I = 1, HowManyAqOrgChems
		MakeBulk%AqOrgChems(I) = 0.0
	END DO

	
	DO I = 3+HowManyAqChems+HowManyOrgChems+1,3+HowManyAqChems+HowManyOrgChems+NumbInsolTypes
		IF (AerosolModes(1,I) .GT. 0) &
		CALL ERROR("MakeBulk() Failed.  You cannot specify insoluble cores to aerosol and then try to run bulk mode.  "// &
		           "This program is not set up to do that.  The simplest system with that includes insoluble cores is "// &
				   "the single-bin sectional representation.")
	END DO

	!! As a token amount of water to avoid initial numerical problems, set the 
	!! water mixing ratio at 0.5
	MakeBulk%AqChems(1) = II*1.
        MakeBulk%Dry = .FALSE.

	NULLIFY(MakeBulk%Next)
	!MakeBulk%Next => Null()

	RETURN

END FUNCTION MakeBulk

!! This is a routine for simple debugging purposes.  Only include it if absolutely need to.
INCLUDE "SCTailorMakeParticle.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Add each of the individual modes to the Particles() array     !!
!! NOTE: Creates a single, internally mixed size distribution    !!
!! as the sum of the input modes.				 !!       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PopulateParticlesSectionsRightAway (Boundary)

	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, Transcript, &
                                        SetParticleID, ERROR

	USE ModelParameters,     ONLY : lengthscale,			&
					DomainX, DomainY, DomainZ,	&
					LagrSectSizeCutoff,		&
					ThermoBulkMode,			&
					HowManyBins,			&
					BinEdges,			&
					micron,				&
					SmallestAerosolPossible,        &
					Monodisperse,                   &
                                        Pi, AverageAerosolDensity

	USE GridPointFields,     ONLY : GetTemp, GetRelativeHumidity

	USE Chemistry,		ONLY :	HowManyAqChems,		        &
					AqMolecularMass,	        &
					howmanyaqanions,                &
                                        howmanyaqcations,               &
					HowManyOrgChems,	        &
					OrgMolecularMass,               &
                                        AqPhaseChemicalNames,           &
                                        HowManyAqOrgChems

	IMPLICIT NONE

        !!Input Variables
        LOGICAL :: Boundary !TRUE if boundary distribution, FALSE if main distribution

	!! Local Variables
	INTEGER :: I, J, K, X, Y, Z, Q,	&
	           NumbParticles,	&
	           SectionalParticles,	&
		   LagrangianParticles

	REAL*8 :: II, Temperature, ln_sigma_sq, vtot, mtot, d_mass_mean
	REAL*8 :: del_diam, vrat, diam, mass_bin, store_num, d_num_mean
	REAL*8 :: sum_num, new_num, sum_mass, num_all_modes, rescale
	REAL*8 :: RescalingFactors(HowManyAerosolModes,HowManyAqChems+HowManyOrgChems)
        REAL*8 :: RH, InorgWaterContent

	CHARACTER (len= 1024) :: str
	LOGICAL :: Scaffolding, Transcribe, Jacobson
	TYPE(Particle), POINTER :: NewSection
	TYPE(Particle), POINTER :: Current, Trailer, SectionRunner

        !!This holds the AerosolModes or EnvAerosolModes arrays, depending on the value of Boundary
        REAL*8, ALLOCATABLE :: ModesArray(:,:)
        INTEGER :: NumModes, allocation_error

	Scaffolding= .FALSE.
	Transcribe = .TRUE.
        Jacobson = .TRUE. !Use Jacobson init or MELAM init


	IF(Scaffolding .OR. Transcribe) CALL Transcript("")
	IF(Scaffolding.OR. Transcribe)	CALL Transcript(">>Entering PopulateParticlesSectionsRightAway()<<")
	IF(Transcribe)	CALL Transcript("_Assigning_Particles_to_Particle_Linked_List_Structure_")

        IF (.NOT.(BOUNDARY)) THEN
              ALLOCATE (ModesArray(HowManyAerosolModes,3+HowManyAqChems+HowManyOrgChems) , stat = allocation_error)
        ELSE
              ALLOCATE (ModesArray(HowManyEnvAerosolModes,3+HowManyAqChems+HowManyOrgChems) , stat = allocation_error)
        ENDIF

	if (allocation_error > 0) CALL ERROR("Allocation of ModesArray Failed in PopulateParticlesSectionsRightAway()")

        IF (.NOT.(BOUNDARY)) THEN
              ModesArray = AerosolModes
              NumModes = HowManyAerosolModes
        ELSE
              ModesArray = EnvAerosolModes
              NumModes = HowManyEnvAerosolModes              
        ENDIF  

	!! Start the aerosol ID numbering at 1
	IF (.NOT.(BOUNDARY)) CALL SetParticleID (1)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Each of the modes has its own sectional representation !!
	!! for particles below the certain size cutoff.			  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	IF (HowManyBins .GT. 0. .OR. ThermoBulkMode) THEN
		IF (.NOT.ThermoBulkMode .AND. .NOT.Monodisperse) THEN
			DO J = 1, HowManyBins
				NewSection => MakeSection(I,J) 
                                IF (.NOT.(Boundary)) THEN
				    !!Budge the new mode in at the head of the list
				    IF (ASSOCIATED(Particles%First)) THEN	
                                        NewSection%Next	=> Particles%First
                                        Particles%First => NewSection         
				    ELSE
					Particles%First => NewSection
					NULLIFY(NewSection%Next)
					!NewSection%Next=> Null()
				    END IF
                                ELSE
				    !!Budge the new mode in at the head of the list
				    IF (ASSOCIATED(BoundaryParticles%First)) THEN	
                                        NewSection%Next	=> BoundaryParticles%First
                                        BoundaryParticles%First => NewSection         
				    ELSE
					BoundaryParticles%First => NewSection
					NULLIFY(NewSection%Next)
					!NewSection%Next=> Null()
				    END IF
                                END IF    
                           
			END DO

		ELSE IF (Monodisperse) THEN
			!Monodisperse Distribution
			!WRITE(*,*) "Check 0"
			HowManyBins = 1
                        IF (.NOT.(Boundary)) THEN
			   NewSection => MakeMonodisperse(.FALSE.)
			   Particles%First => NewSection
		        ELSE
			   NewSection => MakeMonodisperse(.TRUE.)
			   BoundaryParticles%First => NewSection
                        ENDIF
		ELSE
			!Bulk Mode
                        IF (.NOT.(Boundary)) THEN
          		   NewSection => MakeBulk(.FALSE.)
			   Particles%First => NewSection
                        ELSE
          		   NewSection => MakeBulk(.TRUE.)
			   BoundaryParticles%First => NewSection 
                        ENDIF                       
		END IF

	END IF

	IF(Scaffolding) WRITE(*,*) "Check 1", ThermoBulkMode 
	IF (.NOT.ThermoBulkMode .AND. .NOT.Monodisperse) THEN

!! Allocate all of the aerosol and place them in the appropriate structures:
                sum_num = 0.0
                sum_mass = 0.0
                num_all_modes = 0.0
	        DO I = 1, MAX(1,NumModes)
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Scale the particle concentration by the number !!
		!! concentration and the domain volume.	          !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			NumbParticles = ANINT(ModesArray(I,1) * DomainX/LengthScale * DomainY/LengthScale * DomainZ/LengthScale)
			SectionalParticles  = 0
			LagrangianParticles = 0
			IF(ModesArray(I,3) .NE. 0.0) THEN !Log-normal mode

                                !MJA 05-29-2012 Old method, only integer numbers of particles
			    IF (.NOT.(Jacobson)) THEN
                                   DO J = 1, NumbParticles
					!Matt changed this to hopefully put particles in sections right away
					IF (.NOT.(BOUNDARY)) THEN
                                            CALL MakeParticle(I, K, .FALSE., (DFLOAT(J)-0.5)/NumbParticles, .FALSE.)
                                        ELSE
                                            CALL MakeParticle(I, K, .FALSE., (DFLOAT(J)-0.5)/NumbParticles, .TRUE.)
                                        ENDIF                    
					SectionalParticles  = SectionalParticles + 1 - K
					LagrangianParticles = LagrangianParticles + K
				  END DO
                            ELSE
                                !MJA 05-29-2012 New method
			        !The idea here is to populate the sections with mass and number based on Jacobson,
                                !Fundamentals of Atmospheric Modeling, 2005 
                                !Eq. 13.21 and 13.24 - MJA, 05-29-2012

                                !Calculate parameters out of loop if possible
                                ln_sigma_sq = LOG(ModesArray(I,3))**2
                                d_num_mean = ModesArray(I,2)*micron !Convert from um to cm 
                                !Jacobson, Eq. 13.25
                                vtot = (pi/6.0)*(d_num_mean**3.0)*exp(4.5*ln_sigma_sq)*ModesArray(I,1)
                                mtot = vtot*AverageAerosolDensity

                                num_all_modes = num_all_modes +  ModesArray(I,1)            
 
                                !Seinfeld and Pandis, 1998, Eq. 7.52
                                d_mass_mean = exp(log(d_num_mean)+3*ln_sigma_sq)

                                !Jacobson, Eq 13.6
                                vrat = (BinEdges(2,2)/BinEdges(2,1))**3.0


                                !Loop over sections
                                IF (.NOT.(BOUNDARY)) THEN
				   Current => PARTICLES%First
				ELSE
				   Current => BOUNDARYPARTICLES%First
				ENDIF
 
                ! CMB(AER): Rewrite this conditional to play nice with gfortran
 99             if ((current%sectional) .and. current%edges(2) .lt. 1.0E4) then    
! 99			        IF (Current%Sectional .EQ. .TRUE. .AND. Current%Edges(2) .LT. 1.0E4 ) THEN
                           
                                    !Jacobson, Eq 13.9
                                    del_diam = 2.0*(Current%Edges(2)-Current%Edges(1)) 

                                    !Jacobson, Eq. 13.9
                                    diam = del_diam/(2.0**(1.0/3.0))
                                    diam = diam*(1.0+vrat)**(1.0/3.0)
                                    diam = diam/(vrat**(1.0/3.0)-1.0)
                                    !Write(*,*) vrat, del_diam, diam, 2*Current%Edges(2), 2*Current%Edges(1), d_num_mean
                                    !STOP
                                    !Jacobson, Eq. 13.24: Number of particles
                                    store_num = Current%NumberOfParticles

                                    if (Current%Edges(1) .NE. 0.0) then !Skip first bin for now
                                         new_num = (ModesArray(I,1)*(del_diam/diam) &
                                                       *exp(-0.5*(log(diam/d_num_mean)**2)/ln_sigma_sq) &
                                                       /sqrt(2*Pi*ln_sigma_sq))
                                         if(new_num .gt. 1.0e-6) then !Skip bins with less than 1e-6 particles/cm3
                                           sum_num = sum_num + new_num
                                           Current%NumberOfParticles = Current%NumberOfParticles + new_num
                                         endif
                                    else
                                         Current%NumberOfParticles = 0.0 !Skip first bin for now
                                    endif

                                    !Jacobson, Eq. 13.21: Mass of particles
                                    if (Current%Edges(1) .NE. 0.0) then !Skip first bin for now
                                       mass_bin = mtot*(del_diam/diam) &
                                                      *exp(-0.5*(log(diam/d_mass_mean)**2)/ln_sigma_sq) &
                                                      /sqrt(2*Pi*ln_sigma_sq)
                                       !write(*,*) "mass check: ", mass_bin, d_mass_mean
                                       if(new_num .gt. 1.0e-6) then !Skip bins with less than 1e-6 particles/cm3
                                          sum_mass = sum_mass+mass_bin
                                       endif
                                    else
                                       sum_mass = 0.0
                                    endif

                                    if(Current%NumberOfParticles .gt. 1.0e-6 .and. Current%Edges(1) .NE. 0.0) then 
                                     !Skip bins with less than 1e-6 particles/cm3
                                     !and first bin
	                              DO J = 1, HowManyAqChems !MJA, 02-28-2012
                                        !WRITE(*,*) "MakeParticle: ", AverageAerosolDensity
                                         Current%AqChems(J) = (Current%AqChems(J)*store_num &
                                                              + ModesArray(I,3+J) *mass_bin/AqMolecularMass(J))		&
							     /Current%NumberOfParticles
			              END DO
			
			              DO J = 1, HowManyOrgChems
				         Current%OrgChems(J) = (Current%OrgChems(J)*store_num &
                                                              + ModesArray(I,3+HowManyAqChems+J) *mass_bin/OrgMolecularMass(J)) &
							     /Current%NumberOfParticles
			              END DO
                                    else
                                      DO J = 1, HowManyAqChems
                                        !WRITE(*,*) "MakeParticle: ", AverageAerosolDensity
                                         Current%AqChems(J) = 0.0
			              END DO
			
			              DO J = 1, HowManyOrgChems
				         Current%OrgChems(J) = 0.0
			              END DO
                                    endif

			            !NOTE: We initialize with all of the organic compounds in 
			            !the organic phase, and none in the hydrophilic phase
			            DO J = 1, HowManyAqOrgChems
				         Current%AqOrgChems(J) = 0.0
		 	            END DO
                                ENDIF

                                IF(ASSOCIATED(Current%Next)) THEN
				     Current => Current%Next
				     GOTO 99
                                ENDIF
                        
                            ENDIF
		            !! Report the results and allocations of the population
			    IF (TRANSCRIBE) CALL TRANSCRIPT ("")
			    IF (TRANSCRIBE) CALL TRANSCRIPT ("For the '"//TRIM(AerosolModeNames(I))//"' aerosol mode")
			    IF (TRANSCRIBE) CALL TRANSCRIPT ("Allocated "//TRIM(INT2STR(LagrangianParticles))//" Lagrangian particles")
			    IF (TRANSCRIBE) CALL TRANSCRIPT ("and "//TRIM(INT2STR(SectionalParticles))//" Sectional particles")
			ELSE !Exponential Mode
				IF (.NOT.(BOUNDARY)) THEN
                                   CALL MakeParticleExponential(I, .FALSE.)
                                ELSE
                                   CALL MakeParticleExponential(I, .TRUE.)
                                ENDIF   
			END IF
		END DO
	

	END IF


	IF(Scaffolding) WRITE(*,*) "Check 1.5"

	!! If the user specified that the contain a certain mass of each 
	!! species then rescale everything
	IF (InputFlagRatioOrMass .EQ. 1) THEN
		
		!!This should only be called for monodisperse or bulk aerosol
		IF(.NOT.(Monodisperse) .AND. .NOT.(ThermoBulkMode)) &
			CALL ERROR ("You selected a sectional aerosol.  In this mode, you must specify the abundance of chemicals by "// &
						"the relative mass percentage using the second flag in AerosolModes.in (0); you have selected another "// &
						"option.  Please remedy.")	

		DO I = 1, NumModes
		DO J = 1, HowManyAqChems+HowManyOrgChems
			RescalingFactors(I,J) = 0.
		END DO ; END DO

		IF (.NOT.(BOUNDARY)) THEN
		    Current => PARTICLES%First
                ELSE
		    Current => BOUNDARYPARTICLES%First
                ENDIF

                ! CMB (AER): Cleaned syntax a bit to play well with gfortran
	10	IF (ASSOCIATED(Current)) THEN

			DO I = 2,HowManyAqChems
				RescalingFactors(Current%ParticleDistribution,I) = &
                                    RescalingFactors(Current%ParticleDistribution,I) + &
				Current%AqChems(I)*Current%NumberOfParticles
			END DO

			DO I = 1,HowManyOrgChems
				RescalingFactors(Current%ParticleDistribution,HowManyAqChems+I) = &
                                    RescalingFactors(Current%ParticleDistribution,HowManyAqChems+I) + &
				Current%OrgChems(I)*Current%NumberOfParticles
			END DO

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
		END IF 

		IF(Scaffolding) WRITE(*,*) "Check 2"


		!! Then Calculate the Scaling
		DO I = 1, MAX(1,NumModes-1)
			DO J = 2, HowManyAqChems
				IF (RescalingFactors(I,J) .GT. 0.) &
					RescalingFactors(I,J) = ModesArray(I,J+3) /		&
											RescalingFactors(I,J) / 1.e12 / AqMolecularMass(J)
			END DO 
	
			DO J = 1, HowManyOrgChems
				IF (RescalingFactors(I,HowManyAqChems+J) .GT. 0.) &
					RescalingFactors(I,HowManyAqChems+J) = ModesArray(I,HowManyAqChems+J+3) /		&
											RescalingFactors(I,HowManyAqChems+J) / 1.e12 / OrgMolecularMass(J)
			END DO 
		
		
		END DO

		IF(Scaffolding) WRITE(*,*) "Check 3"
		IF (.NOT.(BOUNDARY)) THEN
		    Current => PARTICLES%First
                ELSE
		    Current => BOUNDARYPARTICLES%First
                ENDIF

	20	IF (ASSOCIATED(Current)) THEN

			DO I = 2,HowManyAqChems
				Current%AqChems(I) = Current%AqChems(I) * RescalingFactors(Current%ParticleDistribution,I)
			END DO

			DO I = 1,HowManyOrgChems
				Current%OrgChems(I) = Current%OrgChems(I) * RescalingFactors(Current%ParticleDistribution,HowManyAqChems+I)
			END DO
				
			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 20
			END IF
		END IF

        ELSE
		!!This should not be called for monodisperse or bulk aerosol
		IF(Monodisperse .OR. ThermoBulkMode) &
		   CALL ERROR ("You selected a monodisperse or bulk aerosol.  In this mode, you must specify the abundance of chemicals by "// &
					"the total mass conc. using the second flag in AerosolModes.in (1); you have selected another "// &
					"option.  Please remedy.")
	END IF


	IF(Scaffolding) WRITE(*,*) "Check 4"


	!! Now that all of the consitutents have been added,
	!! walk through the lists and equilibrate everything
	IF (.NOT.(BOUNDARY)) THEN
	    Current => PARTICLES%First
        ELSE
            Current => BOUNDARYPARTICLES%First
        ENDIF
	NULLIFY(Trailer)
	!Trailer => Null()

30	IF (ASSOCIATED(Current)) THEN

		IF (Current%NumberOfParticles .GT. 0.0) THEN
			II = 0.
			DO I = 2, HowManyAqChems
				!IF (Current%AqChems(I) .NE. 0.0) WRITE(*,*) AqPhaseChemicalNames(I)
				II = II + Current%AqChems(I)
			END DO
			DO I = 2, HowManyOrgChems
				!IF (Current%AqChems(I) .NE. 0.0) WRITE(*,*) AqPhaseChemicalNames(I)
				II = II + Current%OrgChems(I)
			END DO

                        !WRITE(*,*) "H2O: ", Current%AqChems(1)
			!! Add a token amount of water to avoid initial numerical problems
			IF (II .GT. 0.0) THEN
                           Current%AqChems(1) = II*1.
                        ELSE
                           !If only BC specified, particle is dry
                           Current%Dry = .TRUE.
                        ENDIF
          


			IF (.NOT.(Current%Dry)) THEN 
                           CALL KusikMeissner (Current)
			   IF(Scaffolding) WRITE(*,*) "Before Equilibrate"

                           CALL FindElectrolyteEquilibrium(Current, .TRUE., FirstEquilibration=.TRUE., ReturnType = Q)
			   IF(Scaffolding) WRITE(*,*) "After Equilibrate"

		           RH =  GetRelativeHumidity ()

                           CALL EquilibriumWaterContentAmount (Current, RH, Q, .TRUE., InorgWaterContent)
				!write(*,*) 'InorgWaterContent=',InorgWaterContent
                           Current%AqChems(1) = InorgWaterContent
			   IF(Scaffolding) WRITE(*,*) "After Water EQ"
			ENDIF
                         !WRITE(*,*) "H2O Second: ", Current%AqChems(1)
                        CALL RecalculateRadius (Current)
				IF(Scaffolding) WRITE(*,*) "After RecalculateRadius"

                        !CALL ShellRefIndAndRad(Current)
				IF(Scaffolding) WRITE(*,*) "After ShellRefIndAndRad"
		
		END IF

		Trailer => Current
		Current => Current%Next		
                GOTO 30
	END IF

	!! Report what ended up with
	IF(.NOT.(Monodisperse) .AND. .NOT.(ThermoBulkMode)) THEN
           CALL TRANSCRIPT("")
	   CALL TRANSCRIPT("Reallocated "//TRIM(INT2STR(SectionalParticles))//" particles to sections and left "//		&
	       TRIM(INT2STR(LagrangianParticles))//" as lagrangian particles.")
	   CALL TRANSCRIPT("")
        ELSEIF(Monodisperse) THEN
           CALL TRANSCRIPT("")
	   CALL TRANSCRIPT("Monodisperse distribution.") 
	   CALL TRANSCRIPT("")
        ELSE
           CALL TRANSCRIPT("")
	   CALL TRANSCRIPT("Bulk Thermodynamics Mode.") 
	   CALL TRANSCRIPT("")
        ENDIF
	RETURN
END SUBROUTINE PopulateParticlesSectionsRightAway 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine reads in the information on the       !!
!! inital aerosol modes.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadDistributionsOrganic ()

        ! CB: Added striptoken
	USE InfrastructuralCode, ONLY : GetFileHandle,		        &
					ReturnFileHandle,	        &
					GetLine,			&
					GetToken,			&
					IsReal,				&
					IsInteger,			&
					INT2STR,			&
					STR2INT,			&
					STR2REAL,			&
					EOF,				&
					ERROR,				&
					WARN,				&
					Transcript,			&
					REAL2STR,                       &
                                        stripToken

	USE ModelParameters,	 ONLY : lengthscale,		&
					InputDeckSubDir,	&
					micron,	grams, cm,	&
					ThermoBulkMode,		&
					SetSectionalParameters,	&
					Monodisperse

	USE Chemistry,		 ONLY : FindChem, & 
                                        HowManyAqChems, HowManyOrgChems

	IMPLICIT NONE

	!! Internal Variables
	INTEGER				:: I, J, K, L, FH, allocation_error
	REAL*8				:: II,JJ,KK, ChemPercentage
	LOGICAL				:: Scaffolding, Transcribe
	CHARACTER (len=256) :: CurrLine, Token, TempLine
	LOGICAL             :: LagrangianOrNot

	Scaffolding = .TRUE.	!! Blurt Notes to Standard Out and Transcript
	Transcribe  = .TRUE.	!! Record Activities for Posterity

	!! SCSCSCSC -- Announce Arrival
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript(">>Entering ReadDistributionsOrganic()<<")

	!! Begin by presuming not using Bulk Mode or Monodisperse Mode
	ThermoBulkMode = .FALSE.

	!! Open the Input Deck
	FH = GetFileHandle()
        OPEN(UNIT=FH, FILE=TRIM(InputDeckSubDir)//'AerosolModes.in', STATUS='OLD')
   	WRITE(*,*) "AerosolModes.in Read from: ", InputDeckSubDir

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Parse the INPUT FLAGS ahead of !!
	!! the distribution descriptions  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	!! Get whether or not we should allow lagrangian particles
	CurrLine = GetLine(FH)
      
	IF (TRIM(CurrLine) .EQ. "F") THEN
		LagrangianOrNot = .FALSE. !No Lagrangian particles
	ELSE IF (TRIM(CurrLine) .EQ. "T") THEN
		LagrangianOrNot = .TRUE.  !Yes Lagrangian Particles
	ELSE
		CALL ERROR("Did not find the Lagrangian-particles-or-not flag for the aerosol distribution.  The first non-comment "// &
	           "line of the input file AerosolModes.in should be a flag telling the model if you are allowing "//     &
			   "Lagrangian particles (T) or forcing a secitonal distribution (F).")
	END IF


	!! Get the input mass-or-number type flag here
	CurrLine = GetLine(FH)
	IF (.NOT.IsInteger(CurrLine))  &
	CALL ERROR("Did not find the input-type flag for the aerosol distribution.  The second non-comment line of the input file "// &
	           "AerosolModes.in should be a flag telling the model if you are inputting the constituents based on their "//       &
			   "relative molecular abundance (0) or the total mass burden (in mg / m3) in the mode (1).  It should be an integer.")
	InputFlagRatioOrMass = STR2INT(TRIM(CurrLine))

	IF (InputFlagRatioOrMass .NE. 0 .AND. InputFlagRatioOrMass .NE. 1) &
	CALL ERROR("The input-type flag for the aerosol distribution is out of range.  The second non-comment line of the input "// &
			   "file AerosolModes.in should be a flag telling the model if you are inputting the constituents based on their "// &
			   " relative molecular abundance (0) or the total mass burden (in mg / m3) in the mode (1).")


	!! Get the SECTIONAL representation parameters

	!! The lagrangian / sectional cutoff radius (in microns)
	CurrLine = GetLine(FH)
	IF (.NOT.IsReal(CurrLine)) &
	CALL ERROR("Did not find the Sectional / Lagrangian cutoff for the aerosol distribution.  The third non-comment line of "// &
           "the input file AerosolModes.in.  It should be a positive real number (or zero).  I found: >>"//TRIM(CurrLine)//"<<")

	II = STR2REAL(TRIM(CurrLine))*micron
	IF (II .LT. 0.) &
	CALL ERROR("Did not find the Sectional / Lagrangian cutoff for the aerosol distribution.  The third non-comment line of "// &
	           "the input file AerosolModes.in.  It should be a positive real number.  I found: >>"//TRIM(CurrLine)//"<<")

	!! The Number of Bins
	CurrLine = GetLine(FH)
	IF (.NOT.IsInteger(CurrLine)) &
	CALL ERROR("Did not find the Number of Bins input for the aerosol distribution.  The fourth non-comment line of the input "// &
	           "file AerosolModes.in.  It should be a positive integer or zero.  I found: >>"//TRIM(CurrLine)//"<<")

	I = STR2INT(TRIM(CurrLine))
	IF (I .LE. 0) CALL ERROR("Did not find the Number of Bins input for the aerosol distribution (input out of range) in "// &
	                         "ReadDistributions().  The fifth non-comment line of the input file AerosolModes.in.  It "//    &
							 "should be a positive integer.  Instead of citing zero bins, have one with a very small "//     &
							 "threshold.  I found: >>"//TRIM(CurrLine)//"<<")

	!! The Upper Cutoff of the Lowest Bin Size
	CurrLine = GetLine(FH)
	IF (.NOT.IsReal(CurrLine)) &
	CALL ERROR("Did not find the cutoff of the lowest Sectional bin for the aerosol distribution.  The sixth non-comment "// &
	           "line of the input file AerosolModes.in.  It should be a positive real number.  I found: >>"//TRIM(CurrLine)//"<<")
	JJ = STR2REAL(TRIM(CurrLine))*micron
	IF (JJ .LT. 0 .OR. JJ .GT. II) &
	CALL ERROR("The cutoff of the lowest Sectional bin for the aerosol distribution is out of bounds.  The sixth non-comment "// &
	           "line of the input file AerosolModes.in.  It should be a positive real number.  I found: >>"//TRIM(CurrLine)//"<<")
	IF (I .EQ. 1 .AND. LagrangianOrNot .AND. JJ .NE. II) &
	CALL ERROR("The cutoff of the lowest Sectional bin for the aerosol distribution is out of bounds.  The sixth "// &
	           "non-comment line of the input file AerosolModes.in.  If the number of bins is set to one and we are "// &
			   "tracking lagrangian particles this must equal the Sectional / Lagrangian cutoff value.  I found: >>"// &
			   TRIM(CurrLine)//"<<")

	IF (I .EQ. 2 .AND. .NOT.LagrangianOrNot  .AND. JJ .NE. II) &
	CALL ERROR("The cutoff of the lowest Sectional bin for the aerosol distribution is out of bounds.  The sixth non-comment "// &
	           "line of the input file AerosolModes.in.  If the number of bins is set to two and we are not allowing "// &
			   "lagrangian particles this must equal the Sectional / Lagrangian cutoff value.  I found: >>"//TRIM(CurrLine)//"<<")

	CALL SetSectionalParameters (II,JJ,I,LagrangianOrNot)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Parse the Distribution Descriptions !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Make a first pass through the file to determine 
	!! the number of Modes Specified.  You need this 
	!! number early to allocate all of the arrays.
	HowManyAerosolModes = 0
	K = 0
10	CurrLine = GetLine(FH)
	IF (CurrLine .NE. EOF) THEN
		CALL GetToken(CurrLine,";",Token)
		IF (ISREAL(Token)) HowManyAerosolModes = HowManyAerosolModes+1
		GOTO 10
	END IF

	!! Allocate the AerosolModes Array
	ALLOCATE (AerosolModes(HowManyAerosolModes,3+HowManyAqChems+HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of AerosolModes Failed in ReadDistributions()")
	ALLOCATE (AerosolModeNames(HowManyAerosolModes), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of AerosolModeNames Failed in ReadDistributions()")

	IF (Transcribe) THEN
		CALL TRANSCRIPT ("")
		CALL TRANSCRIPT ("_While_Parsing_AerosolModes.in,_found "//TRIM(INT2STR(MAX(1,HowManyAerosolModes-1)))// &
						 "_Aerosol_Modes_to_Sample_From:_")
	END IF

	!! Back up for a second pass to parse the modal traits
	REWIND(FH)


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! skip the input flag line !!
	!! Change this if you add   !!
	!! more flags in.           !!
	DO I = 1, 5 !!!!!!!!!!!!!!!!!!
		CurrLine = GetLine(FH)
	END DO


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Parse the distributions into their appropriate subcomponents !!
	I = 0
20	CurrLine = GetLine(FH)

	CALL GetToken(CurrLine,";",Token)

	!! Normalize the Chemicals Specified if specifying by ratio
	IF (InputFlagRatioOrMass .EQ. 0) THEN
	  IF (TRIM(Token) .EQ. EOF .OR. (IsReal(Token) .AND. I .NE. 0)) THEN
		IF (ChemPercentage .EQ. 0) CALL ERROR("Must specify some composition for aerosol mode "//TRIM(AerosolModeNames(I)))
                DO J = 4, HowManyAqChems+HowManyOrgChems+3
			AerosolModes(I,J) = AerosolModes(I,J) / ChemPercentage
		END DO
	  END IF 
	END IF

	!! This is a new mode and will accept the log-normal parameters
	IF (TRIM(Token) .NE. EOF) THEN
	  
	  IF (IsReal(TRIM(Token))) THEN

		I = I + 1
		IF (I > MAX(1,HowManyAerosolModes)) &
		CALL ERROR ("In ReadDistributions(), while allocating aerosol modes, parsed the input deck incorrectly.  At line: >>"// &
				    TRIM(CurrLine)//"<<")

		IF (Transcribe) CALL Transcript("")
		IF (Transcribe) CALL Transcript(TRIM(INT2STR(I))//"st Mode:")

		!! Nt
		AerosolModes(I,1) = STR2REAL(Token)

		!! Check for bulk mode
		IF (AerosolModes(I,1) .EQ. 0) THEN ! yes bulk mode!!

			IF (HowManyAerosolModes  .NE. 1) &
			CALL ERROR ("By entering 0 for the number of aerosol in one of the modes in AerosolModes.in, you are "// &
					    "setting the model into BULK THERMODYNAMIC MODE.  In that mode, you may not specify more than one "// &
						"aerosol mode (it assumes that your aerosol solution is like a cup of water, not a distribution of "// &
						"particles).  Please remedy.  If you want the model to ignore a mode you've specifies in the file, "// &
						"you must comment it out...")

			IF (InputFlagRatioOrMass .NE. 1) &
			CALL ERROR ("By entering 0 for the number of aerosol in one of the modes in AerosolModes.in, you are setting the "// &
						"model into BULK THERMODYNAMIC MODE.  In that mode, you must specify the abundance of chemicals by "// &
						"the total mass concentration using the second flag in AerosolModes.in; you have selected another "// &
						"option.  Please remedy.")

			IF (Monodisperse) &
			CALL ERROR ("By entering 0 for the number of aerosol in one of the modes in AerosolModes.in, you are setting the "// &
						"model into BULK THERMODYNAMIC MODE.  However, you selected a monodisperse aerosol in the main input deck."// &
						" This is incorrect. Please pick one.")

			CALL TRANSCRIPT ("The model will be run in BULK THERMODYNAMICS MODE, in which size distributions are "// &
							 "completely unresolved.")

			ThermoBulkMode = .TRUE.

			AerosolModes(I,2) = 1
			AerosolModes(I,3) = 1

		ELSE

			IF (Monodisperse .AND. HowManyAerosolModes  .NE. 1) then
                            CALL ERROR ("You selected a monodisperse aerosol in the main input deck. "//&
                                    "In that mode, you may not specify more than one "// &
				    "aerosol mode. Please remedy.  If you want the model to ignore "//&
                                    "a mode you've specified in the file, "// &
				    "you must comment it out...")
                        endif
			
			IF (Monodisperse .AND. InputFlagRatioOrMass .NE. 1) then
                            CALL WARN ("You selected a monodisperse aerosol in the main input deck.  "//&
                                    "In that mode, you must specify the abundance of chemicals by "// &
				    "the total mass concentration using the second flag in "//&
                                    "AerosolModes.in (1); you have selected another "// &
				    "option.  Please remedy.")
                        endif

			IF (Monodisperse) then
                            CALL TRANSCRIPT ("The model will be run in MONODISPERSE MODE, in which "//&
                                            "only one particle  is considered.")
			endif
			IF(.NOT.(Monodisperse) .AND. InputFlagRatioOrMass .NE. 0) &
			CALL WARN ("You selected a sectional aerosol.  In this mode, you must specify "//&
                                "the abundance of chemicals by "// &
				"the relative mass percentage using the second flag in "//&
                                "AerosolModes.in (0); you have selected another "// &
				"option.  Please remedy.")
			
			IF (Transcribe) CALL Transcript(" Nt:       "//Trim(Token))

			!! Dbar_g_N 
			IF (.NOT.IsReal(Token)) &
			CALL ERROR ("In ReadDistributions(), while allocating aerosol modes, expected "//&
                                "a real number and got "//TRIM(Token)// &
				" in line: >>"//TRIM(CurrLine)//"<<")

			CALL GetToken(CurrLine,";",Token)
			AerosolModes(I,2) = STR2REAL(Token)

			IF (Transcribe) CALL Transcript(" Dbar_g_N: "//Trim(Token))

			!! Rho_g  
			CALL GetToken(CurrLine,";",Token)
			IF (.NOT.IsReal(Token)) &
			CALL ERROR ("In ReadDistributionsOrganic(), while allocating aerosol modes, "//&
                                "expected a real number and got "//TRIM(Token)//&
                                " in line: >>"//TRIM(CurrLine)//"<<")

			AerosolModes(I,3) = STR2REAL(Token)

			IF (Transcribe) CALL Transcript(" Rho_g:    "//Trim(Token))

			IF (AerosolModes(I,3) .LT. 1 .AND. AerosolModes(I,3) .NE. 0) &
			CALL ERROR("In AerosolModes.in, Rho_g, the geometric standard deviation "//&
                                "of aerosol mode "//TRIM(INT2STR(I))//" is less than one.  It cannot "//&
                                "be, mathematically, as it is ill defined.  Please correct.")
			IF (AerosolModes(I,3) .EQ. 1.) &
			CALL ERROR("In AerosolModes.in, Rho_g, the geometric standard deviation "//&
                                "of aerosol mode "//TRIM(INT2STR(I))// &
				" is equal to one, which is ill defined (an infinite concentration "//&
                                "in zero radius space.  If you "//    &
				"want a mono-dispersion set that flag in MELAMinputdeck.in. ")
			IF (AerosolModes(I,3) .LT. 1.001 .AND. AerosolModes(I,3) .NE. 0) &
			CALL WARN("In AerosolModes.in, Rho_g, the geometric standard deviation "//&
                                "of aerosol mode "//TRIM(INT2STR(I))// &
				 " is very near to one, which may cause problems.  If you want a "//&
                                 "mono-dispersion set that flag in "//        &
				  "MELAMinputdeck.in. ")
            ! CMB (AER): Edit this conditional for gfortran
            IF (AerosolModes(I,3) .EQ. 0.0 .and. (.not. Monodisperse)) &
!			IF (AerosolModes(I,3) .EQ. 0.0 .AND. Monodisperse .EQ. .FALSE.) &
			CALL WARN("In AerosolModes.in, Rho_g, the geometric standard deviation of "//&
                                "aerosol mode #"//TRIM(INT2STR(I))// &
				" is zero. This is the flag for an exponential distribution, "//&
                                "and will be treated as such. ")
		END IF

		!WRITE(*,*) "Hi Matt!"
		
		!! Initialize the Aqueous Arrays
                DO J = 4, HowManyAqChems+HowManyOrgChems+3
			AerosolModes(I,J) = 0.
		END DO
		ChemPercentage = 0

		!! Get the Name of the Mode
                !TempLine = GetLine(FH)
                !do i = 1, 256
                !    num = ichar(tempLine(i:i))
                !    if (num .eq. 13) then
                !        tempLine(i:i) = ' '
            !    endif
                !enddo
                !call GetToken(GetLine(FH), ';', tempLine)  
		AerosolModeNames(I) = TRIM(GetLine(FH))
		!AerosolModeNames(I) = stripToken(tempLine)
                IF (AerosolModeNames(I) .eq. EOF) &
		CALL ERROR ("Reached the end of the Aerosol Input Deck while looking for the name of the "//TRIM(INT2STR(I))//  &
					"th aerosol mode.")


		IF (Transcribe) CALL Transcript("The Mode is Called: '"//Trim(AerosolModeNames(I))//"'")

	!! It is a chemical component and a continuation of the specification of the same mode
	ELSE

		!! Electrolyte or Org Compound

			!! Determine Which Chemical We're Talking About
			J = FindChem(TRIM(Token), 1, SuppressError=.TRUE.) 
			IF (J .EQ. 0) THEN
				
				!Check if it is an organic compound (Phase=2)
				L = FindChem(TRIM(Token), 2, SuppressError=.TRUE.)
				IF (L .EQ. 0) CALL ERROR("Couldn't Figure Out Which Chemical Was Intended by "//TRIM(Token)//	&
									 " in the AerosolModes.in input deck while reading from ReadDistributions()")
				!Put organic compounds after aqueous ones
				J = L + HowManyAqChems
			END IF

			!! Error check about double specification
			!WRITE(*,*) I, 3+J
			IF (AerosolModes(I,3+J) > 0) &
			CALL ERROR("It looks like you double specified the chemical concentration of "// TRIM(Token)//	&
					   " in the "//TRIM(INT2STR(I))//"th mode specied in the AerosolModes.in input deck while "// &
					   "reading from ReadDistributions()")

			IF (Transcribe) CALL Transcript(" containing "//Trim(Token))

			!! Get the Proportional Value of the Chemical
			CALL GetToken(CurrLine,";",Token)
			IF (.NOT.IsReal(Token)) CALL ERROR("The following token is supposed to be a real number >>"//TRIM(Token)//	&
									 "<< in the AerosolModes.in while reading from ReadDistributions()")

			II = STR2REAL(Token)
			ChemPercentage = II + ChemPercentage

			AerosolModes(I,3+J) = II
	END IF

	!! Get the next line if not the end of the file
	GOTO 20
	END IF

	CLOSE(FH)
	CALL ReturnFileHandle(FH)

	InsolCoresYesOrNo = .FALSE.
     
	IF (Scaffolding) CALL Transcript(">>Exiting ReadDistributionsOrganic()<<")
	
	RETURN
	WRITE(*,*) "Number In Dist: ",AerosolModes(1,1)
	WRITE(*,*) "Diam In Dist: ", AerosolModes(1,2)
	WRITE(*,*) "Sigma In Dist: ", AerosolModes(1,3)
END SUBROUTINE ReadDistributionsOrganic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initializes and returns a single sectional pseudo-particle    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION MakeMonodisperse(Boundary)

	USE GridPointFields,	 ONLY : GetTemp, GetPress

	USE Chemistry,	 	 ONLY : HowManyAqChems,			&
					HowManyAqCations,		&
					HowManyAqAnions,		&
					HowManyAqEqReactions,           &
					HowManyOrgChems,                &
					HowManyAqOrgChems

	USE InfrastructuralCode, ONLY : ERROR

	USE ModelParameters,     ONLY : micron,NumbInsolTypes, Lengthscale, &
					DomainX, DomainY, DomainZ	

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle), POINTER :: MakeMonodisperse
        LOGICAL :: Boundary !TRUE if boundary distribution, FALSE if main distribution

	!! Local Variables
	INTEGER :: allocation_error, I,J
	REAL*8  ::  II

	!! If particle is large enough to be tracked, 
	!! create the New Particle and allocate its arrays
	ALLOCATE (MakeMonodisperse, stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse Failed in MakeMonodisperse()")
	ALLOCATE (MakeMonodisperse%AqChems(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse%AqChems Failed in MakeMonodisperse()")
	ALLOCATE (MakeMonodisperse%OrgChems(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse%OrgChems Failed in MakeMonodisperse()")
	ALLOCATE (MakeMonodisperse%AqOrgChems(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse%AqOrgChems Failed in MakeMonodisperse()")



	!! Allocate the Appropriate Electrolytic Arrays
	ALLOCATE (MakeMonodisperse%GammaMixed(HowManyAqEqReactions), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse%GammaMixed Failed in MakeMonodisperse()")



	!! Initialize
	DO I = 1, HowManyAqEqReactions
		MakeMonodisperse%GammaMixed(I) = 0
	END DO



	!Allocate the hydrophobic acitivty coefficients arrays
	ALLOCATE (MakeMonodisperse%HydrophobicActivityCoeffs(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse%HydrophobicActivityCoeffs Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyOrgChems
		MakeMonodisperse%HydrophobicActivityCoeffs(I) = 1.0
	END DO



	!Allocate the hydrophobic acitivty coefficients arrays
	ALLOCATE (MakeMonodisperse%HydrophilicActivityCoeffs(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of MakeMonodisperse%HydrophilicActivityCoeffs Failed in MakeSection()")

	!! Initialize
	DO I = 1, HowManyAqOrgChems
		MakeMonodisperse%HydrophilicActivityCoeffs(I) = 1.0
	END DO

	

	MakeMonodisperse%IonicStr = 1.
	MakeMonodisperse%EmbryoRadius	 = 0.1*micron
	MakeMonodisperse%InsolubleRadius = 0.
	MakeMonodisperse%EffectiveRadius = 0.1*micron


	MakeMonodisperse%ParticleID           = 0.
	MakeMonodisperse%ParticleDistribution = 1

        IF (.NOT.(BOUNDARY)) THEN
	     MakeMonodisperse%NumberOfParticles = &
                ANINT(AerosolModes(1,1) * DomainX/LengthScale * DomainY/LengthScale * DomainZ/LengthScale)
        ELSE
	     MakeMonodisperse%NumberOfParticles = &
                ANINT(EnvAerosolModes(1,1) * DomainX/LengthScale * DomainY/LengthScale * DomainZ/LengthScale)
        ENDIF

	MakeMonodisperse%Sectional	      = .TRUE.

	!! For Now, the Initial Temperature of the Particle is that of the 
	!! Local Environment (1/2002)
	MakeMonodisperse%Temperature = GetTemp()
	MakeMonodisperse%OriginPressureLevel  = GetPress()


	!! Use a fake initial water activity:
	MakeMonodisperse%WaterActivity = 1.

	DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		MakeMonodisperse%AqChems(I) = 0.
	END DO

	II = 0.
	DO I = 2, HowManyAqChems
            IF (.NOT.(BOUNDARY)) THEN
		MakeMonodisperse%AqChems(I) = AerosolModes(1,3+I)
            ELSE
		MakeMonodisperse%AqChems(I) = EnvAerosolModes(1,3+I)
            ENDIF
		II = II + MakeMonodisperse%AqChems(I)
	END DO

	DO I = 1, HowManyOrgChems
	    IF (.NOT.(BOUNDARY)) THEN
                MakeMonodisperse%OrgChems(I) = AerosolModes(1,3+HowManyAqChems+I)
            ELSE
                MakeMonodisperse%OrgChems(I) = EnvAerosolModes(1,3+HowManyAqChems+I)
            ENDIF
	END DO
	
	!Initialize all organic species in hydrophobic phase
	DO I = 1, HowManyAqOrgChems
		MakeMonodisperse%AqOrgChems(I) = 0.0
	END DO

	
	DO I = 3+HowManyAqChems+HowManyOrgChems+1,3+HowManyAqChems+HowManyOrgChems+NumbInsolTypes
		IF (AerosolModes(1,I) .GT. 0) &
		CALL ERROR("MakeMonodisperse() Failed.  You cannot specify insoluble cores "//&
                        "to aerosol and then try to run in monodisperse mode.  "// &
		           "This program is not set up to do that.")
	END DO

	!! As a token amount of water to avoid initial numerical problems, set the 
	!! water mixing ratio at 0.5
	MakeMonodisperse%AqChems(1) = II*1.
        MakeMonodisperse%Dry = .FALSE.

	NULLIFY(MakeMonodisperse%Next)
	!MakeMonodisperse%Next => Null()

	RETURN

END FUNCTION MakeMonodisperse

!!This subroutine adds in an exponential distribution to the sectional
!!size distribution
SUBROUTINE MakeParticleExponential(WhichPopulation, Boundary)

	USE GridPointFields,	 ONLY :  GetTemp, GetPress

	USE ModelParameters,	 ONLY : LagrSectSizeCutoff,		&
					lengthscale,			&
					Pi, EulerMass,			&
					micron,	moles,			&
					AverageAerosolDensity,		&
					DomainX, DomainY, DomainZ,	&
					Monodisperse,			&
					BinEdges, HowManyBins,		&
					NumbInsolTypes,			&
					InsolRadii,			&
					InsolDensity,			&
					InsolContactAngle

	USE Chemistry,			 ONLY : HowManyAqChems,		&
						HowManyAqCations,	&
					        HowManyAqAnions,	&
					        AqMolecularMass,	&
						HowManyAqEqReactions,   &
						HowManyOrgChems,        &
						OrgMolecularMass,       &
						HowManyAqOrgChems,      &
						AqOrgMolecularMass

	USE InfrastructuralCode, ONLY : REAL2STR, INT2STR, Transcript,  &
                                        Warn, Error, GetParticleID

	IMPLICIT NONE

	!! Input Variables
	INTEGER :: WhichPopulation
        LOGICAL :: Boundary !TRUE if boundary distribution, FALSE if main distribution

	!! Local Variables
	REAL*8 :: JJ, N12, M12, Exp1, Exp2, InitNum, VolSum1, VolSum2
	TYPE(Particle), POINTER :: CurrentParticle
	INTEGER				    :: I, allocation_error
	LOGICAL, PARAMETER		:: Scaffolding = .FALSE. !.TRUE.

        !!This holds the AerosolModes or EnvAerosolModes arrays, depending on the value of Boundary
        REAL*8, ALLOCATABLE :: ModesArray(:,:)

	IF (SCAFFOLDING) CALL Transcript(">>>Entering MakeParticleExponential()<<<")
	IF (SCAFFOLDING) CALL Transcript("")

        IF (.NOT.(BOUNDARY)) THEN
              ALLOCATE (ModesArray(HowManyAerosolModes,3+HowManyAqChems+HowManyOrgChems) , stat = allocation_error)
        ELSE
              ALLOCATE (ModesArray(HowManyEnvAerosolModes,3+HowManyAqChems+HowManyOrgChems) , stat = allocation_error)
        ENDIF

	if (allocation_error > 0) CALL ERROR("Allocation of ModesArray Failed in MakeParticle()")

        IF (.NOT.(BOUNDARY)) THEN
              ModesArray = AerosolModes
        ELSE
              ModesArray = EnvAerosolModes
        ENDIF          

	IF (SCAFFOLDING) CALL Transcript(">>>Entering MakeParticleExponential()<<<")
	IF (SCAFFOLDING) CALL Transcript("")

	IF (.NOT.(BOUNDARY)) THEN
            CurrentParticle => Particles%First
        ELSE
            CurrentParticle => BoundaryParticles%First
        ENDIF

	!! Find the appropriate section to put it in, or err
	IF (ASSOCIATED(CurrentParticle)) THEN		
		
	10  InitNum = CurrentParticle%NumberOfParticles
	
		Exp1 = EXP(-8.0*((CurrentParticle%Edges(1)/micron)**3.0)/(ModesArray(WhichPopulation,2)**3))
		Exp2 = EXP(-8.0*((CurrentParticle%Edges(2)/micron)**3.0)/(ModesArray(WhichPopulation,2)**3))
		
		!Number added to this section from Exponential Mode
		N12 = ModesArray(WhichPopulation,1)*(Exp1-Exp2) !particles/cm3
		VolSum1 = (4.0*Pi/3.0)*((CurrentParticle%Edges(1)/micron)**3.0) &
				 + (Pi/6.0)*(ModesArray(WhichPopulation,2)**3)
		
		VolSum2 = (4.0*Pi/3.0)*((CurrentParticle%Edges(2)/micron)**3.0) &
				 + (Pi/6.0)*(ModesArray(WhichPopulation,2)**3)
	
		!Total mass added to this section from Exponential Mode (g/cm3 air)
		M12 = 1.0e-12*AverageAerosolDensity*ModesArray(WhichPopulation,1) &
		      *(VolSum1*Exp1 - VolSum2*Exp2)
			  	
		!Add the appropriate number of particles to each section
		IF(N12 .GT. 1e-6) THEN
			CurrentParticle%NumberOfParticles = CurrentParticle%NumberOfParticles  &
		                                    + N12
			DO I = 2, HowManyAqChems
			
				CurrentParticle%AqChems(I) = (InitNum*CurrentParticle%AqChems(I)  &
						+ ModesArray(WhichPopulation,3+I)*M12/AqMolecularMass(I))/CurrentParticle%NumberOfParticles
			END DO
			
			DO I = 1, HowManyOrgChems
				CurrentParticle%OrgChems(I) = (InitNum*CurrentParticle%OrgChems(I)  &
						+ ModesArray(WhichPopulation,3+HowManyAqChems+I)*M12/OrgMolecularMass(I))/CurrentParticle%NumberOfParticles
			END DO
			
			!NOTE: We initialize with all of the organic compounds in 
			!the organic phase, and none in the hydrophilic phase
			DO I = 1, HowManyAqOrgChems
				CurrentParticle%AqOrgChems(I) = CurrentParticle%AqOrgChems(I)+ 0.0
			END DO
		ELSE
			Currentparticle%numberofparticles = CurrentParticle%NumberOfParticles
			DO I = 2, HowManyAqChems
			
				CurrentParticle%AqChems(I) = CurrentParticle%AqChems(I) 
			END DO
			
			DO I = 1, HowManyOrgChems
				CurrentParticle%OrgChems(I) = CurrentParticle%OrgChems(I) 
			END DO
			
			!NOTE: We initialize with all of the organic compounds in 
			!the organic phase, and none in the hydrophilic phase
			DO I = 1, HowManyAqOrgChems
				CurrentParticle%AqOrgChems(I) = CurrentParticle%AqOrgChems(I)+ 0.0
			END DO
		END IF
		
		IF(ASSOCIATED(CurrentParticle%Next)) THEN
			CurrentParticle => CurrentParticle%Next
			GOTO 10
		END IF

		IF (SCAFFOLDING) CALL Transcript(">>>Exiting MakeParticleExponential()<<<")
		RETURN
	END IF

END SUBROUTINE MakeParticleExponential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine reads in the information on the	 !!
!! boundary aerosol modes.                               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadBoundaryDistributionsOrganic ()

	USE InfrastructuralCode, ONLY : GetFileHandle,		&
					ReturnFileHandle,	&
					GetLine,			&
					GetToken,			&
					IsReal,				&
					IsInteger,			&
					INT2STR,			&
					STR2INT,			&
					STR2REAL,			&
					EOF,				&
					ERROR,				&
					WARN,				&
					Transcript,			&
					REAL2STR

	USE ModelParameters,	 ONLY : lengthscale,		&
					InputDeckSubDir,	&
					micron,	grams, cm,	&
					ThermoBulkMode,		&
					SetSectionalParameters,	&
					Monodisperse

	USE Chemistry,		 ONLY : FindChem, &
                                        HowManyAqChems, HowManyOrgChems

	IMPLICIT NONE

	!! Internal Variables
	INTEGER				:: I, J, K, L, FH, allocation_error
	REAL*8				:: II,JJ,KK, ChemPercentage
	LOGICAL				:: Scaffolding, Transcribe
	CHARACTER (len=256) :: CurrLine, Token
	LOGICAL             :: LagrangianOrNot

	Scaffolding = .TRUE.	!! Blurt Notes to Standard Out and Transcript
	Transcribe  = .TRUE.	!! Record Activities for Posterity

	!! SCSCSCSC -- Announce Arrival
	IF (Scaffolding) CALL Transcript("")
	IF (Scaffolding) CALL Transcript(">>Entering ReadBoundaryDistributionsOrganic()<<")

	!! Open the Input Deck
	FH = GetFileHandle()
    OPEN(UNIT=FH, FILE=TRIM(InputDeckSubDir)//'EnvAerosolModes.in', STATUS='OLD')

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Parse the Distribution Descriptions !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Make a first pass through the file to determine 
	!! the number of Modes Specified.  You need this 
	!! number early to allocate all of the arrays.
	HowManyEnvAerosolModes = 0
	K = 0
10	CurrLine = GetLine(FH)
	IF (CurrLine .NE. EOF) THEN
		CALL GetToken(CurrLine,";",Token)
		IF (ISREAL(Token)) HowManyEnvAerosolModes = HowManyEnvAerosolModes+1
		GOTO 10
	END IF

	!! Allocate the AerosolModes Array
	ALLOCATE (EnvAerosolModes(HowManyEnvAerosolModes,3+HowManyAqChems+HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of EnvAerosolModes Failed in ReadBoundaryDistributionsOrganic()")
	ALLOCATE (EnvAerosolModeNames(HowManyEnvAerosolModes), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of EnvAerosolModeNames Failed in ReadBoundaryDistributionsOrganic()")

	IF (Transcribe) THEN
		CALL TRANSCRIPT ("")
		CALL TRANSCRIPT ("_While_Parsing_EnvAerosolModes.in,_found "//TRIM(INT2STR(MAX(1,HowManyEnvAerosolModes-1)))// &
						 "_Aerosol_Modes_to_Sample_From:_")
	END IF

	!! Back up for a second pass to parse the modal traits
	REWIND(FH)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Parse the distributions into their appropriate subcomponents !!
	I = 0
20	CurrLine = GetLine(FH)

	CALL GetToken(CurrLine,";",Token)

	!! Normalize the Chemicals Specified if specifying by ratio
	IF (InputFlagRatioOrMass .EQ. 0) THEN
	  IF (TRIM(Token) .EQ. EOF .OR. (IsReal(Token) .AND. I .NE. 0)) THEN
		IF (ChemPercentage .EQ. 0) CALL ERROR("Must specify some composition for aerosol mode "//TRIM(AerosolModeNames(I)))
		DO J = 4, HowManyAqChems+HowManyOrgChems+3
			EnvAerosolModes(I,J) = EnvAerosolModes(I,J) / ChemPercentage
		END DO
	  END IF 
	END IF

	

	!! This is a new mode and will accept the log-normal parameters
	IF (TRIM(Token) .NE. EOF) THEN
	  
	  IF (IsReal(TRIM(Token))) THEN

		I = I + 1
		IF (I > MAX(1,HowManyEnvAerosolModes)) &
		CALL ERROR ("In ReadBoundaryDistributionsOrganic(), while allocating aerosol modes, "//&
                        "parsed the input deck incorrectly.  At line: >>"//&
				    TRIM(CurrLine)//"<<")

		IF (Transcribe) CALL Transcript("")
		IF (Transcribe) CALL Transcript(TRIM(INT2STR(I))//"st Mode:")

		!! Nt
		EnvAerosolModes(I,1) = STR2REAL(Token)

		!! Check for bulk mode
		IF (ThermoBulkMode) THEN ! yes bulk mode!!

			IF (HowManyEnvAerosolModes  .NE. 1) &
			CALL ERROR ("By entering 0 for the number of aerosol in one of the modes in AerosolModes.in, you are "// &
					    "setting the model into BULK THERMODYNAMIC MODE.  In that mode, you may not specify more than one "// &
						"aerosol mode (it assumes that your aerosol solution is like a cup of water, not a distribution of "// &
						"particles).  Please remedy.  If you want the model to ignore a mode you've specifies in the file, "// &
						"you must comment it out...")

		    IF (Monodisperse) &
			CALL ERROR ("By entering 0 for the number of aerosol in one of the modes in AerosolModes.in, you are setting the "// &
						"model into BULK THERMODYNAMIC MODE.  However, you selected a monodisperse aerosol in the main input deck."// &
						" This is incorrect. Please pick one.")

			
			EnvAerosolModes(I,2) = 1
			EnvAerosolModes(I,3) = 1

		ELSE

			IF (Monodisperse .AND. HowManyEnvAerosolModes  .NE. 1) &
			CALL ERROR ("You selected a monodisperse aerosol in the main input deck. In that mode, you may not specify more than one "// &
						"aerosol mode. Please remedy.  If you want the model to ignore a mode you've specified in the file, "// &
						"you must comment it out...")
			
			IF (Monodisperse) CALL TRANSCRIPT ("The model will be run in MONODISPERSE MODE, in which only one particle  is considered.")
			
			IF (Transcribe) CALL Transcript(" Nt:       "//Trim(Token))

			!! Dbar_g_N 
			IF (.NOT.IsReal(Token)) &
			CALL ERROR ("In ReadBoundaryDistributionsOrganic(), while allocating "//&
                                "aerosol modes, expected a real number and got "//TRIM(Token)// &
						" in line: >>"//TRIM(CurrLine)//"<<")

			CALL GetToken(CurrLine,";",Token)
			EnvAerosolModes(I,2) = STR2REAL(Token)

			IF (Transcribe) CALL Transcript(" Dbar_g_N: "//Trim(Token))

			!! Rho_g  
			CALL GetToken(CurrLine,";",Token)
			IF (.NOT.IsReal(Token)) &
			CALL ERROR ("In ReadBoundaryDistributionsOrganic(), while allocating aerosol "//&
                                "modes, expected a real number and got "//TRIM(Token)// &
						" in line: >>"//TRIM(CurrLine)//"<<")

			EnvAerosolModes(I,3) = STR2REAL(Token)

			IF (Transcribe) CALL Transcript(" Rho_g:    "//Trim(Token))

			IF (EnvAerosolModes(I,3) .LT. 1 .AND. EnvAerosolModes(I,3) .NE. 0) &
			CALL ERROR("In EnvEnvAerosolModes.in, Rho_g, the geometric standard deviation of aerosol mode "//TRIM(INT2STR(I)) &
					   //" is less than one.  It cannot be, mathematically, as it is ill defined.  Please correct.")
			IF (EnvAerosolModes(I,3) .EQ. 1.) &
			CALL ERROR("In EnvAerosolModes.in, Rho_g, the geometric standard deviation of aerosol mode "//TRIM(INT2STR(I)) &
					   //" is equal to one, which is ill defined (an infinite concentration in zero radius space.  If you "//    &
					   "want a mono-dispersion set that flag in MELAMinputdeck.in. ")
			IF (EnvAerosolModes(I,3) .LT. 1.001 .AND. EnvAerosolModes(I,3) .NE. 0) &
			CALL WARN("In EnvAerosolModes.in, Rho_g, the geometric standard deviation of aerosol mode "//TRIM(INT2STR(I))// &
					  " is very near to one, which may cause problems.  If you want a mono-dispersion set that flag in "//        &
					  "MELAMinputdeck.in. ")
			IF (EnvAerosolModes(I,3) .EQ. 0.0) &
			CALL WARN("In EnvAerosolModes.in, Rho_g, the geometric standard deviation of aerosol mode #"//TRIM(INT2STR(I))// &
					  " is zero. This is the flag for an exponential distribution, and will be treated as such. ")
		END IF

		!! Initialize the Aqueous Arrays
		DO J = 4, HowManyAqChems+HowManyOrgChems+3
			EnvAerosolModes(I,J) = 0.
		END DO
		ChemPercentage = 0

		!! Get the Name of the Mode
		EnvAerosolModeNames(I) = TRIM(GetLine(FH))
		IF (EnvAerosolModeNames(I) .eq. EOF) &
		CALL ERROR ("Reached the end of EnvAerosolModes.in while looking for the name of the "//TRIM(INT2STR(I))//  &
					"th aerosol mode.")


		IF (Transcribe) CALL Transcript("The Mode is Called: '"//Trim(EnvAerosolModeNames(I))//"'")

	!! It is a chemical component and a continuation of the specification of the same mode
	ELSE

		!! Electrolyte or Org Compound

			!! Determine Which Chemical We're Talking About
			J = FindChem(TRIM(Token), 1, SuppressError=.TRUE.) 
			IF (J .EQ. 0) THEN
				
				!Check if it is an organic compound (Phase=2)
				L = FindChem(TRIM(Token), 2, SuppressError=.TRUE.)
				IF (L .EQ. 0) CALL ERROR("Couldn't Figure Out Which Chemical Was Intended by "//TRIM(Token)//	&
									 " in the EnvAerosolModes.in input deck while reading from ReadDistributions()")
				!Put organic compounds after aqueous ones
				J = L + HowManyAqChems
			END IF

			!! Error check about double specification
			IF (EnvAerosolModes(I,3+J) > 0) &
			CALL ERROR("It looks like you double specified the chemical concentration of "// TRIM(Token)//	&
					   " in the "//TRIM(INT2STR(I))//"th mode specied in the EnvAerosolModes.in input deck while "// &
					   "reading from ReadDistributions()")

			IF (Transcribe) CALL Transcript(" containing "//Trim(Token))

			!! Get the Proportional Value of the Chemical
			CALL GetToken(CurrLine,";",Token)
			IF (.NOT.IsReal(Token)) CALL ERROR("The following token is supposed to be a real number >>"//TRIM(Token)//	&
									 "<< in the EnvAerosolModes.in while reading from ReadBoundaryDistributionsOrganic()")

			II = STR2REAL(Token)
			ChemPercentage = II + ChemPercentage

			EnvAerosolModes(I,3+J) = II
	END IF

	!! Get the next line if not the end of the file
	GOTO 20
	END IF

	CLOSE(FH)
	CALL ReturnFileHandle(FH)

	InsolCoresYesOrNo = .FALSE.
	
	IF (Scaffolding) CALL Transcript(">>Exiting ReadBoundaryDistributionsOrganic()<<")

	RETURN
END SUBROUTINE ReadBoundaryDistributionsOrganic

