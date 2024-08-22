!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! SortAerosol.h
!! This file contains the sorting routines for condensation, 
!! coagulation, and the moving-center size distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY					                     !!
!!									     !!
!! Month  Year   Name              Description	                             !!
!! 07     2006   Matt Alvarado     Began Update History	                     !!
!! 08/30  2006   Matt Alvarado     Edited RegridAerosol		             !!
!!					to prevent unnecessary equilibration !!
!! 09/07  2006   Matt Alvarado     Created SortBoundaryPart                  !!
!! 07/26  2007   Matt Alvarado     Added line to RegridAerosol               !!
!!			 to recalculate particle radius and density          !!
!!			  	every time it is called.                     !!
!! 10/05  2007   Matt Alvarado     Redid RegridAerosol to prevent            !!
!!                                 zeros in size bins                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	    !!	
!! 1. SUBROUTINE SortAerosolAtGridPointForCondensation ()           !!
!! 2. SUBROUTINE SortAerosolAtGridPointForCoagulation ()            !!
!! 3. SUBROUTINE RegridAerosol ()			            !! 
!! 4. SUBROUTINE SortBoundaryPart				    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For Condensation, you want to minimize   !!
!! array size by pushing all of the empty   !!
!! sections to the end.  They're then       !!
!! ignored.  There is a different critereon !!
!! for coagulation.                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SortAerosolAtGridPointForCondensation ()

	IMPLICIT NONE
	
	TYPE(Particle), POINTER :: FirstParticle, CurrentParticle, PlaceHolder

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!   SORTING -- SORTING -- SORTING -- SORTING     !!
		!! If there are empty sections, we don't want to  !!
		!! waste time on them, move them to the end       !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		FirstParticle   => Particles%First
		CurrentParticle => FirstParticle
		PlaceHolder     => Null()

7		FirstParticle => Particles%First

		!! Look at the head of the list first.
		IF(((FirstParticle%Sectional .AND. FirstParticle%NumberOfParticles .EQ. 0) .OR. &
			FirstParticle%Dry) .AND. ASSOCIATED(FirstParticle%Next)) THEN
			Particles%First => FirstParticle%Next

			IF (.NOT.ASSOCIATED(PlaceHolder)) THEN
				PlaceHolder        => FirstParticle
				FirstParticle%Next => Null()
			ELSE
				FirstParticle%Next => PlaceHolder
				PlaceHolder        => FirstParticle
			END IF
			GOTO 7
		END IF

		!! Back up to the head of the aerosol list
		FirstParticle   => Particles%First
		CurrentParticle => FirstParticle

		!! Take one step with FirstParticle
		IF (ASSOCIATED(FirstParticle%Next)) FirstParticle => FirstParticle%Next

10		IF(((FirstParticle%Sectional .AND. FirstParticle%NumberOfParticles .EQ. 0) .OR. &
			FirstParticle%Dry) .AND. ASSOCIATED(FirstParticle%Next)) THEN

			CurrentParticle%Next => FirstParticle%Next

			IF (.NOT.ASSOCIATED(PlaceHolder)) THEN
				PlaceHolder        => FirstParticle
				FirstParticle%Next => Null()
			ELSE
				FirstParticle%Next => PlaceHolder
				PlaceHolder        => FirstParticle
			END IF
			FirstParticle => CurrentParticle%Next
			GOTO 10
		END IF

		IF (ASSOCIATED(FirstParticle%Next)) THEN
			CurrentParticle => CurrentParticle %Next
			FirstParticle   => FirstParticle%Next
			GOTO 10
		END IF

		FirstParticle => Particles%First
	    IF (ASSOCIATED(PlaceHolder)) THEN
20          IF(ASSOCIATED(FirstParticle%Next)) THEN
				FirstParticle => FirstParticle%Next
				GOTO 20
			END IF
			FirstParticle%Next => PlaceHolder
		END IF

		RETURN
END SUBROUTINE SortAerosolAtGridPointForCondensation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For Coagulation, Sort so the sections are first !!
!! and in order (1 is smallest, up to biggest).    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SortAerosolAtGridPointForCoagulation ()

        use infrastructuralcode, only : int2str
	
        IMPLICIT NONE


	
	!! Internal Variables	
	TYPE(Particle), POINTER :: Particle1, Particle2, Particle3, PlaceHolder
	INTEGER :: I, J
	LOGICAL :: LocalFlag

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!   SORTING -- SORTING -- SORTING -- SORTING     !!
	!! Push all lagrangian particles to the end       !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Return if the chain is only one item long 
	IF (.NOT. ASSOCIATED(Particles%First%Next)) RETURN


	!! First, count the number of sectional terms
	I = 0
	LocalFlag = .TRUE. ! True until first out of order or non-sectional bin

	!! Return if the gridpoint is empty
	IF (.NOT. ASSOCIATED(Particles%First)) RETURN


	Particle1 => Particles%First


	IF (Particle1%Sectional) I = I+1  	!! The first one gets missed by the loop

	!! Check if they are in order and count the number of sections
	DO WHILE (ASSOCIATED(Particle1%Next)) 

		!! If the two links are not in order
		IF (Particle1%Next%Sectional .AND. .NOT. (Particle1%Sectional .AND. Particle1%Next%Sectional .AND.	 &
		    (Particle1%Edges(2)  .EQ.  Particle1%Next%Edges(1)  .AND.				       &
		     Particle1%ParticleDistribution .EQ. Particle1%Next%ParticleDistribution) .OR. &
		   (Particle1%ParticleDistribution .EQ. Particle1%Next%ParticleDistribution - 1))) &
		LocalFlag = .FALSE.

		Particle1 => Particle1%Next

			
		!! Keep track of how many sectional items
		IF (Particle1%Sectional) I = I+1

	END DO


	!! Return if already in order or Only One
	IF (LocalFlag) RETURN


	!2. double walk (staggered) and cram all sections at beginning
	Particle1 => Particles%First
	Particle2 => Particle1%Next
	J         = I
	IF (Particle1%Sectional) THEN
		LocalFlag = .TRUE. ! True until first non-sectional bin
	ELSE
		LocalFlag = .FALSE.
	END IF


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Count down on Sectional Pieces while   !!
	!! Grabbing all Sections to the front of  !!
	!! the list (successfully)                !!
10	IF (Particle2%Sectional) THEN !!!!!!!!!!!!!!
		J = J - 1

		!! Excise the link if after at least one non-sectional link
		IF (.NOT.LocalFlag) THEN

			Particle1%Next => Particle2%Next
			Particle2%Next => Particles%First 
			Particles%First => Particle2
			Particle2      => Particle1%Next
			GOTO 10
		END IF

	ELSE
		LocalFlag = .FALSE.
	END IF
		
	IF (ASSOCIATED(Particle2%Next) .AND. J .NE. 0) THEN
		Particle1 => Particle2
		Particle2 => Particle2%Next
		GOTO 10
	END IF

	!! We start checking and counting after the first link, which is 
	!! sectional.  So permanently dunk one of these.
	I = I - 1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Sort to make sure each link is the minimum at its spot !!
	
	LocalFlag = .TRUE. ! this signals that we're working on the first slot

	!! Particle1 sits at the head, Particle2 Trails, and Particle3 Leads
	Particle1   => Particles%First
20	Particle2   => Particle1
	Particle3   => Particle1%Next   ! this won't err because if only one link it've exited already
	J           =  I

	!! If Particle3 belongs in front of Particle1, switch them
30	IF ((Particle1%Edges(2) .GT. Particle3%Edges(1) .AND.							&
		 Particle1%ParticleDistribution .EQ. Particle3%ParticleDistribution)		&
		.OR. (Particle1%ParticleDistribution .GT. Particle3%ParticleDistribution)) THEN 

		!! Budge Particle3 at head of line
		Particle2%Next => Particle3%Next
		Particle3%Next => Particle1
		IF (LocalFlag) THEN         !! This step is different for the first and later particles
			Particles%First => Particle3
		ELSE
			PlaceHolder%Next => Particle3
		END IF

		!! Regain Composure
		Particle1 => Particle3
		Particle3 => Particle2%Next

		J = J - 1
		IF (J .GT. 0) GOTO 30

	END IF

	!! Step forward if still more to sort
	IF (J .GT. 1) THEN
		Particle2 => Particle3
		Particle3 => Particle3%Next
		J = J - 1
		GOTO 30
	END IF


	!! We should have the correct particle in slot PlaceHolder now
	!! Now recursively move back
	IF (I .GT. 1) THEN
		IF (LocalFlag) THEN
			PlaceHolder => Particles%First 
			LocalFlag = .FALSE.
		ELSE
			PlaceHolder => PlaceHolder%Next
		END IF
		I = I - 1
		Particle1   => PlaceHolder%Next

		GOTO 20
	END IF
		RETURN
END SUBROUTINE SortAerosolAtGridPointForCoagulation

SUBROUTINE RegridAerosol ()

!! Copyright Matt Alvarado, 2005
!! This subroutine reallocates sectional particles when the
!! "moving center" of a bin crosses grid boundaries.
!! Calling this after growth & coagulation will make MELAM's
!! sectional distribution behave like a moving center model.
!! NOTE: Also updates temperature, radius, density and 
!! calls inorganic aqueous eq
!! NOTE: Particles are only moved one bin forward or backward, with
!! no additional check to make sure this is the correct bin. Take small
!! time steps!

	!! Include modules
	USE InfrastructuralCode, ONLY : Transcript, &
					ERROR, &
					REAL2STR, &
					INT2STR

        ! CB: Need to comment out first two objects; ifort is complaining
	USE Chemistry,		ONLY : &!HowManyAqCations, &
					!HowManyAqAnions, &
					AqPhaseChemicalNames, &	
					AqCationNames, &
					AqAnionNames

	USE GridPointFields,	 ONLY : GetTempfromGrid

	USE ModelParameters
	!USE OutputRoutines
	IMPLICIT NONE
	
	!! Define the local variables
	integer :: i, j, k, l, seedarray(2),q, r, NumBins, NumAqChems, NumOrgChems, NumAqOrgChems, GG
	real*8  :: ii, jj, kk, y(200)
	real*8  :: Mgas, Temp, Dum1, Dum2, MassFactor, Store, SUM

	REAL*8, ALLOCATABLE	::	NumConc(:)
	TYPE(Particle), POINTER :: CurrentParticle, BiggerParticle, SmallerParticle
	LOGICAL :: IonReset
				
	!Reorder particles, smallest to biggest 
	!WRITE(*,*) " Just in Regrid Aerosol()"
	CALL SortAerosolAtGridPointForCoagulation ()
	!WRITE(*,*) "After SortAerosolAtGridPointForCoagulation ()"		
	!Count # of bins, # of Aq. Chems, 
	!call recalculate radius to get radius and density
	CurrentParticle => particles%first
	NumBins = 0
  8 if (associated(currentparticle)) then
		NumBins = NumBins + 1
		NumAqChems = size(currentparticle%AqChems)
		NumOrgChems = size(currentparticle%OrgChems)
		NumAqOrgChems = size(currentparticle%AqOrgChems)
		!WRITE(*,*) 'currentparticle%EffectiveRadius = ',currentparticle%EffectiveRadius,NumBins
		!WRITE(*,*) "Before RecalculateRadius()"
		!WRITE(*,*) 'currentparticle%EffectiveRadius = ',currentparticle%EffectiveRadius,NumBins
		CALL RecalculateRadius(currentparticle)	
		!WRITE(*,*) "After RecalculateRadius()"	
		!WRITE(*,*) 'currentparticle%EffectiveRadius = ',currentparticle%EffectiveRadius,NumBins
		if (associated(currentparticle%next)) then
			currentparticle => currentparticle%next
			goto 8
		end if
	end if
	ALLOCATE(NumConc(NumBins))

	!WRITE(*,*) "Before loop in Regrid Aerosol()"

	CurrentParticle => particles%first
	J = 1
 40 if (associated(currentparticle)) then    								
					
		!If the particle has no particles, ignore it
		IF(CurrentParticle%numberofparticles .GT. 0.0) THEN
						
			!If particle has grown past upper boundary, and a bigger bin exists,
		    !add its number and mass to the next largest bin
			!Note that this check makes sure this routine does nothing
			!when the model is run in bulk mode.
			
			!MJA 100507: The goal of the new scheme is to leave some number
			!and mass in each bin for numerical stability, while conserving mass and number
			!The MassFactor places the remaining particles in the middle of the
			!size bin (in a volume sense).

			!write(*,*) "In loop CurrentParticle%EffectiveRadius = ",CurrentParticle%EffectiveRadius,J

			IF(CurrentParticle%EffectiveRadius .GT.  CurrentParticle%Edges(2) &
				.AND. associated(currentparticle%next) .AND. CurrentParticle%numberofparticles .GT. 2.e-4) THEN
						
				MassFactor = (CurrentParticle%Edges(1)**3 + CurrentParticle%Edges(2)**3) &
				              /(2.0*CurrentParticle%EffectiveRadius**3)
				
				!WRITE(*,*) "Growth", J, MassFactor
				BiggerParticle => currentparticle%next
							
				!Add particle mass to larger bin
				!Remember AqChems values are in units of mol/particle
				I=1
				DO I=1, NumAqChems
					Store = BiggerParticle%AqChems(I)
					BiggerParticle%AqChems(I) = (BiggerParticle%AqChems(I)*BiggerParticle%numberofparticles &
					                             + CurrentParticle%AqChems(I)*CurrentParticle%numberofparticles &
												 -MassFactor*CurrentParticle%AqChems(I)*1.0e-4) &
												 /(BiggerParticle%numberofparticles + CurrentParticle%numberofparticles - 1.0e-4)
					IF(BiggerParticle%AqChems(I) .LT. 0.0 .OR. IonReset) THEN
							!WRITE(*,*) "Below 0!", I
							BiggerParticle%AqChems(I) = Store
					END IF
					CurrentParticle%AqChems(I) = MassFactor*CurrentParticle%AqChems(I)
				END DO				
				DO I=1, NumOrgChems
					Store = BiggerParticle%OrgChems(I)
					BiggerParticle%OrgChems(I) = (BiggerParticle%OrgChems(I)*BiggerParticle%numberofparticles &
					                             + CurrentParticle%OrgChems(I)*CurrentParticle%numberofparticles &
												 -MassFactor*CurrentParticle%OrgChems(I)*1.0e-4)&
												 /(BiggerParticle%numberofparticles + CurrentParticle%numberofparticles - 1.0e-4)
					IF(BiggerParticle%OrgChems(I) .LT. 0.0) BiggerParticle%OrgChems(I) = Store
					CurrentParticle%OrgChems(I) = MassFactor*CurrentParticle%OrgChems(I)
				END DO
				DO I=1, NumAqOrgChems
					Store = BiggerParticle%AqOrgChems(I)
					BiggerParticle%AqOrgChems(I) = (BiggerParticle%AqOrgChems(I)*BiggerParticle%numberofparticles &
					                             + CurrentParticle%AqOrgChems(I)*CurrentParticle%numberofparticles &
												 -MassFactor*CurrentParticle%AqOrgChems(I)*1.0e-4) &
												 /(BiggerParticle%numberofparticles + CurrentParticle%numberofparticles - 1.0e-4)
					IF(BiggerParticle%AqOrgChems(I) .LT. 0.0) BiggerParticle%AqOrgChems(I) = Store
					CurrentParticle%AqOrgChems(I) = MassFactor*CurrentParticle%AqOrgChems(I)
				END DO		

								
				!Add particle number to larger bin
				Store = BiggerParticle%numberofparticles	
				BiggerParticle%numberofparticles = &
					BiggerParticle%numberofparticles + CurrentParticle%numberofparticles - 1.0e-4
				IF(BiggerParticle%numberofparticles .LT. 0.0) BiggerParticle%numberofparticles = Store
				CurrentParticle%numberofparticles = 1.0e-4

				! consider add check for if there is no mass, reset num to zero - CRL

				!WRITE(*,*) "Water: ", BiggerParticle%AqChems(1)
				!Calculate density to reset radius
				!WRITE(*,*) "Density 1CurrentParticle%EffectiveRadius = ",CurrentParticle%EffectiveRadius,J
				Dum1 = ParticleDensity(BiggerParticle)
				!WRITE(*,*) "Density 2CurrentParticle%EffectiveRadius = ",CurrentParticle%EffectiveRadius,J
				Dum2 = ParticleDensity(CurrentParticle)
				!WRITE(*,*) "After Density 2CurrentParticle%EffectiveRadius = ",CurrentParticle%EffectiveRadius,J
				!Reequilibrate particles internally (inorganic only)
				!WRITE(*,*) "Before EqCurrentParticle%EffectiveRadius = ",CurrentParticle%EffectiveRadius,J
				CALL FindElectrolyteEquilibrium (BiggerParticle, UpdateThermo=.TRUE., FirstEquilibration=.FALSE., ReturnType = GG)
				CALL FindElectrolyteEquilibrium (CurrentParticle, UpdateThermo=.TRUE., FirstEquilibration=.FALSE., ReturnType = GG)
				!WRITE(*,*) "After Eq CurrentParticle%EffectiveRadius = ",CurrentParticle%EffectiveRadius,J
			!This accounts for shrinking particles, which are rare in clouds and plumes
			!The second condition says do nothing to smallest particle or bulk
			
			ELSE IF (CurrentParticle%EffectiveRadius .LT.  CurrentParticle%Edges(1) &
					.AND. J .GT. 1) THEN

			    IF(J .EQ. NumBins) THEN
					MassFactor = 2.0
				ELSE
					MassFactor = (CurrentParticle%Edges(1)**3 + CurrentParticle%Edges(2)**3) &
				              /(2.0*CurrentParticle%EffectiveRadius**3)
					WRITE(*,*) "In SortAerosol.h CurrentParticle%EffectiveRadius = ", CurrentParticle%EffectiveRadius
					WRITE(*,*) "In SortAerosol.h (CurrentParticle%Edges(1)= ", CurrentParticle%Edges(1)
					WRITE(*,*) "In SortAerosol.h (CurrentParticle%Edges(2)= ", CurrentParticle%Edges(2)
					WRITE(*,*) "In SortAerosol.h MassFactor = ", MassFactor
					!WRITE(*,*) "Calling Spillbeans()"
					!call spillbeans()
				END IF
				!IF (CurrentParticle%EffectiveRadius .EQ. 0.) THEN
				!	MassFactor = 0.
				!	Write(*,*) "**Caution: CurrentParticle%EffectiveRadius equals zero, so making MassFactor 0."
				!END IF				
				WRITE(*,*) "Shrink", J, CurrentParticle%EffectiveRadius, CurrentParticle%Edges(1),CurrentParticle%Edges(2), MassFactor
				!Scan back through particle list to find one just below
				SmallerParticle => particles%first
				K = 1
			 50 if (associated(smallerparticle)) then
					IF (K .EQ. J-1) THEN
						!Add particle mass to smaller bin
						I=1
						!IonReset makes sure that if one ion is reset
						!to previous values, they all are
						IonReset = .FALSE.
						DO I=1, NumAqChems
							!IF(SmallerParticle%numberofparticles .LT. 1.01e-4) IonReset = .TRUE.
							Store = SmallerParticle%AqChems(I)
							SmallerParticle%AqChems(I) = (SmallerParticle%AqChems(I)*SmallerParticle%numberofparticles &
					                             + CurrentParticle%AqChems(I)*CurrentParticle%numberofparticles) &
												 /(SmallerParticle%numberofparticles + CurrentParticle%numberofparticles)
							IF(SmallerParticle%AqChems(I) .LT. 0.0 .OR. IonReset) THEN
								WRITE(*,*) "Below 0!", I, SmallerParticle%numberofparticles
								SmallerParticle%AqChems(I) = Store
								IF(I .GT. HowManyAqChems) IonReset = .TRUE.
							END IF
							CurrentParticle%AqChems(I) = MassFactor*CurrentParticle%AqChems(I)
						END DO				
						DO I=1, NumOrgChems
							Store = SmallerParticle%OrgChems(I)
							SmallerParticle%OrgChems(I) = (SmallerParticle%OrgChems(I)*SmallerParticle%numberofparticles &
					                             + CurrentParticle%OrgChems(I)*CurrentParticle%numberofparticles) &
												 /(SmallerParticle%numberofparticles + CurrentParticle%numberofparticles)
							IF(SmallerParticle%OrgChems(I) .LT. 0.0) SmallerParticle%OrgChems(I) = Store
							CurrentParticle%OrgChems(I) = MassFactor*CurrentParticle%OrgChems(I)
						END DO
						DO I=1, NumAqOrgChems
							Store = SmallerParticle%AqOrgChems(I)
							SmallerParticle%AqOrgChems(I) = (SmallerParticle%AqOrgChems(I)*SmallerParticle%numberofparticles &
					                             + CurrentParticle%AqOrgChems(I)*CurrentParticle%numberofparticles) &
												 /(SmallerParticle%numberofparticles + CurrentParticle%numberofparticles)
							IF(SmallerParticle%AqOrgChems(I) .LT. 0.0) SmallerParticle%AqOrgChems(I) = Store
							CurrentParticle%AqOrgChems(I) = MassFactor*CurrentParticle%AqOrgChems(I)
						END DO	

								
						!Add particle number to smaller bin
						Store = SmallerParticle%numberofparticles
						SmallerParticle%numberofparticles = &
										SmallerParticle%numberofparticles + CurrentParticle%numberofparticles
						IF(SmallerParticle%numberofparticles .LT. 0.0) SmallerParticle%numberofparticles = Store
						CurrentParticle%numberofparticles = 1.0e-4
							
						!Calculate density to reset radius
						!WRITE(*,*) "Density 1"
						Dum1 = ParticleDensity(SmallerParticle)
						!WRITE(*,*) "Density 2"
						Dum2 = ParticleDensity(CurrentParticle)
						!WRITE(*,*) "After Density 2"
						!Reequilibrate particles internally (inorganic only)
						!WRITE(*,*) "Before Eq 1"
						CALL FindElectrolyteEquilibrium (SmallerParticle, UpdateThermo=.TRUE., FirstEquilibration=.FALSE., ReturnType=GG)
						!WRITE(*,*) "Before Eq 2"
						CALL FindElectrolyteEquilibrium (CurrentParticle, UpdateThermo=.TRUE., FirstEquilibration=.FALSE., ReturnType=GG)
						!WRITE(*,*) "After Eq 2"
					END IF

					if (associated(smallerparticle%next)) then
						smallerparticle => smallerparticle%next
						K = K+1
						goto 50
					end if			
				END IF !Smaller particle Search loop
			END IF !Grow or Shrink
		END IF !Check for Non-zero

		!Update Temperature
		Currentparticle%Temperature = GetTempfromGrid ()
		
		!Go to next particle
		if (associated(currentparticle%next)) then

			currentparticle => currentparticle%next
			J = J+1
			goto 40
		end if	 
	
	END IF !Main particle loop
	
	! Doing a final check to see if there is mass when numberofparticles exists- CRL
	CurrentParticle => particles%first

	if (associated(CurrentParticle)) then    								
	   if(CurrentParticle%numberofparticles .GT. 0.) THEN
		SUM = 0.
		DO I=1, NumAqChems
			SUM = SUM + CurrentParticle%AqChems(I)
		END DO

		DO I=1, NumOrgChems
			SUM = SUM + CurrentParticle%OrgChems(I)
		END DO

		DO I=1, NumAqOrgChems
			SUM = SUM + CurrentParticle%AqOrgChems(I)
		END DO	
		
		if(SUM .EQ. 0. .AND. CurrentParticle%numberofparticles .GT. 0.) THEN
			write(*,*) '**Warning numberofparticles exists without mass!!!'
			CurrentParticle%numberofparticles = 0.
		end if

	   end if

	   CurrentParticle => CurrentParticle%next

	end if	

		
	RETURN
END SUBROUTINE RegridAerosol
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine is similar to			    !!
!! SortAerosolAtGridPointForCoagulation, but for the!!
!! boundary particles.				    !!
!! Sort so the sections are first		    !!
!! and in order (1 is smallest, up to biggest).	    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SortBoundaryPart()

use infrastructuralcode, only : int2str
	IMPLICIT NONE
	
	!! Internal Variables	
	TYPE(Particle), POINTER :: Particle1, Particle2, Particle3, PlaceHolder
	INTEGER :: I, J
	LOGICAL :: LocalFlag



	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!   SORTING -- SORTING -- SORTING -- SORTING     !!
	!! Push all lagrangian particles to the end       !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Return if the chain is only one item long 
	IF (.NOT. ASSOCIATED(BoundaryParticles%First%Next)) RETURN


	!! First, count the number of sectional terms
	I = 0
	LocalFlag = .TRUE. ! True until first out of order or non-sectional bin

	!! Return if the gridpoint is empty
	IF (.NOT. ASSOCIATED(BoundaryParticles%First)) RETURN


	Particle1 => BoundaryParticles%First


	IF (Particle1%Sectional) I = I+1  	!! The first one gets missed by the loop

	!! Check if they are in order and count the number of sections
	DO WHILE (ASSOCIATED(Particle1%Next)) 

		!! If the two links are not in order
		IF (Particle1%Next%Sectional .AND. .NOT. (Particle1%Sectional .AND. Particle1%Next%Sectional .AND.	 &
		    (Particle1%Edges(2)  .EQ.  Particle1%Next%Edges(1)  .AND.				       &
		     Particle1%ParticleDistribution .EQ. Particle1%Next%ParticleDistribution) .OR. &
		   (Particle1%ParticleDistribution .EQ. Particle1%Next%ParticleDistribution - 1))) &
		LocalFlag = .FALSE.

		Particle1 => Particle1%Next

			
		!! Keep track of how many sectional items
		IF (Particle1%Sectional) I = I+1

	END DO


	!! Return if already in order or Only One
	IF (LocalFlag) RETURN


	!2. double walk (staggered) and cram all sections at beginning
	Particle1 => BoundaryParticles%First
	Particle2 => Particle1%Next
	J         = I
	IF (Particle1%Sectional) THEN
		LocalFlag = .TRUE. ! True until first non-sectional bin
	ELSE
		LocalFlag = .FALSE.
	END IF


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Count down on Sectional Pieces while   !!
	!! Grabbing all Sections to the front of  !!
	!! the list (successfully)                !!
10	IF (Particle2%Sectional) THEN !!!!!!!!!!!!!!
		J = J - 1

		!! Excise the link if after at least one non-sectional link
		IF (.NOT.LocalFlag) THEN

			Particle1%Next => Particle2%Next
			Particle2%Next => BoundaryParticles%First 
			BoundaryParticles%First => Particle2
			Particle2      => Particle1%Next
			GOTO 10
		END IF

	ELSE
		LocalFlag = .FALSE.
	END IF
		
	IF (ASSOCIATED(Particle2%Next) .AND. J .NE. 0) THEN
		Particle1 => Particle2
		Particle2 => Particle2%Next
		GOTO 10
	END IF

	!! We start checking and counting after the first link, which is 
	!! sectional.  So permanently dunk one of these.
	I = I - 1

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Sort to make sure each link is the minimum at its spot !!
	
	LocalFlag = .TRUE. ! this signals that we're working on the first slot

	!! Particle1 sits at the head, Particle2 Trails, and Particle3 Leads
	Particle1   => BoundaryParticles%First
20	Particle2   => Particle1
	Particle3   => Particle1%Next   ! this won't err because if only one link it've exited already
	J           =  I

	!! If Particle3 belongs in front of Particle1, switch them
30	IF ((Particle1%Edges(2) .GT. Particle3%Edges(1) .AND.							&
		 Particle1%ParticleDistribution .EQ. Particle3%ParticleDistribution)		&
		.OR. (Particle1%ParticleDistribution .GT. Particle3%ParticleDistribution)) THEN 

		!! Budge Particle3 at head of line
		Particle2%Next => Particle3%Next
		Particle3%Next => Particle1
		IF (LocalFlag) THEN         !! This step is different for the first and later particles
			BoundaryParticles%First => Particle3
		ELSE
			PlaceHolder%Next => Particle3
		END IF

		!! Regain Composure
		Particle1 => Particle3
		Particle3 => Particle2%Next

		J = J - 1
		IF (J .GT. 0) GOTO 30

	END IF

	!! Step forward if still more to sort
	IF (J .GT. 1) THEN
		Particle2 => Particle3
		Particle3 => Particle3%Next
		J = J - 1
		GOTO 30
	END IF


	!! We should have the correct particle in slot PlaceHolder now
	!! Now recursively move back
	IF (I .GT. 1) THEN
		IF (LocalFlag) THEN
			PlaceHolder => BoundaryParticles%First 
			LocalFlag = .FALSE.
		ELSE
			PlaceHolder => PlaceHolder%Next
		END IF
		I = I - 1
		Particle1   => PlaceHolder%Next

		GOTO 20
	END IF

		RETURN
END SUBROUTINE SortBoundaryPart
