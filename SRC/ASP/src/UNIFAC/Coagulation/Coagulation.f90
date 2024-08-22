!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Coagulation.f90
!! This file contains the coagulation numerical integration routine 
!! for the internally mixed, moving-center size distribution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 08/26  2006   Matt Alvarado     Fixed bug in CollisionEfficiency          !!
!! 07/26  2007   Matt Alvarado     Added RecalculateRadius to                !!
!!                                  monodisperse coag	                     !!
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES				    !!
!! 1. ModelParameters			    !!
!! 2. InfrastructuralCode		    !!
!! 3. GridPointFields			    !!
!! 4. Chemistry				    !!
!! 5. Aerosols				    !!
!! 6. OutputRoutines			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES						    !!
!! 1. CoagulationKernels.h					    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. REAL*8 FUNCTION CollisionEfficiency (Particle1, Particle2)	!!
!! 2. SUBROUTINE StepSectionalCoagulationJacobson (TimeStep)            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


MODULE Coagulation

	PRIVATE

	REAL*8, PARAMETER :: CE = 1.0 !Coagulation Efficiency
					!(fraction of collisions where
					!particles stick together)
	
	PUBLIC :: StepSectionalCoagulationJacobson
			  
CONTAINS

INCLUDE "CoagulationKernels.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Collision Efficiency 
!! using formula from Jacobson, "Fundamentals of Atmospheric Modeling"
!! , 2nd ed., p. 510 (2005)
REAL*8 FUNCTION CollisionEfficiency (Particle1, Particle2)

	USE Aerosols,	ONLY : Particle, TerminalVelocity, ReynoldsNumber

	USE ModelParameters, ONLY : g
	USE InfrastructuralCode, ONLY : Error
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle), POINTER :: Particle1, Particle2

	!! Internal Variables
	REAL*8  :: Stokes, Reynolds, Ev, Ea

	!From Jacobson, "Fundamentals of Atmospheric Modeling", p. 510
	IF(Particle2%EffectiveRadius .GE. Particle1%EffectiveRadius) THEN
		Stokes = TerminalVelocity(Particle1) &
                         *ABS(TerminalVelocity(Particle2) - TerminalVelocity(Particle1)) &
                         /Particle2%EffectiveRadius/g
		Reynolds = ReynoldsNumber(Particle2)
		Ea = (Stokes**2.0)/((Stokes+0.5)**2.0)
		IF(Stokes .LE. 1.214) THEN
			Ev = 0.
		ELSE
			Ev = (1 + (0.75*LOG(2*Stokes))/(Stokes-1.214) )**(-2.0)
		END IF

		CollisionEfficiency = (60*Ev + Ea*Reynolds)/(60 + Reynolds)

	ELSE
		Stokes = TerminalVelocity(Particle2)*ABS(TerminalVelocity(Particle1) - TerminalVelocity(Particle2))/Particle1%EffectiveRadius/g
		Reynolds = ReynoldsNumber(Particle1)
		Ea = (Stokes**2.0)/((Stokes+0.5)**2.0)
		IF(Stokes .LE. 1.214) THEN
			Ev = 0.
		ELSE
			Ev = (1 + (0.75*LOG(2*Stokes))/(Stokes-1.214) )**(-2.0)
		END IF

		CollisionEfficiency = (60*Ev + Ea*Reynolds)/(60 + Reynolds)
	
	
	END IF

	IF (CollisionEfficiency .GT. 1.0 .OR. CollisionEfficiency .LT. 0.0) THEN
		WRITE(*,*) CollisionEfficiency
		CALL ERROR("CollisionEfficiency in Coagulation.f90 is greater than 1 or less than 0.")
	END IF
END FUNCTION CollisionEfficiency 

SUBROUTINE StepSectionalCoagulationJacobson (TimeStep)  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Matt Alvarado created this by simplifying StepCoagulation to only 
!! work for a sectional distribution and to always use mean values.
!! This follows the semi-implicit method of Jacobson - see
!! Jacobson, "Fundamentals of Atmospheric Modeling"
!! , 2nd ed., Section 15.2, (2005)

	USE InfrastructuralCode, ONLY : INT2STR, ERROR

	USE Aerosols,            ONLY : Particles, Particle,		&
				SortAerosolAtGridPointForCoagulation,	&
					RecalculateRadius,		&
					RegridAerosol

	USE ModelParameters,	  ONLY : ThermoBulkMode,		  &
					 Pi, &
					 ConstantKernel, &
					 DoCoagulationFlag, &
					 Monodisperse, &
					 CoagKernel

	USE Chemistry,	          ONLY : HowManyAqChems,		&
					 HowManyAqCations,	&
				         HowManyAqAnions, &
					 AqMolecularMass, &
					 HowManyOrgChems, &
					 OrgMolecularMass, &
					 OrgDensity, &
					 HowManyAqOrgChems, &
					 AqOrgMolecularMass, &
					 AqOrgDensity, &
					 SolidSaltDensity, &
					 AqPhaseChemicalNames
	USE OutputRoutines

	IMPLICIT NONE

	!! External Variables
	REAL*8, INTENT(IN)  :: TimeStep

	!! Internal Variables
	INTEGER	:: HowManyLinks, I, J, K, L, M, Q, Count, JJ

	REAL*8	:: AA, BB, II, LL, Kfac, &
		   BigV, LowerSum, UpperSum, InnerSum, TotalAqOrgVolConc, & 
                   TotalOrgVolConc, TotalAqVolConc, &
		   InVolConc, OutVolConc, SolnSum, AvgSolutionDensity, & 
                   OutSumVolumeConc, InSumVolumeConc, StoreDensity, Monokernel

	REAL*8, ALLOCATABLE	:: CoagArray(:,:), SingleParticleVolume(:), &
                                   VolumeConc(:), &
				   VolumeFractions(:,:,:), NumConc(:), &
                                   AqVolConc(:,:), &
				   OrgVolConc(:,:), AqOrgVolConc(:,:), &
                                   InitNumConc(:), SumVolConc(:)

	TYPE(Particle), POINTER :: Particle1, Particle2, ParticleK, & 
                                   Runner, ReplicaPlaceholder, &
		DifferenceParticle1, DifferenceParticle2, DifferenceRunner , &
			           ParticleI, ParticleJ

	!IF main input deck said to skip coagulation, skip it
	IF(.NOT.DoCoagulationFlag) RETURN
	
	!! If we're in bulk mode, coagulation doesn't do anything
	IF (ThermoBulkMode) RETURN

	!! Sort the gridcell, sections in order and first
	CALL SortAerosolAtGridPointForCoagulation ()

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If there are no aerosol at the gridcell, !!
	!! then there is nothing to do.  Return.    !!
	Particle1 => Particles%First !!!!!!!





5	IF (Particle1%NumberOfParticles .EQ. 0.0D0) THEN
		IF (ASSOCIATED(Particle1%Next)) THEN
			Particle1 => Particle1%Next
			GOTO 5
		END IF
		RETURN
	END IF

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Everything should now be sorted           !!
	!! Back up to the head of the aerosol list   !!
	Particle1 => Particles%First !!!!!!!!
	Particle2 => Particle1
	Runner    => Particle1 !! In the second pass, this may need an initial position


	!! COUNT the number of LINKS in the gridcell (total number of sections)
	HowManyLinks =  0
15	HowManyLinks = HowManyLinks + 1

	IF(ASSOCIATED(Particle2%Next)) THEN
		Particle2 => Particle2%Next
		GOTO 15
	END IF
	Particle2 => Particle1


	!! BUILD N VECTOR for Single particle volumes
	ALLOCATE (SingleParticleVolume(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate SingleParticleVolume at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Volume Concentrations
	ALLOCATE (VolumeConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate VolumeConc at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Number Concentrations
	ALLOCATE (NumConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate NumConc at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Initial Number Concentrations
	ALLOCATE (InitNumConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate NumConc at size "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N VECTOR for Summed Volume Concentrations
	ALLOCATE (SumVolConc(HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate NumConc at size "//TRIM(INT2STR(HowManyLinks)))

	
	!! BUILD N by N ARRAY for COAGULATION KERNELS
	ALLOCATE (CoagArray(HowManyLinks,HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate CoagArray at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyLinks)))

	!! BUILD N by M ARRAY for Species Volume Concentrations
	ALLOCATE (AqVolConc(HowManyLinks, HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate AqVolConc at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))

	!! BUILD N by M ARRAY for Species Volume Concentrations
	ALLOCATE (OrgVolConc(HowManyLinks, HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate OrgVolConc at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyOrgChems)))

	!! BUILD N by M ARRAY for Species Volume Concentrations
	ALLOCATE (AqOrgVolConc(HowManyLinks, HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In StepSectionalCoagulationJacobson(), couldn't allocate AqOrgVolConc at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyAqOrgChems)))

	
	!! BUILD N by N by N ARRAY for VolumeFractions
	ALLOCATE (VolumeFractions(HowManyLinks, HowManyLinks, HowManyLinks), STAT = I)
	IF (I > 0) CALL ERROR("In StepCoagulation(), couldn't allocate VolumeFractions at size "//TRIM(INT2STR(HowManyLinks))//	&
						  " by "//TRIM(INT2STR(HowManyLinks))//" by "//TRIM(INT2STR(HowManyLinks)))

	

	!! If we're in monodisperse mode, do monodisperse coag
	IF (Monodisperse) THEN
		Particle1 => Particles%First
		Particle2 => Particle1
		IF(ConstantKernel) THEN
			MonoKernel = CoagKernel	
		ELSE
			MonoKernel = (CoagKernBrownianDiff   (Particle1, Particle2) &
					                  + CoagKernGravCollection (Particle1, Particle2)) !! Add Brownian and Convective Contributions
		END IF
		
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(1,I) = Particle1%AqChems(I)*Particle1%NumberofParticles !mol/cm3 air
		END DO
		
		DO I = 1, HowManyOrgChems
			OrgVolConc(1,I) = Particle1%OrgChems(I)*Particle1%NumberofParticles !mol/cm3 air
		END DO

		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(1,I) = Particle1%AqOrgChems(I)*Particle1%NumberofParticles !mol/cm3 air
		END DO

		!Explicit coagulation integration
		Particle1%NumberofParticles = Particle1%NumberofParticles - Timestep*MonoKernel*Particle1%NumberofParticles**2
		
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Particle1%AqChems(I) = AqVolConc(1,I)/Particle1%NumberofParticles !mol/particle
		END DO
		
		DO I = 1, HowManyOrgChems
			Particle1%OrgChems(I) = OrgVolConc(1,I)/Particle1%NumberofParticles !mol/particle
		END DO

		DO I = 1, HowManyAqOrgChems
			Particle1%AqOrgChems(I) = AqOrgVolConc(1,I)/Particle1%NumberofParticles !mol/particle
		END DO
		
		!Does radius and density
		CALL RecalculateRadius(Particle1)		
			
		RETURN		
	END IF


	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! LOOP AND FILL CoagArray with Kernel Values   !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	DO I = 1, HowManyLinks
		DO J = 1, HowManyLinks

			!! If either of these is zero we don't have to actually calculate anything
			IF (Particle1%NumberOfParticles .EQ. 0.0D0  .OR.  Particle2%NumberOfParticles .EQ. 0.0D0) THEN
				CoagArray(I,J) = 0.
				GOTO 40
			END IF

			IF (Particle1%EffectiveRadius   .EQ. 0. .OR. Particle2%EffectiveRadius   .EQ. 0. .OR. &
				Particle1%NumberOfParticles .EQ. 0.0D0 .OR. Particle2%NumberOfParticles .EQ. 0.0D0) THEN

				CoagArray(I,J) = 0.

			ELSE

				
				IF(ConstantKernel) THEN
					CoagArray(I,J) = CoagKernel

				ELSE
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					!! Add together compenents of the Coagulation Kernel (units cm3 part.-1 s-1) !!
					CoagArray(I,J) = (CoagKernBrownianDiff   (Particle1, Particle2) &
					                  + CoagKernGravCollection (Particle1, Particle2)) !! Add Brownian and Convective Contributions
					!write(*,*) 'CRL: CoagKernBrownianDiff   (Particle1, Particle2) =',CoagKernBrownianDiff   (Particle1, Particle2)
					!write(*,*) 'CRL: CoagKernGravCollection (Particle1, Particle2) =',CoagKernGravCollection (Particle1, Particle2)
					!write(*,*) ''

				END IF
			END IF
			
			IF(CoagArray(I,J) .LT. 0.0) WRITE(*,*) "CoagArray", I, J, CoagArray(I,J)
				
			!! Step the second particle (always either smaller than the larger or with a lower distribution id)
40			Particle2 => Particle2%Next
		END DO

		!! Advance the larger particle and reset the smaller
		Particle1 => Particle1%Next
		Particle2 => Particles%First

	END DO
	!write(*,*) 'CRL: After CoagArray =',CoagArray

!	WRITE(*,*) "Coag Check 1"


	!! Get average solution density
	!! Use this so that mass is conserved in transfers between particles
	Particle1 => Particles%First
	SolnSum =  0
	Count = 0
22	IF(Particle1%Numberofparticles .NE. 0.0) THEN
		SolnSum = SolnSum + Particle1%SolutionDensity
		Count = Count + 1
	END IF
	IF(ASSOCIATED(Particle1%Next)) THEN
		Particle1 => Particle1%Next
		GOTO 22
	END IF

	AvgSolutionDensity = SolnSum/Count
	
	!! Calculate single particle volume
	!!and volume concentrations for each component
	Particle1 => Particles%First
	JJ = 1
	InSumVolumeConc = 0.0
45	IF (Particle1%numberofparticles .GT. 0.0D0) THEN
		SingleParticleVolume(JJ) = (4.0*Pi/3.0)*(Particle1%effectiveradius)**3 !cm3
		
		NumConc(JJ) = Particle1%numberofparticles
		
		InitNumConc(JJ) = Particle1%numberofparticles 

		VolumeConc(JJ) = SingleParticleVolume(JJ)*Particle1%numberofparticles !cm3/cm3


		IF(VolumeConc(JJ) .LT. 0.0) WRITE(*,*) "Negative Volume Conc!"
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		
		!Water
		AqVolConc(JJ,1)  = Particle1%AqChems(1) * Particle1%NumberOfParticles*AqMolecularMass(1)/AvgSolutionDensity !Particle1%SolutionDensity

		!Solid Salts
		DO I = 2, HowManyAqChems
			AqVolConc(JJ,I)  = Particle1%AqChems(I) * Particle1%NumberOfParticles*AqMolecularMass(I)/SolidSaltDensity(I)
			IF(AqVolConc(JJ,I) .LT. 0.0) WRITE(*,*) "Negative AqVolume Conc!"	
			IF(AqVolConc(JJ,I) .GT. 1.0e30) WRITE(*,*) "Error", AqPhaseChemicalNames(I), AqVolConc(JJ,I)
		END DO

		!Inorganic Ions
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(JJ,I)  = Particle1%AqChems(I) * Particle1%NumberOfParticles*AqMolecularMass(I)/AvgSolutionDensity !Particle1%SolutionDensity
			IF(AqVolConc(JJ,I) .LT. 0.0) WRITE(*,*) "Negative AqVolume Conc!"	
		END DO

		!Convert org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(JJ,I)  = Particle1%OrgChems(I) * Particle1%NumberOfParticles*OrgMolecularMass(I)/OrgDensity(I)
			IF(OrgVolConc(JJ,I) .LT. 0.0) WRITE(*,*) "Negative OrgVolume Conc!"	
		END DO
	
		!Convert hydrophilic org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(JJ,I)  = Particle1%AqOrgChems(I) * Particle1%NumberOfParticles*AqOrgMolecularMass(I)/AqOrgDensity(I)
			IF(AqOrgVolConc(JJ,I) .LT. 0.0) THEN
				WRITE(*,*) "Negative AqOrgVolume Conc!", JJ, I, Particle1%NumberOfParticles, Particle1%AqOrgChems(I)
			END IF
		END DO
					

		TotalAqVolConc = 0.0
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalAqVolConc = TotalAqVolConc + AqVolConc(JJ,I)
		END DO

		TotalOrgVolConc = 0.0
		DO I = 1, HowManyOrgChems
			TotalOrgVolConc = TotalOrgVolConc + OrgVolConc(JJ,I)
		END DO

		TotalAqOrgVolConc = 0.0
		DO I = 1, HowManyAqOrgChems
			TotalAqOrgVolConc = TotalAqOrgVolConc + AqOrgVolConc(JJ,I)
		END DO

		InSumVolumeConc = InSumVolumeConc + (TotalAqVolConc+TotalOrgVolConc+TotalAqOrgVolConc)
	
	ELSE
		
		SingleParticleVolume(JJ) = (4.0*Pi/3.0)*(Particle1%edges(1)**3.0)
		
		NumConc(JJ) = 0.0
		InitNumConc(JJ) = 0.0
		
		VolumeConc(JJ) = 0.0
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(JJ,I)  = 0.0
		END DO
		!Convert organic concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(JJ,I)  = 0.0
		END DO
		!Convert aqueous organic concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(JJ,I)  = 0.0
		END DO
	END IF

	!WRITE(*,*) VolumeConc(JJ), NumConc(JJ)

	!! Go to the next particle
	IF (ASSOCIATED(Particle1%Next)) THEN
		Particle1 => Particle1%Next
		JJ = JJ+1
		GOTO 45
	END IF

!	WRITE(*,*) "Coag Check 2"

	InVolConc = 0.0
	DO JJ = 1, HowManyLinks
		InVolConc = InVolConc + VolumeConc(JJ)
	END DO
	!WRITE(*,*) "In: ", InVolConc
	
	!!Calculate the volume fractions 

	DO I = 1, HowManyLinks
		DO J = 1, HowManyLinks
			!Calculate volume of the combined particle
			BigV = SingleParticleVolume(I) + SingleParticleVolume(J)
		
			DO K = 1, HowManyLinks
				
				IF(K .EQ. HowManyLinks .AND. BigV .GE. SingleParticleVolume(K)) THEN
					VolumeFractions(I,J,K) = 1.0
				ELSE IF(K .LT. HowManyLinks .AND. BigV .GE. SingleParticleVolume(K) .AND. BigV .LT. SingleParticleVolume(K+1)) THEN
					VolumeFractions(I,J,K) = (SingleParticleVolume(K)/BigV)*((SingleParticleVolume(K+1)-BigV)/(SingleParticleVolume(K+1)-SingleParticleVolume(K)))
					IF(VolumeFractions(I,J,K) .GT. 1.0) THEN
						WRITE(*,*) "Volume Fraction out of bounds"
						STOP
					END IF
				ELSE IF(K .GT. 1 .AND. BigV .GE. SingleParticleVolume(K-1) .AND. BigV .LT. SingleParticleVolume(K)) THEN
					VolumeFractions(I,J,K) = 1-VolumeFractions(I,J,K-1)
				ELSE
					VolumeFractions(I,J,K) = 0.0
				END IF
			
			END DO
		END DO
	END DO

!	WRITE(*,*) "Coag Check 3"
	!Calculate change in total volume conc. and use that to get change in number
	DO K = 1, HowManyLinks    
        LowerSum=0
		DO J=1, HowManyLinks
			IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(K,J) .NE. 0.0 .AND. VolumeFractions(K,J,K) .NE. 1.0) THEN 
				LowerSum = LowerSum + (1-VolumeFractions(K,J,K))*CoagArray(K,J)*InitNumConc(J)
			END IF
		END DO
    
        UpperSum=0
		IF(K .EQ. 1) THEN
			UpperSum=0
		ELSE
			DO J=1,K
				InnerSum=0
				DO I=1, K-1
					IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(I,J) .NE. 0.0 .AND. VolumeFractions(I,J,K) .NE. 0.0) THEN 

						InnerSum = InnerSum + VolumeFractions(i,j,k)*CoagArray(i,j)*VolumeConc(i)*InitNumConc(J)


					END IF
				end DO 
				UpperSum = UpperSum + InnerSum
			end DO
		END IF    
		
		IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine, early")


		VolumeConc(K) = (VolumeConc(K)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)

		IF (VolumeConc(K) .GT. 0.0 .AND. SingleParticleVolume(K) .GT. 0.0) THEN
			!Number concentration estimated from change in volume concentration
			NumConc(K) = VolumeConc(K)/SingleParticleVolume(K)
		END IF

	END DO       

!	WRITE(*,*) "Coag Check 4"

	!Calculate change in volume conc. of individual Aq. species
	DO q = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		DO K = 1, HowManyLinks    
			LowerSum=0
			DO J=1, HowManyLinks
				IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(K,J) .NE. 0.0 .AND. VolumeFractions(K,J,K) .NE. 1.0) THEN 
					!WRITE(*,*) K, J, VolumeFractions(K,J,K), CoagArray(K,J), InitNumConc(J)
					LowerSum = LowerSum + (1.0-VolumeFractions(K,J,K))*CoagArray(K,J)*InitNumConc(J)
				END IF
			END DO
    
			UpperSum=0
			IF(K .EQ. 1) THEN
				UpperSum=0
			ELSE
				DO J=1,K
					InnerSum=0
					DO I=1, K-1
						IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(I,J) .NE. 0.0 .AND. VolumeFractions(I,J,K) .NE. 0.0) THEN 
							InnerSum = InnerSum + VolumeFractions(i,j,k)*CoagArray(i,j)*AqVolConc(i,q)*InitNumConc(J)
	
						END IF
					END DO 
					
					IF (InnerSum .LT. 0.0 ) CALL ERROR("Problem in InnerSum in Coag routine")
					UpperSum = UpperSum + InnerSum
				END DO
			END IF    
				
			IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine")
			IF (UpperSum .LT. 0.0 ) CALL ERROR("Problem in UpperSum in Coag routine")
	
			AqVolConc(K,q) = (AqVolConc(K,q)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)

			IF(AqVolConc(K,Q) .LT. 0.0) CALL ERROR("Negative inorganic concentrations from Coag Routine.")
			
		END DO  
	END DO

	!Calculate change in volume conc. of individual hydrophobic organic species
	DO q = 1, HowManyOrgChems
		DO K = 1, HowManyLinks    
				LowerSum=0
				DO J=1, HowManyLinks
					IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(K,J) .NE. 0.0 .AND. VolumeFractions(K,J,K) .NE. 1.0) THEN 
						!WRITE(*,*) K, J, VolumeFractions(K,J,K), CoagArray(K,J), NumConc(J)
						LowerSum = LowerSum + (1-VolumeFractions(K,J,K))*CoagArray(K,J)*InitNumConc(J)
					END IF
				END DO
    
				UpperSum=0
				IF(K .EQ. 1) THEN
					UpperSum=0
				ELSE
					DO J=1,K
						InnerSum=0
						DO I=1, K-1
							IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(I,J) .NE. 0.0 .AND. VolumeFractions(I,J,K) .NE. 0.0) THEN 
								InnerSum = InnerSum + VolumeFractions(i,j,k)*CoagArray(i,j)*OrgVolConc(i,q)*InitNumConc(J)
							END IF
						END DO 
						UpperSum = UpperSum + InnerSum
					END DO
				END IF    
				
				IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine")
				OrgVolConc(K,q) = (OrgVolConc(K,q)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)
				IF(OrgVolConc(K,Q) .LT. 0.0) CALL ERROR("Negative hydrophobic organic concentrations from Coag Routine.")
		END DO  
	END DO

	!Calculate change in volume conc. of individual hydrophilic organic species
	DO q = 1, HowManyAqOrgChems
		DO K = 1, HowManyLinks    
				LowerSum=0
				DO J=1, HowManyLinks
					IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(K,J) .NE. 0.0 .AND. VolumeFractions(K,J,K) .NE. 1.0) THEN 
						!WRITE(*,*) K, J, VolumeFractions(K,J,K), CoagArray(K,J), NumConc(J)
						LowerSum = LowerSum + (1-VolumeFractions(K,J,K))*CoagArray(K,J)*InitNumConc(J)
					END IF
				END DO
    
				UpperSum=0
				IF(K .EQ. 1) THEN
					UpperSum=0
				ELSE
					DO J=1,K
						InnerSum=0
						DO I=1, K-1
							IF (InitNumConc(J) .GT. 0.0D0 .AND. CoagArray(I,J) .NE. 0.0 .AND. VolumeFractions(I,J,K) .NE. 0.0) THEN 
								InnerSum = InnerSum + VolumeFractions(i,j,k)*CoagArray(i,j)*AqOrgVolConc(i,q)*InitNumConc(J)
							END IF
						END DO 
						UpperSum = UpperSum + InnerSum
					END DO
				END IF    
				
				IF (LowerSum .LT. 0.0 ) CALL ERROR("Problem in LowerSum in Coag routine")
				AqOrgVolConc(K,q) = (AqOrgVolConc(K,q)+TimeStep*UpperSum)/(1+TimeStep*LowerSum)
				IF(AqOrgVolConc(K,Q) .LT. 0.0) CALL ERROR("Negative aqueous organic concentrations from Coag Routine.")
			
		END DO  
	END DO

	OutSumVolumeConc = 0.0
	DO JJ = 1, HowManyLinks
		!Volume Conc Check
		TotalAqVolConc = 0.0
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalAqVolConc = TotalAqVolConc + AqVolConc(JJ,I)
			!IF(AqVolConc(JJ,I) .NE. 0.0) WRITE(*,*) JJ, AqPhaseChemicalNames(I), AqVolConc(JJ,I) 
		END DO

		TotalOrgVolConc = 0.0
		DO I = 1, HowManyOrgChems
			TotalOrgVolConc = TotalOrgVolConc + OrgVolConc(JJ,I)
		END DO

		TotalAqOrgVolConc = 0.0
		DO I = 1, HowManyAqOrgChems
			TotalAqOrgVolConc = TotalAqOrgVolConc + AqOrgVolConc(JJ,I)
		END DO

		SumVolConc(JJ) = TotalAqVolConc+TotalOrgVolConc+TotalAqOrgVolConc
		!WRITE(*,*) JJ, SumVolConc(JJ)
		OutSumVolumeConc = OutSumVolumeConc + (TotalAqVolConc+TotalOrgVolConc+TotalAqOrgVolConc)

	END DO

	OutVolConc = 0.0
	DO JJ = 1, HowManyLinks
		OutVolConc = OutVolConc + VolumeConc(JJ)
	END DO
	!WRITE(*,*) "Volume Ratio", OutVolConc/InVolConc
	IF(ABS((OutVolConc/InVolConc)-1.0) .GT. 0.0001) CALL ERROR("Error in StepSectionalCoagulationJacobson(). Total Volume concentration" &
			           //" is lost during step.")
	
	!WRITE(*,*) "Sum Volume Ratio", OutSumVolumeConc/InSumVolumeConc
	IF(ABS((OutSumVolumeConc/InSumVolumeConc)-1.0) .GT. 0.0001) CALL ERROR("Error in StepSectionalCoagulationJacobson(). Summed Total Volume concentration" &
			           //" is lost during step.")


!	WRITE(*,*) "Coag Check 5"

	!Save adjustments back to particle list
	!!Calculate single particle volume
	!!and volume concentrations for each component
	Particle1 => Particles%First
	JJ = 1
55	IF (NumConc(JJ) .GT. 0.0D0) THEN
		!WRITE(*,*) 'In Gt 0 Coauglation.f90 Particle1%AqChem(1) = ',Particle1%AqChems(1),NumConc(JJ),Particle1%NumberOfParticles

		Particle1%numberofparticles = NumConc(JJ)
		!Convert aqueous concentrations back to mol/particle
		IF(Particle1%SolutionDensity .EQ. 0.0) THEN
			Particle1%SolutionDensity = StoreDensity
		END IF
		StoreDensity = Particle1%SolutionDensity
		Particle1%AqChems(1) = AqVolConc(JJ,1)*AvgSolutionDensity/AqMolecularMass(1)/ Particle1%NumberOfParticles	
		!write(*,*) 'AqVolConc(JJ,1) = ',AqVolConc(JJ,1)
		!write(*,*) 'AvgSolutionDensity = ',AvgSolutionDensity 
		!write(*,*) 'AqMolecularMass(1) = ',AqMolecularMass(1) 
		!write(*,*) 'Particle1%NumberOfParticles = ',Particle1%NumberOfParticles 






		DO I = 2, HowManyAqChems
			Particle1%AqChems(I) = AqVolConc(JJ,I)*SolidSaltDensity(I)/AqMolecularMass(I)/ Particle1%NumberOfParticles		
		END DO
		
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Particle1%AqChems(I) = AqVolConc(JJ,I)*AvgSolutionDensity/AqMolecularMass(I)/ Particle1%NumberOfParticles			
		END DO
		
		!Convert organic concentrations back to mol/particle
		DO I = 1, HowManyOrgChems
			Particle1%OrgChems(I) = OrgVolConc(JJ,I)*OrgDensity(I)/OrgMolecularMass(I)/ Particle1%NumberOfParticles		
		END DO

		!Convert hydrophilic organic concentrations back to mol/particle
		DO I = 1, HowManyAqOrgChems
			Particle1%AqOrgChems(I) = AqOrgVolConc(JJ,I)*AqOrgDensity(I)/AqOrgMolecularMass(I)/ Particle1%NumberOfParticles		
		END DO

		!! This will do RADIUS, EffectiveRadius, & DENSITY
		!write(*,*) 'Before RecalculateRadius(Particle1) Particle1%AqChems(1) = ',Particle1%AqChems(1),NumConc(JJ)
		CALL RecalculateRadius(Particle1) !Calls ParticleDensity!
		!write(*,*) 'After RecalculateRadius(Particle1)Particle1%AqChems(1) = ',Particle1%AqChems(1),NumConc(JJ)

	ELSE
		!WRITE(*,*) 'In Else Coauglation.f90 Particle1%AqChem(1) = ',Particle1%AqChems(1),NumConc(JJ)
		Particle1%numberofparticles = 0.0D0
		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			Particle1%AqChems(I) = 0.0
		END DO
		
		DO I = 1, HowManyOrgChems
			Particle1%OrgChems(I) = 0.0
		END DO
		DO I = 1, HowManyAqOrgChems
			Particle1%AqOrgChems(I) = 0.0
		END DO
		!WRITE(*,*) 'In Coauglation.f90 Particle1%AqChem(1) = ',Particle1%AqChems(1),NumConc(JJ)
	END IF



	!! Go to the next particle
	IF (ASSOCIATED(Particle1%Next)) THEN
		Particle1 => Particle1%Next
		JJ = JJ+1
		GOTO 55
	END IF

!	WRITE(*,*) "Coag Check 6"

	!Regrid the aerosol to make sure sizes work out correctly
	!(Also calls aqueous equilibrium, and resets temperature)
	!write(*,*) 'Before RegridAerosol in Coagulation.f90 calling spillbeans...'
	!call spillbeans()
	CALL RegridAerosol()
	!write(*,*) 'After RegridAerosol in Coagulation.f90 calling spillbeans...'
	!call spillbeans()	
	!! Houseclean memory
	DEALLOCATE(CoagArray)
	DEALLOCATE(SingleParticleVolume)
	DEALLOCATE(VolumeFractions)
	DEALLOCATE(VolumeConc)
	DEALLOCATE(NumConc)
	DEALLOCATE(AqVolConc)
	DEALLOCATE(OrgVolConc)
	DEALLOCATE(AqOrgVolConc)
			
!	WRITE(*,*) "Coag Check 7"
	RETURN


   	!Particle1 => particles%first
	!DO WHILE(associated(Particle1))
	!	write(*,*) 'CRL: In StepSectionalCoagulation - current%AqChems = ',Particle1%AqChems
	!	write(*,*) 'CRL: In StepSectionalCoagulation - current%numberofparticles = ',Particle1%numberofparticles
	!	Particle1 => Particle1%next
	!END DO !END STEP PARTICLE DILUTION



END SUBROUTINE StepSectionalCoagulationJacobson



END MODULE Coagulation
