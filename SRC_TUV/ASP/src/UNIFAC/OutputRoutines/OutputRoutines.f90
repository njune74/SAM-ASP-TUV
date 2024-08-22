!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! OutputRoutines.f90
!! This file contains many useful ways for summarizing aerosol
!! data for output to ASCII text files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							    !!
!!									    !!
!! Month  Year   Name              Description				    !!
!! 07     2006   Matt Alvarado     Began Update History		            !!
!! 08/29  2006   Matt Alvarado     1. Changed TotalPOA and TotalSOA         !!
!!				     calculation in SizeDistOutput	    !!
!!				   2. Changed currentparticle to current    !!
!!				      in SizeDistOutput			    !!
!! 09/01  2006   Matt Alvarado     1. Added TotalOrgCar calculation to	    !!
!!				SizeDistOutput		                    !!
!! 09/07  2006   Matt Alvarado Updated SpillBeans to spill BoundaryParticles!!
!! 11/01  2006   Matt Alvarado     Added Optical properties to SpillBeans   !!
!! 10/03  2007   Matt Alvarado     Added lines to SpillBeans to help with   !!
!!						initializing CRM6           !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP       !!
!!                                 a one-box model or subroutine.           !!
!! 07/08  2013   Matt Alvarado     Created Subroutines GasOutput,           !!
!!                                 GasEnhanceOutput, and AerosolOutput      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES		                    !!
!! 1. ModelParameters			    !!
!! 2. InfrastructuralCode		    !!
!! 3. GridPointFields			    !!
!! 4. Chemistry				    !!
!! 5. Aerosols				    !!
!! 6. Condensation			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE EachAerosolOutput (OutputFileName)                     !!
!! 2. SUBROUTINE BinEntireDomainForPrinting (OutputFileName, inHowManyBins)
!! 3. SUBROUTINE SummarizeDomainMassForOutput (FileName, FirstCall)     !!
!! 4. SUBROUTINE BinLagrangianParticlesForOutput (XP,YP,ZP,NumberOfBins,!!
!!               SmallEdge,LargeEdge, FH, FirstComment)                 !!
!! 5. SUBROUTINE SpillBeans(XP, YP, ZP)                                 !!
!! 6. SUBROUTINE SizeDistOutput(XP, YP, ZP, t, Filename)                !!
!! 7. SUBROUTINE GasOutput(FH1, FH2, FH3,t,First, Thermo)               !!
!! 8. SUBROUTINE GasEnhanceOutput(FH, t, First)                         !!
!! 9. SUBROUTINE AerosolOutput(InParticle, FH, t, First,, Thermo)       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!
!! Modifications made by CMB to play nice with gfortran
!!

MODULE OutputRoutines

PRIVATE
PUBLIC ::         SpillBeans,                           &
		  SizeDistOutput,                       &
                  GasOutput,                            &
                  GasEnhanceOutput,                     &
                  AerosolOutput 

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Produces a data file that is appropriate      !!
!! for printing the aerosol distribution of the  !!
!! entire domain ensemble aerosol.               !!
SUBROUTINE EachAerosolOutput (OutputFileName)

	USE ModelParameters,	 ONLY : OutputDeckSubDir, micron, grams, moles

	USE Aerosols,            ONLY : Particles, Particle, AerosolPH
	
	USE Chemistry,           ONLY : HowManyAqChems,			&
					HowManyAqCations,		&
					HowManyAqAnions,		&
					AqPhaseChemicalNames,		&
					AqCationNames,AqAnionNames

	USE InfrastructuralCode, ONLY : WARN, ERROR, INT2STR, REAL2STR,	 &
					GetFileHandle, ReturnFileHandle

	IMPLICIT NONE

	CHARACTER *(*) :: OutputFileName
	TYPE(Particle), POINTER :: Current
	LOGICAL :: ChemicalPresent(HowManyAqChems+HowManyAqCations+HowManyAqAnions)
	INTEGER :: I, FH

	!! Fill the Logical Array "ChemicalPresent()" with trues if the
	!! chemical is found in at least one aerosol particle.  This will
	!! be used in a second pass through the particles for outputting 
	!! their contents.
	DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		ChemicalPresent(I) = .FALSE.
	END DO
	ChemicalPresent(1) = .TRUE. !! Assume water is present

		Current => Particles%First

10		IF (ASSOCIATED(Current)) THEN
			DO I = 2, HowManyAqChems+HowManyAqCations+HowManyAqAnions
				IF (Current%AqChems(I) .GT. 0.) ChemicalPresent(I) = .TRUE.
			END DO
		END IF

		IF (ASSOCIATED(Current%Next)) THEN
			Current => Current%Next
			GOTO 10
		END IF

!DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
!	Print *, I,ChemicalPresent(I) 
!END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Output the traits of each particle !!
	!! to a file                          !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FH = GetFileHandle()
    OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFileName)//".csv")
	WRITE (FH, '(a)', advance='no') "ParticleID, "//		&
		"Particle Distribution #,"// &
		"Number Of Particles,"//	 &
		"Effective Radius,"//		 &
		"Solution Embryo Radius,"//	 &
		"InsolCoreRadius,"//		 &
		"CosContactAngle,"//		 &
		"InsolCoreMass,"//		&
		"Solution Density,"//		 &
		"Particle Density,"//		 &
		"Surface Tension,"//		 &
		"PH,"//				&
		"Origin Pressure Level"

	DO I = 1, HowManyAqChems
		IF (ChemicalPresent(I)) &
		WRITE (FH, '(a)', advance='no') ","//TRIM(AqPhaseChemicalNames(I))
	END DO
	DO I = 1, HowManyAqCations
		IF (ChemicalPresent(I+HowManyAqChems)) &
		WRITE (FH, '(a)', advance='no') ","//TRIM(AqCationNames(I))
	END DO
	DO I = 1, HowManyAqAnions
		IF (ChemicalPresent(I+HowManyAqChems+HowManyAqCations)) &
		WRITE (FH, '(a)', advance='no') ","//TRIM(AqAnionNames(I))
	END DO

	WRITE (FH, '(a)') ""

		Current => Particles%First

20		IF (ASSOCIATED(Current)) THEN

		    WRITE (FH, '(a)', advance='no') &
                      TRIM(INT2STR(Current%ParticleID))//", "// &
	              TRIM(INT2STR(Current%ParticleDistribution))//", "// &
		      TRIM(REAL2STR(Current%NumberOfParticles,8))//", "

			
	            WRITE (FH, '(a)', advance='no') &
                      TRIM(REAL2STR(Current%EffectiveRadius/micron,8))//", "//&
		      TRIM(REAL2STR(Current%EmbryoRadius/micron,8))//", "// &
		      TRIM(REAL2STR(Current%SolutionDensity,8))//", "//	&
		      TRIM(REAL2STR(Current%ParticleDensity,8))//", "//	&
		      TRIM(REAL2STR(Current%SurfaceTension,8))//", "//	&
		      TRIM(REAL2STR(AerosolPH(Current)))//", "//	&
		      TRIM(REAL2STR(Current%OriginPressureLevel))

			DO I = 1, HowManyAqChems
				IF (ChemicalPresent(I)) &
				WRITE (FH, '(a)', advance='no') ","//TRIM(REAL2STR(Current%AqChems(I)/moles))
			END DO
			DO I = 1, HowManyAqCations
				IF (ChemicalPresent(I+HowManyAqChems)) &
				WRITE (FH, '(a)', advance='no') ","//TRIM(REAL2STR(Current%AqChems(I+HowManyAqChems)/moles))
			END DO
			DO I = 1, HowManyAqAnions
				IF (ChemicalPresent(I+HowManyAqChems+HowManyAqCations)) &
				WRITE (FH, '(a)', advance='no') ","//TRIM(REAL2STR(Current%AqChems(I+HowManyAqChems+HowManyAqCations)/moles))
			END DO

			WRITE (FH, '(a)') ""

		END IF


		IF (ASSOCIATED(Current%Next)) THEN
			Current => Current%Next
			GOTO 20
		END IF

	CLOSE(FH)
	CALL ReturnFileHandle(FH)

	RETURN
END SUBROUTINE EachAerosolOutput 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Produces a data file that is appropriate      !!
!! for printing the aerosol distribution of the  !!
!! entire domain ensemble aerosol.               !!
SUBROUTINE BinEntireDomainForPrinting (OutputFileName, inHowManyBins)

	USE Aerosols,            ONLY : Particles, Particle,		&
					ParticleArray,			&
					RecalculateRadius,		&
					HypotheticalElectrolyteConcentrations 

	USE Chemistry,		ONLY : HowManyAqChems,			&
					HowManyAqCations,		&
					HowManyAqAnions,		&
					AqMolecularMass,		&
					AqCationMass, AqAnionMass,     &
					AqPhaseChemicalNames,		&
					AqCationNames,AqAnionNames


	USE InfrastructuralCode, ONLY : WARN, ERROR, INT2STR, REAL2STR,	 &
					GetFileHandle, ReturnFileHandle
	
	USE ModelParameters,	 ONLY : OutputDeckSubDir, micron

	IMPLICIT NONE

	!! External Variables
	INTEGER	       :: inHowManyBins
	CHARACTER *(*) :: OutputFileName

	!! Internal Variables
	INTEGER ::  I, X,Y,Z, Bin,		&
		HowManyBins,			&
		HowManyAerosol,			&
		HowManyElectroEquiv,		&
		HowManyContents

	REAL*8  ::	SmallestRadius, LargestRadius,		&
				BinBoundaries(inHowManyBins+1),		&
				RadiusRatio,				&
				NumbHistogram(inHowManyBins),		&
			ElectrolyteEquiv(inHowManyBins,HowManyAqChems), &
MassHistogram(inHowManyBins,HowManyAqChems+HowManyAqCations+HowManyAqAnions), &
			HypotheticalElectroConcs(HowManyAqChems),	&
				MaxNum, MinNum,				&
				MaxEquiv, MinEquiv,			&
				MaxContent, MinContent

	LOGICAL  :: HypotheticalContent(HowManyAqChems), &
               ActualContent(HowManyAqChems+HowManyAqCations+HowManyAqAnions)

	TYPE(Particle), POINTER :: Current


	HowManyBins = inHowManyBins

	!! Make a first pass to determine the largest and smallest radius
	Current => Particles%First
	CALL RecalculateRadius(Current)
	SmallestRadius = Current%EffectiveRadius
	LargestRadius = Current%EffectiveRadius

		Current => Particles%First
		
10		IF (ASSOCIATED(Current)) THEN

			HowManyAerosol = HowManyAerosol + Current%NumberOfParticles

			IF (Current%EffectiveRadius .LT. SmallestRadius) SmallestRadius = Current%EffectiveRadius
			IF (Current%EffectiveRadius .GT. LargestRadius)  LargestRadius  = Current%EffectiveRadius

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
		END IF


	!! Check to see if the specified number of bins makes sense with the number of particles
	IF (HowManyAerosol/10. .LT. HowManyBins) THEN
		X = MAX (1, FLOOR(HowManyAerosol / 10.))
		HowManyBins = X
	END IF

	!! Get the volume ratio for the bins
	RadiusRatio = (LargestRadius / SmallestRadius) ** (1. / HowManyBins)

	!! And set the boundaries of the radius bins
	BinBoundaries(1) = SmallestRadius
	DO I = 1, HowManyBins
		BinBoundaries(I+1) = BinBoundaries(I) * RadiusRatio 
	END DO
	DO I = 1, HowManyBins+1
		BinBoundaries(I)   = BinBoundaries(I) / micron
	END DO

	!! Initialize histogram outputs
	DO I = 1, HowManyBins

		NumbHistogram(I) = 0
		
		DO X = 1, HowManyAqChems
			ElectrolyteEquiv(I,X) = 0.
		END DO

		DO X = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			MassHistogram(I,X) = 0.
		END DO
	END DO

	!! Tally everything

		Current => Particles%First
		
20		IF (ASSOCIATED(Current)) THEN

			!! The largest radius particle get places out by the inversion, so correct
			Bin = FLOOR (DLOG(Current%EffectiveRadius / SmallestRadius) / DLOG(RadiusRatio)) + 1
			IF (Bin .EQ. HowManyBins + 1) Bin = HowManyBins

			NumbHistogram(Bin) = NumbHistogram(Bin)+ANINT(Current%NumberOfParticles)
			
			!! Hypothetical electrolyte concentrations
			!! (translate to mg / m3)
			HypotheticalElectroConcs = HypotheticalElectrolyteConcentrations (Current)
			DO I = 1, HowManyAqChems
				ElectrolyteEquiv(Bin,I) = ElectrolyteEquiv(Bin,I) + Current%NumberOfParticles*(Current%AqChems(I)	&
										  + HypotheticalElectroConcs(I)) * &
										  AqMolecularMass(I) * 1.e12
			END DO

			!! Actual concentrations including ions
			!! (translate to mg / m3)
			DO I = 1, HowManyAqChems
				MassHistogram(Bin,I) = MassHistogram(Bin,I) + Current%AqChems(I)	&
									   * AqMolecularMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqCations
				MassHistogram(Bin,HowManyAqChems+I) =								&
										MassHistogram(Bin,HowManyAqChems+I) +		&
										Current%AqChems(HowManyAqChems+I) *			&
										AqCationMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqAnions
				MassHistogram(Bin,HowManyAqChems+HowManyAqCations+I) =							&
										MassHistogram(Bin,HowManyAqChems+HowManyAqCations+I) +	&
										Current%AqChems(HowManyAqChems+HowManyAqCations+I) *	&
										AqAnionMass(I) * 1.e12 * Current%NumberOfParticles
			END DO

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 20
			END IF
		END IF

	MaxNum		= NumbHistogram(1)
	MinNum		= NumbHistogram(1)
	MaxEquiv	= ElectrolyteEquiv(1,1)
	MinEquiv	= ElectrolyteEquiv(1,1)
	MaxContent	= MassHistogram(1,1)
	MinContent	= MassHistogram(1,1)

	!! Sum everything to stack the plots
	DO X = 1, HowManyBins

		DO Y = 1, HowManyAqChems-1
			DO Z = Y+1, HowManyAqChems
				IF(ElectrolyteEquiv(X,Y) .GT. 0.) ElectrolyteEquiv(X,Y) = ElectrolyteEquiv(X,Y) + ElectrolyteEquiv(X,Z)
		END DO ; END DO 

		!! And general rescaleing
		DO Y = 1, HowManyAqChems
			ElectrolyteEquiv(X,Y) = ElectrolyteEquiv(X,Y) / (DLOG10(BinBoundaries(X+1)) - DLOG10(BinBoundaries(X)))
		END DO

		DO Y = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions - 1
			DO Z = Y+1, HowManyAqChems + HowManyAqCations + HowManyAqAnions 
				IF(MassHistogram(X,Y) .GT. 0.) MassHistogram(X,Y) = MassHistogram(X,Y) + MassHistogram(X,Z)
		END DO ; END DO 

		!! And general rescaleing
		DO Y = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions - 1
			IF (MassHistogram(X,Y) .GT. 0) &
				MassHistogram(X,Y) = MassHistogram(X,Y) / (DLOG10(BinBoundaries(X+1)) - DLOG10(BinBoundaries(X)))
		END DO

		!! Find High and Lows for each type
		IF (MaxNum	   .LT. NumbHistogram(X))		MaxNum		= NumbHistogram(X)
		IF (MinNum	   .GT. NumbHistogram(X)											&
				.AND. NumbHistogram(X) .GT. 0.)		MinNum		= NumbHistogram(X)
		IF (MaxEquiv   .LT. ElectrolyteEquiv(X,1))	MaxEquiv	= ElectrolyteEquiv(X,1)
		IF (MaxContent .LT. MassHistogram(X,1))		MaxContent	= MassHistogram(X,1)

		DO Y = 1, HowManyAqChems
			IF (MinEquiv .GT. ElectrolyteEquiv(X,Y)										&
				.AND. ElectrolyteEquiv(X,Y) .GT. 0.) MinEquiv	= ElectrolyteEquiv(X,Y)
		END DO

		DO Y = 1, HowManyAqChems + HowManyAqCations + HowManyAqAnions 
			IF (MinContent .GT. MassHistogram(X,Y)										&
				.AND. MassHistogram(X,Y) .GT. 0.) MinContent	= MassHistogram(X,Y)
		END DO
	END DO

	!! Double the Radii into Diameters
	DO X = 1, HowManyBins+1
		BinBoundaries(X) = BinBoundaries(X)*2.
	END DO

	MaxNum			= 10.**CEILING(DLOG10(MaxNum))
	MinNum			= 10.**FLOOR  (DLOG10(MinNum))
	MaxEquiv		= 10.**CEILING(DLOG10(MaxEquiv))
	MinEquiv		= 10.**FLOOR  (DLOG10(MinEquiv))
	MaxContent		= 10.**CEILING(DLOG10(MaxContent))
	MinContent		= 10.**FLOOR  (DLOG10(MinContent))
	LargestRadius	= 10.**CEILING(DLOG10(LargestRadius*2./micron))
	SmallestRadius	= 10.**FLOOR  (DLOG10(SmallestRadius*2./micron))

	!! Figure out which content is there (don't assume one  
	!! bin is representative, although probably could)
	DO X = 1, HowManyAqChems
		HypotheticalContent(X) = .FALSE.
	END DO

	!! Figure out which species are present
	DO X = 1, HowManyBins
		DO Y = 1, HowManyAqChems
			IF (ElectrolyteEquiv(X,Y) .GT. 0.) HypotheticalContent(Y) = .TRUE.
	END DO ; END DO

	HowManyElectroEquiv = 0
	DO Y = 1, HowManyAqChems
		IF (HypotheticalContent(Y)) HowManyElectroEquiv = HowManyElectroEquiv + 1
	END DO

	DO X = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		ActualContent(X) = .FALSE.
	END DO

	DO X = 1, HowManyBins
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (MassHistogram(X,Y) .GT. 0.) ActualContent(Y) = .TRUE.
	END DO ; END DO

	HowManyContents     = 0
	DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		IF (ActualContent(Y)) HowManyContents = HowManyContents+1
	END DO

	!! Output the tallied data as a MATLAB file
	I = GetFileHandle()
    OPEN(UNIT=I, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFileName)//".m")

	!! The Number Histogram
	WRITE (I,'(a)') "% Output from MELAM produced by SUBROUTINE BinEntireDomainForPrinting()"
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') "clear MassHistogram ElectrolyteEquiv Radius NumHistogram ;"
	WRITE (I,'(a)') "clear C D Content ;"
	WRITE (I,'(a)') "clear ActualContentNames HypotheticalContentNames;"
	WRITE (I,'(a)') "clf;"
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') ""

	WRITE (I,'(a)') "% Number Histogram Information"
	WRITE (I,'(a)') "% bin number, radius, number conc"

	WRITE (I,'(a)') "NumHistogram = ["

	DO X = 1, HowManyBins

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !//","
		WRITE (I,'(a)')              TRIM(REAL2STR(DBLE(1.e-50)))//"; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !//","
		WRITE (I,'(a)')              TRIM(REAL2STR(MAX(1.e-50,NumbHistogram(X))))//"; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !", "
		WRITE (I,'(a)')				 TRIM(REAL2STR(MAX(1.e-50,NumbHistogram(X))))//";"

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !//","
		WRITE (I,'(a)')              TRIM(REAL2STR(DBLE(1.e-50)))//"; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !//","
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(DBLE(1.e-50)))

		IF (X .LT. HowManyBins) THEN
			WRITE (I,'(a)') "; "
		ELSE
			WRITE (I,'(a)') "]; "
		END IF
	END DO

	WRITE (I,'(a)') "subplot(3,1,1)"
	WRITE (I,'(a)') "loglog(NumHistogram(:,2),NumHistogram(:,3),'k');"
	WRITE (I,'(a)') "hold on;"
	WRITE (I,'(a)') "fill(NumHistogram(:,2),NumHistogram(:,3),1);"
	WRITE (I,'(a)') "AXIS(["//TRIM(REAL2STR(SmallestRadius))//" "//TRIM(REAL2STR(LargestRadius))//" "//	&
					TRIM(REAL2STR(MinNum))//" "//TRIM(REAL2STR(MaxNum))//" ])"
	WRITE (I,'(a)') "xlabel('Particle Diameter (microns)');"
	WRITE (I,'(a)') "ylabel('d N / d log D');"
	WRITE (I,'(a)') "title('Number Concentration');"

	WRITE (I,'(a)') "hold off;"

	!! The Hypothetical Electrolyte Content
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') ""

	WRITE (I,'(a)') "% Hypothetical Electrolyte Concentration Histogram Information"
	WRITE (I,'(a)') "% First cite the names:"
	WRITE (I,'(a)') ""


	WRITE (I,'(a)') ""
	WRITE (I,'(a)') "% bin number, radius, Mass Conc of Each (Hypothetical for Electrolytes) Content..."

	WRITE (I,'(a)') "ElectrolyteEquiv = ["

	DO X = 1, HowManyBins

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50"
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,ElectrolyteEquiv(X,Y)),8)) 
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,ElectrolyteEquiv(X,Y)),8)) 
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50"
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " ! ", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems
			IF (HypotheticalContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50"
		END DO

		IF (X .LT. HowManyBins) THEN
			WRITE (I,'(a)') "; "
		ELSE
			WRITE (I,'(a)') "]; "
		END IF
	END DO

	WRITE (I,'(a)') "subplot(3,1,2)"
	WRITE(I,'(a)',advance='no') "HypotheticalContentNames = {"
	DO X = 1, HowManyAqChems
		IF (HypotheticalContent(X)) THEN
			WRITE(I,'(a)',advance='no') "'"//TRIM(AqPhaseChemicalNames(X))//"' ; "
		END IF
	END DO
	WRITE(I,'(a)') "};"
	WRITE (I,'(a)') "C = [1:"//TRIM(INT2STR(HowManyElectroEquiv))//"];"
	WRITE (I,'(a)') "loglog(ElectrolyteEquiv(:,2),ElectrolyteEquiv(:,3:"//TRIM(INT2STR(2+HowManyElectroEquiv))//"),'k');"
	WRITE (I,'(a)') "AXIS(["//TRIM(REAL2STR(SmallestRadius))//" "//TRIM(REAL2STR(LargestRadius))//" "//TRIM(REAL2STR(MinEquiv))	&
					//" "//TRIM(REAL2STR(MaxEquiv))//" ])"
	WRITE (I,'(a)') "xlabel('Particle Diameter (microns)');"
	WRITE (I,'(a)') "ylabel('d Mass / d log D');"
	WRITE (I,'(a)') "title('Particle Mass Content as Reconstituted Electrolytes');"
	WRITE (I,'(a)') "hold on;"
	WRITE (I,'(a)') "H=fill(ElectrolyteEquiv(:,2),ElectrolyteEquiv(:,3:"//TRIM(INT2STR(2+HowManyElectroEquiv))//"),C);"
	WRITE (I,'(a)') "legend(H,HypotheticalContentNames);"
	WRITE (I,'(a)') "hold off;"

	!! The Actual Content
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') ""

	WRITE (I,'(a)') "% Electrolyte & Other Concentration Histogram Information"
	WRITE (I,'(a)') "% First cite the names:"
	WRITE (I,'(a)') ""

	!! Then the content
	WRITE (I,'(a)') ""
	WRITE (I,'(a)') "% bin number, radius, Mass Conc of Each (Hypothetical for Electrolytes) Content..."

	WRITE (I,'(a)') "Content = ["

	DO X = 1, HowManyBins
		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,MassHistogram(X,Y)),8)) !", "
		END DO
		WRITE (I,'(a)') "; "


		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " "//TRIM(REAL2STR(MAX(1.e-50,MassHistogram(X,Y)),8)) !", "
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
		END DO
		WRITE (I,'(a)') "; "

		WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
		WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X)))//" " !", "
		DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
		END DO

		IF (X .LT. HowManyBins) THEN
			WRITE (I,'(a)') "; "
		ELSE
			WRITE (I,'(a)') ""
			WRITE (I,'(a)',advance='no') TRIM(INT2STR(X))//" " !", "
			WRITE (I,'(a)',advance='no') TRIM(REAL2STR(BinBoundaries(X+1)))
			DO Y = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
				IF (ActualContent(Y)) WRITE (I,'(a)',advance='no') " 1.e-50 "
			END DO
			WRITE (I,'(a)') "]; "
		END IF
	END DO

	WRITE (I,'(a)') "subplot(3,1,3)"
	!! Put out names
	WRITE(I,'(a)',advance='no') "ActualContentNames = {"
	DO X = 1, HowManyAqChems
		IF (ActualContent(X)) WRITE(I,'(a)',advance='no') "'"//TRIM(AqPhaseChemicalNames(X))//"'; "
	END DO
	DO X = 1, HowManyAqCations
		IF (ActualContent(HowManyAqChems+X)) WRITE(I,'(a)',advance='no') "'"//TRIM(AqCationNames(X))//"'; "
	END DO
	DO X = 1, HowManyAqAnions
		IF (ActualContent(HowManyAqChems+HowManyAqCations+X)) WRITE(I,'(a)',advance='no') "'"//TRIM(AqAnionNames(X))//"'; "
	END DO
	WRITE (I,'(a)') "};"
	WRITE (I,'(a)') "D = [1:"//TRIM(INT2STR(HowManyContents))//"];"
	WRITE (I,'(a)') "loglog(Content(:,2),Content(:,3:"//TRIM(INT2STR(2+HowManyContents))//"),'k');"
	WRITE (I,'(a)') "AXIS(["//TRIM(REAL2STR(SmallestRadius))//" "//TRIM(REAL2STR(LargestRadius))//" "//		&
					TRIM(REAL2STR(MinContent))//" "//TRIM(REAL2STR(MaxContent))//" ])"
	WRITE (I,'(a)') "xlabel('Particle Diameter (microns)');"
	WRITE (I,'(a)') "ylabel('d Mass / d log D');"
	WRITE (I,'(a)') "title('Particle Mass Content');"
	WRITE (I,'(a)') "hold on;"
	WRITE (I,'(a)') "H = fill(Content(:,2),Content(:,3:"//TRIM(INT2STR(2+HowManyContents))//"),D);"
	WRITE (I,'(a)') "legend(H,ActualContentNames);"
	WRITE (I,'(a)') "hold off;"

	!! Close the file
	CLOSE(I)
	CALL ReturnFileHandle(I)
	
	RETURN
END SUBROUTINE BinEntireDomainForPrinting

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Splash everything into one bin by constituent !!
!! and print in text file, vaguely mirroring     !!
!! ISORROPIA's output type.                      !!
!!												 !!!!!!
!! FOR REPEATED CALLS, ONLY CALL AFTER EQUILIBRATION !!
!! OR IT WILL INCREASE THE NUMBER OF COLUMNS FROM ONE!!
!! CALL TO ANOTHER.									 !!
SUBROUTINE SummarizeDomainMassForOutput (FileName, FirstCall)

	USE GridPointFields,	ONLY :	GetGasChemVector,		&
					GetRelativeHumidity,	&
					GetTemp

	USE Aerosols,           ONLY :	Particles, Particle, &
									ParticleArray

	USE Chemistry,			ONLY :  HowManyAqChems,					&
									HowManyAqCations,				&
									HowManyAqAnions,				&
									AqMolecularMass,				&
									GasMolecularMass,				&
									AqCationMass, AqAnionMass,		&
									GasPhaseChemicalNames,			&
									AqPhaseChemicalNames,			&
									AqCationNames,AqAnionNames,		&
									HowManyGasChems

	USE ModelParameters,    ONLY :  Avogadro, OutputDeckSubDir

	USE InfrastructuralCode, ONLY : GetFileHandle,					&
									ReturnFileHandle,				&
									ERROR,							&
									REAL2STR

	IMPLICIT NONE

	!! External Variables
	!INTEGER,OPTIONAL         :: FileHandle
	CHARACTER *(*), OPTIONAL :: FileName
	LOGICAL, OPTIONAL :: FirstCall

	!! Internal Variables
	INTEGER :: I,  FH
	REAL*8  :: AerosolMass(HowManyAqChems+HowManyAqCations+HowManyAqAnions), &
			   GasMass(HowManyGasChems), GridPointGasVec(HowManyGasChems),	 &
			   RH, Temperature, NumbAerosol

	TYPE(Particle), POINTER :: Current


	!! Error check the inputs
	!IF (.NOT.PRESENT(FileName) .AND. .NOT.PRESENT(FileHandle) .OR.  &
	!         PRESENT(FileName) .AND.      PRESENT(FileHandle))      &
	!	CALL ERROR("Error in call to SummarizeDomainMassForOutput(), you must supply either a file handle or a file name "//	&
	!			   "and supplied either neither or both.  Please check the syntax.")

	!IF (PRESENT(FileHandle)) FH = FileHandle
	IF (PRESENT(FileName)) THEN
		FH = GetFileHandle()
		OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(FileName)//".txt")
	END IF

	!! Initialize
	DO I = 1, HowManyGasChems
		GasMass(I) = 0.
	END DO

	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		AerosolMass(I) = 0.
	END DO

	NumbAerosol = 0.
	RH          = 0.
	Temperature = 0.

	!! Tally everything

		GridPointGasVec = GetGasChemVector ()
		RH          = RH + GetRelativeHumidity()
		Temperature = Temperature + GetTemp()

		DO I = 1, HowManyGasChems
			GasMass(I) = GasMass(I) + GridPointGasVec(I) &
					     * GasMolecularMass(I) &
						 * 1.e12 / Avogadro
		END DO

		Current => Particles%First
10		IF (ASSOCIATED(Current)) THEN

			NumbAerosol = NumbAerosol + Current%NumberOfParticles
			DO I = 1, HowManyAqChems
				AerosolMass(I) =	AerosolMass(I) + Current%AqChems(I)	&
									* AqMolecularMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqCations
				AerosolMass(HowManyAqChems+I) =							&
									AerosolMass(HowManyAqChems+I) +		&
									Current%AqChems(HowManyAqChems+I) *	&
									AqCationMass(I) * 1.e12 * Current%NumberOfParticles
			END DO
			DO I = 1, HowManyAqAnions
				AerosolMass(HowManyAqChems+HowManyAqCations+I) =						 &
									AerosolMass(HowManyAqChems+HowManyAqCations+I) +	 &
									Current%AqChems(HowManyAqChems+HowManyAqCations+I) * &
									AqAnionMass(I) * 1.e12 * Current%NumberOfParticles
			END DO

			IF (ASSOCIATED(Current%Next)) THEN
				Current => Current%Next
				GOTO 10
			END IF
		END IF

	!! Scale RH, Temp, and Number Conc
	RH          = RH 
	Temperature = Temperature 
	NumbAerosol = NumbAerosol 

	!! Now You Have all of the data.  
	!! Output the tallied data as an EXCEL type comma delimited file.
	IF (PRESENT(FirstCall) .AND. FirstCall) THEN

		!! RH 
		WRITE (FH,'(a)',advance='no') "Mean RH, "

		!! Temperature
		WRITE (FH,'(a)',advance='no') "Mean Temperature, "

		!! Number content
		WRITE (FH,'(a)',advance='no') "Aerosol Numb (n/cm3), "
		
		!! Print names of content
		DO I = 1, HowManyGasChems
			IF (GasMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(GasPhaseChemicalNames(I))//" (g) (microg / m3), "
		END DO

		IF (AerosolMass(1) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(AqPhaseChemicalNames(1))//" (l) (microg / m3), "

		DO I = 2, HowManyAqChems
			IF (AerosolMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(AqPhaseChemicalNames(I))//" (aq) (microg / m3), "
		END DO
		DO I = 1, HowManyAqCations
			IF (AerosolMass(HowManyAqChems+I) .GT. 0.) &
				WRITE (FH,'(a)',advance='no') TRIM(AqCationNames(I))//" (aq) (microg / m3), "
		END DO
		DO I = 1, HowManyAqAnions
			IF (AerosolMass(HowManyAqChems+HowManyAqCations+I) .GT. 0.) &
				WRITE (FH,'(a)',advance='no') TRIM(AqAnionNames(I))//" (aq) (microg / m3), "
		END DO

		WRITE (FH,'(a)') ""
	END IF

	!! Report the environmental variables
	WRITE (FH,'(a)',advance='no') &
		TRIM(REAL2STR(RH,8))//", "//TRIM(REAL2STR(Temperature,8))//", "//TRIM(REAL2STR(NumbAerosol,8))//", "

	!! Report the current concentrations / mass loadings
	DO I = 1, HowManyGasChems
		IF (GasMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(REAL2STR(GasMass(I),8))//", "
	END DO
	DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		IF (AerosolMass(I) .GT. 0.) WRITE (FH,'(a)',advance='no') TRIM(REAL2STR(AerosolMass(I),8))//", "
	END DO
	WRITE (FH,'(a)') ""

	!! Close the file if we're controlling it locally
	IF (PRESENT(FileName)) THEN
		CLOSE(FH)
		CALL ReturnFileHandle(FH)
	END IF

END SUBROUTINE SummarizeDomainMassForOutput


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Take Lagrangian Particles and put them in set of !!
!! bins for outputting                              !!
!!                                                  !!
!! Endpoints are treated just as they are for the   !!
!! sectional definitions: the lowest bin is from    !!
!! zero to an endpoint and the largest is from an   !!
!! endpoint to infinity.                            !!
SUBROUTINE BinLagrangianParticlesForOutput (NumberOfBins,SmallEdge,LargeEdge, FH, FirstComment)

	USE InfrastructuralCode, ONLY : ERROR, REAL2STR
	USE Aerosols,            ONLY : Particles, Particle, RecalculateRadius
	USE ModelParameters, ONLY : LagrangianOrNot,km

	IMPLICIT NONE

	INTEGER :: NumberOfBins,FH
	REAL*8  :: SmallEdge, LargeEdge, BinRadiusRatio
	REAL*8  :: BinEdges(NumberOfBins,2),BinVol(NumberOfBins),BinNumb(NumberOfBins)
	INTEGER :: I
	TYPE(Particle), POINTER :: Runner
	CHARACTER *(*) :: FirstComment

	!! Exception Handle on Settings
	IF (.NOT.LagrangianOrNot) &
		CALL ERROR ("Called SCBinLagrangianParticlesForOutput() even though we are not tracking lagrangian " // &
					"particles on this run.  This is not allowed.")

	IF (NumberOfBins .LE. 1) &
		CALL ERROR ("SCBinLagrangianParticlesForOutput() is only to be used for number distributions of aerosol " // &
				    "(at least two bins).  You've asked for one or fewer bins to be calculated.  We don't accommodate that.")

	IF (NumberOfBins .EQ. 2 .AND. SmallEdge .NE. LargeEdge) &
		CALL ERROR ("If you want to use only two bins in SCBinLagrangianParticlesForOutput(), you have to set the small " // &
					"bin boundary and the large bin boundary equal.  The structure will go from zero to this dividor and " // &
					"then from that dividor to infinity (radius).")

	DO I = 1, NumberOfBins
		BinVol(I) = 0.
		BinNumb(I) = 0.
	END DO


	BinEdges(1,1) = 0.
	BinEdges(1,2) = SmallEdge

	IF (NumberOfBins .EQ. 2) THEN
		!! Calculate the Bin Edges
		BinEdges(2,1) = SmallEdge
		BinEdges(2,2) = 1.*km
	ELSE
		!! Calculate the radius ratio
		BinRadiusRatio = (LargeEdge/SmallEdge)**(1./(NumberOfBins-2.))
		DO I = 2, NumberOfBins-1
			BinEdges(I,1) = BinEdges(I-1,2)
			BinEdges(I,2) = BinEdges(I,1) * BinRadiusRatio
		END DO

		BinEdges(NumberOfBins-1,2) = LargeEdge !! just so this is exactly as asked with no rounding error
		BinEdges(NumberOfBins,1)   = LargeEdge
		BinEdges(NumberOfBins,2)   = 1.*km
	END IF


	!! Count up the bins
	Runner => Particles%First

100	IF (ASSOCIATED(Runner)) THEN

		CALL RecalculateRadius(Runner)

		DO I = 1, NumberOfBins
			IF (Runner%EffectiveRadius .GE. BinEdges(I,1) .AND. &
				Runner%EffectiveRadius .LT. BinEdges(I,2)) EXIT
		END DO

		BinVol(I)  = BinVol(I)  + Runner%EffectiveRadius**3.
		BinNumb(I) = BinNumb(I) + Runner%NumberOfParticles

	END IF
	
	IF (ASSOCIATED(Runner%Next)) THEN
		Runner => Runner%Next
		GOTO 100
	END IF

	WRITE (FH,'(a)',advance='no') TRIM(FirstComment)

	DO I = 1, NumberOfBins
		WRITE (FH,'(a)',advance='no') ","//TRIM(REAL2STR(BinNumb(I)))
	END DO

	DO I = 1, NumberOfBins
		IF (BinNumb(I) .GT. 0.) THEN
			BinVol(I) = (BinVol(I)/BinNumb(I))**(1./3.)
		ELSE
			BinVol(I) = (BinEdges(I,1)+BinEdges(I,2))/2.
		END IF

		WRITE (FH,'(a)',advance='no') ","//TRIM(REAL2STR(BinVol(I)))

	END DO

	WRITE (FH,'(a)') ""


	RETURN
END SUBROUTINE BinLagrangianParticlesForOutput

! CMB: Modified to accept an output unit
SUBROUTINE SpillBeans(Env, tout_unit)
	
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, INT2STR
	USE Aerosols,            ONLY : Particles, &
					Particle, &
					RecalculateRadius, &
					AerosolPh, &
					CurvatureCorrection, &
					OrgCurvatureCorrection, &
					BoundaryParticles, &
					HypotheticalElectrolyteConcentrations
	USE Chemistry,		 ONLY : AqPhaseChemicalNames, &
					OrgPhaseChemicalNames, &
					HowManyOrgChems, &
					HowManyAqOrgChems, &
					AqOrgPhaseChemicalNames, &
					AqMolecularMass, &
					OrgMolecularMass, &
					AqOrgMolecularMass, &
					HowManyAqChems, FindChem, &
					GasPhaseBackground
	USE ModelParameters,	 ONLY : micron, Avogadro, ProtonIndex
	USE GridPointFields,     ONLY : GetAirDensity, GetM, &
				        GridGasChem,GetRelativeHumidity

	Implicit None

	!External
	LOGICAL, OPTIONAL :: Env !True if we want the boundary distribution
        integer, optional :: tout_unit
    
	!Internal
	INTEGER		:: I, NumAqChems, NumOrgChems, A
	REAL :: Dens, ppbmCO, Mgas, RH
	REAL*8 :: SaltConc(HowManyAqChems)
	integer :: out_unit
	
	TYPE(Particle), POINTER :: CurrentParticle
	
	! cmb 
	if (present(tout_unit)) then
	    out_unit = tout_unit
	else
	    out_unit = 6    ! default to stdout
	endif
	! end cmb
	
			
	! spill the beans...
	IF (.NOT.(PRESENT(Env))) THEN
		CurrentParticle => particles%first
	ELSE
		IF (Env) THEN
			CurrentParticle => BoundaryParticles%first
		ELSE
			CurrentParticle => particles%first
		END IF
	END IF

	
        Dens = GetAirDensity() !g/cm3
	Mgas = GetM()
	RH = GetRelativeHumidity()
	WRITE(out_unit,*) "RH: ", RH
	WRITE(out_unit,*) "Dens: ", Dens
	ppbmCO = (GridGasChem(FindChem("CO",0)) - &
	          GasPhaseBackground(FindChem("CO",0)))*(1e9/Mgas)*(28/28.973)
	WRITE(out_unit,*) "ppbm CO: ", ppbmCO
			
	NumAqChems = size(currentparticle%AqChems)
	NumOrgChems = HowManyOrgChems
	I = 1
	99 if (associated(currentparticle)) then
		write(out_unit, *)""
		write(out_unit, *)"Number ",trim(INT2STR(I)),",  ",currentparticle%numberofparticles,currentparticle%sectional
		write(out_unit, *)"Effective Radius: ",trim(real2str(CurrentParticle%EffectiveRadius/micron)),",  (",trim(real2str(CurrentParticle%Edges(1)/micron))," to ",trim(real2str(CurrentParticle%Edges(2)/micron)),")"
		write(out_unit, *)"Embryo Radius: ",trim(real2str(CurrentParticle%EmbryoRadius/micron))
		write(out_unit, *)"Insoluble Radius: ",trim(real2str(CurrentParticle%InsolubleRadius/micron))
		write(out_unit, *)"Absorbing Core Radius: ",trim(real2str(CurrentParticle%AbsCoreRad/micron))
		write(out_unit, *)"Shell Real Ref. Ind. 550 nm",trim(real2str(currentparticle%ShellRealRefracInd(1)))
		write(out_unit, *)"Shell Imag. Ref. Ind. 550 nm",trim(real2str(currentparticle%ShellImagRefracInd(1)))
		write(out_unit, *)"Extinction Efficiency: 550 nm",trim(real2str(CurrentParticle%ExtEff(1)))
		write(out_unit, *)"Scattering Efficiency: 550 nm",trim(real2str(currentparticle%ScaEff(1)))
		write(out_unit, *)"BackScattering Efficiency: 550 nm",trim(real2str(currentparticle%BackScaEff(1)))
		write(out_unit, *)"Assymetry Parameter: 550 nm",trim(real2str(currentparticle%AssymParam(1)))

	
		write(out_unit, *)"Temperature ",trim(real2str(currentparticle%Temperature))
		write(out_unit, *)"Ionic Strength ",trim(real2str(currentparticle%IonicStr))
		write(out_unit, *)"Water Activity ",trim(real2str(currentparticle%WaterActivity))
		write(out_unit, *)"SolutionDensity ",trim(real2str(currentparticle%SolutionDensity))
		write(out_unit, *)"InsolubleDensity ",trim(real2str(currentparticle%InsolubleDensity))
		write(out_unit, *)"ParticleDensity ",trim(real2str(currentparticle%ParticleDensity))
		write(out_unit, *)"Ph ", trim(real2str(AerosolPh(currentparticle)))
		write(out_unit, *)"Curvature Correction ", trim(real2str(CurvatureCorrection(currentparticle)))
		write(out_unit, *)"OrgCurvature Correction ", trim(real2str(OrgCurvatureCorrection(currentparticle)))	
		write(out_unit, *)"Solution Surface Tension", trim(real2str(currentparticle%SurfaceTension))
		
		write(out_unit, *)"Temperature ", trim(real2str(currentparticle%temperature))
		
		if (currentparticle%numberofparticles .GT. 0.0) then
                    SaltConc = HypotheticalElectrolyteConcentrations (currentparticle)
                endif   
         
        write(out_unit, *) 'AqChems:'
        write(out_unit, *) '--------'
		DO A = 1, NumAqChems		
			!IF(currentparticle%aqchems(A) .GT. 1e-40) THEN
				write(out_unit,*) trim(int2str(A))," : ",AqPhaseChemicalNames(A) ," (mol/cm3 air): ", trim(real2str(currentparticle%aqchems(A)*currentparticle%numberofparticles))
                                !WRITE(*,*) currentparticle%aqchems(A)
				!IF(A .LE. HowManyAqChems) THEN
				!write(*,*) trim(int2str(A))," : ",AqPhaseChemicalNames(A) ," (ppbm): ", trim(real2str((SaltConc(A)+currentparticle%aqchems(A))*currentparticle%numberofparticles* &
				!																						AqMolecularMass(A)*1.0e9/Dens))
				!END IF
				!IF(A .LE. HowManyAqChems) THEN
				!write(*,*) trim(int2str(A))," : ",AqPhaseChemicalNames(A) ," (ppbm/ppbm CO): ", trim(real2str((SaltConc(A)+currentparticle%aqchems(A))*currentparticle%numberofparticles* &
				!																						AqMolecularMass(A)*1.0e9/Dens/ppbmCO))
				!END IF
			!END IF
		END DO
		write(out_unit, *) 'OrgChems:'
        write(out_unit, *) '---------'
		DO A = 1, NumOrgChems		
			!IF(currentparticle%Orgchems(A) .GT. 1e-40) THEN
				write(out_unit,*) trim(int2str(A))," : ",OrgPhaseChemicalNames(A) ," (mol/cm3 air): ", trim(real2str(currentparticle%Orgchems(A)*currentparticle%numberofparticles))
                                !WRITE(*,*) trim(int2str(A))," : ",currentparticle%Orgchems(A)," (molec/particle): ", trim(real2str(currentparticle%Orgchems(A)*Avogadro))
				!write(*,*) trim(int2str(A))," : ",OrgPhaseChemicalNames(A) ," (ppbm): ", trim(real2str(currentparticle%Orgchems(A)*currentparticle%numberofparticles* &
				!																							OrgMolecularMass(A)*1.0e9/Dens))
				!write(*,*) trim(int2str(A))," : ",OrgPhaseChemicalNames(A) ," (ppbm/ppbm CO): ", trim(real2str(currentparticle%Orgchems(A)*currentparticle%numberofparticles* &
				!																							OrgMolecularMass(A)*1.0e9/Dens/ppbmCO))

			!END IF
		END DO
		write(out_unit, *) 'AqOrgChems:'
        write(out_unit, *) '-----------'
		DO A = 1, HowManyAqOrgChems		
			!IF(currentparticle%AqOrgchems(A) .GT. 1e-40) THEN
				write(out_unit,*) trim(int2str(A))," : ",AqOrgPhaseChemicalNames(A) ," (mol/cm3 air): ", trim(real2str(currentparticle%AqOrgchems(A)*currentparticle%numberofparticles))
                !                WRITE(*,*) currentparticle%AqOrgchems(A)
				!write(*,*) trim(int2str(A))," : ",AqOrgPhaseChemicalNames(A) ," (ppbm): ", trim(real2str(currentparticle%AqOrgchems(A)*currentparticle%numberofparticles* &
				!																							AqOrgMolecularMass(A)*1.0e9/Dens))
				!write(*,*) trim(int2str(A))," : ",AqOrgPhaseChemicalNames(A) ," (ppbm/ppbm CO): ", trim(real2str(currentparticle%AqOrgchems(A)*currentparticle%numberofparticles* &
				!																							AqOrgMolecularMass(A)*1.0e9/Dens/ppbmCO))
			!END IF
		END DO
		I = I + 1

		if (associated(currentparticle%next)) then
			currentparticle => currentparticle%next
			goto 99
		end if
	end if !Spill Beans

END SUBROUTINE SpillBeans

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine outputs a file at a given time step that                  !!
!! Contains the information necessary to plot the mass and                   !!
!! number distributions.                                                     !!
!! First written by Matt Alvarado, 10/19/05.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE SizeDistOutput(t, Filename)
	
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, INT2STR, GetFileHandle

	USE Aerosols,            ONLY : Particles, &
					Particle, &
					RecalculateRadius, &
					AerosolPh

	USE Chemistry,		 ONLY : AqPhaseChemicalNames, &
					OrgPhaseChemicalNames, &
					GasPhaseChemicalNames, &
					HowManyOrgChems, &
					FindChem, &
					OrgMolecularMass, &
					AqMolecularMass, &
					HowManyAqOrgChems, &
					AqOrgMolecularMass, &
					AqOrgNumCarbon, &
                                        GasMolecularMass

	USE ModelParameters,	 ONLY : micron, OutputDeckSubDir, &
					Monodisperse, ThermoBulkMode, &
                                        Avogadro

	USE GridPointFields,	 ONLY : GridGasChem
	
	Implicit None

	!External
	REAL*8      :: t
	CHARACTER *(*), OPTIONAL :: FileName

	!Internal
	INTEGER		:: I, NumAqChems, NumOrgChems, A, J

	TYPE(Particle), POINTER :: Current
	Integer :: 	KNO3i, KCli, KHSO4i, &
			Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
			NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
			Cli, NO3i, SO4i, HSO4i, K2SO4i, PotassiumIndex, &
			K2SO4Index, FH, NumBins	, COi, Traceri, Cai, &
                        KOHi, NaOHi, Mgi, MgNO32i, MgCl2i, MgHSO42i, &
                        MgSO4i, MgOH2i, CaNO32i, CaCl2i, CaHSO42i, CaSO4i, &
                        CaOH2i, LEVi, CO2i
	REAL*8 :: TotalPOA, TotalSOA, TotalK, TotalNa, TotalNH4, TotalCl, &
	          TotalNO3, TotalSO4, DeltaLogRad, Mgas, TotalOrgCar, &
                  TotalMg, TotalCa		
					
	
	!Find chemical indices for later
	PotassiumIndex = FindChem("K+", 1)
	KNO3i = FindChem("KNO3", 1)
	KCli = FindChem("KCl", 1)
	KHSO4i = FindChem("KHSO4", 1)
	K2SO4Index = FindChem("K2SO4", 1)
	K2SO4i = K2SO4Index
        KOHi = FindChem("KOH", 1)
		
	Nai = FindChem("Na+", 1)
	NaNO3i = FindChem("NaNO3", 1)
	NaCli = FindChem("NaCl", 1)
	NaHSO4i = FindChem("NaHSO4", 1)
	Na2SO4i = FindChem("Na2SO4", 1)
        NaOHi = FindChem("NaOH", 1)

	Mgi = FindChem("Mg++", 1)
	MgNO32i = FindChem("Mg(NO3)2", 1)
	MgCl2i = FindChem("MgCl2", 1)
	MgHSO42i = FindChem("Mg(HSO4)2", 1)
	MgSO4i = FindChem("MgSO4", 1)
        MgOH2i = FindChem("Mg(OH)2", 1)

	Cai = FindChem("Ca++", 1)
	CaNO32i = FindChem("Ca(NO3)2", 1)
	CaCl2i = FindChem("CaCl2", 1)
	CaHSO42i = FindChem("Ca(HSO4)2", 1)
	CaSO4i = FindChem("CaSO4*2H2O", 1)
        CaOH2i = FindChem("Ca(OH)2", 1)
		
	NH3i = FindChem("NH3", 1)
	NH4i = FindChem("NH4+", 1)
	NH4NO3i = FindChem("NH4NO3", 1)
	NH4Cli = FindChem("NH4Cl", 1)
	NH4HSO4i = FindChem("NH4HSO4", 1)
	NH42SO4i = FindChem("(NH4)2SO4", 1)
	LEVi = FindChem("(NH4)3H(SO4)2", 1)

	Cli = FindChem("Cl-", 1)
	NO3i = FindChem("NO3-", 1)
	SO4i = FindChem("SO4--", 1)
	HSO4i = FindChem("HSO4-", 1)

	Traceri = FindChem("Trace", 0)
	COi = FindChem("CO", 0)
	CO2i = FindChem("CO2", 0)

	!Count # of bins, # of Aq. Chems
	Current => particles%first
	NumBins = 0
	I = 1
	7 if (associated(Current)) then
		NumBins = NumBins + 1
		NumAqChems = size(Current%AqChems)
		NumOrgChems = size(Current%OrgChems)
		IF (I .EQ. 2) THEN
			DeltaLogRad = log10(Current%Edges(2)*1.0e4) - log10(Current%Edges(1)*1.0e4)
		END IF
			
		if (associated(Current%next)) then
			Current => Current%next
			I = I + 1
			goto 7
		end if
	end if

	FH = GetFileHandle()
	
	OPEN(UNIT = FH, FILE=TRIM(OutputDeckSubDir)//TRIM(Filename), STATUS="REPLACE")
	WRITE(FH, FMT='(A4,A11)', ADVANCE='NO') "Time", " "
	WRITE(FH, FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(Traceri)
	WRITE(FH, FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(COi)
	WRITE(FH, FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(CO2i)

	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Edge 1", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Edge 2", " "
	
	WRITE(FH,FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Number", " "

	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "OrgCar", " "
			
	WRITE(FH,FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotOA", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
	WRITE(FH,FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
	WRITE(FH,FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
	WRITE(FH,FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotMg", " "
	WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "TotCa", " "
			
	WRITE(FH, FMT='(A1)') ""
		
	Current => particles%first
		DO I = 1, NumBins
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') t/60

                        !Do gas concentrations in ug/m3 to match aerosol
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') GridGasChem(Traceri)* &
                                              GasMolecularMass(Traceri)*1.0e12/Avogadro	
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') GridGasChem(COi)* &
                                              GasMolecularMass(COi)*1.0e12/Avogadro
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') GridGasChem(CO2i)* &
                                              GasMolecularMass(CO2i)*1.0e12/Avogadro
			
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Edges(1)*1.0e4
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Edges(2)*1.0e4
			
			IF(Current%effectiveradius*1.0e4 .GT. 0.0) THEN
				WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%effectiveradius*1.0e4
			ELSE
				IF (I .EQ. NumBins) THEN			
					IF((.not. Monodisperse) .AND. (.not. ThermoBulkMode)) THEN
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 10**(log10(Current%Edges(1)*1.0e4)+0.5*DeltaLogRad)
					ELSE
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
					END IF
				ELSE
					IF((.not. Monodisperse) .AND. (.not. ThermoBulkMode)) THEN
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 10**(log10(Current%Edges(2)*1.0e4)-0.5*DeltaLogRad)
					ELSE
						WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
					END IF
				END IF
			END IF
			
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%numberofparticles
						
			TotalOrgCar = 0.0
			DO J = 2, NumOrgChems
				TotalOrgCar = TotalOrgCar + Current%AqOrgChems(J-1)*AqOrgNumCarbon(J-1)*Current%numberofparticles*1.0e12*12.0
				TotalOrgCar = TotalOrgCar + Current%OrgChems(J)*AqOrgNumCarbon(J-1)*Current%numberofparticles*1.0e12*12.0
			END DO
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalOrgCar

			!BC Mass Conc (ug/m3)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Orgchems(1)*Current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
			!OA Mass Conc (ug/m3)
			TotalPOA = 0.
			DO J = 2, NumOrgChems
				TotalPOA = TotalPOA + Current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + Current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA
	
			!Water Mass Conc. (ug/m3)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') Current%Aqchems(1)*Current%numberofparticles*1.0e12*AqMolecularMass(1)
					
			!K+ Mass Conc (ug/m3)
			TotalK = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i .OR. J .EQ. KOHi) THEN
					TotalK = TotalK + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. K2SO4Index) THEN
					TotalK = TotalK + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalK = TotalK*AqMolecularMass(PotassiumIndex)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
			!Na+ Mass Conc (ug/m3)
			TotalNa = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i .OR. J .EQ. NaOHi) THEN
					TotalNa = TotalNa + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. Na2SO4i) THEN
					TotalNa = TotalNa + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNa = TotalNa*AqMolecularMass(Nai)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
			!NH4+ Mass Conc (ug/m3)
			TotalNH4 = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
					TotalNH4 = TotalNH4 + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. NH42SO4i) THEN
					TotalNH4 = TotalNH4 + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
                                ELSE IF (J .EQ. LEVi) THEN
                                        TotalNH4 = TotalNH4 + 3*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

			!Cl- Mass Conc (ug/m3)
			TotalCl = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. MgCl2i .OR. J .EQ. CaCl2i) THEN
						TotalCl = TotalCl + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
                                END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalCl = TotalCl*AqMolecularMass(Cli)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

			!NO3- Mass Conc (ug/m3)
			TotalNO3 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
					TotalNO3 = TotalNO3 + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				ELSE IF (J .EQ. MgNO32i .OR. J .EQ. CaNO32i) THEN
					TotalNO3 = TotalNO3 + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
                                END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

			!SO4- Mass Conc (ug/m3)
			TotalSO4 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
				    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i &
                                    .OR. J .EQ. CaSO4i .OR. J .EQ. MgSO4i) THEN
					TotalSO4 = TotalSO4 + Current%Aqchems(J)*Current%numberofparticles*1.0e12
                                ELSE IF (J .EQ. LEVi .OR. J .EQ. MgHSO42i .OR. J .EQ. CaHSO42i) THEN
					TotalSO4 = TotalSO4 + 2*Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

			
			!Mg++ Mass Conc (ug/m3)
			TotalMg = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Mgi .OR. J .EQ. MgNO32i .OR. J .EQ. MgCl2i .OR. J .EQ. MgHSO42i .OR. J .EQ. MgOH2i &
                                     .OR. J .EQ. MgSO4i) THEN
					TotalMg = TotalMg + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalMg = TotalMg*AqMolecularMass(Mgi)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalMg

			!Ca++ Mass Conc (ug/m3)
			TotalCa = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Cai .OR. J .EQ. CaNO32i .OR. J .EQ. CaCl2i .OR. J .EQ. CaHSO42i .OR. J .EQ. CaOH2i &
                                     .OR. J .EQ. CaSO4i) THEN
					TotalCa = TotalCa + Current%Aqchems(J)*Current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalCa = TotalCa*AqMolecularMass(Cai)
			WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalCa

			WRITE(FH, FMT='(A1)') ""
					
			IF (associated(Current%next)) THEN
				Current => Current%next
			ELSE
				CYCLE
			END IF

		END DO

		CLOSE(UNIT = FH)


END SUBROUTINE SizeDistOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine writes gas phase output to text files in units of ppbv.   !!
!! First written by Matt Alvarado, 20130709.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GasOutput(FH1, FH2, FH3,t,First, Thermo)

	USE GridPointFields,	 ONLY : GridGasChem, &
					GetM
	USE Chemistry,	         ONLY : GasPhaseChemicalNames, &
                                        HowManyGasChems

	Implicit None

        !EXTERNAL Variables
        INTEGER, INTENT(IN) :: FH1, FH2, FH3
        REAL,    INTENT(IN) :: t
        LOGICAL, INTENT(IN) :: First, Thermo

        !INTERNAL Variables
        REAL                :: Mgas
        INTEGER             :: q

	IF (First) THEN !On first call, write header line
 
		IF (Thermo) THEN
		            !write(*, *)'--------------------------- DERPPPPPPPPPPPPPPPPP --------'
                   WRITE(FH1, FMT='(A4,A12)', ADVANCE='NO') "RH  ", " "
                   WRITE(FH2, FMT='(A4,A12)', ADVANCE='NO') "RH  ", " "
                   WRITE(FH3, FMT='(A4,A12)', ADVANCE='NO') "RH  ", " "
                ELSE
                   WRITE(FH1, FMT='(A4,A12)', ADVANCE='NO') "Time", " "
                   WRITE(FH2, FMT='(A4,A12)', ADVANCE='NO') "Time", " "
                   WRITE(FH3, FMT='(A4,A12)', ADVANCE='NO') "Time", " "
                ENDIF
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(FH1,FMT='(A16)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(FH2,FMT='(A16)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(FH3,FMT='(A16)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO                

                WRITE(FH1, FMT='(A1)') ""
		WRITE(FH2, FMT='(A1)') ""
		WRITE(FH3, FMT='(A1)') ""
         ELSE
		WRITE(FH1,FMT='(ES16.5E2)', ADVANCE='NO') t
		WRITE(FH2,FMT='(ES16.5E2)', ADVANCE='NO') t
		WRITE(FH3,FMT='(ES16.5E2)', ADVANCE='NO') t

		Mgas = GetM()
                Print *,'time = ',t
		DO q = 1, HowManyGasChems		
                   IF (q .LE. 70) THEN
                      WRITE(FH1,FMT='(ES16.5E2)', ADVANCE='NO') GridGasChem(q)*1.0e9/Mgas !from molec/cm3 to ppb
                      Print *, GasPhaseChemicalNames(q),' = ', GridGasChem(q)*1.0e9/Mgas
                   ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
                      WRITE(FH2,FMT='(ES16.5E2)', ADVANCE='NO') GridGasChem(q)*1.0e9/Mgas
                      Print *, GasPhaseChemicalNames(q),' = ', GridGasChem(q)*1.0e9/Mgas
                   ELSE 
                      WRITE(FH3,FMT='(ES16.5E2)', ADVANCE='NO') GridGasChem(q)*1.0e9/Mgas
                      Print *, GasPhaseChemicalNames(q),' = ', GridGasChem(q)*1.0e9/Mgas
                   END IF
                END DO
		
		WRITE(FH1, FMT='(A1)') ""
		WRITE(FH2, FMT='(A1)') ""
		WRITE(FH3, FMT='(A1)') ""
          ENDIF

END SUBROUTINE GasOutput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine writes gas phase molar enhancement ratios of selected     !!
!! species to a text files  in units of ppbv.                                !!
!! First written by Matt Alvarado, 20130709.                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GasEnhanceOutput(FH, t, First)      

	USE GridPointFields,	 ONLY : GridGasChem, &
					GetM
	USE Chemistry,	         ONLY : GasPhaseChemicalNames, &
                                        HowManyGasChems,       &
                                        GasPhaseBackground,    &
                                        FindChem
	Implicit None

        !EXTERNAL Variables
        INTEGER, INTENT(IN) :: FH
        REAL,    INTENT(IN) :: t
        LOGICAL, INTENT(IN) :: First

        !INTERNAL Variables
        REAL                :: Mgas, del_CO, del_PAN, del_AN
        INTEGER             :: q, COi, O3i, HNO3i, NOi, NO2i, &
                               PANi, N2O5i, HONOi, ANi, NO3i

	IF (First) THEN !On first call, write header line
		WRITE(FH, FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dO3/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dNOx/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dHNO3/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dPAN/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dTotPNs/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dTotANs/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dNO3/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dN2O5/dCO"
		WRITE(FH, FMT='(A13)', ADVANCE='NO') "dHONO/dCO"
                WRITE(FH, FMT='(A1)') ""
         ELSE
                COi = FindChem("CO", 0)
                O3i = FindChem("O3", 0)
                NOi = FindChem("NO", 0)
                NO2i = FindChem("NO2", 0)
                HNO3i = FindChem("HNO3", 0)
                PANi =  FindChem("PAN", 0)              
                ANi =  FindChem("CH3NO3", 0)
                NO3i = FindChem("NO3", 0)
                N2O5i = FindChem("N2O5", 0)
                HONOi = FindChem("HONO", 0)

                Mgas = GetM()

                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') t
                Print *,'time = ',t

                del_CO = (GridGasChem(COi)*1e9/Mgas)-(GasPhaseBackground(COi)*1e9/Mgas)
                !Print *,'del_CO_convertedbefore = ', del_CO
                
                del_CO = GridGasChem(COi)-GasPhaseBackground(COi)
                !Print *, 'del_CO_unconverted = ', del_CO
                !Print *, 'del_CO_convereted after = ',del_CO*1e9/Mgas


                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(O3i)-GasPhaseBackground(O3i))/del_CO
                !Print *, 'GridGasChem(O3i)= ',GridGasChem(O3i)*1e9/Mgas
                !Print *, 'GasPhaseBackground(O3i) =',GasPhaseBackground(O3i)*1e9/Mgas

                !Print *, 'del_O3_after = ', (GridGasChem(O3i)-GasPhaseBackground(O3i))*1e9/Mgas                
                !Print *, 'del_O3_before = ', ((GridGasChem(O3i)*1e9/Mgas)-(GasPhaseBackground(O3i)*1e9/Mgas))

                !Print *, 'del_O3/del_CO_converted = ', ((GridGasChem(O3i)*1e9/Mgas)-(GasPhaseBackground(O3i)*1e9/Mgas))/(del_CO*1e9/Mgas)
                !Print *, 'del_O3/del_CO_RATIO = ', (GridGasChem(O3i)-GasPhaseBackground(O3i))/del_CO


                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(NOi)-GasPhaseBackground(NOi) &
                                                           +GridGasChem(NO2i)-GasPhaseBackground(NO2i))/del_CO
                !Print *, '(GridGasChem(NOi) = ',GridGasChem(NOi)*1e9/Mgas
                !Print *, '-GasPhaseBackground(NOi) = ',GasPhaseBackground(NOi)*1e9/Mgas
                !Print *, 'GridGasChem(NO2i) = ', GridGasChem(NO2i)*1e9/Mgas
                !Print *, 'GasPhaseBackground(NO2i)) = ', GasPhaseBackground(NO2i)*1e9/Mgas
                !Print *, 'del_NOx_after = ',(GridGasChem(NOi)-GasPhaseBackground(NOi) &
                !                                           +GridGasChem(NO2i)-GasPhaseBackground(NO2i))*1e9/Mgas

                !Print *, 'del_NOx_before = ',((GridGasChem(NOi)*1e9/Mgas)-(GasPhaseBackground(NOi)*1e9/Mgas) &
                !                                           +(GridGasChem(NO2i)*1e9/Mgas)-(GasPhaseBackground(NO2i)*1e9/Mgas))



                !Print *, 'del_NOx_del_CO_convert = ',(((GridGasChem(NOi)*1e9/Mgas)-(GasPhaseBackground(NOi)*1e9/Mgas)) &
                !                                           +((GridGasChem(NO2i)*1e9/Mgas)-(GasPhaseBackground(NO2i)*1e9/Mgas)))/(del_CO*1e9/Mgas)


                !Print *, 'del_NOx_del_CO_ratio = ', (GridGasChem(NOi)-GasPhaseBackground(NOi) &
                !                                           +GridGasChem(NO2i)-GasPhaseBackground(NO2i))/del_CO

                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(HNO3i)-GasPhaseBackground(HNO3i))/del_CO
                !Print *, 'GridGasChem(HNO3i) = ',GridGasChem(HNO3i)*1e9/Mgas
                !Print *, 'GasPhaseBackground(HNO3i) = ',GasPhaseBackground(HNO3i)*1e9/Mgas
                !Print *, 'del_HNO3i = ', ((GridGasChem(HNO3i)*1e9/Mgas)-(GasPhaseBackground(HNO3i)*1e9/Mgas))
                !Print *, 'del_HNO3_del_CO_converted = ',((GridGasChem(HNO3i)*1e9/Mgas)-(GasPhaseBackground(HNO3i)*1e9/Mgas))/(del_CO*1e9/Mgas)
                !Print *, 'del_HNO3_del_CO_ratio = ',(GridGasChem(HNO3i)-GasPhaseBackground(HNO3i))/del_CO


                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(PANi)-GasPhaseBackground(PANi))/del_CO
                !Print *,'(GridGasChem(PANi) = ',GridGasChem(PANi)*1e9/Mgas
                !Print *,'GasPhaseBackground(PANi) = ',GasPhaseBackground(PANi)*1e9/Mgas
                !Print *, 'del_PAN = ', (GridGasChem(PANi)*1e9/Mgas)-(GasPhaseBackground(PANi)*1e9/Mgas)
                !Print *, 'del_PAN/del_CO_converted = ', ((GridGasChem(PANi)*1e9/Mgas)-(GasPhaseBackground(PANi)*1e9/Mgas))/(del_CO*1e9/Mgas)
                !Print *, 'del_PAN/del_CO_ratio = ', (GridGasChem(PANi)-GasPhaseBackground(PANi))/del_CO


                del_PAN = 0.0
                DO q = 0,32
                     del_PAN = del_PAN+GridGasChem(PANi+q)-GasPhaseBackground(PANi+q)
                END DO 
                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') del_PAN/del_CO
                !Print *, 'del_PN = ', del_PAN*1e9/Mgas
                !Print *, 'dal_PN/del_CO_converted', (del_PAN*1e9/Mgas)/(del_CO*1e9/Mgas)
                !Print *, 'del_PN/del_CO = ',del_PAN/del_CO


                del_AN = 0.0
                DO q = 0,76
                     del_AN = del_AN+GridGasChem(ANi+q)-GasPhaseBackground(ANi+q)
                END DO 
                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') del_AN/del_CO
                !Print *, 'del_AN = ',del_AN*1e9/Mgas
                !Print *, 'del_AN/del_CO_converted = ',(del_AN*1e9/Mgas)/(del_CO*1e9/Mgas)
                !Print *, 'del_AN/del_CO_ratio = ',del_AN/del_CO

                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(NO3i)-GasPhaseBackground(NO3i))/del_CO
                !Print *, 'GridGasChem(NO3i) = ',GridGasChem(NO3i)*1e9/Mgas
                !Print *, 'GasPhaseBackground(NO3i) = ',GasPhaseBackground(NO3i)*1e9/Mgas
                !Print *, 'del_NO3i = ', ((GridGasChem(NO3i)*1e9/Mgas)-(GasPhaseBackground(NO3i)*1e9/Mgas))
                !Print *, 'del_NO3i/del_CO_converted = ', ((GridGasChem(NO3i)*1e9/Mgas)-(GasPhaseBackground(NO3i)*1e9/Mgas))/(del_CO*1e9/Mgas)
                !Print *, 'del_NO3i/del_CO_ratio = ', (GridGasChem(NO3i)-GasPhaseBackground(NO3i))/del_CO


                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(N2O5i)-GasPhaseBackground(N2O5i))/del_CO
                !Print *, 'GridGasChem(N2O5i) = ',GridGasChem(N2O5i)*1e9/Mgas
                !Print *, 'GasPhaseBackground(N2O5i) = ',GasPhaseBackground(N2O5i)*1e9/Mgas
                !Print *, 'del_N2O5i = ', ((GridGasChem(N2O5i)*1e9/Mgas)-(GasPhaseBackground(N2O5i)*1e9/Mgas))
                !Print *, 'del_N2O5i/del_CO_converted = ', ((GridGasChem(N2O5i)*1e9/Mgas)-(GasPhaseBackground(N2O5i)*1e9/Mgas))/(del_CO*1e9/Mgas)
                !Print *, 'del_N2O5i/del_CO_ratio = ', (GridGasChem(N2O5i)-GasPhaseBackground(N2O5i))/del_CO


                WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') (GridGasChem(HONOi)-GasPhaseBackground(HONOi))/del_CO
                !Print *, 'GridGasChem(HONOi) = ',GridGasChem(HONOi)*1e9/Mgas
                !Print *, 'GasPhaseBackground(HONOi) = ',GasPhaseBackground(HONOi)*1e9/Mgas
                !Print *, 'del_HONOi = ',(GridGasChem(HONOi)*1e9/Mgas)-(GasPhaseBackground(HONOi)*1e9/Mgas)
                !Print *, 'del_HONOi/del_CO_converted = ',((GridGasChem(HONOi)*1e9/Mgas)-(GasPhaseBackground(HONOi)*1e9/Mgas))/(del_CO*1e9/Mgas)
                !Print *, 'del_HONOi/del_CO_ratio = ',(GridGasChem(HONOi)-GasPhaseBackground(HONOi))/del_CO
		WRITE(FH, FMT='(A1)') ""
         ENDIF
END SUBROUTINE GasEnhanceOutput    


SUBROUTINE AerosolOutput(InParticle, FH, t, First, Thermo)

	USE Chemistry,		 ONLY : AqPhaseChemicalNames, &
					OrgPhaseChemicalNames, &
					AqOrgPhaseChemicalNames, &
					FindChem, &
					OrgMolecularMass, &
					AqMolecularMass, &
					AqOrgMolecularMass

	USE Aerosols,            ONLY : Particle, &
					AerosolPh
	Implicit None

	!External
	TYPE(Particle), POINTER  :: InParticle
	REAL*8,  INTENT(IN)      :: t
	INTEGER, INTENT(IN)      :: FH
        LOGICAL, INTENT(IN)      :: First, Thermo

	!Internal
	INTEGER	:: J, NumAqChems, NumOrgChems, NumAqOrgChems, &
                   KNO3i, KCli, KHSO4i, &
		   Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
		   NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
		   Cli, NO3i, SO4i, HSO4i, K2SO4i, PotassiumIndex, &
		   Cai, KOHi, NaOHi, Mgi, MgNO32i, MgCl2i, MgHSO42i, &
                   MgSO4i, MgOH2i, CaNO32i, CaCl2i, CaHSO42i, CaSO4i, &
                   CaOH2i, LEVi

	REAL*8  :: TotalPOA, TotalSOA, TotalK, TotalNa, TotalNH4, TotalCl, &
	           TotalNO3, TotalSO4, TotalMg, TotalCa

        NumAqChems = size(InParticle%AqChems)
        NumOrgChems = size(InParticle%OrgChems)
        NumAqOrgChems = size(InParticle%AqOrgChems)

	IF (First) THEN !On first call, write header line
           
           IF (.NOT.(Thermo)) THEN
               WRITE(FH, FMT='(A4,A7)', ADVANCE='NO') "Time", " "
               !WRITE(FH, FMT='(A4,A11)', ADVANCE='NO') "Time", " "
           ELSE
               WRITE(FH, FMT='(A2,A9)', ADVANCE='NO') "RH",""
               !WRITE(FH, FMT='(A4,A11)', ADVANCE='NO') "Time", " "
           ENDIF              
           !WRITE(FH,FMT='(A6, A7)', ADVANCE='NO') "Number", " "
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "Number",""
           !WRITE(FH,FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
           WRITE(FH,FMT='(A9,A2)', ADVANCE='NO') "EffRadius",""
           !WRITE(FH,FMT='(A10, A3)', ADVANCE='NO') "Ionic.Str.", " "
           WRITE(FH,FMT='(A8,A3)', ADVANCE='NO') "IonicStr",""
           !WRITE(FH,FMT='(A8, A5)', ADVANCE='NO') "H2O.Act.", " "
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "H2OAct",""
           WRITE(FH,FMT='(A7,A4)', ADVANCE='NO') "SolDens",""
           WRITE(FH,FMT='(A8,A3)', ADVANCE='NO') "PartDens",""
           WRITE(FH,FMT='(A2,A9)', ADVANCE='NO') "Ph",""
           WRITE(FH,FMT='(A4,A7)', ADVANCE='NO') "Temp"	,""		
           WRITE(FH,FMT='(A5,A6)', ADVANCE='NO') "TotBC",""
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "TotPOA",""
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "TotSOA",""
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "TotH2O",""
           WRITE(FH,FMT='(A4,A7)', ADVANCE='NO') "TotK",""
           WRITE(FH,FMT='(A5,A6)', ADVANCE='NO') "TotNa",""
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "TotNH4",""
           WRITE(FH,FMT='(A5,A6)', ADVANCE='NO') "TotCl",""
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "TotNO3",""
           WRITE(FH,FMT='(A6,A5)', ADVANCE='NO') "TotSO4",""
           WRITE(FH,FMT='(A5,A6)', ADVANCE='NO') "TotMg",""
           WRITE(FH,FMT='(A5,A6)', ADVANCE='NO') "TotCa",""
           DO J = 1, NumAqChems
              WRITE(FH,FMT='(A15)', ADVANCE='NO') AqPhaseChemicalNames(J)
           END DO
           DO J = 1, NumOrgChems		
              WRITE(FH,FMT='(A15)', ADVANCE='NO') OrgPhaseChemicalNames(J)
           END DO
           DO J = 1, NumAqOrgChems		
              WRITE(FH,FMT='(A15)', ADVANCE='NO') AqOrgPhaseChemicalNames(J)
           END DO
           WRITE(FH, FMT='(A1)') ""
        ELSE
           PotassiumIndex = FindChem("K+", 1)
           KNO3i = FindChem("KNO3", 1)
           KCli = FindChem("KCl", 1)
           KHSO4i = FindChem("KHSO4", 1)
           K2SO4i = FindChem("K2SO4", 1)
		
           Nai = FindChem("Na+", 1)
           NaNO3i = FindChem("NaNO3", 1)
           NaCli = FindChem("NaCl", 1)
           NaHSO4i = FindChem("NaHSO4", 1)
           Na2SO4i = FindChem("Na2SO4", 1)
		
           NH3i = FindChem("NH3", 1)
           NH4i = FindChem("NH4+", 1)
           NH4NO3i = FindChem("NH4NO3", 1)
           NH4Cli = FindChem("NH4Cl", 1)
           NH4HSO4i = FindChem("NH4HSO4", 1)
           NH42SO4i = FindChem("(NH4)2SO4", 1)

	   Mgi = FindChem("Mg++", 1)
	   MgNO32i = FindChem("Mg(NO3)2", 1)
	   MgCl2i = FindChem("MgCl2", 1)
	   MgHSO42i = FindChem("Mg(HSO4)2", 1)
	   MgSO4i = FindChem("MgSO4", 1)
           MgOH2i = FindChem("Mg(OH)2", 1)

	   Cai = FindChem("Ca++", 1)
	   CaNO32i = FindChem("Ca(NO3)2", 1)
	   CaCl2i = FindChem("CaCl2", 1)
	   CaHSO42i = FindChem("Ca(HSO4)2", 1)
	   CaSO4i = FindChem("CaSO4*2H2O", 1)
           CaOH2i = FindChem("Ca(OH)2", 1)
		
           Cli = FindChem("Cl-", 1)
           NO3i = FindChem("NO3-", 1)
           SO4i = FindChem("SO4--", 1)
           HSO4i = FindChem("HSO4-", 1)

           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') t			
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%numberofparticles
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%effectiveradius
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%IonicStr
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%WaterActivity
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%SolutionDensity
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%ParticleDensity
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(InParticle)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%Temperature
			
           !BC Mass Conc (ug/m3)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%Orgchems(1)*InParticle%numberofparticles &
                                                      *1.0e12*OrgMolecularMass(1)
					
           !POA Mass Conc (ug/m3)
           TotalPOA = 0.
           DO J = 2, 9
              TotalPOA = TotalPOA + InParticle%Orgchems(J)*InParticle%numberofparticles &
                                    *1.0e12*OrgMolecularMass(J)
              TotalPOA = TotalPOA + InParticle%AqOrgchems(J-1)*InParticle%numberofparticles &
                                    *1.0e12*AqOrgMolecularMass(J-1)
           END DO
           DO J = 20, NumOrgChems
              TotalPOA = TotalPOA + InParticle%Orgchems(J)*InParticle%numberofparticles &
                                    *1.0e12*OrgMolecularMass(J)
              TotalPOA = TotalPOA + InParticle%AqOrgchems(J-1)*InParticle%numberofparticles &
                                    *1.0e12*AqOrgMolecularMass(J-1)
           END DO
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

           !SOA Mass Conc (ug/m3)
           TotalSOA = 0.
           DO J = 10, 19
              TotalSOA = TotalSOA + InParticle%Orgchems(J)*InParticle%numberofparticles &
                                    *1.0e12*OrgMolecularMass(J)
              TotalSOA = TotalSOA + InParticle%AqOrgchems(J-1)*InParticle%numberofparticles &
                                    *1.0e12*AqOrgMolecularMass(J-1)
           END DO
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
           !Water Mass Conc. (ug/m3)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%Aqchems(1)*InParticle%numberofparticles &
                                                     *1.0e12*AqMolecularMass(1)
					
           !K+ Mass Conc (ug/m3)
           TotalK = 0.		
           DO J = 1, NumAqChems
              IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli &
                  .OR. J .EQ. KHSO4i .OR. J .EQ. KOHi) THEN
                 TotalK = TotalK + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. K2SO4i) THEN
                 TotalK = TotalK + 2*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalK = TotalK*AqMolecularMass(PotassiumIndex)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
           !Na+ Mass Conc (ug/m3)
           TotalNa = 0.		
           DO J = 1, NumAqChems
              IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i .OR. J .EQ. NaOHi) THEN
                 TotalNa = TotalNa + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. Na2SO4i) THEN
                 TotalNa = TotalNa + 2*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalNa = TotalNa*AqMolecularMass(Nai)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
           !NH4+ Mass Conc (ug/m3)
           TotalNH4 = 0.		
           DO J = 1, NumAqChems
              IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i &
                  .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
                 TotalNH4 = TotalNH4 + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. NH42SO4i) THEN
                 TotalNH4 = TotalNH4 + 2*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. LEVi) THEN
                 TotalNH4 = TotalNH4 + 3*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

           !Cl- Mass Conc (ug/m3)
           TotalCl = 0.		
           DO J = 1, NumAqChems
              IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
                 TotalCl = TotalCl + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. MgCl2i .OR. J .EQ. CaCl2i) THEN
                 TotalCl = TotalCl + 2*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalCl = TotalCl*AqMolecularMass(Cli)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

           !NO3- Mass Conc (ug/m3)
           TotalNO3 = 0.
           DO J = 1, NumAqChems
              IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
                 TotalNO3 = TotalNO3 + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. MgNO32i .OR. J .EQ. CaNO32i) THEN
                 TotalNO3 = TotalNO3 + 2*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

           !SO4- Mass Conc (ug/m3)
           TotalSO4 = 0.
           DO J = 1, NumAqChems
              IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
                   .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i &
                   .OR. J .EQ. CaSO4i .OR. J .EQ. MgSO4i) THEN
                 TotalSO4 = TotalSO4 + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              ELSE IF (J .EQ. LEVi .OR. J .EQ. MgHSO42i .OR. J .EQ. CaHSO42i) THEN
                 TotalSO4 = TotalSO4 + 2*InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

			
           !Mg++ Mass Conc (ug/m3)
           TotalMg = 0.		
           DO J = 1, NumAqChems
              IF (J .EQ. Mgi .OR. J .EQ. MgNO32i .OR. J .EQ. MgCl2i .OR. J .EQ. MgHSO42i .OR. J .EQ. MgOH2i &
                  .OR. J .EQ. MgSO4i) THEN
                 TotalMg = TotalMg + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalMg = TotalMg*AqMolecularMass(Mgi)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalMg

           !Ca++ Mass Conc (ug/m3)
           TotalCa = 0.		
           DO J = 1, NumAqChems
              IF (J .EQ. Cai .OR. J .EQ. CaNO32i .OR. J .EQ. CaCl2i .OR. J .EQ. CaHSO42i .OR. J .EQ. CaOH2i &
                  .OR. J .EQ. CaSO4i) THEN
                 TotalCa = TotalCa + InParticle%Aqchems(J)*InParticle%numberofparticles*1.0e12
              END IF
           END DO
           !Convert from umol/m3 to ug/m3
           TotalCa = TotalCa*AqMolecularMass(Cai)
           WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') TotalCa

           !mol/particle
           DO J = 1, NumAqChems
              WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%AqChems(J)
           END DO

           DO J = 1, NumOrgChems		
              WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%Orgchems(J)
           END DO

           DO J = 1, NumAqOrgChems		
              WRITE(FH,FMT='(7ES13.5E2)', ADVANCE='NO') InParticle%AqOrgchems(J)
           END DO

           WRITE(FH, FMT='(A1)') ""
	ENDIF
END SUBROUTINE AerosolOutput

END MODULE OutputRoutines

