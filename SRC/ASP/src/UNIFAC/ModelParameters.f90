!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ModelParameters.f90
!! This file holds parameters and constants germaine to the aerosol definitions
!! and internal processes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 07/19  2007   Matt Alvarado     Changed MaxIonicStrKM from 100 to 50      !!
!! 08/01  2007   Matt Alvarado     Changed MinimumWater from 1.0e-21 to 
!!                                   5.0e-21
!! 08/29  2007   Matt Alvarado     Fixed value of epsilon to match 
!!                                   Jacobson value
!! 02/15  2012   Matt Alvarado     Removed Eulerian coordinates
!! 03/01  2012   Matt Alvarado     Made it so AverageAerosolDensity can be   !!
!!                                  set and reset automatically.             !!
!! 11/06  2012   Matt Alvarado     Updated main input deck to give model more!!
!!                                 flexibility without recompiling.          !!
!! 12/28  2012   Matt Alvarado     Fixed units conversion bug in SetHetero   !!
!! 01/29  2013   Matt Alvarado     Changed AqThermoNumbAllowedIterations     !!
!!                                    from 50000 to 100.                     !!
!! 07/03  2013   Matt Alvarado     Added subroutine SetPhotoFlag and         !!
!!                                    SetPhotoEndTime                        !!
!! 07/03  2013   Matt Alvarado     Updated SetEulerianBox to include Fire CO !!
!!                                   emission rate.                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES							    !!
!! 1. NONE							    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	
!! 1. SUBROUTINE SetDomainSize ()
!! 2. SUBROUTINE SetBackgroundRH (BRH)
!! 3. SUBROUTINE SetDissolutionEquilibriumOrFlux(Input)
!! 4. SUBROUTINE SetSectionalParameters (inLagrSectSizeCutoff,inLowestBinCutoff,inHowManyBins, inStochasticOrNot, inLagrangianOrNot)
!! 5. SUBROUTINE SetArraysForInsolubleAerosolCores (NumbDists)
!! 6. SUBROUTINE DeallocateDataForInsolubleAerosolCores ()
!! 7. SUBROUTINE SpecifyNewInsolType(Radius, Density, CosContactAngle, DistNumb
!! 8. SUBROUTINE SetCondensationFlag(InFlag)
!! 9. SUBROUTINE SetDissolutionFlag(InFlag)
!!10. SUBROUTINE SetHydrophobicOrgDissolutionFlag(InFlag)
!!11. SUBROUTINE SetHydrophilicOrgDissolutionFlag(InFlag)
!!12. SUBROUTINE SetCoagulationFlag(InFlag)
!!13. SUBROUTINE SetGasChemistryFlag(InFlag)
!!14. SUBROUTINE SetAqChemistryFlag(InFlag)
!!15. SUBROUTINE SetWaterEquilibriumIndex(WEI)
!!16. SUBROUTINE SetTranscriptFH(FH)
!!17. SUBROUTINE SetInputDeckSubDir(iInputDeckSubDir)
!!18. SUBROUTINE SetOutputDeckSubDir(iOutputDeckSubDir)
!!19. SUBROUTINE SetHowManyVariousAqChems(iHowManyAqChems, iHowManyAqCations, iHowManyAqAnions)
!!20. SUBROUTINE SetProtonIndex(iProtonIndex)
!!21. SUBROUTINE SetHydroxyIndex(iHydroxyIndex)
!!22. SUBROUTINE SetMonodisperse(iMonodisperse)
!!23. SUBROUTINE SetConstantKernel
!!24. SUBROUTINE SetMetastable(iMetastable)
!!25. SUBROUTINE SetAmmoniaFlux(iAmmoniaFlux)
!!26. SUBROUTINE SetIgnoreDiffuseCorrect(iDiffuseCorrect)
!!27. SUBROUTINE SetInitialDens(iDens)
!!28. SUBROUTINE SetLagrangianParcelModel(iMod)
!!29. SUBROUTINE SetEulerianBoxModel(iMod)
!!30. SUBROUTINE SetSmogChamberModel(iMod)
!!31. SUBROUTINE SetOpticalPropertiesModel(iMod)
!!32. SUBROUTINE SetThermoTest(iMod)
!!33. SUBROUTINE SetCoagTest(iMod)
!!34. SUBROUTINE SetCondTest(iMod)
!!35. SUBROUTINE SetSubroutineModel(iMod)
!!36. SUBROUTINE SetStartTime(UTC)
!!37. SUBROUTINE SetRunLength(time)
!!38. SUBROUTINE SetChemStep(step)
!!39. SUBROUTINE SetMixStep(step)
!!40. SUBROUTINE SetCondStep(step)
!!41. SUBROUTINE SetCoagStep(step)
!!42. SUBROUTINE SetInversionHeight(height)
!!43. SUBROUTINE SetHetero(num, rad)
!!44. SUBROUTINE SetLagrangianMix(mix)
!!45. SUBROUTINE SetEulerianBox(wind,length,width, fireCO)
!!46. SUBROUTINE SetIncludeAerosols(flag)
!!47  SUBROUTINE SetPhotoFlag(flag)
!!48. SUBROUTINE SetPhotoEndTime(UTC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This one holds all of the physical constants that are used in the code
MODULE ModelParameters

	!! These parameters set the scaling of internal representations of
	!! Length, Mass, and the like
	!!
	!! BUT DON'T CHANGE THEM!!! DEBUGGING THESE SUCKS.  I KNOW THAT THE CHEMSCALES 
	!! WORK BUT I HAVEN'T YET (02/13/04) DEBUGGED THE OTHER BUGGERS.
	!!
	real*8, parameter :: LengthScale = 1.			 ! The lengthscale defines the scale at which we make internal
	real*8, parameter :: MassScale   = 1.			 ! Convert from grams to g/MassScale
	real*8, parameter :: ChemScale   = 1.			 ! Convert from ppb to ppb/MassScale
	real*8, parameter :: AqChemScale = ChemScale	 ! Convert to scaled moles DON"T CHANGE -- I SCREWED UP BY NOT 
													 ! SCALING MOLALITY AND BY CONFUSING CHEMSCALE AND AQCHEMSCALE
	real*8, parameter :: Pressurescale = MassScale / LengthScale
!! Do not provide this service for pressure.  It is not written into the code and it will only screw things
!! up to convert everything to non-milibars.  mb = 1., absolutely
!!!	real*8, parameter :: PressureScale = 1.	 ! Convert from milibars to mb / PressureScale

	!! These parameters may be used to specify which length scale is specified, and 
	!! converts those units to the working scale, defined in centimeters multiplied by the 
	!! LengthScale specified above.
	real*8, parameter :: kilograms	         = MassScale * 1.d3
	real*8, parameter :: grams		 = MassScale * 1.		! Base Mass is Grams

	real*8, parameter :: km			 = LengthScale * 1.d5
	real*8, parameter :: m			 = LengthScale * 100.
	real*8, parameter :: cm			 = LengthScale * 1.		! Base Length is Centimeters
	real*8, parameter :: mm			 = LengthScale * 0.1
	real*8, parameter :: micron		 = LengthScale * 1.d-4
	real*8, parameter :: nm			 = LengthScale * 1.d-7
	real*8, parameter :: Angstrom	         = LengthScale * 1.d-8
	real*8, parameter :: fermi		 = LengthScale * 1.d-13

	real*8, parameter :: mile		 = LengthScale * 160938.
	real*8, parameter :: yard		 = LengthScale * 91.44
	real*8, parameter :: foot		 = LengthScale * 30.48
	real*8, parameter :: inch		 = LengthScale * 2.540

	!! Convert Gas Phase Chemistry variables to ppb
	real*8, parameter :: ppb		 = ChemScale * 1.
	real*8, parameter :: ppt		 = ChemScale * 0.001
	real*8, parameter :: ppm		 = ChemScale * 1000.
	real*8, parameter :: ppbscale    = ChemScale * 1.d-9

	!! Convert Aqueous Phase Chemistry variables to Moles
	real*8, parameter :: moles		 = AqChemScale * 1.

	real*8, parameter :: seconds = 1.
	real*8, parameter :: minutes = seconds * 60.
	real*8, parameter :: hours   = minutes * 60.
	real*8, parameter :: days    = hours * 24.

	!! Don't scale these with a pressure scale.  It's too much to convert back and forth, Since it's often used in 
	!! the conversion routines of other quantities and would add extra overhead.  The internal calculation is based
	!! on milibars and nothing else.
	real*8, parameter :: mbar		 = 1.*Pressurescale		 ! Base Pressure is Milibars  
	real*8, parameter :: mbars		 = 1.*Pressurescale		 ! synonymns
	real*8, parameter :: Atmospheres = 1013.25*Pressurescale !  (all factors from Iribarne and Godson)
	real*8, parameter :: Pa			 = 0.01*Pressurescale    ! kg / m / s^2
	real*8, parameter :: DynCm2		 = Pressurescale	     !	grams / cm			! Dynes per CM^2
	real*8, parameter :: DynCm		 = grams 
	real*8, parameter :: kPa		 = 10.*PressureScale

	!! Don't scale energies either.  Base energy is JOULES
	real*8, parameter :: joules		  = 1.*kilograms*m*m	 ! Base energy is JOULES (J = Kg m2 s-2)
	real*8, parameter :: ergs		  = 1.d-7*joules		 ! One Joule = 1.d7 ergs


	!!!!!!!!!!!!!!!!!!!!!!!!
	!! Physical constants !!
	!!!!!!!!!!!!!!!!!!!!!!!!
	real*8, parameter  :: pi=3.1415926535898				  ! Pi
	real*8, parameter  :: EulerMass=.577215664901532860606512 ! Euler-Mascheroni Constant
	real*8, parameter  :: kB=1.3807d-16 * grams * cm * cm	  ! Boltzmann's Constant (g cm^2 / s^2 K molec)
	real*8, parameter  :: AirMolecMass=4.8096d-26 * kilograms ! average mass of an air molecule (g)
	real*8, parameter  :: WaterMolecMass=18.016*grams/moles
	real*8, parameter  :: eps=0.622							  ! thermodynamic epsilon (see Jacobson, eqn 2.28 (2nd ed. )
	real*8, parameter  :: Rdry=2870.4 * cm * cm * cm * mbar & ! Gas Constant for Dry air (cm^3 * mb / g K)
							/ grams
	real*8, parameter  :: g=981. * cm 						  ! gravity g in cm / s^2
	real*8, parameter  :: Avogadro=6.0221367d23				  ! Avogardo's Number
	real*8, parameter  :: Cpd=1004.67 * joules / kilograms	  ! specific heat of dry air at constant P (J / kg K) at 298 K
	real*8, parameter  :: CpV = 1865.1 * joules / kilograms	  ! specific heat of water vapor at constant P (J / kg K) at 298 K
	real*8, parameter  :: Rstar=8.3145 * kilograms * m * m / moles ! Universal Gas Constant (Kg m2 / s2 mole K)
	real*8, parameter  :: RstarMB=8.3145*10000. * cm * cm * cm * mbar / moles ! Universal Gas Constant (cm3 mb / mole / K)
	real*8, parameter  :: ThmVelFac = 1.87508d-08 !(8*kB / Pi )**0.5 ! The constant leading factor of the thermal velocity of air
															  ! (retreival function in "Temperature.f90" file in Eulerian
															  ! Environment folder).
	REAL*8, PARAMETER  :: ThermalAccommodationCoeff = 0.96


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Should we do the different processes? !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	LOGICAL :: DoDissolutionFlag
	LOGICAL :: DoHydrophobicOrgDissolutionFlag
	LOGICAL :: DoHydrophilicOrgDissolutionFlag
	LOGICAL :: DoCondensationFlag
	LOGICAL :: DoCoagulationFlag
	LOGICAL :: DoGasChemistryFlag
	LOGICAL :: DoAqChemistryFlag
	LOGICAL :: ConstantKernel
	LOGICAL :: Metastable

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! These are parameters really stored other places,		!!
	!! but are here simply to provide a way around circular !!
	!! dependencies											!!
	INTEGER :: HowManyAqChems, HowManyAqCations, HowManyAqAnions
	INTEGER :: ProtonIndex, HydroxyIndex
	INTEGER :: HowManyOrgChems
	INTEGER :: HowManyAqOrgChems


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Locate the Output Decks !!
	CHARACTER (len=512) :: InputDeckSubDir
	CHARACTER (len=512) :: OutputDeckSubDir

	!! This is set to -1 if no transcript is to be kept
	INTEGER :: TranscriptFH

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! THERMODYNAMIC Variables !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8, PARAMETER  :: MaxIonicStrKM = 50.0   !Matt put in new value 8/25/05
					             !! The maximum ionic strength allowed by the Kusik-Meissner routine 
						     !! (larger will just be done using this)
	REAL*8, PARAMETER  :: ThermoTref  = 298.15	        !! The reference temperature for the Kusik-Meissner Calculations
	REAL*8, PARAMETER  :: EqCoeffTref = 298			!! The reference temperature for the Equilibrium 
								!! Coefficient Calculations
	REAL*8, PARAMETER  :: AqThermoEquilibriumError = 0.0005 !! This is an error value for how close a dissociation reaction must be 
								!! to convergence to return.  It's in a Kcalc / Ksupplied value.	 If 
								!! it's too small, then no aq equilibrium reactions will be demoted and
								!! conflicts may occur.
	INTEGER, PARAMETER :: AqThermoNumbAllowedIterations = 100
	INTEGER, PARAMETER :: AqThermoNumbAllowedIterationsForHypLoop = 20
	REAL*8, PARAMETER  :: AerosolWaterEquilibriumRHThreshold = 0.98 !! If RH is above this value, simply don't equilibrate the 
									!! water.  Trust what is there
	LOGICAL            :: ThermoBulkMode                  !! Resolve the particle distributions or not?
	
	INTEGER :: WaterEquilibriumIndex
	LOGICAL :: CalculateWaterGamma			!! If cannot figure out how to calculate the water act coeff, just use one!

	!! If WaterContentPrecision number is too large, then equilibrium water content may oscillate grandly, 
	!! preventing AqThermoEquilibriumError from being satisfied 

	REAL*8 :: WaterContentPrecision = 0.0001   !! The Relative Humidity and Water Activity must be within this amount, or LSODES 
											   !! uses this as the error tolerance...
	LOGICAL, PARAMETER :: OpenSystemForCondensation = .FALSE.
	LOGICAL			   :: DissolutionEquilibriumOrFlux		! False is equilibrium, True is flux
	REAL*8            :: AverageAerosolDensity  != 1.5 * grams / cm / cm / cm
	REAL*8, PARAMETER :: SmallestAerosolPossible = 0.0015 * micron   ! Any electrolyte sphere below this size doesn't take up water
																	! If this is too small, say O(0.001 microns), there will be 
																	! only O(100) molecules of water and the condensation routine 
																	! will make this go negative...
	REAL*8, PARAMETER :: MinimumWater  = 5.0e-21   ! If particle water content is this small (mol/particle), ignore aqueous equilibrium
	REAL*8, Parameter :: CosContactAngle = 0.      !Cosine of contact angle of water embryo with insoluble sphere														
	REAL*8, Parameter :: CoagKernel = 1.0e-9   !Kernel in cm3/s for constant kernel cases
	REAL*8 :: InversionHeight ! = 1.0*km !Height of inversion layer (stored as cm)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The extent of the physical DOMAIN is the same for all grids	!!
	!! in the model and is taken to be rectangular.  This set of	!!
	!! parameters define that domain size.				!!
	!!								!!
	!! The origin is at 0,0,0, and the domain stretches this far in	!!
	!! each direction.						!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 :: DomainX, DomainY, DomainZ


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The SECTIONAL REPRESENTATION has	    !!
	!! several parameters, allocated here	    !!
	!! and defined in SetSectionalParameters()  !!
	REAL*8  :: LagrSectSizeCutoff,	&
	           LowestBinCutoff
	INTEGER :: HowManyBins
	LOGICAL :: StochasticOrNot
	REAL*8, ALLOCATABLE  :: BinEdges(:,:)
	REAL*8  :: BinRadiusRatio
	LOGICAL :: LagrangianOrNot !! Are we tracking lagrangian aerosol or does everything go into a section?
	LOGICAL :: MONODISPERSE !! Are the aerosol monodisperse?
	LOGICAL :: AmmoniaFlux !!Is ammonia dissolution a flux (1) or eq(0) process?
	LOGICAL, PARAMETER :: UseDonnanSurfaceTension = .TRUE. !If true, use Donnan's formula, else use Jacobson
	LOGICAL :: IgnoreCorrectionstoDiffusivity
	LOGICAL, PARAMETER :: CombinedInorgDiss = .TRUE.
	LOGICAL, PARAMETER :: ZSROrgWater = .TRUE.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The INSOLUBLE CORES have	a number of characteristics, !!
	!! allocated here.                                       !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER              :: NumbInsolTypes
	REAL*8, ALLOCATABLE  :: InsolRadii(:)
	REAL*8, ALLOCATABLE  :: InsolDensity(:)
	REAL*8, ALLOCATABLE  :: InsolContactAngle(:)
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!Added to allow user to input Background RH in main input deck!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 :: Back_RH

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Set Model Run Type using logical variables (MJA, 11/06/2012) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        LOGICAL :: LagrangianParcelModel
        LOGICAL :: EulerianBoxModel
        LOGICAL :: SmogChamberModel      
        LOGICAL :: OpticalPropertiesModel
        LOGICAL :: ThermoTest
        LOGICAL :: CoagTest
        LOGICAL :: CondTest
        LOGICAL :: SubroutineModel 
        LOGICAL :: IncludeAerosols
        LOGICAL :: PhotoFlag

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Set start time, run length, and time steps (MJA, 11/06/2012) !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
        REAL*8  :: UTCStartTime, RunLength, ChemStep, MixStep, CondStep, CoagStep, PhotoEndTime

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Set Lagrangian Parecl Parameters (MJA, 11/06/2012)           !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL*8  :: b_mix, hetero_num, hetero_rad  

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!Set Eulerian Box Parameters (MJA, 11/06/2012)                !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        REAL*8  :: windspeed, boxlength, boxwidth, fire_co_emis    

	!!!!!!! DEFINE SUBROUTINES !!!!!!!
	CONTAINS

! 23 September 2016 CMB: Add these functions so that InitializeASP() can work properly
    real*8 function getHowManyBins()
        getHowManyBins = HowManyBins
        return
    end function
    real*8 function getHowManyAqChems()
        getHowManyAqChems = HowManyAqChems
        return
    end function
    real*8 function getHowManyAqCations()
        getHowManyAqCations = HowManyAqCations
        return
    end function     
    real*8 function getHowManyAqAnions()
        getHowManyAqAnions = HowManyAqAnions
        return
    end function
    real*8 function getHowManyOrgChems()
        getHowManyOrgChems = HowManyOrgChems
        return
    end function
    real*8 function getHowManyAqOrgChems()
        getHowManyAqOrgChems = HowManyAqOrgChems
        return
    end function
! end CMB 

	!! Always call SetGridPointVolumes() with this!  
        !! (Can't from here because of heirarchy)
        !! MJA, 02-15-2012 Should always be 1 cm x 1 cm x 1 cm
        !! so that concentrations are per cm^3
	SUBROUTINE SetDomainSize ()

		IMPLICIT NONE

		DomainX = 1*cm
		DomainY = 1*cm
		DomainZ = 1*cm

		RETURN
	END SUBROUTINE SetDomainSize

	SUBROUTINE SetBackgroundRH (BRH)

		IMPLICIT NONE
		REAL*8 :: BRH

		Back_RH = BRH

		RETURN
	END SUBROUTINE SetBackgroundRH

	SUBROUTINE SetDissolutionEquilibriumOrFlux(Input)

	  IMPLICIT NONE
	  INTEGER :: Input

	  IF (Input .EQ. 0) THEN
		DissolutionEquilibriumOrFlux = .FALSE. !Flag to equilibrate
	  ELSE IF (Input .EQ. 1) THEN
		DissolutionEquilibriumOrFlux = .TRUE. !Flag to do growth
	  END IF

	  RETURN
	END SUBROUTINE SetDissolutionEquilibriumOrFlux

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set everything to do with the Representation of Sections       !!
        !! in the Model !!
	SUBROUTINE SetSectionalParameters (inLagrSectSizeCutoff,inLowestBinCutoff,inHowManyBins, inLagrangianOrNot)

		IMPLICIT NONE
		REAL*8  :: inLagrSectSizeCutoff,inLowestBinCutoff
		INTEGER :: inHowManyBins
		LOGICAL :: inLagrangianOrNot

		INTEGER :: I

		LagrSectSizeCutoff     = inLagrSectSizeCutoff
		LowestBinCutoff        = inLowestBinCutoff
		HowManyBins            = inHowManyBins

		LagrangianOrNot = inLagrangianOrNot

		IF (HowManyBins .NE. 0) THEN

			!! Allocate the array for bin edges
			ALLOCATE (BinEdges(HowManyBins,2))


			IF (LagrangianOrNot .AND. HowManyBins .EQ. 2) THEN
				!! Calculate the Bin Edges
				BinEdges(1,1) = 0.
				BinEdges(1,2) = LowestBinCutoff
				BinEdges(2,1) = LowestBinCutoff
				BinEdges(2,2) = LagrSectSizeCutoff
				RETURN
			END IF

			IF (.NOT.LagrangianOrNot .AND. HowManyBins .EQ. 1) THEN
				BinEdges(1,1) = 0.
				BinEdges(1,2) = 1.*km  !! this is an artificially large number equal enough to "infinity"
			END IF

			!! Calculate the Bin Edges
			BinEdges(1,1) = 0.
			BinEdges(1,2) = LowestBinCutoff

			IF (.NOT.LagrangianOrNot) THEN
				BinEdges(HowManyBins,1) = LagrSectSizeCutoff
				BinEdges(HowManyBins,2) = 1.*km
			END IF

			IF (HowManyBins .GT. 1 .AND. LagrangianOrNot) THEN

				!! Calculate the radius ratio
				BinRadiusRatio = (LagrSectSizeCutoff/LowestBinCutoff)**(1./(HowManyBins-1.))

				DO I = 2, HowManyBins
					BinEdges(I,1) = BinEdges(I-1,2)
					BinEdges(I,2) = BinEdges(I,1) * BinRadiusRatio
				END DO

				BinEdges(HowManyBins,2) = LagrSectSizeCutoff
			END IF

			IF (HowManyBins .GT. 2 .AND. .NOT.LagrangianOrNot) THEN

				!! Calculate the radius ratio
				BinRadiusRatio = (LagrSectSizeCutoff/LowestBinCutoff)**(1./(HowManyBins-2.))

				DO I = 2, HowManyBins-1
					BinEdges(I,1) = BinEdges(I-1,2)
					BinEdges(I,2) = BinEdges(I,1) * BinRadiusRatio
				END DO
				BinEdges(HowManyBins-1,2) = LagrSectSizeCutoff

			END IF

		END IF

		RETURN
	END SUBROUTINE SetSectionalParameters

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set flags that control performance    !!
	!! of the various microphysical routines !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetCondensationFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoCondensationFlag = InFlag
	END SUBROUTINE SetCondensationFlag

	SUBROUTINE SetDissolutionFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoDissolutionFlag = InFlag
	END SUBROUTINE SetDissolutionFlag

	SUBROUTINE SetHydrophobicOrgDissolutionFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoHydrophobicOrgDissolutionFlag = InFlag
	END SUBROUTINE SetHydrophobicOrgDissolutionFlag

	SUBROUTINE SetHydrophilicOrgDissolutionFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoHydrophilicOrgDissolutionFlag = InFlag
	END SUBROUTINE SetHydrophilicOrgDissolutionFlag
	
	SUBROUTINE SetCoagulationFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoCoagulationFlag = InFlag
	END SUBROUTINE SetCoagulationFlag

	SUBROUTINE SetGasChemistryFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoGasChemistryFlag = InFlag
	END SUBROUTINE SetGasChemistryFlag

	SUBROUTINE SetAqChemistryFlag(InFlag)
		IMPLICIT NONE
		LOGICAL :: InFlag
		DoAqChemistryFlag = InFlag
	END SUBROUTINE SetAqChemistryFlag





	SUBROUTINE SetWaterEquilibriumIndex(WEI)

		IMPLICIT NONE
		INTEGER :: WEI

		WaterEquilibriumIndex = WEI

	END SUBROUTINE SetWaterEquilibriumIndex


	SUBROUTINE SetTranscriptFH(FH)

		IMPLICIT NONE
		INTEGER :: FH

		TranscriptFH = FH

	END SUBROUTINE SetTranscriptFH

	SUBROUTINE SetInputDeckSubDir(iInputDeckSubDir)
		
		IMPLICIT NONE

		CHARACTER*(*) :: iInputDeckSubDir

		InputDeckSubDir = TRIM(iInputDeckSubDir)

	END SUBROUTINE SetInputDeckSubDir

	SUBROUTINE SetOutputDeckSubDir(iOutputDeckSubDir)
		
		IMPLICIT NONE

		CHARACTER*(*) :: iOutputDeckSubDir

		OutputDeckSubDir = TRIM(iOutputDeckSubDir)
	END SUBROUTINE SetOutputDeckSubDir
	
	SUBROUTINE SetHowManyVariousAqChems(iHowManyAqChems, iHowManyAqCations, iHowManyAqAnions)
		
		IMPLICIT NONE

		INTEGER :: iHowManyAqChems, iHowManyAqCations, iHowManyAqAnions

		HowManyAqChems	  = iHowManyAqChems
		HowManyAqCations  = iHowManyAqCations
		HowManyAqAnions   = iHowManyAqAnions

	END SUBROUTINE SetHowManyVariousAqChems
	
	SUBROUTINE SetHowManyOrgChems(iHowManyOrgChems)
		
		IMPLICIT NONE

		INTEGER :: iHowManyOrgChems

		HowManyOrgChems	  = iHowManyOrgChems

	END SUBROUTINE SetHowManyOrgChems


	SUBROUTINE SetHowManyAqOrgChems(iHowManyAqOrgChems)
		
		IMPLICIT NONE

		INTEGER :: iHowManyAqOrgChems

		HowManyAqOrgChems	  = iHowManyAqOrgChems

	END SUBROUTINE SetHowManyAqOrgChems

	SUBROUTINE SetProtonIndex(iProtonIndex)
		
		IMPLICIT NONE

		INTEGER :: iProtonIndex

		ProtonIndex = iProtonIndex

	END SUBROUTINE SetProtonIndex

	SUBROUTINE SetHydroxyIndex(iHydroxyIndex)
		
		IMPLICIT NONE

		INTEGER :: iHydroxyIndex

		HydroxyIndex = iHydroxyIndex

	END SUBROUTINE SetHydroxyIndex

	SUBROUTINE SetMonodisperse(iMonodisperse)
		
		IMPLICIT NONE

		LOGICAL :: iMonodisperse

		MONODISPERSE = iMonodisperse

	END SUBROUTINE SetMonodisperse

	SUBROUTINE SetConstantKernel(iConstantKernel)
		
		IMPLICIT NONE

		LOGICAL :: iConstantKernel

		ConstantKernel = iConstantKernel

	END SUBROUTINE SetConstantKernel

	SUBROUTINE SetMetastable(iMetastable)
		
		IMPLICIT NONE

		LOGICAL :: iMetastable

		Metastable = iMetastable

	END SUBROUTINE SetMetastable

	SUBROUTINE SetAmmoniaFlux(iAmmoniaFlux)
		
		IMPLICIT NONE

		LOGICAL :: iAmmoniaFlux

		AmmoniaFlux = iAmmoniaFlux

	END SUBROUTINE SetAmmoniaFlux

	SUBROUTINE SetIgnoreDiffuseCorrect(iDiffuseCorrect)
		
		IMPLICIT NONE

		LOGICAL :: iDiffuseCorrect

		IgnoreCorrectionstoDiffusivity = iDiffuseCorrect

	END SUBROUTINE SetIgnoreDiffuseCorrect

	SUBROUTINE SetInitialDens(iDens)
		
		IMPLICIT NONE

		REAL*8 :: iDens

		AverageAerosolDensity = iDens * grams / cm / cm / cm

	END SUBROUTINE SetInitialDens

	SUBROUTINE SetLagrangianParcelModel(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                    LagrangianParcelModel = .TRUE.
                ELSE 
                    LagrangianParcelModel = .FALSE.
                ENDIF

	END SUBROUTINE SetLagrangianParcelModel

	SUBROUTINE SetEulerianBoxModel(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                    EulerianBoxModel = .TRUE.
                ELSE 
                    EulerianBoxModel = .FALSE.
                ENDIF

	END SUBROUTINE SetEulerianBoxModel

	SUBROUTINE SetSmogChamberModel(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                    SmogChamberModel = .TRUE.
                ELSE 
                    SmogChamberModel = .FALSE.
                ENDIF

	END SUBROUTINE SetSmogChamberModel

	SUBROUTINE SetOpticalPropertiesModel(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                   OpticalPropertiesModel = .TRUE.
                ELSE 
                   OpticalPropertiesModel = .FALSE.
                ENDIF

	END SUBROUTINE SetOpticalPropertiesModel

	SUBROUTINE SetThermoTest(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                   ThermoTest = .TRUE.
                ELSE 
                   ThermoTest = .FALSE.
                ENDIF

	END SUBROUTINE SetThermoTest

	SUBROUTINE SetCoagTest(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                    CoagTest = .TRUE.
                ELSE 
                    CoagTest = .FALSE.
                ENDIF

	END SUBROUTINE SetCoagTest

	SUBROUTINE SetCondTest(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                    CondTest = .TRUE.
                ELSE 
                    CondTest = .FALSE.
                ENDIF

	END SUBROUTINE SetCondTest

	SUBROUTINE SetSubroutineModel(iMod)
		
		IMPLICIT NONE

		INTEGER :: iMod

		IF (iMod .EQ. 1) THEN 
                      SubroutineModel = .TRUE.
                ELSE 
                     SubroutineModel = .FALSE.
                ENDIF
	END SUBROUTINE SetSubroutineModel

	SUBROUTINE SetStartTime(UTC)
		
		IMPLICIT NONE

		REAL*8 :: UTC

		UTCStartTime = UTC
	END SUBROUTINE SetStartTime

	SUBROUTINE SetRunLength(time)
		
		IMPLICIT NONE

		REAL*8 :: time !!In min

		RunLength = time
	END SUBROUTINE SetRunLength

	SUBROUTINE SetChemStep(step)
		
		IMPLICIT NONE

		REAL*8 :: step !!in sec

		ChemStep = step
	END SUBROUTINE SetChemStep

	SUBROUTINE SetMixStep(step)
		
		IMPLICIT NONE

		REAL*8 :: step !!in sec

		MixStep = step
	END SUBROUTINE SetMixStep

	SUBROUTINE SetCondStep(step)
		
		IMPLICIT NONE

		REAL*8 :: step !!in sec

		CondStep = step
	END SUBROUTINE SetCondStep

	SUBROUTINE SetCoagStep(step)
		
		IMPLICIT NONE

		REAL*8 :: step !!in sec

		CoagStep = step
	END SUBROUTINE SetCoagStep

	SUBROUTINE SetInversionHeight(height)
		
		IMPLICIT NONE

		REAL*8 :: height !!in km

		InversionHeight = height
	END SUBROUTINE SetInversionHeight

	SUBROUTINE SetHetero(num, rad)
		
		IMPLICIT NONE

		REAL*8 :: num !!in cm-3
		REAL*8 :: rad !!in um
		hetero_num = num
                hetero_rad = rad/1.0e6 !Convert to m
	END SUBROUTINE SetHetero

	SUBROUTINE SetLagrangianMix(mix)
		
		IMPLICIT NONE

		REAL*8 :: mix !!in s

		b_mix = mix
	END SUBROUTINE SetLagrangianMix

	SUBROUTINE SetEulerianBox(wind,length,width, fireCO)
		
		IMPLICIT NONE

		REAL*8 :: wind !!in m/s
		REAL*8 :: length, width !!in km
                REAL*8 :: fireCO !!in kg/s

		windspeed = wind
		boxlength = length
                boxwidth = width  
                fire_CO_emis = fireCO              
	END SUBROUTINE SetEulerianBox

	SUBROUTINE SetIncludeAerosols(flag)
		
		IMPLICIT NONE

		INTEGER :: flag 
	        IF (flag .EQ. 1) THEN 
                    IncludeAerosols = .TRUE.
                ELSE 
                    IncludeAerosols = .FALSE.
                ENDIF
	END SUBROUTINE SetIncludeAerosols

	SUBROUTINE SetPhotoFlag(flag)
		
		IMPLICIT NONE

		INTEGER :: flag 
	        IF (flag .EQ. 1) THEN 
                    PhotoFlag = .TRUE.
                ELSE 
                    PhotoFlag = .FALSE.
                ENDIF
	END SUBROUTINE SetPhotoFlag

	SUBROUTINE SetPhotoEndTime(UTC)
		
		IMPLICIT NONE

		REAL*8 :: UTC

		PhotoEndTime = UTC
	END SUBROUTINE SetPhotoEndTime

end module ModelParameters

