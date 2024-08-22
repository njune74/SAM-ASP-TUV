!! ASP (c), 2004-2012, Matt Alvarado (malvarad@ser.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ParticleStructure.h
!! This file is the definition of the structure used for aerosol 
!! particles in ASP.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 11/01  2006   Matt Alvarado     Added Optical Properties                  !!
!! 05/16  2007   Matt Alvarado     Revised to match optical                  !!
!!				       properties to CRM6                    !!
!! 08/30  2010   Matt Alvarado     Changed dimensions of optical properties  !!
!!                                 to allow 60 bins between 290 nm and 410 nm!!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/16  2012   Matt Alvarado     Changed dimensions of optical properties !!
!!                                 to allow 451 bins from 250 nm to 700 nm.  !!
!! 08/17  2012   Matt Alvarado     Add extinction and scattering efficiency 
!!                                  and asymm param for a particle
!!                                  with the same effective size but 
!!                                  (a) all BC or (b) all shell.
!! 11/08  2012   Matt Alvarado     Add extinction and scattering efficiency 
!!                                  and asymm param for a particle using
!!                                  MaxWell -Garnett mixing rule
!! 11/08  2012   Matt Alvarado     Add extinction and scattering efficiency 
!!                                  and asymm param for a particle using
!!                                  volume average dielectric constant mixing
!!                                  rule.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file defines the following data types:	!!
!! 1. TYPE Particle				!!
!! 2. TYPE ParticleArray                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! This is the primary structure to be used
TYPE Particle

	!! Probably we'll want to redo this to have a linked list for 
	!! aerosol that combine
	INTEGER :: ParticleID
	INTEGER :: ParticleDistribution  
        !! Corresponds to the initial section or -1 if from coagulation, 
	!!if mixed in will have fractional component as well.
	REAL*8  :: OriginPressureLevel

	!! When the system is running as a sectional model each of these
	!! particles will represent more than one particle
	LOGICAL :: Sectional
	REAL*8    :: NumberOfParticles
	REAL*8    :: Edges(2)

	!! velocity is a 3-space (x,y,z) expression 
        !! (in Model Length Scale and Seconds)
	REAL*8 :: velocity(3)

	LOGICAL :: Dry         !! True if the particle never picks up water 
                               !! (i.e., contact angle = 2pi or no electrolyte)

	!! These are physical parameters that each particle is tagged with
	REAL*8 :: Temperature
	
	REAL*8 :: EmbryoRadius    !! This is the radius of the solution embryo
	REAL*8 :: InsolubleRadius 
        !! This is the radius of the insoluble part of the aerosol
	REAL*8 :: EffectiveRadius 
        !! This is the radius of the sphere with this particles' volume
		
	!! Thermodynamic Properties
	REAL*8 :: IonicStr
	REAL*8 :: WaterActivity
	REAL*8 :: InorgSolnDensity   
        !! Only includes water and ions - for optical property calcs
	REAL*8 :: SolutionDensity    
        !! This is the solution density, including water, ions, 
        !! and aqueous organics
	REAL*8 :: InsolubleDensity   
        !! This is the average density of the combined solution 
        !! and solid salts, not including insoluble compounds
	REAL*8 :: ParticleDensity    !! This includes insoluble compounds
	REAL*8 :: SurfaceTension

	!! This gives the mean activity coefficients for each electrolyte combo
	REAL*8, POINTER :: GammaMixed(:)   
        !! Mixed solution mean activity coefficients
        !! These are indexed by AqRxn number, *not* electrolyte index!

	!! These are the composition vectors
	!! 
	!! Aqueous Composition is tracked in moles of the particular chemical
	!! Functions are provided to change into molality or other variables,
	!! but since the water content varies quickly recalculation of these 
	!! vectors would become a problem for conservation, etc.
	!!
	!! AqChems lists first the Chemicals 
        !!  (Electrolytes and Non-Electrolytes),
	!!  then the Cations and then the Anions.  
        !! AqChems(1) is hardcoded to be water.
        !!
	!! OrgChems lists the chemicals in the organic phase, 
        !! with BC hard coded as the first compound
	REAL*8, POINTER :: AqChems(:)
	REAL*8, POINTER :: OrgChems(:)
	REAL*8, POINTER :: AqOrgChems(:)

	REAL*8, POINTER :: HydrophobicActivityCoeffs(:) 
        !These are indexed by org compound number
	REAL*8, POINTER :: HydrophilicActivityCoeffs(:) 
        !These are indexed by aq org compound number

	REAL*8, DIMENSION(451) :: ShellRealRefracInd 
        !Solar (550 nm) and IR (10 um)
	REAL*8, DIMENSION(451) :: ShellImagRefracInd 
        !Solar (550 nm) and IR (10 um)
	REAL*8, DIMENSION(451) :: MGRealRefracInd 
        !Solar (550 nm) and IR (10 um)
	REAL*8, DIMENSION(451) :: MGImagRefracInd 
        !Solar (550 nm) and IR (10 um)
	REAL*8, DIMENSION(451) :: VARealRefracInd 
        !Solar (550 nm) and IR (10 um)
	REAL*8, DIMENSION(451) :: VAImagRefracInd 
        !Solar (550 nm) and IR (10 um)
	REAL*8				  :: AbsCoreRad
!In the CRM6 radiation scheme, six  and  12 bands are selected for solar 
! and thermal IR regions, respectively. The spectral division is below: 
! Band Range              Refractive Index  Proxy Wavelength
! 0.2 - 0.7 um              Solar (550 nm)      550 nm
! 0.7 - 1.3 um              Solar (550 nm)      1.0 um
! 1.3 - 1.9 um	            Solar (550 nm)      1.6 um
! 1.9 - 2.5 um              Solar (550 nm)      2.2 um
! 2.5 -3.5 um               Solar (550 nm)      3.0 um
! 3.5 - 4.0 um              Solar (550 nm)      3.75 um  
! 2200 - 1900 cm**-1,           IR (10 um)	4.79 um
! 1900 - 1700 cm**-1,		IR (10 um)	5.45 um 
! 1700 -1400  cm**-1,		IR (10 um)	6.51 um
!  1400 - 1250 cm**-1,		IR (10 um)	7.57 um
!  1250 - 1100 cm**-1,		IR (10 um)	8.55 um
! 1100 - 980 cm**-1,		IR (10 um)	9.65 um
! 980 - 800 cm**-1,		IR (10 um)	11.35 um
!  800 - 670 cm**-1,		IR (10 um)	13.71 um
!  670 - 540 cm**-1,		IR (10 um)	16.72 um
! 540 - 400 cm**-1,		IR (10 um)      21.76 um
!  400 - 280 cm**-1,		IR (10 um) 	30.36 um
!  280 - 0 cm**-1,		IR (10 um)	71.43 um
	
!For photolysis, the first 60 are from 290 nm to 410 nm in 2 nm increments, 
! then 550 nm and 700 nm for comparison with aerosol optical property 
! measurements.

	REAL*8, DIMENSION(451) :: ExtEff !Extinction Efficiency
	REAL*8, DIMENSION(451) :: ScaEff !Scattering Efficiency
	REAL*8, DIMENSION(451) :: BackScaEff !Back-Scattering Efficiency
	REAL*8, DIMENSION(451) :: AssymParam !Assymetry Parameter

        !MJA 08-17-2012 Below for a hypothetical particle of
        !same effective size but all BC
	REAL*8, DIMENSION(451) :: ExtEffBC !Extinction Efficiency
	REAL*8, DIMENSION(451) :: ScaEffBC !Scattering Efficiency
	REAL*8, DIMENSION(451) :: BackScaEffBC !Back-Scattering Efficiency
	REAL*8, DIMENSION(451) :: AssymParamBC !Assymetry Parameter

        !MJA 08-17-2012 Below for a hypothetical particle of
        !same effective size and same ref. ind. as the shell
	REAL*8, DIMENSION(451) :: ExtEffShell !Extinction Efficiency
	REAL*8, DIMENSION(451) :: ScaEffShell !Scattering Efficiency
	REAL*8, DIMENSION(451) :: BackScaEffShell !Back-Scattering Efficiency
	REAL*8, DIMENSION(451) :: AssymParamShell !Assymetry Parameter

        !MJA 11-08-2012 Below for Maxwell-Garnett Mixing (MG)
	REAL*8, DIMENSION(451) :: ExtEffMG !Extinction Efficiency
	REAL*8, DIMENSION(451) :: ScaEffMG !Scattering Efficiency
	REAL*8, DIMENSION(451) :: BackScaEffMG !Back-Scattering Efficiency
	REAL*8, DIMENSION(451) :: AssymParamMG !Assymetry Parameter

        !MJA 11-08-2012 Below for Volume average dielectric constant mixing (VA)
	REAL*8, DIMENSION(451) :: ExtEffVA !Extinction Efficiency
	REAL*8, DIMENSION(451) :: ScaEffVA !Scattering Efficiency
	REAL*8, DIMENSION(451) :: BackScaEffVA !Back-Scattering Efficiency
	REAL*8, DIMENSION(451) :: AssymParamVA !Assymetry Parameter
	
	!! These area stored as linked lists
	TYPE (Particle), POINTER :: Next

END TYPE Particle

!! The Particles Are Stored in an Array of these
TYPE ParticleArray
	TYPE (Particle), POINTER :: First
END TYPE ParticleArray
