!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ASP Model (c) 2004-2015 Matt Alvarado (malvarad@aer.com)		     !!
!!                                                                           !!
!!  Originally developed by Matt Alvarado for his thesis                     !!
!!  with Prof. R.G. Prinn (rprinn@mit.edu) at MIT with funding from          !!
!!  an NSF Graduate Fellowship, a Martin Family Fellowship in                !!
!!  Sustainability, an MIT Norman B. Leventhal Presidential Fellowship,      !!
!!  NSF Grant ATM-0120468, DOE grant DE-FG02-94ER61937, and the industrial   !!
!!  and foundation sponsors of the MIT Joint Program on the Science and      !!
!!  Policy of Global Change.(2004-2008)		                             !!
!!									     !!
!!  Recent updates to ASP (2011-present) were funded through NSF Grant       !!
!!  AGS-1144165 and NASA Grants NNX11AN72G and NNX14AP45G.                   !!
!!                                                                           !!
!!  ASP is an updated, expanded version of the				     !!
!!  Mixed Eulerian-Lagrangian Aerosol Model (MELAM) developed by 	     !!
!!  Donnan Steele  for his thesis project		                     !!
!!  with Prof. R.G.Prinn at MIT. (2000-2004)				     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						             !!
!!								             !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 08/25  2006   Matt Alvarado     Changed GOTO loops to DO WHILE loops      !!
!! 09/06  2006   Matt Alvarado     Added particle deposition to size dist.   !!
!!                                   tests                                   !!
!! 09/07  2006   Matt Alvarado     Added gas deposition to size dist. test   !!
!! 10/03  2007   Matt Alvarado     Minor Changes                             !!
!! 08/27  2010   Matt Alvarado     Minor Changes                             !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/28  2012   Matt Alvarado     Added section to run optical property     !!
!!                                  tests.                                   !!
!! 08/01  2012   Tara Soni         Minor Changes                             !!
!! 11/06  2012   Matt Alvarado     Updated main input deck to give model more!!
!!                                 flexibility without recompiling.          !!
!! 11/08  2012   Matt Alvarado     Updated call to AerosolOptProps to include!!
!!                                  mixing flag                              !!
!! 12/28  2012   Matt Alvarado     Set Lagrangian Parcel to always use 1 s   !!
!!                                  time step in first 10 minutes of run     !!
!! 07/08  2013   Matt Alvarado     Moved a lot of output code to             !!
!!                                  OutputRputines.f90                       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM ASP

	!! Include modules
	USE ModelParameters,      ONLY : SetDomainSize,                 &
                                         SetInitialDens,                &
                                         SetSectionalParameters,        &
                                         OutputDeckSubDir,              &
                                         Pi,                            &
                                         micron,                        &
                                         ThermoBulkMode,                &
                                         Monodisperse,                  &
                                         DissolutionEquilibriumorFlux,  &
                                         IncludeAerosols,               &
                                         LagrangianParcelModel,         &
                                         EulerianBoxModel,              &
                                         OpticalPropertiesModel,        &
                                         ThermoTest,                    &
                                         CoagTest,                      &
                                         CondTest,                      &
                                         SmogChamberModel,              &
                                         MixStep,                       &
                                         ChemStep,                      &
                                         CoagStep,                      &
                                         CondStep,                      &
                                         b_mix,                         &
                                         InversionHeight,               &
                                         UTCStartTime,                  &
                                         RunLength,                     &
                                         Hetero_Num, Hetero_Rad
                                         

	USE InfrastructuralCode,  ONLY : Transcript,                    &
                                         Error,                         &
                                         Real2Str,                      &
                                         Int2Str,                       &
                                         SetFileHandleCounter,          & 
                                         GetFileHandle,                 &
                                         ReturnFileHandle     

	USE GridPointFields,      ONLY : GetM,                          &
                                         GetTemp,                       &
                                         SetRelativeHumidity,           &
                                         GridGasChem,                   &
                                         UpdateChemicalConcentrations

	USE Chemistry,            ONLY : FindChem,                      &
                                         IfGas, IfAq,                   &
                                         AqMolecularMass,               &
                                         OrgMolecularMass,              &
                                         HowManyGasChems,               &
                                         HowManyEvolveGasChems,         &
                                         GasPhaseBackground,            &
                                         GasPhaseDepVel,                &
                                         SetChemistryParams,            &
                                         ReadPhotolysisRates,           &
                                         RecalculateReactionRates,      &
                                         RecalculatePhotoRates,         &
                                         RecalculateHeteroRates,        &
                                         StepGasChemistry

	USE Aerosols,             ONLY : Particles,                     &
                                         Particle,                      &
                                         AerosolModes,                  &
                                         BoundaryParticles,             &
                                         TerminalVelocity,              &
                                         ParticleMass,                  &
                                         InputFlagRatioOrMass,          &
                                         RecalculateRadius,             &
                                         RegridAerosol,                 &
                                         InitializeParticlesSectional,  &
                                         ReadDistributionsOrganic,      &
                               ReadBoundaryDistributionsOrganic,        &
                               PopulateParticlesSectionsRightAway,      &
                               SortBoundaryPart,                        &
                               SortAerosolatGridPointForCoagulation,    &
                                         ShellRefIndandRad,             &
                                         AerosolOptProp,                &
                               FindAqElectrolyteEquilibriumForGridPoint
                                         		
	USE Condensation,         ONLY : SetAllDissolution,             &
                                         OrganicDissolutionData,        &
                                         StepCondensationAll,           &
                                         EqAllWater,                    &
                               EquilibrateInternallyAtGridPoint,        &
                                         EquilibrateGridPoint

	USE OutputRoutines,       ONLY : SpillBeans,                    &
                                         SizeDistOutput,                &
                                         GasOutput,                     &
                                         GasEnhanceOutput,              &
                                         AerosolOutput

	USE Coagulation,          ONLY : StepSectionalCoagulationJacobson
	
	IMPLICIT NONE

	!! Define the local variables
	INTEGER :: i, j, q, r, NumBins, HeteroBins, &
	           HNO3GasIndex, NumAqChems, NH3GasIndex, PotassiumIndex, &
                   HClGasIndex, SulfateIndex, BisulfateIndex, &
                   K2SO4Index, SOA2GasIndex, NumOrgChems, &
                   KNO3i, KCli, KHSO4i, K2SO4i, UR21i, N, HNO3Index, &
                   NH3Gas, HNO3Gas, HClGas, H2SO4Gas, NumAqOrgChems, &
                   ReturnType, Hi, NumModes, clock, BCi, status, &
                   NaHSO4i, Na2SO4i, NH4HSO4i, NH42SO4i, &
                   SO4i, HSO4i,  LEVi, FH1, EnhFH, NH3i, NH4i, NH4NO3i, &
                   NH4Cli, NO2i, ETHEi, HCHOi, Aceti, NOi

	REAL*8  :: ii, T_mix, Tau_Dep, slope, b, ChemTimestep, MixTimestep, &
		   CoagTimeStep, CondTimeStep, t, Mgas, Temp, &
                   NewRH, HNO3Conc, NH3Conc, ut, &
                   Potassium, HClConc, TotalNumber, TotalVolume, &
                   TotalPotassiumMass, RH, &
                   K2SO4Mass, Kmass, TotalMass, TimeStepSize, PressEq, &
                   P1, EnvMass, InitialMass, NewMass, ChemPercentage, sum, &
                   Total, Yield, Store_H2O, Store_Dens,InvHgt, BC, Sulf, &
                   OA, TotalNH4, TermV

	REAL*4  :: etimearray(2), etime
	
	TYPE(Particle),    POINTER :: Current, Env
	
	LOGICAL :: GADS

	INTEGER :: time(8)

	REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSingScat
	REAL*8	:: Term(250)

        INTEGER, ALLOCATABLE :: outFH(:)
	REAL*8, ALLOCATABLE  :: Chem(:), NumConc(:), &
				HeteroNumberConc(:), HeteroRadii(:), &
                                InitNumConc(:), EnvNumConc(:)

	
	!! Initialize the random number generator 
	CALL Random_Seed()

	!! Send a Welcome Message
	CALL Transcript ("**********************************************************")
	CALL Transcript ("** Welcome to the ASP Model developed by Matt Alvarado  **")
	CALL Transcript ("** of AER (malvarad@aer.com).                           **")
	CALL Transcript ("**********************************************************")

	!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SET the DOMAIN SIZE !!
        !! MJA 02-15-2012 always 1 cm^3
	CALL SetDomainSize ()
	
	!! Call subroutines to set module data
	CALL SetFileHandleCounter(12)
	CALL ReadMainInputDeck()

	IF(EulerianBoxModel) THEN
           CALL TRANSCRIPT("")
	   CALL TRANSCRIPT("Eulerian Box Model Selected.")
           CALL EulerianBox()

        ELSE IF(OpticalPropertiesModel) THEN !Optical Properties Only
         !MJA, 02-28-2012: The general idea here is to:
         ! 1. Read the size distribution parameters, including multiple modes
         ! 2. Read the non-H2O species as relative mass fractions
         ! 3. For each mode, create a monodisperse aerosol and calculate the 
         !    H2O content in equilibrium with the non-H2O species
         !    and the particle density.
         ! 4. Use the calculated relative mass fraction of H2O and 
         !    density from step 3 to initialize a sectional distribution for this mode.
         ! 5. Sum all modes..
         ! 6. Calculate aerosol optical properties.
         ! WARNING: This works for multiple modes, but calculates the density and 
         ! water mass fraction based on the first mode only. 
         ! So if the composition of the modes is radically different, you
         ! could have problems.
          
             CALL TRANSCRIPT("")
	     CALL TRANSCRIPT("Optical Properties Model Selected.")
             GADS = .FALSE. !Flag for GADS comparison

             ! INITIALIZE the CHEMISTRY 
	     CALL SetChemistryParams

	     !! INITIALIZE the DISSOLUTION routine
	     CALL SetAllDissolution
             
             ii = etime(etimearray)
	     CALL TRANSCRIPT("")
	     CALL TRANSCRIPT("CPU time after chem initialization "//trim(real2str(ii)))

             Particles%First => Null()
             CALL SetInitialDens(1.500)
	     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     !! Open the Particle Modes Input Deck and Read It In !!
	     !! And then use it to populate the AerosolModes Array!!
	     CALL ReadDistributionsOrganic() !!!!!!!!!!!!!!!!!!!!!!!

             !!Make and equilibrate a Monodisperse aerosol
             Monodisperse = .TRUE.
             InputFlagRatioOrMass = 1
             Call PopulateParticlesSectionsRightAway(.FALSE.)
	     !Output size distribution
	     CALL SizeDistOutput(0.0,"SizeDist_Opt_Mono.txt")
             
             ii = etime(etimearray)
	     CALL TRANSCRIPT("")
	     CALL TRANSCRIPT("CPU time after mono particle initialization "//trim(real2str(ii)))

             !!WARNING Force initial water equilibrium between gas and aerosol
	     !!(Uses open system for water eq.)
	     CALL EquilibrateInternallyAtGridPoint(EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
	     CALL SpillBeans() 

             ii = etime(etimearray)
	     CALL TRANSCRIPT("")
	     CALL TRANSCRIPT("CPU time after mono particle eq "//trim(real2str(ii)))

             !Save water amount and density to initialize size distribution
             current => particles%first
             Store_H2O = current%AqChems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1) !18.015 !ug/m3
             Store_Dens = current%ParticleDensity
             !NumAqChems = size(current%AqChems)
             !Need to subtract out ions
             Hi = FindChem("H+", 1)
             NumAqChems = Hi-1
	     NumOrgChems = size(current%OrgChems)
             NumModes = ubound(AerosolModes, 1)
             !WRITE(*,*) "NumModes: ", NumModes
             !WRITE(*,*) "Store_Dens: ", Store_Dens
             !WRITE(*,*) "NumAqChems: ", NumAqChems
             !WRITE(*,*) "NumOrgChems: ", NumOrgChems

             !Set Sectional Distribution
             Monodisperse = .FALSE.
             Particles%First => Null()
             DO I = 1, NumModes
               CALL SetInitialDens(Store_Dens)

               IF (.NOT.(GADS)) THEN
                  AerosolModes(I,4) = Store_H2O 
               ELSE 
                  AerosolModes(I,4) = 0.0
               ENDIF
               InputFlagRatioOrMass = 0
               !WRITE(*,*) "Store_H2O: ", Store_H2O
               !!Renormalize mass fractions including water
               ChemPercentage = 0
               Sum = 0
               DO J = 4, NumAqChems+NumOrgChems+3
			!WRITE(*,*) J, AerosolModes(I,J)
                        ChemPercentage = AerosolModes(I,J) + ChemPercentage
               END DO
               DO J = 4, NumAqChems+NumOrgChems+3
			AerosolModes(I,J) = AerosolModes(I,J)/ChemPercentage
                        Sum = Sum + AerosolModes(I,J)
               END DO
               !WRITE(*,*) "After Renorm: ", Sum, I, AerosolModes(I,4), ChemPercentage
             ENDDO

             !!Reset Sectional Parameters, currently hardwired distribution
             CALL SetSectionalParameters (10.0*micron,0.005*micron,40,.FALSE.)
             !CALL SetSectionalParameters (20.0*micron,0.005*micron,100,.FALSE.)
 
             !!CALL ReadDistributionsOrganic()
             InputFlagRatioOrMass = 0
             Call PopulateParticlesSectionsRightAway(.FALSE.)
             ii = etime(etimearray)
	     CALL TRANSCRIPT("")
	     CALL TRANSCRIPT("CPU time after sectional particle initialization "//trim(real2str(ii)))
             CALL SizeDistOutput(0.0,"SizeDist_Opt_Sect_init.txt")

             CALL EquilibrateInternallyAtGridPoint(EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
             CALL SpillBeans()
	     !Output size distribution
	     CALL SizeDistOutput(0.0,"SizeDist_Opt_Sect.txt")

             ii = etime(etimearray)
	     CALL TRANSCRIPT("")
	     CALL TRANSCRIPT("CPU time after sectional particle eq "//trim(real2str(ii)))

             !Calculate optical properties for each particle
	     current => particles%first
	     I = 1
	     DO WHILE(associated(current))
		CALL RecalculateRadius(current)
			
		!Recalculate optical parameters
		CALL ShellRefIndAndRad(Current)
			
		I = I + 1
		
		current => current%next
				
	      END DO !END STEP Particle Optical Properties

	      !Core-in-shell
              CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSingScat, 0)

              ALLOCATE(OutFH(4), STAT = status)
              IF (status > 0) THEN
                CALL ERROR("Allocation error for OutFH in Optical Property Model of DevelopmentDriver.f90.")
              ENDIF

	      OutFH(1) = GetFileHandle()
              OPEN(UNIT=OutFH(1), FILE=TRIM(OutputDeckSubDir)//"AeroOptCIS.txt", STATUS="REPLACE")
              WRITE(OutFH(1), FMT='(7A17)') "Wavelength", "ExtCoeff (Mm^-1)", "SSA", &
                                            "<g>", "BackScat (Mm^-1)", "SubExtCoeff", "SubSSA"
              DO I = 1, 451 
                 WRITE(OutFH(1), FMT='(F17.1,6E17.6)') 250.0+(I-1), ExtCoeff(I), SingScat(I), &
                                       Assym(I), BackScatCoeff(I), SubExtCoeff(I), SubSingScat(I)
              ENDDO

	      !External Mixture
              CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSingScat, 1)

	      OutFH(2) = GetFileHandle()
              OPEN(UNIT=OutFH(2), FILE=TRIM(OutputDeckSubDir)//"AeroOptEXT.txt", STATUS="REPLACE")
              WRITE(OutFH(2), FMT='(7A17)') "Wavelength", "ExtCoeff (Mm^-1)", "SSA", "<g>", &
                                             "BackScat (Mm^-1)", "SubExtCoeff", "SubSSA"
              DO I = 1, 451 
                 WRITE(OutFH(2), FMT='(F17.1,6E17.6)') 250.0+(I-1), ExtCoeff(I), SingScat(I), &
                                       Assym(I), BackScatCoeff(I), SubExtCoeff(I), SubSingScat(I)
              ENDDO	

	      !Volume average dielectric constant mixing rule
              CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSingScat, 2)

	      OutFH(3) = GetFileHandle()
              OPEN(UNIT=OutFH(3), FILE=TRIM(OutputDeckSubDir)//"AeroOptVA.txt", STATUS="REPLACE")
              WRITE(OutFH(3), FMT='(7A17)') "Wavelength", "ExtCoeff (Mm^-1)", "SSA", &
                                            "<g>", "BackScat (Mm^-1)", "SubExtCoeff", "SubSSA"
              DO I = 1, 451 
                 WRITE(OutFH(3), FMT='(F17.1,6E17.6)') 250.0+(I-1), ExtCoeff(I), SingScat(I), &
                                      Assym(I), BackScatCoeff(I), SubExtCoeff(I), SubSingScat(I)
              ENDDO

	      !Maxwell-Garnett mixing rule
              CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSingScat, 3)

	      OutFH(4) = GetFileHandle()
              OPEN(UNIT=OutFH(4), FILE=TRIM(OutputDeckSubDir)//"AeroOptMG.txt", STATUS="REPLACE")
              WRITE(OutFH(4), FMT='(7A17)') "Wavelength", "ExtCoeff (Mm^-1)", "SSA", &
                                            "<g>", "BackScat (Mm^-1)", "SubExtCoeff", "SubSSA"
              DO I = 1, 451 
                 WRITE(OutFH(4), FMT='(F17.1,6E17.6)') 250.0+(I-1), ExtCoeff(I), SingScat(I), &
                                            Assym(I), BackScatCoeff(I), SubExtCoeff(I), SubSingScat(I)
              ENDDO

              DEALLOCATE(OutFH, STAT = status)

        ELSE IF(LagrangianParcelModel) THEN 
           CALL TRANSCRIPT("")
	   CALL TRANSCRIPT("Lagrangian Parcel Model Selected.")

           !Convert Inversion Height from km to cm
           InvHgt = InversionHeight*1.0e5

           IF(.NOT.(IncludeAerosols)) THEN !Gas Chemistry Only
	        CALL TRANSCRIPT("")
	        CALL TRANSCRIPT("Gas Chemistry Only Selected.")
		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
	
		!Calculate reaction rates for non-photolysis and non-hetero gas rxns
		CALL RecalculateReactionRates(IFGAS,IFAQ)
		
		!Calculate initial hetero. Rxn rates
		HeteroBins = 1
		ALLOCATE(HeteroNumberConc(HeteroBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for HeteroNumberConc in " &
                             //"Gas Only Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		ALLOCATE(HeteroRadii(HeteroBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for HeteroRadii in " &
                             //"Gas Only Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF
		HeteroNumberConc(1) = hetero_num !For Williams fire, see file 
                !/nas/project/p1704/asp_rundir/Yokelson/README_Williams_aero_init
                !1.5e5 !Guess for Williams fire
                !3.0e4 !Guess for Alaska
		!3.1e4 !Otavi (Haywood et al, 2003)  
		!1.433e5 !Timbavati (Trentmann et al, 2005)
		HeteroRadii(1) = hetero_rad ! in m
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		
		!Calculate Initial Photolysis Rates
		ut = UTCStartTime !Yokelson Williams fire, 11 AM local time
                !0.217 !MILAGRO Yucatan fire: in the fire measuring photolysis rates  at 73007, 
                       !smoke age starts at 72989
                !3.0 !Yokelson Williams fire, 11 AM local time
                !19.0 !Yokelson Lake McKay fire, 11 AM local time
                !2.51 !ARCTAS Lake McKay Late Pass time in hours 
                ! 0.87 !Alaska Universal time in hours
		!11.00 !Otavi Initial Universal Time in hours
		!8.75 !Timbavati Initial Universal Time in hours
	
                CALL RecalculatePhotoRates (ut, .TRUE.)
      
		ALLOCATE(Chem(HowManyGasChems), STAT=status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for Chem in " &
                             //"Gas Only Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		TEMP = GetTemp()
		Mgas = GetM   ()
         
		ALLOCATE(OutFH(4), STAT=status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for OutFH in " &
                             //"Gas Only Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		OutFH(1) = GetFileHandle()
		OutFH(2) = GetFileHandle()
		OutFH(3) = GetFileHandle()
                
		OPEN(UNIT=OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		OPEN(UNIT=OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		OPEN(UNIT=OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")

                CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),0.0,.TRUE.,.FALSE.)	
	
		!MJA, 12-28-2012
                !New output file for enhancement ratios
                !Of important species for GEOS-Chem study	
                OutFH(4) = GetFileHandle()
	        OPEN(UNIT=OutFH(4), & 
                     FILE=TRIM(OutputDeckSubDir)//"EnhanceRatios.txt", &
                     STATUS="REPLACE")
                CALL GasEnhanceOutput(OutFH(4),0.0,.TRUE.)

		t=0
		ChemTimestep = 1
		MixTimestep = 1 !Always start as 1s for first 600 seconds, MJA 12-28-2012
		DO r = 1, INT(600+((RunLength-10)*60/MixStep))+1
		
                   !Write Gas Output to files 
                   IF (MOD(t,60.0)<=0.05) THEN
                      CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),t,.FALSE.,.FALSE.)
                      CALL GasEnhanceOutput(OutFH(4),t,.FALSE.)
                   END IF !Write Output

                   CALL RecalculateReactionRates(IFGAS,IFAQ)
                   !WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
                   !after this step, or those rates will be set to 0!
			
                   !Recalculate Photolysis Rates
                   CALL RecalculatePhotoRates(ut, .FALSE.)
                        
                   !Recalculate Heterogeneous Rates
                   CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
			
                   !Step Gas Chemistry
                   CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep))

                   Slope = 2.00
                   b=b_mix
                   !15  !Yokelson Williams Fire fast 
                   !106.87357 !Yokelson Williams Fire (California, Nov. 17 2009) 
                               !estimated from curvefit of CO data
                   !212.22/3600.0 !Dilution time, Williams, slow, in hr.
                   !4000 ! slow fit for MILAGRO Fire (Yucatan March 23, 2006)
                   !100 !fast fit for MILAGRO Fire (Yucatan March 23, 2006)
                   !529.7 ! MILAGRO Fire (Yucatan March 23, 2006), middle curvefit
                   !7884 !ARCTAS Lake McKay Late Pass from Jason St. Clair peak CO data  
                                 !(Ky = 8632 m2/s, yo = 16.5 km)
    
                   !1387.1  !Alaska Base
                   !1731.0 !Arctas Early pass
                   !1387.1 !Alaska Base
                   !7000.0 !Alaska Slow
                   !350.0 !Alaska Fast
	
                   !156.25 !Otavi
			
                   !142 + (6.0/7.0) !Timbavati

                   !400 !Timbavati slow
                   !66 + (2.0/3.0) !Timbavati Fast  
                   !1.0e50 !No Dilution 	
                   !12000	 
                   T_mix = Slope*t + b
		
                   !Step Lagrangian Dilution of Gas Phase Chems
                   DO q = 1, HowManyGasChems
                      Chem(q) = GridGasChem(q)
                      Chem(q) = (Chem(q) + (GasPhaseBackground(q)*MixTimestep/T_mix)) &
                           /(1+ MixTimestep/T_mix+ MixTimestep*GasPhaseDepVel(q)/InvHgt)
                   END DO

                   !Update chemical concentrations
                   CALL UpdateChemicalConcentrations(Chem(1:HowManyEvolveGasChems))
			
                   !Step Dilution of "particles" for heterogeneous chemistry
                   DO q = 1, HeteroBins
                      HeteroNumberConc(q) = HeteroNumberConc(q)/(1+ MixTimestep/T_mix)
                   END DO

                   !Update time since start and UTC time for photolysis
                   t = t + MixTimestep
                   ut = ut + (MixTimestep/3600.0)
	        
                   IF(ut .GE. 24.0) ut = 0.0
                   IF(t .GE. 600) THEN !Reset to input mix time after 600 s
                      ChemTimestep = ChemStep
                      MixTimestep = MixStep    
                   ENDIF
                   WRITE(*,*) t	
		END DO !Time Loop

                DEALLOCATE(HeteroNumberConc, HeteroRadii, Chem, OutFH, STAT=status)

	    ELSE IF (IncludeAerosols) THEN!Bulk, Mono, or Size dist aerosols
	        CALL TRANSCRIPT("")
	        CALL TRANSCRIPT("Gas and Aerosol Chemisty Selected.")

		! INITIALIZE the CHEMISTRY 
		CALL SetChemistryParams
		CALL ReadPhotolysisRates
		
		CALL RecalculateReactionRates(IFGAS,IFAQ)
			
		!Calculate Initial Photolysis Rates
		ut = UTCStartTime !Timbavati Initial Universal Time in hours
		!11.00 !Otavi Initial Universal Time in hours
		CALL RecalculatePhotoRates (ut, .TRUE.)

		ALLOCATE(Chem(HowManyGasChems), STAT=status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for Chem in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF
	
		!! INITIALIZE the DISSOLUTION routine
		CALL SetAllDissolution

		!! INITIALIZE the PARTICLES	
		CALL InitializeParticlesSectional
		CALL RegridAerosol ()

		!Find chemical indices for later
		PotassiumIndex = FindChem("K+", 1)
		KNO3i = FindChem("KNO3", 1)
		KCli = FindChem("KCl", 1)
		KHSO4i = FindChem("KHSO4", 1)
		K2SO4Index = FindChem("K2SO4", 1)
		K2SO4i = K2SO4Index

                NaHSO4i = FindChem("NaHSO4", 1)
                Na2SO4i = FindChem("Na2SO4", 1)
                NH4HSO4i = FindChem("NH4HSO4", 1)
                NH42SO4i = FindChem("(NH4)2SO4", 1)
                LEVi = FindChem("(NH4)3H(SO4)2", 1)
                SO4i = FindChem("SO4--", 1)
                HSO4i = FindChem("HSO4-", 1)

                NH3i = FindChem("NH3", 1)
                NH4i = FindChem("NH4+", 1)
                NH4NO3i = FindChem("NH4NO3", 1)
                NH4Cli = FindChem("NH4Cl", 1)

		!Count # of bins, # of Aq. Chems
		current => particles%first
		NumBins = 0
		Potassium = 0.0
                BC = 0.0
                Sulf = 0.0
		DO WHILE(associated(current))
                   NumBins = NumBins + 1
                   NumAqChems = size(current%AqChems)
                   NumOrgChems = size(current%OrgChems)
                   NumAqOrgChems = size(current%AqOrgChems)
                   
                   BC        = BC + current%OrgChems(1)*current%numberofparticles

                   Potassium = Potassium + current%AqChems(PotassiumIndex)*current%numberofparticles
                   Potassium = Potassium + current%AqChems(KNO3i)*current%numberofparticles
                   Potassium = Potassium + 2*current%AqChems(K2SO4i)*current%numberofparticles
                   Potassium = Potassium + current%AqChems(KCli)*current%numberofparticles
                   Potassium = Potassium + current%AqChems(KHSO4i)*current%numberofparticles

                   Sulf = Sulf +  current%AqChems(K2SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(KHSO4i)*current%numberofparticles
                   Sulf = Sulf +  current%AqChems(Na2SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(NaHSO4i)*current%numberofparticles
                   Sulf = Sulf +  current%AqChems(NH42SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(NH4HSO4i)*current%numberofparticles
                   Sulf = Sulf +  2*current%AqChems(LEVi)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(HSO4i)*current%numberofparticles
                		
                   current => current%next
		END DO

		!WRITE(*,*) "Potassium (ug/m3): ", Potassium*AqMolecularMass(PotassiumIndex)*1.0E12
		!WRITE(*,*) "BC        (ug/m3): ", BC*OrgMolecularMass(1)*1.0E12
                !WRITE(*,*) "Sulfate   (ug/m3): ",Sulf*AqMolecularMass(SO4i)*1.0E12
		!CALL SpillBeans()
                
		!! Call Boundary distributions
		CALL ReadBoundaryDistributionsOrganic
		CALL PopulateParticlesSectionsRightAway(.TRUE.)
		CALL SortBoundaryPart
		!CALL SpillBeans(Env=.TRUE.)
		
		!!WARNING Force initial water equilibrium between gas and aerosol
		!!(Uses open system for water eq.)
		CALL EquilibrateInternallyAtGridPoint(EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
		CALL RegridAerosol ()
		!CALL SpillBeans()
                	
		TEMP = GetTemp()
		Mgas = GetM   ()

		!Count # of bins, # of Aq. Chems
		current => boundaryparticles%first
		NumBins = 0
		Potassium = 0.0
                BC = 0.0
                Sulf = 0.0
                OA = 0.0
                TotalNH4 = 0.0
		DO WHILE(associated(current))
                   NumBins = NumBins + 1
                   NumAqChems = size(current%AqChems)
                   NumOrgChems = size(current%OrgChems)
                   NumAqOrgChems = size(current%AqOrgChems)

                   BC        = BC + current%OrgChems(1)*current%numberofparticles

                   OA        = OA + current%OrgChems(14)*current%numberofparticles

                   Potassium = Potassium + current%AqChems(PotassiumIndex)*current%numberofparticles
                   Potassium = Potassium + current%AqChems(KNO3i)*current%numberofparticles
                   Potassium = Potassium + 2*current%AqChems(K2SO4i)*current%numberofparticles
                   Potassium = Potassium + current%AqChems(KCli)*current%numberofparticles
                   Potassium = Potassium + current%AqChems(KHSO4i)*current%numberofparticles

                   Sulf = Sulf +  current%AqChems(K2SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(KHSO4i)*current%numberofparticles
                   Sulf = Sulf +  current%AqChems(Na2SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(NaHSO4i)*current%numberofparticles
                   Sulf = Sulf +  current%AqChems(NH42SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(NH4HSO4i)*current%numberofparticles
                   Sulf = Sulf +  2*current%AqChems(LEVi)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(SO4i)*current%numberofparticles  
                   Sulf = Sulf +  current%AqChems(HSO4i)*current%numberofparticles

                   TotalNH4 = TotalNH4 + current%Aqchems(NH3i)*current%numberofparticles
                   TotalNH4 = TotalNH4 + current%Aqchems(NH4i)*current%numberofparticles
                   TotalNH4 = TotalNH4 + current%Aqchems(NH4NO3i)*current%numberofparticles
                   TotalNH4 = TotalNH4 + current%Aqchems(NH4Cli)*current%numberofparticles
                   TotalNH4 = TotalNH4 + current%Aqchems(NH4HSO4i)*current%numberofparticles
                   TotalNH4 = TotalNH4 + 2*current%Aqchems(NH42SO4i)*current%numberofparticles
                   TotalNH4 = TotalNH4 + 3*current%Aqchems(LEVi)*current%numberofparticles

                   current => current%next
		END DO

		!WRITE(*,*) "Env Ammonium  (ug/m3): ", TotalNH4*AqMolecularMass(NH4i)*1.0E12
		!WRITE(*,*) "Env BC        (ug/m3): ", BC*OrgMolecularMass(1)*1.0E12
                !WRITE(*,*) "Env Sulfate   (ug/m3): ",Sulf*AqMolecularMass(SO4i)*1.0E12
                !WRITE(*,*) "Env OA        (ug/m3): ",OA*OrgMolecularMass(14)*1.0E12                

		ALLOCATE(NumConc(NumBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for NumConc in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		ALLOCATE(InitNumConc(NumBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for InitNumConc in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		ALLOCATE(EnvNumConc(NumBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for EnvNumConc in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF
		
		!Set heterogeneous Reactions
		HeteroBins = NumBins
		ALLOCATE(HeteroNumberConc(HeteroBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for HeteroNumberConc in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF
		ALLOCATE(HeteroRadii(HeteroBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for HeteroRadii in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		!Get particle number conc. and radii for hetero. rate calculation
		IF(ThermoBulkMode) THEN
		    HeteroNumberConc(1) = hetero_num 
	  	    HeteroRadii(1) = hetero_rad ! in m
                ELSE
                  current => particles%first
		  I = 1
		  DO WHILE(associated(current))
                     HeteroNumberConc(I) = current%NumberofParticles
                     HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
                     !ResetHeteroRates expects radii in units of m
                     current => current%next
                     I = I+1
		  END DO
		ENDIF
		
		!Calculate Heterogeneous rates from size dist. info
		CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		CALL RecalculatePhotoRates(ut, .TRUE.)
		WRITE(*,*) "After Photo Rates"			
		TEMP = GetTemp()
		Mgas = GetM   ()

		ALLOCATE(OutFH(3+NumBins), STAT = status)
                IF (status > 0) THEN
                  CALL ERROR("Allocation error for OutFH in " &
                             //"Gas&Aero Lagrangian Parcel Model of DevelopmentDriver.f90.")
                ENDIF

		DO I = 1, 3+Numbins
			OutFH(I) = GetFileHandle()
		END DO
	
		OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")

                CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),0.0,.TRUE.,.FALSE.)	

		!MJA, 12-28-2012
                !New output file for enhancement ratios
                !Of important species for GEOS-Chem study	
                EnhFH = GetFileHandle()
	        OPEN(UNIT=EnhFH, & 
                     FILE=TRIM(OutputDeckSubDir)//"EnhanceRatios.txt", &
                     STATUS="REPLACE")
                CALL GasEnhanceOutput(EnhFH,0.0,.TRUE.)

		WRITE(*,*) "After GasOutput1"		
		current => particles%first
		DO I = 4, 3 + NumBins
                   OPEN(UNIT = OutFH(I), FILE=TRIM(OutputDeckSubDir)//"Aerosol"//TRIM(INT2STR(I-3))&
                        //".txt", STATUS="REPLACE")
                   CALL AerosolOutput(current, OutFH(I), 0.0, .TRUE.,.FALSE.)
                   IF (associated(current%next)) THEN
                      current => current%next
                   ELSE
                      CYCLE
                   END IF
		END DO

		t=0	

                CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),t,.FALSE.,.FALSE.)	
                CALL GasEnhanceOutput(EnhFH,t,.FALSE.)
		WRITE(*,*) "After GasOutput2"

		current => particles%first
		DO I = 4, 3+NumBins
                   CALL AerosolOutput(current, OutFH(I), t, .FALSE.,.FALSE.)
					
                   IF (associated(current%next)) THEN
                      current => current%next
                   ELSE
                      CYCLE
                   END IF

		END DO
		
		!Output size distribution
		CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
	
		Slope = 2.0 
		b = b_mix 
			
		ChemTimestep = 1
		MixTimestep = 1 !Always start as 1s for first 600 seconds, MJA 12-28-2012
		CoagTimeStep = CoagStep
		CondTimeStep = CondStep
		
		!Run the main loop with above timesteps (1s timestep to start)
		DO r = 1, INT(600+((RunLength-10)*60/MixStep))
                   WRITE(*,*) "In Loop"
			
                   !Calculate Temperature
                   TEMP = GetTemp()

                   CALL RecalculateReactionRates(IFGAS,IFAQ)
                   !WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
                   !after this step, or those rates will be set to 0!
			
                   !Recalculate Photolysis Rates
                   CALL RecalculatePhotoRates(ut, .FALSE.)
                   WRITE(*,*) "Photo Okay"

                   !Step Dilution of bulk "particles" for heterogeneous chemistry
                   IF (ThermoBulkMode) THEN
                      if (r .NE. 1) HeteroNumberConc(1) = HeteroNumberConc(1)/(1+ MixTimestep/T_mix)
                   ELSE!Get particle number conc. and radii for hetero. rate calculation		  
                      current => particles%first
                      I = 1
                      DO WHILE(associated(current))
                         HeteroNumberConc(I) = current%NumberofParticles
                         HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
                         !ResetHeteroRates expects radii in units of m
                         current => current%next
                         I = I+1
                      END DO
                   ENDIF
	
                   !Calculate Heterogeneous rates from size dist. info
                   CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
                   WRITE(*,*) "Hetero Okay"
						
                   !Step Gas Chemistry
                   CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep))
                   WRITE(*,*) "StepGasChem Okay"
			
                   !Step Condensation
                   IF (MOD(t,CondTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN
                      !Step Condensation, if desired
                      IF(.TRUE.) THEN
                         CALL RegridAerosol ()
                         CALL FindAqElectrolyteEquilibriumForGridPoint ( UpdateThermo = .TRUE.) 
                         WRITE(*,*) "AqEq Okay"
                         !CALL SpillBeans()
                         CALL StepCondensationAll (CondTimeStep,  1)
                         !Note above calls RegridAerosol automatically
                         WRITE(*,*) "Cond Okay"
                      END IF
                   END IF

                   !Step Coagulation
                   IF (MOD(t,CoagTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN

                      !Step Coagulation, if desired
                      CALL StepSectionalCoagulationJacobson(CoagTimeStep)
                      !Note that above already calls RegridAerosol
                      WRITE(*,*) t, "Coag Okay"
					
                   END IF

                   !Calculate T_mix
                   T_mix = Slope*t+b
		
                   !Step Dilution of Chems
                   DO q = 1, HowManyGasChems
                      Chem(q) = GridGasChem(q)
                      Chem(q) = (Chem(q) + (GasPhaseBackground(q)*MixTimestep/T_mix)) &
                           /(1+ MixTimestep/T_mix + MixTimestep*GasPhaseDepVel(q)/InvHgt)
                   END DO

!FOR SC ONLY, TRY PUMPING UP NOX!
!                   if(t .gt. 2695 .and. t .lt. 2705) then
!                    if(t .gt. 1795 .and. t .lt. 1805) then
!                          NO2i = FindChem("NO2", 0)
!                          NOi = FindChem("NO", 0)
!                          ETHEi = FindChem("HCHO", 0)
!                          HCHOi = FindChem("ETHE", 0)
!                          Aceti = FindChem("Acetaldehyde",0)
!                          Mgas = GetM()
!                          Chem(NO2i) = 60.0*1.0e-9*Mgas                          
!                          Chem(NOi) = 20.0*1.0e-9*Mgas  
!                          Chem(ETHEi) = 3.0*1.0e-9*Mgas   
!                          Chem(HCHOi) = 5.5*1.0e-9*Mgas   
!                          Chem(Aceti) = 1.0*1.0e-9*Mgas   
!                   endif


                   !Update chemical concentrations
                   CALL UpdateChemicalConcentrations( Chem(1:HowManyEvolveGasChems))

                   !Step Dilution and Deposition of Particles
                   CALL SortAerosolAtGridPointForCoagulation () 
                   current => particles%first
                   Env => BoundaryParticles%first
                   I = 1
                   DO WHILE(associated(current))
                      
                      IF(Env%Edges(1) .NE. Current%Edges(1)) THEN
                         !WRITE(*,*) Env%Edges(1), Current%Edges(1)
                         !CALL ERROR("Env. Aerosol not lined up!")					
                      END IF
                      CALL RecalculateRadius(current)

                      !Get Particle Terminal Velocity
                      IF (ThermoBulkMode) THEN 
                          TermV = 0.0 
                      ELSE 
                          TermV = TerminalVelocity(Current)
                      ENDIF

                      !Get Number Concentration of particles in particles/cm3
                      InitNumConc(I) = current%numberofparticles 
                      EnvNumConc(I) = Env%numberofparticles
                      
                      !Dilute and Dep number of particles (assume no background for now)
                      current%numberofparticles = (InitNumConc(I) + (EnvNumConc(I)*MixTimestep/T_mix))/ &
                           (1.0 + MixTimestep*(1.0/T_mix + TermV/InvHgt))
                      !WRITE(*,*) I, InitNumConc(I), current%numberofparticles, TerminalVelocity(Current)
                      IF(current%Numberofparticles .GT. 0.) THEN
                         DO J = 1, NumAqChems
                            InitialMass = InitNumConc(I)*current%AqChems(J)
                            EnvMass = EnvNumConc(I)*Env%AqChems(J)
                            NewMass = (InitialMass + EnvMass*MixTimestep/T_mix)/ &
                                 (1+ MixTimestep*(1.0/T_mix + TermV/InvHgt))
                            Current%AqChems(J) = NewMass/current%numberofparticles
					
                         END DO
				
                         DO J = 1, NumOrgChems
                            InitialMass = InitNumConc(I)*current%OrgChems(J)
                            EnvMass = EnvNumConc(I)*Env%OrgChems(J)
                            NewMass = (InitialMass + EnvMass*MixTimestep/T_mix)/ &
                                 (1+ MixTimestep*(1.0/T_mix + TermV/InvHgt))
                            Current%OrgChems(J) = NewMass/current%numberofparticles
                         END DO
				
                         DO J = 1, NumAqOrgChems
                            InitialMass = InitNumConc(I)*current%AqOrgChems(J)
                            EnvMass = EnvNumConc(I)*Env%AqOrgChems(J)
                            NewMass = (InitialMass + EnvMass*MixTimestep/T_mix)/ &
                                 (1+ MixTimestep*(1.0/T_mix + TermV/InvHgt))
                            Current%AqOrgChems(J) = NewMass/current%numberofparticles
                         END DO
                      ELSE !No particles
                         DO J = 1, NumAqChems
                            Current%AqChems(J) = 0.
                            
                         END DO
				
                         DO J = 1, NumOrgChems
                            Current%OrgChems(J) = 0.
                         END DO
				
                         DO J = 1, NumAqOrgChems
                            Current%AqOrgChems(J) = 0.
                         END DO
                      END IF
				
                      !Force particle temperature to be env. Temperature
                      !WARNING: This is a kludge, since the program is calculating
                      !a huge temperature for the largest particles, and I'm not sure why
                      current%Temperature = TEMP

                      CALL RecalculateRadius(current)
				
                      I = I + 1
		
                      current => current%next
                      Env => Env%next
				
                   END DO !END STEP PARTICLE DILUTION

                        
                   WRITE(*,*) "Before AOD"
                   !Call SpillBeans
		   !WRITE(*,*) "After SpillBeans"
                   !Only call AOD every 10 min
                   IF (MOD(t,600.0) .LT. 0.01) THEN
                      current => particles%first
                      I = 1
                      DO WHILE(associated(current))
                         !WRITE(*,*) "Particle Number ", I
                         CALL ShellRefIndAndRad(Current)
                         I = I + 1
                         current => current%next
                      END DO
                      !WRITE(*,*) "Before AerosolOptProp"
                      CALL AerosolOptProp (ExtCoeff, SingScat, Assym, &
                           BackScatCoeff, SubExtCoeff, SubSingScat, 0)
                      FH1 = GetFileHandle()
                      OPEN(UNIT=FH1, FILE=TRIM(OutputDeckSubDir)//"AeroOptCIS" &
                           //TRIM(INT2STR(INT(t/60)))//"min.txt", STATUS="REPLACE")
                      WRITE(FH1, FMT='(7A17)') "Wavelength", "ExtCoeff (Mm^-1)", "SSA", &
                           "<g>", "BackScat (Mm^-1)", "SubExtCoeff", "SubSSA"
                      DO I = 1, 451 
                         WRITE(FH1, FMT='(F17.1,6E17.6)') 250.0+(I-1), ExtCoeff(I), SingScat(I), &
                              Assym(I), BackScatCoeff(I), SubExtCoeff(I), SubSingScat(I)
                      ENDDO
                      CALL ReturnFileHandle(FH1)
                   ENDIF
			
                   !Write Output to file every minute 
                   IF (MOD(t,60.0) .LT. 0.01) THEN
                      WRITE(*,*) "Before GasOutput3"
                      CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),t,.FALSE.,.FALSE.)
                      CALL GasEnhanceOutput(EnhFH,t,.FALSE.)    
                  
                      current => particles%first
                      DO I = 4, 3+NumBins
                         WRITE(*,*) "Before AerosolOutput, Bin #", I-3
                         CALL AerosolOutput(current, OutFH(I), t, .FALSE.,.FALSE.)
                         
                         IF (associated(current%next)) THEN
                            current => current%next
                         ELSE
                            CYCLE
                         END IF

                      END DO
                   END IF !Write Output

                   IF (MOD(t,600.0)<=0.05 .OR. ABS(t-540) <=0.05 .OR. ABS(t-2820) <=0.05) THEN
                      IF(INT(t) .NE. 0.) THEN
                         CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
                      END IF
                   END IF
			
                   !Update time
                   t = t + MixTimestep
                   ut = ut + MixTimeStep/3600.

                   IF(ut .GE. 24.0) ut = 0.0
                   IF(t .GE. 600) THEN !Reset to input mix time after 600 s
                      ChemTimestep = ChemStep
                      MixTimestep = MixStep    
                   ENDIF

                   WRITE(*,*), trim(real2str(t)), " seconds done"
                   IF (MOD(t,CoagTimeStep)<=0.01) THEN
                      CALL TRANSCRIPT(trim(real2str(t/60))//" minutes done")
                   END IF

		END DO

	    END IF !Aerosol Tests
                
            DEALLOCATE(NumConc, InitNumConc, EnvNumConc, HeteroNumberConc, HeteroRadii, &
                       Chem, OutFH, STAT=status)

	ELSE IF (ThermoTest) THEN !THERMO TESTS (Compare to ISORROPIA)
           CALL TRANSCRIPT("")
           CALL TRANSCRIPT("Bulk Thermodynamics Model Selected.")		
	
           ! INITIALIZE the CHEMISTRY 
           CALL SetChemistryParams
           CALL RecalculateReactionRates(IFGAS,IFAQ)
           
           !Set heterogeneous Reactions
           HeteroBins = 1
           ALLOCATE(HeteroNumberConc(HeteroBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for HeteroNumberConc in " &
                   //"Bulk Thermodynamics Model of DevelopmentDriver.f90.")
           ENDIF

           ALLOCATE(HeteroRadii(HeteroBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for HeteroRadii in " &
                   //"Bulk Thermodynamics Model of DevelopmentDriver.f90.")
           ENDIF

           HeteroNumberConc(1) = 1.433e5 !Set for Timbavati Conditions
           HeteroRadii(1) = 1.0e-7
           CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
		
           ALLOCATE(Chem(HowManyGasChems), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for Chem in " &
                   //"Bulk Thermodynamics Model of DevelopmentDriver.f90.")
           ENDIF

           !! INITIALIZE the DISSOLUTION routine
           CALL SetAllDissolution

           !WRITE(*,*) "Before Initialization"
		
           !! INITIALIZE the PARTICLES
           CALL InitializeParticlesSectional	

           IF(.NOT.(ThermoBulkMode)) THEN
              CALL ERROR("Cannot call Bulk Theronamics mode unless you also selected"// &
                   "bulk aerosol in AerosolModes.in")
           ENDIF
		
           !Spill Beans
           CALL SpillBeans()
                 
           NumBins = 1
		
           ALLOCATE(NumConc(NumBins), STAT=status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for NumConc in " &
                   //"Bulk Thermodynamics Model of DevelopmentDriver.f90.")
           ENDIF
		
           TEMP = GetTemp()
           Mgas = GetM   ()
                
           ALLOCATE(OutFH(4), STAT=status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for OutFH in " &
                   //"Bulk Thermodynamics Model of DevelopmentDriver.f90.")
           ENDIF

           OutFH(1) = GetFileHandle()
           OutFH(2) = GetFileHandle()
           OutFH(3) = GetFileHandle()
           OutFH(4) = GetFileHandle()

           OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOutThermo1.txt", STATUS="REPLACE")
           OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOutThermo2.txt", STATUS="REPLACE")
           OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOutThermo3.txt", STATUS="REPLACE")
           CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),0.0,First=.TRUE.,Thermo=.TRUE.)

                
           IF (NumBins .GT. 1) CALL ERROR("This model assummes a bulk aerosol - " &
                //"please fix your AerosolModes.in file.")

           OPEN(UNIT = OutFH(4), FILE=TRIM(OutputDeckSubDir)//"AerosolOutThermo.txt", STATUS="REPLACE")
           current => particles%first
           CALL AerosolOutput(current, OutFH(4), 0.0, First=.TRUE., Thermo=.TRUE.)

           !Spill Beans
           CALL SpillBeans()
           !WRITE(*,*) "After Spill Beans"
                

           NewRH = 0.92
           DO J = 1, 36
              NewRH = NewRH - 0.02		
              CALL SetRelativeHumidity (NewRH)
              
              WRITE(*,*) "Before Equilibrate"	, NewRH
              CALL EquilibrateGridPoint(EquilibrateWater=.TRUE., &
                   WaterOpenSystem=.TRUE., InorgOnly=.FALSE.)
              WRITE(*,*) "Equilibrate Okay"

              CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),NewRH,.FALSE.,.TRUE.)
		
              CALL AerosolOutput(current, OutFH(4), NewRH, .FALSE., .TRUE.)

              !Spill Beans
              CALL SpillBeans()
 
              !Pause the code
              !READ(*,*)
           END DO
           
           DEALLOCATE(NumConc, HeteroNumberConc, HeteroRadii, Chem, OutFH, STAT=status)
	
	ELSE IF (CoagTest) THEN !COAG TESTS 
           CALL TRANSCRIPT("")
           CALL TRANSCRIPT("Coagulation Test Selected.")
	
           ! INITIALIZE the CHEMISTRY 
           CALL SetChemistryParams
		
           !! INITIALIZE the PARTICLES
           CALL InitializeParticlesSectional

           !! Sort the gridcell, sections in order and first
           CALL SortAerosolAtGridPointForCoagulation ()
           CALL SpillBeans()

           !Count # of bins
           current => particles%first
           NumBins = 0
           DO WHILE(associated(current))
              NumBins = NumBins + 1
              current => current%next
           END DO
           ALLOCATE(NumConc(NumBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for NumConc in " &
                   //"Coagulation Test of DevelopmentDriver.f90.")
           ENDIF

           !Calculate Temperature
           TEMP = GetTemp()
			
           ALLOCATE(OutFH(4), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for OutFH in " &
                   //"Coagulation Test of DevelopmentDriver.f90.")
           ENDIF
           OutFH(4) = GetFileHandle()

           OPEN(UNIT = OutFH(4), FILE=TRIM(OutputDeckSubDir)//"CoagTestOut.txt", STATUS="REPLACE")
           WRITE(OutFH(4), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
           DO Q = 1, NumBins		
              WRITE(OutFH(4),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
              WRITE(OutFH(4),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
           END DO

           WRITE(OutFH(4), FMT='(A1)') ""
			
           !Spill Beans
           CALL RegridAerosol()
           CALL SpillBeans()

           t=0
           CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
           MixTimestep = 60
           CoagTimeStep = 60 
		
           current => particles%first
           I = 1
           TotalNumber = 0.0
           TotalVolume = 0.0
           TotalMass = 0.0
           TotalPotassiumMass = 0.0
           PotassiumIndex = FindChem("K+", 1)
           K2SO4Index = FindChem("K2SO4", 1)
           DO WHILE(associated(current))
              !current%Temperature = TEMP
              TotalNumber = TotalNumber + current%numberofparticles
              TotalVolume = TotalVolume + (4.0/3.0)*Pi*(current%effectiveradius**3)&
                   *current%numberofparticles
              IF(current%numberofparticles .GT. 0.0) THEN
                 TotalMass = TotalMass + ParticleMass(current)*current%numberofparticles
                 TotalPotassiumMass = TotalPotassiumMass + &
                      (current%aqchems(PotassiumIndex) + 2*current%aqchems(K2SO4Index)) &
                      *AqMolecularMass(PotassiumIndex)*current%numberofparticles
              END IF
              I = I + 1
              current => current%next
           END DO
	
           !WRITE(*,*) "Total Initial Number Conc: ", TotalNumber
           !WRITE(*,*) "Total Initial Volume Conc: ", TotalVolume
           !WRITE(*,*) "Total Initial Particle Mass: ", TotalMass
           !WRITE(*,*) "Total Initial Potassium Mass: ", TotalPotassiumMass
           
           !Run the main loop with above timesteps
           DO r = 1, INT(10*60*60/MixTimestep)
              
              IF (MOD(t,60.0)<=0.05) THEN !Write every 10 minutes
                 !Write current time in minutes
                 WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
                 
                 CALL SortAerosolAtGridPointForCoagulation ()
				
                 current => particles%first
                 DO WHILE(associated(current))
                    WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
                    WRITE(OutFH(4),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius/micron
                    current => current%next
                 END DO
				
                 !Advance to next line
                 WRITE(OutFH(4), FMT='(A1)') ""
              END IF
			
              current => particles%first
              I = 1
              TotalNumber = 0.0
              TotalVolume = 0.0
              TotalPotassiumMass = 0.0
              K2SO4mass = 0.0 
              Kmass = 0.0
              DO WHILE(associated(current))
                 IF (current%numberofparticles .GT. 0.0) THEN
                    TotalNumber = TotalNumber + current%numberofparticles
                    TotalVolume = TotalVolume + (4.0/3.0)*Pi*(current%effectiveradius**3) &
                         *current%numberofparticles
                    !WRITE(*,*) I, current%radius/micron, current%aqchems(K2SO4Index)
                    K2SO4mass = K2SO4mass + 2*current%aqchems(K2SO4Index)* &
                         AqMolecularMass(PotassiumIndex)*current%numberofparticles
                    Kmass = Kmass + current%aqchems(PotassiumIndex)* &
                         AqMolecularMass(PotassiumIndex)*current%numberofparticles
                    TotalPotassiumMass = TotalPotassiumMass + &
                         (current%aqchems(PotassiumIndex) + 2*current%aqchems(K2SO4Index)) &
                         *AqMolecularMass(PotassiumIndex)*current%numberofparticles
                 END IF
                 I = I + 1
                 current => current%next
              END DO

              !IF(r .EQ. 2) THEN
                 !WRITE(*,*) "Total Number Conc: ", TotalNumber
                 !WRITE(*,*) "Total Volume Conc: ", TotalVolume
                 !WRITE(*,*) "K2SO4 potassium mass", K2SO4Mass
                 !WRITE(*,*) "K+ mass", Kmass
                 !WRITE(*,*) "Total Potassium Mass: ", TotalPotassiumMass
              !END IF

              !Coagulate
              CALL StepSectionalCoagulationJacobson(CoagTimeStep)
			
              !CALL SpillBeans()
              !Update time
              t = t + MixTimestep

              IF (MOD(t,3600.0)<=0.05) THEN
                 CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
              END IF

              IF (MOD(t,60.0)<=0.05) THEN
                 !WRITE(*,*) trim(real2str(t/60)), " minutes done"
              END IF

           END DO

           current => particles%first
           I = 1
           TotalNumber = 0.0
           TotalVolume = 0.0
           TotalPotassiumMass = 0.0
           TotalMass = 0.0

           DO WHILE(associated(current))
              IF(current%numberofparticles .GT. 0.0) THEN
                 TotalNumber = TotalNumber + current%numberofparticles
                 TotalVolume = TotalVolume + (4.0/3.0)*Pi*(current%effectiveradius**3) &
                      *current%numberofparticles
                 TotalPotassiumMass = TotalPotassiumMass + &
                      (current%aqchems(PotassiumIndex) + 2*current%aqchems(K2SO4Index)) &
                      *AqMolecularMass(PotassiumIndex)*current%numberofparticles
                 TotalMass = TotalMass + ParticleMass(current)*current%numberofparticles
              END IF

              I = I + 1
              current => current%next
           END DO

           !Spill Beans
           !CALL SpillBeans()
           !WRITE(*,*) "Total Final Number Conc: ", TotalNumber
           !WRITE(*,*) "Total Final Volume Conc: ", TotalVolume
           !WRITE(*,*) "Total Final Particle Mass: ", TotalMass
           !WRITE(*,*) "Total Final Potassium Mass: ", TotalPotassiumMass

           DEALLOCATE(NumConc, OutFH, STAT=status)

	ELSE IF (CondTest) THEN !Condensation (flux) TESTS
           CALL TRANSCRIPT("")
           CALL TRANSCRIPT("Condensation Test Selected.")
                
           ! INITIALIZE the CHEMISTRY 
           CALL SetChemistryParams
           CALL ReadPhotolysisRates
		
           CALL RecalculateReactionRates(IFGAS,IFAQ)
			
           !Calculate Initial Photolysis Rates
           ut = 8.75 !Timbavati Initial Universal Time in hours
           !11.00 !Otavi Initial Universal Time in hours
           CALL RecalculatePhotoRates (ut, .TRUE.)

           ALLOCATE(Chem(HowManyGasChems), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for Chem in " &
                   //"Condensation Test of DevelopmentDriver.f90.")
           ENDIF
           !! INITIALIZE the DISSOLUTION routine
           CALL SetAllDissolution

           !! INITIALIZE the PARTICLES	
           CALL InitializeParticlesSectional
		
           !!WARNING Force initial water equilibrium between gas and aerosol
           !!(Uses open system for water eq.)
           CALL EquilibrateInternallyAtGridPoint(EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)
           CALL RegridAerosol ()
           CALL SpillBeans()

           TEMP = GetTemp()
           Mgas = GetM   ()
		
           !PressEq = SaturationVaporPressure (TEMP, 1)
           !WRITE(*,*) "Press Eq: ", PressEq
           !WRITE(*,*) "Diff: " , MolecularDiffusionCoefficient(,INT(OrganicDissolutionData(1,1)))
           !STOP

           !Find chemical indices for later
           NH3Gas = FindChem("NH3", 0)
           HNO3Gas = FindChem("HNO3", 0)
           HClGas = FindChem("HCl", 0)
           H2SO4Gas = FindChem("H2SO4", 0)

           !Count # of bins, # of Aq. Chems, and equilibrate water
           current => particles%first
           NumBins = 0
           DO WHILE(associated(current))
              NumBins = NumBins + 1
              NumAqChems = size(current%AqChems)
              NumOrgChems = size(current%OrgChems)
              NumAqOrgChems = size(current%AqOrgChems)
              CALL EqAllWater(Current, ReturnType, .TRUE.) !Open system!
              current => current%next
           END DO
           ALLOCATE(NumConc(NumBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for NumConc in " &
                   //"Coagulation Test of DevelopmentDriver.f90.")
           ENDIF

           !Set heterogeneous Reactions
           HeteroBins = NumBins
           ALLOCATE(HeteroNumberConc(HeteroBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for HeteroNumberConc in " &
                   //"Coagulation Test of DevelopmentDriver.f90.")
           ENDIF

           ALLOCATE(HeteroRadii(HeteroBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for HeteroRadii in " &
                   //"Coagulation Test of DevelopmentDriver.f90.")
           ENDIF

           !Get particle number conc. and radii for hetero. rate calculation
           current => particles%first
           I = 1
           DO WHILE(associated(current))
              HeteroNumberConc(I) = current%NumberofParticles
              HeteroRadii(I) = current%EffectiveRadius/100.0 !convert from cm to m
              !ResetHeteroRates expects radii in units of m
              current => current%next
              I = I+1
           END DO
				
           !Calculate Heterogeneous rates from size dist. info
           CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
           CALL RecalculatePhotoRates(ut, .TRUE.)
					
           TEMP = GetTemp()
           Mgas = GetM   ()

           ALLOCATE(OutFH(3+NumBins), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for OutFH in " &
                   //"Coagulation Test of DevelopmentDriver.f90.")
           ENDIF
           DO I = 1, 3+Numbins
              OutFH(I) = GetFileHandle()
           END DO
	
           OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
           OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
           OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")

           CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),0.0,.TRUE.,.FALSE.)
           
           current => particles%first
           DO I = 4, 3 + NumBins
              OPEN(UNIT = OutFH(I), FILE=TRIM(OutputDeckSubDir)//"Aerosol"//TRIM(INT2STR(I-3))//".txt", &
                   STATUS="REPLACE")
              current => particles%first
              CALL AerosolOutput(current, OutFH(I), 0.0, .TRUE., .FALSE.)
              IF (associated(current%next)) THEN
                 current => current%next
              ELSE
                 CYCLE
              END IF
           END DO

           t=0	
           CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),t,.FALSE.,.FALSE.)
	
           current => particles%first
           DO I = 4, 3+NumBins
              CALL AerosolOutput(current, OutFH(I), t, .FALSE., .FALSE.)
              IF (associated(current%next)) THEN
                 current => current%next
              ELSE
                 CYCLE
              END IF
              
           END DO
		
           !Output size distribution
           CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
			
           !WRITE(*,*) "Initial NH3 : ", trim(real2str(GridGasChem(NH3Gas)*1e9/Mgas)), " ppb"
           !WRITE(*,*) "Initial HNO3 : ", trim(real2str(GridGasChem(HNO3Gas)*1e9/Mgas)), " ppb"
           !WRITE(*,*) "Initial HCl : ", trim(real2str(GridGasChem(HClGas)*1e9/Mgas)), " ppb"
           !WRITE(*,*) "Initial H2SO4 : ", trim(real2str(GridGasChem(H2SO4Gas)*1e9/Mgas)), " ppb"
				
           !Run the main loop with above timesteps
           TimeStepSize = 60.0
           DO r = 1, INT(30*60/TimeStepSize)
              
              !Call Condensation
              !WRITE(*,*) "Before StepCondensationAll"
              !CALL FindAqElectrolyteEquilibriumForGridPoint ( UpdateThermo = .TRUE.) 
              CALL StepCondensationAll (TimeStepSize,  1)
	
              !PressEq = SaturationVaporPressure (TEMP, 1)
              !P1 = GetGridCellChemBurden ( INT(OrganicDissolutionData(1,1)))*RstarMB*Temp/Avogadro
              !WRITE(*,*) "Press Diff (mbar): ", P1-PressEq
              !WRITE(*,*) "Molecular Weight: ", GasMolecularMass(INT(OrganicDissolutionData(1,1)))
              !WRITE(*,*) "Diff: " , MolecularDiffusionCoefficient(,INT(OrganicDissolutionData(1,1)))
              !STOP
			
			
              !Decrease Temperature
              !TEMP = GetTemp()
              !Temp = Temp -0.01*TimeStepSize
              !CALL SetTempField ( Temp)
              !WRITE(*,*) "Temperature: ", GetTemp()
              !WRITE(*,*) "RH :", GetRelativeHumidity ()
              !Mgas = GetM ()

              !WRITE(*,*) "Before regrid"
			
              !Move particles across bin boundaries, if necessary!
              !(Also updates temperature)
              CALL RegridAerosol ()

              !Update time
              t = t + TimeStepSize

              !WRITE(*,*) "Before output"
					
              !Write Output to file every minute 
              IF (MOD(t,60.0) .LT. 0.01 .OR. DissolutionEquilibriumorFlux .EQV. .FALSE.) THEN
                 
                 CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),t,.FALSE.,.FALSE.)
			
                 current => particles%first
                 DO I = 4, 3+NumBins
                    CALL AerosolOutput(current, OutFH(I), t, .FALSE., .FALSE.)
                    
                    IF (associated(current%next)) THEN
                       current => current%next
                    ELSE
                       CYCLE
                    END IF

                 END DO
              END IF !Write Output
			
              !WRITE(*,*) "After output"
			
              IF (MOD(t,60.0)<=0.05 .OR. .NOT.(DissolutionEquilibriumorFlux)) THEN
                 IF(DissolutionEquilibriumorFlux) THEN
                    CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
                 ELSE
                    CALL SizeDistOutput(t,"SizeDist2.txt")
                 END IF
              END IF

              !WRITE(*,*) "After size dist"
			
              !IF (MOD(t,60.0)<=0.01) THEN
              WRITE(*,*) trim(real2str(t/60)), " minutes done"
              !END IF
              !If just doing equilibrium, exit here
              IF(.NOT.(DissolutionEquilibriumOrFlux)) EXIT
		
           END DO

           CALL SpillBeans()
           WRITE(*,*) "Final NH3 : ", trim(real2str(GridGasChem(NH3Gas)*1e9/Mgas)), " ppb"
           WRITE(*,*) "Final HNO3 : ", trim(real2str(GridGasChem(HNO3Gas)*1e9/Mgas)), " ppb"
           WRITE(*,*) "Final HCl : ", trim(real2str(GridGasChem(HClGas)*1e9/Mgas)), " ppb"
           WRITE(*,*) "Final H2SO4 : ", trim(real2str(GridGasChem(H2SO4Gas)*1e9/Mgas)), " ppb"
           WRITE(*,*) "Final POA1 : ", trim(real2str(GridGasChem(INT( &
                                       OrganicDissolutionData(1,1)))*1e9/Mgas)), " ppb"
	
           DEALLOCATE(Chem, HeteroNumberConc, HeteroRadii, OutFH, STAT=status)

	ELSE IF (SmogChamberModel) THEN !Chamber TESTS
           CALL TRANSCRIPT("")
           CALL TRANSCRIPT("Smog Chamber Model Selected.")

           ! INITIALIZE the CHEMISTRY 
           CALL SetChemistryParams
           CALL ReadPhotolysisRates
		
           CALL RecalculateReactionRates(IFGAS,IFAQ)
           !Recalculate Photolysis Rates
	   CALL RecalculatePhotoRates(ut, .TRUE.)             
           !Recalculate Heterogeneous Rates
           CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .TRUE.)
				
           ALLOCATE(Chem(HowManyGasChems), STAT = status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for Chem in " &
                   //"Smog Chamber Model of DevelopmentDriver.f90.")
           ENDIF
           !! INITIALIZE the DISSOLUTION routine
           WRITE(*,*) IncludeAerosols
           IF (IncludeAerosols) THEN 
              CALL SetAllDissolution
           !WRITE(*,*) "Dissolution Okay"
           !STOP
           !! INITIALIZE the PARTICLES
              CALL InitializeParticlesSectional	

           !WRITE(*,*) "Check 1"
           !!WARNING Force initial equilibrium between gas and aerosol
           !!(Uses open system for water eq.)
              CALL EquilibrateGridPoint(EquilibrateWater=.FALSE., &                   
                WaterOpenSystem=.FALSE., InorgOnly=.FALSE.)
           !WRITE(*,*) "Check 2"
           
              CALL SpillBeans()
		
           !Count # of bins, # of Aq. Chems

              current => particles%first
              NumBins = 0
              DO WHILE(associated(current))
                 NumBins = NumBins + 1
                 NumAqChems = size(current%AqChems)
                 NumOrgChems = size(current%OrgChems)
                 current => current%next
              END DO
              ALLOCATE(NumConc(NumBins), STAT = status)
              IF (status > 0) THEN
                 CALL ERROR("Allocation error for NumConc in " &
                   //"Smog Chamber Model of DevelopmentDriver.f90.")
              ENDIF
           ENDIF
	
           TEMP = GetTemp()
           Mgas = GetM   ()

           ALLOCATE(OutFH(4), STAT =status)
           IF (status > 0) THEN
              CALL ERROR("Allocation error for OutFH in " &
                   //"Smog Chamber Model of DevelopmentDriver.f90.")
           ENDIF
           OutFH(1) = GetFileHandle()
           OutFH(2) = GetFileHandle()
           OutFH(3) = GetFileHandle()
           OutFH(4) = GetFileHandle()

           OPEN(UNIT = OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOutChamber1.txt", STATUS="REPLACE")
           OPEN(UNIT = OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOutChamber2.txt", STATUS="REPLACE")
           OPEN(UNIT = OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOutChamber3.txt", STATUS="REPLACE")
           
           CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),0.0,.TRUE.,.FALSE.)
           IF (IncludeAerosols) THEN !Gas Chemistry Only
              IF (NumBins .GT. 1) CALL ERROR("Smog Chamber Model only works with bulk aerosol.") 
              OPEN(UNIT = OutFH(4), FILE=TRIM(OutputDeckSubDir)//"AerosolOutChamber.txt", STATUS="REPLACE")
              CALL AerosolOutput(current, OutFH(4), 0.0, .TRUE., .FALSE.)
           ENDIF
           t = 0.0
           ChemTimestep = ChemStep
           MixTimestep = Mixstep !Should always be an integer multiple!
           !Run the main loop with above timesteps
           DO r = 1, INT(RunLength*60/MixTimestep+1)

              !Write Output to file every hour in units of ppb

              IF (MOD(t,3600.0)<=0.05) THEN
	      !IF (MOD(t,60.0)<=0.05) THEN			
                 CALL GasOutput(OutFH(1),OutFH(2),OutFH(3),t,.FALSE.,.FALSE.)
                 IF (IncludeAerosols) THEN !Gas Chemistry Only
                    current => particles%first
                    CALL AerosolOutput(current, OutFH(4), t, .FALSE., .FALSE.)	
                 ENDIF

              END IF !Write Output
		
              !Calculate Temperature
		
              TEMP = GetTemp()

			
              CALL RecalculateReactionRates(IFGAS,IFAQ)
              !WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
              !after this step, or those rates will be set to 0!
              CALL RecalculatePhotoRates(ut, .FALSE.)             
              !Recalculate Heterogeneous Rates
              CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
		
              !Step Gas Chemistry
              !UR21i = FindChem("UR21", 0)
              WRITE(*,*) "Before Chem"
              CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep))			
              WRITE(*,*) "After Chem"
              !Only Equilibrate every 60 s (one second before the minute)
              !IF (MOD(t,CondTimeStep) .LT. 59.01 .AND. MOD(t,CondTimeStep) .GT. 58.99) THEN

              !Equilibrate Gases and particles
              !CALL EquilibrateGridPointAll( .TRUE., .TRUE., .TRUE.)
              !WRITE(*,*) "Before Eq", GridGasChem(,UR21i)
              IF (IncludeAerosols) THEN !Gas Chemistry Only
                 CALL EquilibrateGridPoint(EquilibrateWater=.FALSE., &                   
                   WaterOpenSystem=.FALSE., InorgOnly=.FALSE.)	
              ENDIF
              !WRITE(*,*) "After Eq", GridGasChem(,UR21i)
              !END IF
								
              !Update time
              t = t + MixTimestep

              !IF (MOD(t,60.0)<=0.05) THEN
              WRITE(*,*) trim(real2str(t/60)), " minutes done"
              !END IF

           END DO

           !Spill Beans
           IF (IncludeAerosols) THEN !Gas Chemistry Only
              CALL SpillBeans()
           ENDIF
           DEALLOCATE(Chem, NumConc, OutFH, STAT=status)		    
	END IF !Chamber Tests

        ii = etime(etimearray)
	CALL TRANSCRIPT("")
	CALL TRANSCRIPT("Total Execution Time was "//trim(real2str(ii)))

END PROGRAM ASP
