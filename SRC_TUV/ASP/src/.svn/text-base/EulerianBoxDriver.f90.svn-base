!! 11/06  2012  Matt Alvarado Updated intitalization to make parallel to others in
!!                            DevelopmentDriver.f90
!! 07/03  2013  Matt Alvarado Major rewrite to make it work.
SUBROUTINE EulerianBox()
	!! Include modules
	USE ModelParameters
	USE InfrastructuralCode	
	USE GridPointFields
	USE Chemistry
	USE Aerosols		
	USE Condensation
	USE OutputRoutines
	USE Coagulation
	
	implicit none

	!! Define the local variables
	integer :: i, j, k, l, seedarray(2),q, r, NumBins, A, S, HeteroBins, &
	           Gas, Aq, HNO3GasIndex, NumAqChems, NH3GasIndex, PotassiumIndex, &
			   HClGasIndex, SulfateIndex, BisulfateIndex, K2SO4Index, SOA2GasIndex, NumOrgChems
	Integer :: 	KNO3i, KCli, KHSO4i, UR21i, &
				Nai, NaNO3i, NaCli, NaHSO4i, Na2SO4i, &
				NH3i, NH4i, NH4NO3i, NH4Cli, NH4HSO4i, NH42SO4i, &
				Cli, NO3i, SO4i, HSO4i, K2SO4i,LevAeroIndex, AmSulfAeroIndex, H2SO4AeroIndex, N, HNO3Index, &
				NH3Gas, HNO3Gas, HClGas, H2SO4Gas, NumAqOrgChems, ReturnType
	real*8  :: ii, jj, kk, y(200), Dum1
	real*8  :: T_mix, Tau_Dep, slope, b, ChemTimestep, MixTimestep, &
				CoagTimeStep, CondTimeStep, t, Mgas, Temp, &
				X(3), NewRH, HNO3Conc, NH3Conc, ut, &
				Potassium, HClConc, TotalNumber, TotalVolume, TotalSulfateMass, TotalPotassiumMass, SOA2Conc, &
				TotalPOA, TotalSOA, TotalK, TotalNa, TotalNH4, TotalCl, TotalNO3, TotalSO4, &
				NucRate, MolsinNucParticle, GridTotalSOA, RH, tiny, x1, x2, x3, f1, f2, f3, chg, xs, molalbin, &
				K2SO4Mass, Kmass, TotalMass, TimeStepSize, PressEq, P1, EnvMass, InitialMass, NewMass, &
				Xlength, Ylength, Height, FireCO, FluxCO, DeltaCO, spinup

	REAL*4  :: etimearray(2)
	REAL*8, ALLOCATABLE	::	XPASS(:), GAMMA(:), Chem(:), NumConc(:), &
							HeteroNumberConc(:), HeteroRadii(:), &
							InitNumConc(:), EnvNumConc(:), SourceNumConc(:), SourceMassConc(:,:)
        REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym,BackScatCoeff,SubExtCoeff, SubSSA, &
                                          ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                                          SubExtCoeffSep, SubSSASep
        INTEGER, ALLOCATABLE :: OutFH(:), OFH(:)

							
	
	LOGICAL :: DoTestCoagulation, DoHumidityTest
	
	TYPE(Particle), POINTER :: Current, Env

	!Molar emision ratios to CO for Timbavati Fire
	!See Trentman, 2005 
	real*8, ALLOCATABLE :: Gas_ER(:), Flux(:)
        real*8, ALLOCATABLE :: Aero_Num_ER(:), Aero_AqMol_ER(:,:)
        real*8, ALLOCATABLE :: Aero_OrgMol_ER(:,:), Aero_AqOrgMol_ER(:,:)
        real*8, ALLOCATABLE :: Aero_Num_Flux(:), Aero_AqMol_Flux(:,:)
        real*8, ALLOCATABLE :: Aero_OrgMol_Flux(:,:), Aero_AqOrgMol_Flux(:,:)

        spinup = 600 !spinup time in s
	Xlength = boxlength !in km
	Ylength = boxwidth !in km
	Height = inversionheight !in km
		
	FireCO = fire_co_emis !kg/s
	!Want molecules/m2/s
	FluxCO = FireCO/XLength/Ylength !kg/km2/s
	FluxCO = FluxCO*Avogadro*1000/28.0 !molecs/km2/s
	
	! INITIALIZE THE CHEMISTRY 
	CALL SetChemistryParams

        ! INITIALIZE THE REACTION RATES
	CALL RecalculateReactionRates(IFGAS,IFAQ)			
 	CALL ReadPhotolysisRates
	ut = UTCStartTime
	CALL RecalculatePhotoRates (ut, .TRUE.)

	ALLOCATE(CHEM(HowManyGasChems))
	ALLOCATE(Gas_ER(HowManyGasChems)) !molecules gas q/molecules CO
	ALLOCATE(Flux(HowManyGasChems))   !gas fluxes (molecs/km2/s)
	DO q = 1, HowManyGasChems
		Gas_ER(q) = 0.0
		Flux(q) = 0.0
	END DO
        !!This subroutine (a) calculates the ERs for all gasses 
        !!and places them in Gas_ER(:), (b) calculates the Fluxes 
        !!for all gasses and places them in Flux(:), and (c) sets the initial
        !!gas concentrations in the box to the background values
        CALL InitializeGasEulerianBox(HowManyGasChems, FluxCO, Gas_ER, Flux, Chem) 
	
	!! INITIALIZE the DISSOLUTION routine
	CALL SetAllDissolution

	!! INITIALIZE the PARTICLES	
	CALL InitializeParticlesSectional
	CALL RegridAerosol ()
	CALL SpillBeans()

	!! Call Boundary distributions
	CALL ReadBoundaryDistributionsOrganic
	CALL PopulateParticlesSectionsRightAway(.TRUE.)
	CALL SortBoundaryPart
	CALL SpillBeans(Env=.TRUE.)
		
	!!WARNING Force initial water equilibrium between gas and aerosol
	!!(Uses open system for water eq.)
	CALL EquilibrateInternallyAtGridPoint(EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)	
	CALL RegridAerosol ()
	CALL SpillBeans()

	!Count # of bins, # of Aq., Org, and AqOrg Chems
	current => particles%first
	NumBins = 0
	DO WHILE(associated(current))
		NumBins = NumBins + 1
		NumAqChems = size(current%AqChems)
		NumOrgChems = size(current%OrgChems)
		NumAqOrgChems = size(current%AqOrgChems)
		current => current%next
	END DO

	!Set heterogeneous Reactions
	HeteroBins = NumBins
	ALLOCATE(HeteroNumberConc(HeteroBins))
	ALLOCATE(HeteroRadii(HeteroBins))

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
	CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, &
                                     HeteroBins, .TRUE.)

	ALLOCATE(Aero_Num_ER(NumBins)) !particles/molec CO
	ALLOCATE(Aero_AqMol_ER(NumBins,NumAqChems)) !mol/molec CO
	ALLOCATE(Aero_OrgMol_ER(NumBins,NumOrgChems)) !mol/molec CO
	ALLOCATE(Aero_AqOrgMol_ER(NumBins,NumAqOrgChems)) !mol/molec CO

	ALLOCATE(Aero_Num_Flux(NumBins)) !particles/km2/s
        !Note below are *mass* ratios, not molar!
	ALLOCATE(Aero_AqMol_Flux(NumBins,NumAqChems)) !mol/km2/s
	ALLOCATE(Aero_OrgMol_Flux(NumBins,NumOrgChems)) !mol/km2/s
	ALLOCATE(Aero_AqOrgMol_Flux(NumBins,NumAqOrgChems)) !mol/km2/s

	DO q = 1, NumBins
		Aero_Num_ER(q) = 0.0
		Aero_Num_Flux(q) = 0.0
	    DO i = 1, NumAqChems
		Aero_AqMol_ER(q,i) = 0.0
		Aero_AqMol_Flux(q,i) = 0.0
	    END DO
	    DO i = 1, NumOrgChems
		Aero_OrgMol_ER(q,i) = 0.0
		Aero_OrgMol_Flux(q,i) = 0.0
	    END DO
	    DO i = 1, NumAqOrgChems
		Aero_AqOrgMol_ER(q,i) = 0.0
		Aero_AqOrgMol_Flux(q,i) = 0.0
	    END DO
	END DO

        !!This subroutine (a) calculates the ERs for all aerosols 
        !!(b) calculates the Fluxes for all aerosols , and (c) sets the initial
        !!areosol concentrations in the box to the background values
        CALL InitializeAeroEulerianBox(NumBins, NumAqChems, NumOrgChems, &
                                      NumAqOrgChems, FluxCO, 	         &
                                      Aero_Num_ER, Aero_Num_Flux,        &
                                      Aero_AqMol_ER, Aero_AqMol_Flux,    &
                                      Aero_OrgMol_ER, Aero_OrgMol_Flux,  &
                                      Aero_AqOrgMol_ER,                  &
                                      Aero_AqOrgMol_Flux) 
		
	TEMP = GetTemp()
	Mgas = GetM   ()

	!Find chemical indices for later
	PotassiumIndex = FindChem("K+", 1)
	KNO3i = FindChem("KNO3", 1)
	KCli = FindChem("KCl", 1)
	KHSO4i = FindChem("KHSO4", 1)
	K2SO4Index = FindChem("K2SO4", 1)
	K2SO4i = K2SO4Index
		
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

	Cli = FindChem("Cl-", 1)
	NO3i = FindChem("NO3-", 1)
	SO4i = FindChem("SO4--", 1)
	HSO4i = FindChem("HSO4-", 1)

!!!!!!!!!!!!!!OUTPUT STUFF!!!!!!!!!!!!!!!!!!!!!!!!!!!
					
		ALLOCATE(OutFH(3+NumBins))
		DO I = 1, 3+Numbins
			OutFH(I) = GetFileHandle()
		END DO
	
		OPEN(UNIT=OutFH(1), FILE=TRIM(OutputDeckSubDir)//"GasOut1.txt", STATUS="REPLACE")
		WRITE(OutFH(1), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT=OutFH(2), FILE=TRIM(OutputDeckSubDir)//"GasOut2.txt", STATUS="REPLACE")
		WRITE(OutFH(2), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		OPEN(UNIT=OutFH(3), FILE=TRIM(OutputDeckSubDir)//"GasOut3.txt", STATUS="REPLACE")
		WRITE(OutFH(3), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
	
		DO q = 1, HowManyGasChems
			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			ELSE 
				WRITE(OutFH(3),FMT='(A13)', ADVANCE='NO') GasPhaseChemicalNames(q)
			END IF
		END DO
		WRITE(OutFH(3),FMT='(A8,A5)', ADVANCE='NO') "ExtCoeff", " "
		WRITE(OutFH(3),FMT='(A8,A5)', ADVANCE='NO') "SingScat", " "
		WRITE(OutFH(3),FMT='(A5,A8)', ADVANCE='NO') "Assym", " "

		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""

		DO I = 4, 3 + NumBins
			OPEN(UNIT = OutFH(I), FILE=TRIM(OutputDeckSubDir)//"Aerosol"//TRIM(INT2STR(I-3))//".txt", STATUS="REPLACE")
			WRITE(OutFH(I), FMT='(A4,A11)', ADVANCE='NO') "Time", " "
		
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "Number", " "
			WRITE(OutFH(I),FMT='(A11, A2)', ADVANCE='NO') "Eff. Radius", " "
			WRITE(OutFH(I),FMT='(A10, A3)', ADVANCE='NO') "Ionic Str.", " "
			WRITE(OutFH(I),FMT='(A8, A5)', ADVANCE='NO') "H2O Act.", " "
			WRITE(OutFH(I),FMT='(A10, A3)', ADVANCE='NO') "Sol. Dens.", " "
			WRITE(OutFH(I),FMT='(A11, A2)', ADVANCE='NO') "Part. Dens.", " "
			WRITE(OutFH(I),FMT='(A2, A11)', ADVANCE='NO') "Ph", " "
			WRITE(OutFH(I),FMT='(A4, A9)', ADVANCE='NO') "Temp", " "
			DO J = 1, NumAqChems
				WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') AqPhaseChemicalNames(J)
			END DO
			DO J = 1, NumOrgChems		
				!WRITE(OutFH(I),FMT='(A13)', ADVANCE='NO') OrgPhaseChemicalNames(J)
			END DO
			
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotBC", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotPOA", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotSOA", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotH2O", " "
			WRITE(OutFH(I),FMT='(A4, A9)', ADVANCE='NO') "TotK", " "
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotNa", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotNH4", " "
			WRITE(OutFH(I),FMT='(A5, A8)', ADVANCE='NO') "TotCl", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotNO3", " "
			WRITE(OutFH(I),FMT='(A6, A7)', ADVANCE='NO') "TotSO4", " "
			
			WRITE(OutFH(I), FMT='(A1)') ""
		END DO

		t=0	
		!Write current time in minutes
		WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
		WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
					
		!Write gas phase concentrations in ppb
		DO q = 1, HowManyGasChems
		
			Chem(q) = GridGasChem(q) !molecules/cm3

			IF (q .LE. 70) THEN
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			ELSE 
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
			END IF
		END DO
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') 0.0
		WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') 0.0

		!Advance to next line
		WRITE(OutFH(1), FMT='(A1)') ""
		WRITE(OutFH(2), FMT='(A1)') ""
		WRITE(OutFH(3), FMT='(A1)') ""
			
		current => particles%first
		DO I = 4, 3+NumBins
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature
			!mol/particle
			DO J = 1, NumAqChems
				WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(J)
			END DO
			DO J = 1, NumOrgChems		
				!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(J)
			END DO
								
			!BC Mass Conc (ug/m3)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
			!POA Mass Conc (ug/m3)
			TotalPOA = 0.
			DO J = 2, 9
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			DO J = 20, NumOrgChems
				TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

			!SOA Mass Conc (ug/m3)
			TotalSOA = 0.
			DO J = 10, 19
				TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
				TotalSOA = TotalSOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
			END DO
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
			!Water Mass Conc. (ug/m3)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
					
			!K+ Mass Conc (ug/m3)
			TotalK = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
					TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. K2SO4Index) THEN
					TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalK = TotalK*AqMolecularMass(PotassiumIndex)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
			!Na+ Mass Conc (ug/m3)
			TotalNa = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
					TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. Na2SO4i) THEN
					TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNa = TotalNa*AqMolecularMass(Nai)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
			!NH4+ Mass Conc (ug/m3)
			TotalNH4 = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
					TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
				ELSE IF (J .EQ. NH42SO4i) THEN
					TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

			!Cl- Mass Conc (ug/m3)
			TotalCl = 0.		
			DO J = 1, NumAqChems
				IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
						TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalCl = TotalCl*AqMolecularMass(Cli)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

			!NO3- Mass Conc (ug/m3)
			TotalNO3 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
					TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

			!SO4- Mass Conc (ug/m3)
			TotalSO4 = 0.
			DO J = 1, NumAqChems
				IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
				    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
					TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
				END IF
			END DO
			!Convert from umol/m3 to ug/m3
			TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
			WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

			
			WRITE(OutFH(I), FMT='(A1)') ""
					
			IF (associated(current%next)) THEN
				current => current%next
			ELSE
				CYCLE
			END IF

		END DO
		
		!Output size distribution
		CALL SizeDistOutput(t,"SizeDist"//TRIM(INT2STR(INT(t/60)))//"min.txt")
!!!!!END OUTPUT STUFF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			
		ChemTimestep = 1
		MixTimestep = 1 !Always start as 1s for first 600 seconds, MJA 12-28-2012		
 
                CoagTimeStep = CoagStep
		CondTimeStep = CondStep

		DO r = 1, INT(600+((RunLength-10)*60/MixStep))+1

			
			!Calculate Temperature
			TEMP = GetTemp()
			CALL RecalculateReactionRates(IFGAS,IFAQ)
			!WARNING: You MUST call RecalculatePhotoRates and RecalculateHeteroRates
			!after this step, or those rates will be set to 0!
			
			!Recalculate Photolysis Rates
			CALL RecalculatePhotoRates(ut, .FALSE.)
			!WRITE(*,*) "Photo Okay"

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
			CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HeteroBins, .FALSE.)
			!WRITE(*,*) "Hetero Okay"
						
			!Step Gas Chemistry
			CALL StepGasChemistry(ChemTimestep,INT(MixTimestep/ChemTimestep))
			!WRITE(*,*) "StepGasChem Okay"
			
			!Step Condensation
			IF (MOD(t,CondTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN
				!Step Condensation, if desired
				IF(.TRUE.) THEN
					CALL RegridAerosol ()
					CALL FindAqElectrolyteEquilibriumForGridPoint (UpdateThermo = .TRUE.) 
					CALL StepCondensationAll (CondTimeStep, 1)
					!Note above calls RegridAerosol automatically
					!WRITE(*,*) t, "Cond Okay"
				END IF
			END IF

			!Step Coagulation
			IF (MOD(t,CoagTimeStep) .LT. 0.01 .AND. t .NE. 0.0) THEN

				!Step Coagulation, if desired
				DoTestCoagulation = .TRUE.
				IF (DoTestCoagulation) THEN	
					CALL StepSectionalCoagulationJacobson(CoagTimeStep)
					!Note that above already calls RegridAerosol
					!WRITE(*,*) t, "Coag Okay"
				END IF
					
			END IF

		
			!Step Dilution of Chems
			!Flux should have units of molecules/m2/s
                        CALL StepDiluGasEulerian(HowManyGasChems,WindSpeed,& 
                                         Xlength,Height,MixTimestep,Flux, &
                                         t, spinup)
		

		      !Step Dilution and Deposition of Particles
                      CALL StepDiluAeroEulerian(WindSpeed,Xlength,Height,Temp,&
                                                  MixTimestep, t, spinup, &
                                                  Aero_Num_Flux, &
                                                  Aero_AqMol_Flux, &
                                                  Aero_OrgMol_Flux, &
                                                  Aero_AqOrgMol_Flux, &
                                                  NumAqChems, NumOrgChems, &
                                                  NumAqOrgChems)


			!WRITE(*,*) "Before AOD"
			!WRITE(*,*) "AOD: ", AerosolOpticalDepth ()	
			!WRITE(*,*) "Dilution Okay"	
			!CALL SpillBeans()
			CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSSA, 0)
			
			!Write Output to file every minute 
			IF (MOD(t,10.0) .LT. 0.01) THEN
!!!!OUTPUT STUFF!!!!!				
				!Write current time in minutes
				WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') (t)/60
					
				!Write gas phase concentrations in ppb
				DO q = 1, HowManyGasChems
		
					Chem(q) = GridGasChem(q) !molecules/cm3

					IF (q .LE. 70) THEN
						WRITE(OutFH(1),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE IF (q .GT. 70 .AND. q .LE. 140) THEN
						WRITE(OutFH(2),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					ELSE 
						WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Chem(q)*1e9/Mgas
					END IF
				END DO
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') ExtCoeff(1)
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') SingScat(1)
				WRITE(OutFH(3),FMT='(7ES13.5E2)', ADVANCE='NO') Assym(1)

				!Advance to next line
				WRITE(OutFH(1), FMT='(A1)') ""
				WRITE(OutFH(2), FMT='(A1)') ""
				WRITE(OutFH(3), FMT='(A1)') ""
			
				GridTotalSOA = 0.0
				current => particles%first
				DO I = 4, 3+NumBins
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') t/60
				
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%numberofparticles
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%effectiveradius
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%IonicStr
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%WaterActivity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%SolutionDensity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%ParticleDensity
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') AerosolPh(current)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Temperature	
					!mol/particle
					DO J = 1, NumAqChems
						WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%AqChems(J)
					END DO
					DO J = 1, NumOrgChems		
						!WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(J)
					END DO
					
					!BC Mass Conc (ug/m3)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Orgchems(1)*current%numberofparticles*1.0e12*OrgMolecularMass(1)
					
					!POA Mass Conc (ug/m3)
					TotalPOA = 0.
					DO J = 2, 9
						TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
						TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
					END DO
					DO J = 20, NumOrgChems
						TotalPOA = TotalPOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
						TotalPOA = TotalPOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
					END DO
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalPOA

					!SOA Mass Conc (ug/m3)
					TotalSOA = 0.
					DO J = 10, 19
						TotalSOA = TotalSOA + current%Orgchems(J)*current%numberofparticles*1.0e12*OrgMolecularMass(J)
						TotalSOA = TotalSOA + current%AqOrgchems(J-1)*current%numberofparticles*1.0e12*AqOrgMolecularMass(J-1)
					END DO
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSOA
	
					!Water Mass Conc. (ug/m3)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') current%Aqchems(1)*current%numberofparticles*1.0e12*AqMolecularMass(1)
			
					!K+ Mass Conc (ug/m3)
					TotalK = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. PotassiumIndex .OR. J .EQ. KNO3i .OR. J .EQ. KCli .OR. J .EQ. KHSO4i) THEN
							TotalK = TotalK + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. K2SO4Index) THEN
							TotalK = TotalK + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalK = TotalK*AqMolecularMass(PotassiumIndex)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalK
	
					!Na+ Mass Conc (ug/m3)
					TotalNa = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. Nai .OR. J .EQ. NaNO3i .OR. J .EQ. NaCli .OR. J .EQ. NaHSO4i) THEN
							TotalNa = TotalNa + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. Na2SO4i) THEN
							TotalNa = TotalNa + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNa = TotalNa*AqMolecularMass(Nai)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNa
	
					!NH4+ Mass Conc (ug/m3)
					TotalNH4 = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. NH3i .OR. J .EQ. NH4i .OR. J .EQ. NH4NO3i .OR. J .EQ. NH4Cli .OR. J .EQ. NH4HSO4i) THEN
							TotalNH4 = TotalNH4 + current%Aqchems(J)*current%numberofparticles*1.0e12
						ELSE IF (J .EQ. NH42SO4i) THEN
							TotalNH4 = TotalNH4 + 2*current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNH4 = TotalNH4*AqMolecularMass(NH4i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNH4

					!Cl- Mass Conc (ug/m3)
					TotalCl = 0.		
					DO J = 1, NumAqChems
						IF (J .EQ. Cli .OR. J .EQ. NH4Cli .OR. J .EQ. KCli .OR. J .EQ. NaCli) THEN
							TotalCl = TotalCl + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalCl = TotalCl*AqMolecularMass(Cli)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalCl

					!NO3- Mass Conc (ug/m3)
					TotalNO3 = 0.
					DO J = 1, NumAqChems
						IF (J .EQ. NO3i .OR. J .EQ. NH4NO3i .OR. J .EQ. KNO3i .OR. J .EQ. NaNO3i) THEN
							TotalNO3 = TotalNO3 + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalNO3 = TotalNO3*AqMolecularMass(NO3i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalNO3

					!SO4- Mass Conc (ug/m3)
					TotalSO4 = 0.
					DO J = 1, NumAqChems
						IF (J .EQ. SO4i .OR. J .EQ. NH42SO4i .OR. J .EQ. K2SO4i .OR. J .EQ. Na2SO4i &
						    .OR. J .EQ. HSO4i .OR. J .EQ. NH4HSO4i .OR. J .EQ. KHSO4i .OR. J .EQ. NaHSO4i) THEN
							TotalSO4 = TotalSO4 + current%Aqchems(J)*current%numberofparticles*1.0e12
						END IF
					END DO
					!Convert from umol/m3 to ug/m3
					TotalSO4 = TotalSO4*AqMolecularMass(SO4i)
					WRITE(OutFH(I),FMT='(7ES13.5E2)', ADVANCE='NO') TotalSO4

					WRITE(OutFH(I), FMT='(A1)') ""

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
!!!!END OUTPUT STUFF!!!!!	
			
			!Update time
			t = t + MixTimestep
			ut = ut + MixTimeStep/3600.

			IF (MOD(t,CoagTimeStep)<=0.01) THEN
				WRITE(*,*) trim(real2str(t/60)), " minutes done"
			END IF

		END DO


END SUBROUTINE EulerianBox

SUBROUTINE InitializeGasEulerianBox(NumGasChems, FluxCO, Gas_ER, Flux, Chem) 

	!! Include modules
	USE GridPointFields,     ONLY : GridGasChem, &
                                        UpdateChemicalConcentrations
	USE Chemistry,           ONLY : FindChem, &
                                        GasPhaseBackground, &
                                        HowManyEvolveGasChems	
	implicit none

        !External variables
        INTEGER, INTENT(IN)   :: NumGasChems 
        REAL*8, INTENT(IN)    :: FluxCO
        REAL*8, INTENT(INOUT) ::  Gas_ER(:), Flux(:), Chem(:)

        !Internal variables
        INTEGER               :: q, IndCO
        REAL                  :: GasConc

        !Find Index for CO
	IndCO = FindChem("CO", 0)

	DO q = 1, NumGasChems  
             GasConc = GridGasChem(q)
             IF (GasConc .GT. 0.0) THEN
                 !Calculate ERs to CO
                 Gas_ER(q) = (GasConc-GasPhaseBackground(q))/ &
                         (GridGasChem(IndCO)-GasPhaseBackground(IndCO))
                 !Calculate Flux (molecules/km2/s)
                  Flux(q) = Gas_ER(q)*FluxCO
             ENDIF
	     !Set Inital Gas Concentrations to Background
             Chem(q) = GasPhaseBackground(q)
        ENDDO

	!Update chemical concentrations
	CALL UpdateChemicalConcentrations(Chem(1:HowManyEvolveGasChems))

END SUBROUTINE

SUBROUTINE InitializeAeroEulerianBox(NumBins, NumAqChems, NumOrgChems, &
                                      NumAqOrgChems, FluxCO, 	        &
                                      Aero_Num_ER, Aero_Num_Flux, &
                                      Aero_AqMol_ER, Aero_AqMol_Flux, &
                                      Aero_OrgMol_ER, Aero_OrgMol_Flux, &
                                      Aero_AqOrgMol_ER, &
                                      Aero_AqOrgMol_Flux) 
	!! Include modules
	USE GridPointFields,     ONLY : GridGasChem
	USE Chemistry,           ONLY : FindChem, &
                                        GasPhaseBackground	
	USE Aerosols,            ONLY : particle, particles, &
                                        Boundaryparticles, &
                                        SortAerosolAtGridPointForCoagulation, &
                                        RecalculateRadius
	implicit none

        !External variables
        INTEGER, INTENT(IN)     :: NumBins, NumAqChems
        INTEGER, INTENT(IN)     :: NumOrgChems, NumAqOrgChems
        REAL*8,  INTENT(IN)     :: FluxCO
        REAL*8,  INTENT(INOUT)  :: Aero_Num_ER(:), Aero_Num_Flux(:)
        REAL*8,  INTENT(INOUT)  :: Aero_AqMol_ER(:,:), Aero_AqMol_Flux(:,:)
        REAL*8,  INTENT(INOUT)  :: Aero_OrgMol_ER(:,:), Aero_OrgMol_Flux(:,:)
        REAL*8,  INTENT(INOUT)  :: Aero_AqOrgMol_ER(:,:), Aero_AqOrgMol_Flux(:,:) 

        !Internal variables
        INTEGER                 :: q, j, IndCO
        REAL*8                  :: molecCO !Excess CO in molec/cm3
	TYPE(Particle), POINTER :: Current, Env

	CALL SortAerosolAtGridPointForCoagulation () 
        !Find Index for CO
	IndCO = FindChem("CO", 0)
        molecCO = GridGasChem(IndCO)-GasPhaseBackground(IndCO)

	current => particles%first
        Env => BoundaryParticles%first

	DO q = 1, NumBins           

           IF (current%Numberofparticles .GT. 0.0) THEN
               !Calculate ERs to CO (particles/molecules CO)
               Aero_Num_ER(q)= (current%Numberofparticles &
                                    -Env%Numberofparticles) &
                                 /molecCO                              
               !Calculate Flux (particles/km2/s)
               Aero_Num_Flux(q) = Aero_Num_ER(q)*FluxCO
               
               DO j = 1,NumAqChems
                  !Calculate ERs to CO (mols/molecules CO)
                  Aero_AqMol_ER(q,j)=(current%Numberofparticles*current%AqChems(j) &
                                       -Env%Numberofparticles*Env%AqChems(j)) &
                                       /molecCO
                  !Calculate Flux (mols/km2/s)
                  Aero_AqMol_Flux(q,j) = Aero_AqMol_ER(q,j)*FluxCO
	          !Set Initial Mass Concentrations to Background
                  current%AqChems(j)=Env%AqChems(j)
               ENDDO 

               DO j = 1,NumOrgChems
                  !Calculate ERs to CO (mols/molecules CO)
                  Aero_OrgMol_ER(q,j)=(current%Numberofparticles*current%OrgChems(j) &
                                       -Env%Numberofparticles*Env%OrgChems(j)) &
                                       /molecCO
                  !Calculate Flux (mols/km2/s)
                  Aero_OrgMol_Flux(q,j) = Aero_OrgMol_ER(q,j)*FluxCO
	          !Set Initial Mass Concentrations to Background
                  current%OrgChems(j)=Env%OrgChems(j)
               ENDDO 

               DO j = 1,NumAqOrgChems
                  !Calculate ERs to CO (mols/molecules CO)
                  Aero_AqOrgMol_ER(q,j)=(current%Numberofparticles*current%AqOrgChems(j) &
                                       -Env%Numberofparticles*Env%AqOrgChems(j)) &
                                       /molecCO
                  !Calculate Flux (mols/km2/s)
                  Aero_AqOrgMol_Flux(q,j) = Aero_AqOrgMol_ER(q,j)*FluxCO
	          !Set Initial Mass Concentrations to Background
                  current%AqOrgChems(j)=Env%AqOrgChems(j)
               ENDDO 
           ENDIF

	   !Set Initial Number Concentrations to Background
           current%Numberofparticles = Env%Numberofparticles
	   CALL RecalculateRadius(current)
           !Iterate to next particle
           IF(q .LT. NumBins) THEN
	      current => current%next
              Env => Env%next
           ENDIF
	END DO 
END SUBROUTINE InitializeAeroEulerianBox

SUBROUTINE StepDiluGasEulerian(HowManyGasChems,WindSpeed,& 
                               Xlength,Height, MixTimestep,Flux, t, spinup)

	!! Include modules
	USE GridPointFields,     ONLY : GridGasChem, &
                                        UpdateChemicalConcentrations
	USE Chemistry,           ONLY : FindChem, &
                                        GasPhaseBackground, &
                                        HowManyEvolveGasChems, &
                                        gasphasedepvel	

	implicit none

        !External variables
        INTEGER, INTENT(IN) :: HowManyGasChems
        REAL*8,  INTENT(IN) :: WindSpeed !m/s
        REAL*8,  INTENT(IN) :: Xlength, Height !km
        REAL*8,  INTENT(IN) :: MixTimeStep, t, spinup !s
        REAL*8,  INTENT(IN) :: Flux(:) !molecules/km2/s  

        !Internal Variables
        REAL*8, ALLOCATABLE :: Chem(:)
        REAL*8              :: xlen, hi !m
        INTEGER             :: q

        xlen =1000.0*XLength
        hi = 1000.0*Height

        ALLOCATE(Chem(HowManyGasChems))

	DO q = 1, HowManyGasChems
		Chem(q) = GridGasChem(q)
		IF(t .LT. spinup) THEN
                    Chem(q) = (Chem(q) + (WindSpeed*GasPhaseBackground(q)*MixTimestep/xlen)) &
				    /(1+ (MixTimestep*GasPhaseDepVel(q)/(hi*1.0e2)) + (MixTimestep*Windspeed/xlen))
		
		ELSE
					Chem(q) = (Chem(q) + (WindSpeed*GasPhaseBackground(q)*MixTimestep/xlen) + (Flux(q)*MixTimestep/hi/1.0e12)) &
				          /(1+ (MixTimestep*GasPhaseDepVel(q)/(hi*1.0e2)) + (MixTimestep*Windspeed/xlen))
		END IF
	END DO	

	!Update chemical concentrations
	CALL UpdateChemicalConcentrations(Chem(1:HowManyEvolveGasChems))


END SUBROUTINE StepDiluGasEulerian

SUBROUTINE StepDiluAeroEulerian(WindSpeed,Xlength,Height,Temp,&
                                                  MixTimestep, t, spinup, &
                                                  Aero_Num_Flux, &
                                                  Aero_AqMol_Flux, &
                                                  Aero_OrgMol_Flux, &
                                                  Aero_AqOrgMol_Flux, &
                                NumAqChems, NumOrgChems, NumAqOrgChems)
	!! Include modules

	USE Aerosols,            ONLY : particle, particles, &
                                        Boundaryparticles, &
                                        TerminalVelocity, &
                                        SortAerosolAtGridPointForCoagulation, &
                                        ShellRefIndAndRad, &
                                        RecalculateRadius
	implicit none

        !External variables
        REAL*8,  INTENT(IN) :: WindSpeed !m/s
        REAL*8,  INTENT(IN) :: Xlength, Height !km
        REAL*8,  INTENT(IN) :: Temp !K
        REAL*8,  INTENT(IN) :: MixTimeStep, t, spinup !s
        REAL*8,  INTENT(IN) :: Aero_Num_Flux(:) !particles/km2/s  
        REAL*8,  INTENT(IN) :: Aero_AqMol_Flux(:,:) !mols/km2/s  
        REAL*8,  INTENT(IN) :: Aero_OrgMol_Flux(:,:) !mols/km2/s  
        REAL*8,  INTENT(IN) :: Aero_AqOrgMol_Flux(:,:) !mols/km2/s 
        INTEGER, INTENT(IN) :: NumAqChems, NumOrgChems, NumAqOrgChems 

        !Internal Variables
        REAL*8              :: xlen, hi !m
        REAL*8              :: InitNumConc, EnvNumConc
        REAL*8              :: InitialMol, EnvMol, NewMol
        INTEGER             :: I, J
	TYPE(Particle), POINTER :: Current, Env

        xlen =1000.0*XLength
        hi = 1000.0*Height

	CALL SortAerosolAtGridPointForCoagulation () 
	current => particles%first
	Env => BoundaryParticles%first
	I = 1
	DO WHILE(associated(current))
 
            CALL RecalculateRadius(current)

	    !Get Number Concentration of particles in particles/cm3
	    InitNumConc = current%numberofparticles 
	    EnvNumConc = Env%numberofparticles

	    !Dilute and Dep number of particles (assume no background for now)
	    !Flux should have units of particles/km2/s
	    IF(t .LT. spinup) THEN
		current%numberofparticles = (InitNumConc + &
                (WindSpeed*EnvNumConc*MixTimestep/xlen)) &
		/(1.0 + MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                + Windspeed/xlen))
	    ELSE
		current%numberofparticles = (InitNumConc + &
                (WindSpeed*EnvNumConc*MixTimestep/xlen) &
		+ (Aero_Num_Flux(I)*MixTimestep/hi/1.0e12)) &
	        /(1.0 + MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                + Windspeed/xlen))
	    ENDIF
			
	    !Aero_AqMol_Flux in units of mol/km2/s
	    IF(current%Numberofparticles .GT. 0.) THEN
		DO J = 1, NumAqChems
		   InitialMol = InitNumConc*current%AqChems(J)
		   EnvMol = EnvNumConc*Env%AqChems(J)
		   
                   IF(t .LT. spinup) THEN
			NewMol = (InitialMol + &
                         Windspeed*EnvMol*MixTimestep/xlen)/ &
	                 (1.0+MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                         + Windspeed/Xlen))
		   ELSE
			NewMol = (InitialMol + &
                         Windspeed*EnvMol*MixTimestep/xlen &
                         + Aero_AqMol_Flux(I,J)*MixTimestep/hi/1.0e12)/ &
			(1.0+ MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                        + Windspeed/xlen))
		   END IF
			Current%AqChems(J) = NewMol/current%numberofparticles
					
		END DO

		DO J = 1, NumOrgChems
		   InitialMol = InitNumConc*current%OrgChems(J)
		   EnvMol = EnvNumConc*Env%OrgChems(J)
		   
                   IF(t .LT. spinup) THEN
			NewMol = (InitialMol + &
                         Windspeed*EnvMol*MixTimestep/xlen)/ &
	                 (1.0+MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                         + Windspeed/Xlen))
		   ELSE
			NewMol = (InitialMol + &
                         Windspeed*EnvMol*MixTimestep/xlen &
                         + Aero_OrgMol_Flux(I,J)*MixTimestep/hi/1.0e12)/ &
			(1.0+ MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                        + Windspeed/xlen))
		   END IF
			Current%OrgChems(J) = NewMol/current%numberofparticles
					
		END DO

		DO J = 1, NumAqOrgChems
		   InitialMol = InitNumConc*current%AqOrgChems(J)
		   EnvMol = EnvNumConc*Env%AqOrgChems(J)
		   
                   IF(t .LT. spinup) THEN
			NewMol = (InitialMol + &
                         Windspeed*EnvMol*MixTimestep/xlen)/ &
	                 (1.0+MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                         + Windspeed/Xlen))
		   ELSE
			NewMol = (InitialMol + &
                         Windspeed*EnvMol*MixTimestep/xlen &
                         + Aero_AqOrgMol_Flux(I,J)*MixTimestep/hi/1.0e12)/ &
			(1.0+ MixTimestep*(TerminalVelocity(Current)/(hi*100) &
                        + Windspeed/xlen))
		   END IF
			Current%AqOrgChems(J) = NewMol/current%numberofparticles
					
		END DO

		!Force particle temperature to be env. Temperature
		!WARNING: This is a kludge, since the program is calculating
		!a huge temperature for the largest particles, 
                !and I'm not sure why
		current%Temperature = TEMP
				
		!Recalculate optical parameters
		!WRITE(*,*) "Call Optical"
		CALL ShellRefIndAndRad(Current)
		!WRITE(*,*) "Optical Okay"
				
		I = I + 1
		
		current => current%next
		Env => Env%next
	     ENDIF			
	ENDDO
END SUBROUTINE StepDiluAeroEulerian
