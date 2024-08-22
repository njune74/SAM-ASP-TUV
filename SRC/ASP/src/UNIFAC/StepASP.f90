!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! StepASP.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						             !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 06/27  2007   Matt Alvarado     Added InputConcentrations
!!		                   Added OutputConcentrations
!!                                 Added InitializeASP
!! 07/16  2007   Matt Alvarado     Set any aerosol concentration below
!!					1.0e-30 mol/particle to 0.
!! 07/18  2007   Matt Alvarado     Corrected output concentrations for
!!					solid electrolytes
!! 08/31  2007   Matt Alvarado     Added aerosol size to interface
!! 10/09  2007   Matt Alvarado     Changed zero criteria for 
!!                                     aerosol concentrations in 
!!				       InputConcentrations to 
!!                                    1.0e-40 mol/particle
!! 10/12  2007   Matt Alvarado     Removed Scale from photorates: 
!!                                     now in chemtrop.F
!! 02/19  2009   Matt Alvarado     Changed DFPORT to IFPORT 
!!                                     (to use ifort compiler)
!! 08/26  2010   Matt Alvarado     Comment out IFPORT (to use pgf90 compiler)
!! 02/15  2012   Matt Alvarado     Removed Eulerian Coords
!! 05/03  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  BackScatCoeff
!! 08/16  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  SubExtCoeff and SubSSA
!! 08/17  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  SubExtCoeff and SubSSA
!! 11/08  2012   Matt Alvarado     Updated call to AerosolOptProps to include
!!                                  mixing flag
!! 05/18  2016   Matt Alvarado     Made changes to allow aerosols to interface 
!!                                  with STILT-Chem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES	                            !!
!! 1. ModelParameters			    !!
!! 2. GridPointFields			    !!
!! 3. Aerosols				    !!
!! 4. Condensation			    !!
!! 5. Time				    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	    !!	
!! 1. SUBROUTINE StepASPOnce ()					    !!
!! 2. SUBROUTINE InputConcentrations
!! 3. SUBROUTINE OutputConcentrations
!! 4. SUBROUTINE InitializeASP
!! 5. SUBROUTINE ASPInterface
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE StepASP

		IMPLICIT NONE

        ! CMB: aerosol flag is private
		PRIVATE :: isAerosolProc

        ! CMB: add last two methods
		PUBLIC ::	ASPInterface, InitializeASP, setAerosolProc, getAerosolProc

        integer :: isAerosolProc
	CONTAINS

    ! CMB add
    integer function getAerosolProc() result(aeroProc)
        aeroProc = isAerosolProc
    end function

    subroutine setAerosolProc(aeroProc)
        integer, intent(in) :: aeroProc
        isAerosolProc = aeroProc
    end subroutine
    ! end CMB add

        
	SUBROUTINE StepASPOnce (TimeStep)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine recalculates the gas-phase reactions rates,  !!
	!! then steps gas chemistry, condensation, coagulation, and     !!
	!! updates the optical properties of the aerosol.	        !!
	!! All calculations are done in gridbox 1,1,1 MJA 052307	!!
        !! Removed Eulerian Grid Points MJA 021512                      !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


		!USE IFPORT,				 ONLY : RTC, ETIME
		USE ModelParameters
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: TimeStep
		
		!! Internal Variables
		INTEGER :: I, NumBins
		REAL*8 :: TEMP
	
		TYPE(Particle), POINTER :: cur
		TYPE(Particle), POINTER :: Current, Env
		! CMB add
		integer :: doAerosols != getAerosolProc()
		doAerosols = getAerosolProc()
		! end CMB add
		
		!write(*,*) 'Calculate Temperature'
		!Calculate Temperature	
		TEMP = GetTemp()
	
		!write(*,*) 'Before StepGasChemistry calling spillbeans'
		!call spillbeans()
		!Step Gas Chemistry
                print *, 'Stepping gas chemistry'

		CALL StepGasChemistry(TimeStep,1)
		

		!if (doAerosols .gt. 0) then
            !Step Condensation, if desired

            cur => particles%first
    
            IF(.TRUE.) THEN
                !print *, 'before stepcondensationall from StepASP()'

		!write(*,*) 'Calling spillbeans'
		!call spillbeans()
                print *, 'calling stepcondensationall from StepASP()'

		CALL StepCondensationAll (TimeStep, 1)
                !Note above calls RegridAerosol automatically
                !with recalculates radius and density
                !WRITE(*,*) "Cond Okay"
            END IF
    
            !Step Coagulation, if desired
            IF (.TRUE.) THEN
		!write(*,*) 'Calling spillbeans before stepsectionalcoagulationjacobson'
		!call spillbeans()
                print *, 'calling stepsectionalcoagulationjacobson'
                CALL StepSectionalCoagulationJacobson(TimeStep)

                !Note that above already calls RegridAerosol
                !which recalculates radius and density
                !WRITE(*,*) "Coag Okay"
            END IF
    
            !Step Optical Properties Particles
                    !print *, 'calling sortaerosolatgridpointforcoagulation from stepasp'

	    write(*,*) 'calling SortAerosolAtGridPointForCoagulation...'

            CALL SortAerosolAtGridPointForCoagulation () 


            cur => particles%first
            I = 1

            DO WHILE(associated(cur))
                    
                !Force particle temperature to be env. Temperature
                !WARNING: This is a kludge, since the program is 
                            !calculating a huge temperature for the largest 
                            !particles, and I'm not sure why
                            !WRITE (*,*) 'applying a temperature kludge'
                cur%Temperature = TEMP
                 ! WRITE(*,*)    'cur%Temperature = ',TEMP   
                !Recalculate optical parameters
                !WRITE(*,*) "Call Optical, Particle #", I
                ! CMB (AER, Inc): This call significantly slows down
                !                 aerosol processing.  Comment out until
                !                 we figure everything out.
                !CALL ShellRefIndAndRad(cur)
                !WRITE(*,*) "Optical Okay"
                    
                I = I + 1
            
                cur => cur%next
                    
            END DO !END STEP OPTICAL PROPERTIES
        !endif !doAerosols
        
		RETURN
	END SUBROUTINE StepASPOnce

	SUBROUTINE InputConcentrations(Temp, Press, Dens, GasConc, &
				AeroNumConc, AeroMassConc, &
				!PhotoRates, &
				sza, Transmissivity)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine takes in Temperature, pressure,		!!
	!! gas and aerosol concentrations, and photolysis rates from	!!
	!! a 3D model and updates the appropriate concentrations in ASP !!
	!! grid box 1,1,1. MJA 052307					!!
        !! Removed Eulerian coordinates, MJA 021512                     !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
	  !USE IFPORT,				ONLY : RTC, ETIME
	  USE ModelParameters,	ONLY : HowManyBins
	  USE InfrastructuralCode	
	  USE GridPointFields
	  USE Chemistry
	  USE Aerosols		
	  USE Condensation
	  USE OutputRoutines
	  USE Coagulation

	  IMPLICIT NONE

	  !! External Variables
	  REAL*8, INTENT(IN) :: Temp  !(in K)
	  REAL*8, INTENT(IN) :: Press !(in mbar)
	  REAL*8, INTENT(IN) :: Dens  !(in kg/m3 dry air)
	  REAL*8 :: GasConc(HowManyEvolveGasChems) !(in molecules/cm3)
	  REAL*8, INTENT(IN) :: AeroNumConc(HowManyBins) !(in particles/cm3)
	  REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems) !(in ug/m3)
	  REAL*8, INTENT(IN) :: Transmissivity
	  !REAL*8, INTENT(IN) :: PhotoRates(12)

	  !! Internal Variables
	  INTEGER :: I, J, L, RxnIndex, NumAq, GG
	  REAL*8 :: HeteroNumberConc(HowManyBins)
          REAL*8 :: HeteroRadii(HowManyBins)
	  REAL*8 :: Scale

	  TYPE(Particle), POINTER :: cur
	  
	  ! CMB add
	  real*8, intent(in) :: sza
	  integer :: doAerosols
	  ! CRL add
	  double precision minimum_H2O
	  real*8,parameter :: cm3_to_m3 = 1e6 ! convert cm3 to m3
	  real*8,parameter :: umol_to_mol = 1e6 ! convert umol to mol 
	  doAerosols = getAerosolProc()
	  ! end CMB add
	  !write(*,*) 'Entered InputConcenctrations calling spillbeans...'
	  !call spillbeans()  	
	  !Set temperature
          !print *, 'setting temperature = ', Temp
	  CALL SetTempField (Temp)

	  !Set Pressure
          !print *, 'setting pressure = ', press
	  CALL SetPressField (Press)
		
	  !Update chemical concentrations
          !print *, 'updating gas concentrations'
	  DO I = 1, HowManyEvolveGasChems
             IF(GasConc(I) .LT. 0.) GasConc(I) = 0.
             !if (gasconc(i) .lt. 1.0e-29) gasconc(i) = 1.0e-30
	  END DO 
	  CALL UpdateChemicalConcentrations(GasConc(1:HowManyEvolveGasChems))
	
	  !if (doAerosols .gt. 0) then	
          !Update aerosol concentrations
              print *, 'calling sort aerosol at grid point for coag'

	    !write(*,*) 'Before SortAerosolAtGridPointForCoagulation() in InputConcentrations calling spillbeans...'
	    !call spillbeans()
	    !write(*,*) 'Calling SortAerosolAtGridPointForCoagulation() in InputConcentrations'
            CALL SortAerosolAtGridPointForCoagulation () 
	    !write(*,*) 'After SortAerosolAtGridPointForCoagulation() in InputConcentrations calling spillbeans...'
	    !call spillbeans()

            cur => particles%first !current particle
            I = 1
            
            NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
            DO WHILE(associated(cur))
                cur%numberofparticles = AeroNumConc(I)
				
		!Force minimum water
		! CRL add - if mass conc of water is less than radius of bin
                minimum_H2O = 5e-21*AeroNumConc(I)*cm3_to_m3*umol_to_mol*AqMolecularMass(1)
                IF(AeroMassConc(I,1) .LT. minimum_H2O) 	AeroMassConc(I,1) = minimum_H2O		
    
                ! CB: for debugging
		!print *, 'In InputConcentrations'
                !print *, 'Particle ID = ', cur%particleID
                !print *, 'I = ', I
                !print *, 'AeroNumConc(I) = ', AeroNumConc(I)
                !write(*,*) "cur%edges = ", cur%edges
		!write(*,*) 'In loop in InputConcentrations calling spillbeans',I
		!call spillbeans()

		!IF(cur%numberofparticles.GT. 0.) THEN
		!IF(cur%Numberofparticles .GT. 1E-40) THEN
		if(cur%numberofparticles .gt. 1.0E-6) then
                    !print *, '******** number of particles > 0 ********'
		   IF(cur%dry) THEN
		  	cur%dry = .FALSE.
		   END IF
    
 		   DO J = 1, NumAq
			IF(AeroMassConc(I,J) .LT. 0.) AeroMassConc(I,J) = 0. ! 1.0e-30

			cur%AqChems(J) = AeroMassConc(I,J)*1.0e-12 &
                          /(cur%numberofparticles*AqMolecularMass(J))

 			IF(cur%AqChems(J) .LT. 1.0e-40) cur%AqChems(J)=0.  !cur%aqChems(j) = 1.0e-30
 		   END DO
    
 		   DO J = 1, HowManyOrgChems
 			IF(AeroMassConc(I,J+NumAq) .LT. 0.) AeroMassConc(I,J+NumAq) = 0.
              		
			Cur%OrgChems(J) = AeroMassConc(I,J+NumAq)*1.0e-12 &
                            		/(cur%numberofparticles*OrgMolecularMass(J))

 			IF(Cur%OrgChems(J) .LT. 1.0e-40) Cur%OrgChems(J) = 0.
 		   END DO
                
 		   DO J = 1, HowManyAqOrgChems
 			IF(AeroMassConc(I,J+NumAq+HowManyOrgChems) .LT. 0.) &
				AeroMassConc(I,J+NumAq+HowManyOrgChems) = 0.
			
			Cur%AqOrgChems(J) = AeroMassConc(I,J+NumAq+HowManyOrgChems) &
                         	*1.0e-12/(cur%numberofparticles*AqOrgMolecularMass(J))

 			IF(Cur%AqOrgChems(J) .LT. 1.0e-40) Cur%AqOrgChems(J) = 0.
                  END DO
		else !No particles
                    !write(*,*) '**** number of particles not > 0 ********'
              	   	cur%numberofparticles = 0.

                  DO J = 1, NumAq
             		Cur%AqChems(J) = 0.
             		!Cur%AqChems(J) = 1.0e-40
                  END DO
                    
                  DO J = 1, HowManyOrgChems
             		Cur%OrgChems(J) = 0.
             		!Cur%OrgChems(J) = 1.0e-40
                  END DO
                    
                  DO J = 1, HowManyAqOrgChems
             		Cur%AqOrgChems(J) = 0.
             		!Cur%AqOrgChems(j) = 1.0e-40
                  END DO				
		end if
            	!write(*,*) 'Out of loop in InputConcentrations calling spillbeans',I
		!call spillbeans()
                !Force particle temperature to be env. Temperature
                cur%Temperature = Temp
                    
                !Equilibrate, recalculate radius and density
                ! CMB (AER, Inc): Add if condition for this block to not do a find
                !                 here if no water.
                ! Nope, causes some weird errors?
                !if (cur%aqchems(j) .gt. 1.0e-34) then
                !WRITE(*,*) "Before FindElectrolyteEquilibrium"

	    	!write(*,*) 'Before FindElectrolyteEquilibrium() in InputConcentrations calling spillbeans...',I
		!write(*,*) 'AeroNumConc = ',AeroNumConc, I
	    	!call spillbeans()
	    	!write(*,*) 'calling FindElectrolyteEquilibrium() in InputConcentrations'

                CALL FindElectrolyteEquilibrium (cur, UpdateThermo=.TRUE., &
                                                 FirstEquilibration=.TRUE., &
                                                 ReturnType=GG)

            	!write(*,*) 'Done with FindElectrolyteEquilibrium'
	    	!write(*,*) 'After FindElectrolyteEquilibrium() in InputConcentrations calling spillbeans...',I
	    	!call spillbeans()
     
            	!endif

	    	!write(*,*) 'Before RecalculateRadius(cur)in InputConcentrations calling spillbeans...'
	    	!call spillbeans()

                CALL RecalculateRadius(cur)
		
            	!WRITE(*,*) "After Recauculate Radius cur%effectiveradius:", cur%effectiveradius
	    	!write(*,*) 'After RecalculateRadius(cur)in InputConcentrations calling spillbeans...',I
	    	!call spillbeans()

            	I = I + 1
            	cur => cur%next		
            END DO! end cur => particle

	  !endif !end doaerosol
	  !write(*,*) 'End doaerosol loop calling spillbeans'
	  !call spillbeans()  	
      !Recalculate all reaction rates
      !print *, 'CMB: Calculating reaction rates: (ifgas, ifaq) = ', ifgas, ifaq
      call RecalculateReactionRates(IFGAS,IFAQ)

      ! CMB: Recalc photo rates (new)
      call RecalculatePhotoRatesSZA(sza, Transmissivity, .false.)
      
	!write(*,*) 'After RecalculatePhotoRatesSZA InputConcentrations calling spillbeans...',I
	!call spillbeans()
      ! CMB: Debugging...
!      print *, 'Solar Zenith Angle: ', sza
!      do l = 1, GasPhaseReactionRateArraySize(1)
!        print *, 'L = ', L
!        if (int(GasPhaseReactionRates(L, 1)) .eq. 1 .and. &
!                int(GasPhaseReactionRates(L,3)) .eq. 0) then
!            RxnIndex = int(GasPhaseReactionRates(l,4))
!            
!            print *, '  Photolysis Reaction'
!            select case(RxnIndex)
!            case(2)
!                print *, '    RxnIndex = 2, GasPhaseChemicalRates(L) = ', &
!                        GasPhaseChemicalRates(l)
!            case(6)
!                print *, '    RxnIndex = 6, GasPhaseChemicalRates(L) = ', &
!                        GasPhaseChemicalRates(l)
!            end select
!        endif   
!      enddo
!      stop
      ! end CMB
      		
	  !WRITE(*,*) "Before Photorates"
	  !Reset Photolysis Rates
!	  DO L = 1, GasPhaseReactionRateArraySize(1)
!						
!	    !! Only do photolysis reactions
!	    IF (INT(GasPhaseReactionRates(L,1)) .EQ. 1 .AND. &
!                INT(GasPhaseReactionRates(L,3)) .EQ. 0.) THEN	
!								
!		!Find reaction index number
!		RxnIndex = INT(GasPhaseReactionRates(L,4))
!
!		SELECT CASE (RxnIndex)
!		  CASE(4) !NO2 => NO + O
!		    GasPhaseChemicalRates(L) = PhotoRates(1)
!		  CASE(5) !NO3 => NO + O2
!		    GasPhaseChemicalRates(L) = PhotoRates(2)
!		  CASE(6) !NO3 => NO2 + O
!		    GasPhaseChemicalRates(L) = PhotoRates(3)		
!		  CASE(3) !O3 => O + O2
!	            GasPhaseChemicalRates(L) = PhotoRates(4)
!		  CASE(2) !O3 => O1D + O2
!		    GasPhaseChemicalRates(L) = PhotoRates(5)
!		  CASE(12) !HONO => Products
!		    GasPhaseChemicalRates(L) = PhotoRates(6)
!		  CASE(11) !H2O2 => 2OH
!		    GasPhaseChemicalRates(L) = PhotoRates(7)
!		  CASE(15) !HCHO => 2HO2 + CO
!		    GasPhaseChemicalRates(L) = PhotoRates(8)
!		  CASE(16) !HCHO => H2 + CO
!		    GasPhaseChemicalRates(L) = PhotoRates(9)
!                  CASE(17) !ALD2 => Products
!		    GasPhaseChemicalRates(L) = PhotoRates(10)
!		  CASE(23) !KETL => Products
!		    GasPhaseChemicalRates(L) = PhotoRates(11)
!	          CASE(22) !MGLY => Products
!		    GasPhaseChemicalRates(L) = PhotoRates(12)
!		END SELECT
!	     END IF
!          END DO
	  !WRITE(*,*) "After Photorates"

      !if (doAerosols .gt. 0) then
          !Get particle number conc. and radii for hetero. rate calculation
          cur => particles%first
          I = 1
          DO WHILE(associated(cur))
              HeteroNumberConc(I) = cur%NumberofParticles
              HeteroRadii(I) = cur%EffectiveRadius/100.0 !convert from cm to m
              !ResetHeteroRates expects radii in units of m
              cur => cur%next
              I = I+1
          END DO
                
          !Calculate Heterogeneous rates from size dist. info
          CALL RecalculateHeteroRates (HeteroRadii, HeteroNumberConc, HowManyBins, .FALSE.)
      !WRITE(*,*) "Hetero Okay"
     !endif ! doaerosols
     	!write(*,*) 'End of InputConcentrations calling spillbeans...',I
	!call spillbeans()
	END SUBROUTINE InputConcentrations

	SUBROUTINE OutputConcentrations(Temp, Press, Dens, GasConc, &
					AeroNumConc, AeroMassConc, &
					ExtCoeff, SingScat, Assym, Radius, TermVel)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine outputs Temperature, pressure,		!!
	!! gas and aerosol concentrations, and optical properties from	!!
	!! a 3D model and updates the appropriate concentrations in ASP !!
	!! grid box 1,1,1. MJA 052307					!!
        !! Removed Eulerian coordinates, MJA 021512                     !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: Temp, Press, Dens
		REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym,BackScatCoeff,SubExtCoeff, SubSSA, &
                                          ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                                          SubExtCoeffSep, SubSSASep
		REAL*8 :: GasConc(HowManyEvolveGasChems)!molecules/cm3
		REAL*8 :: AeroNumConc(HowManyBins)!particles/cm3
		REAL*8 :: AeroMassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems) !ug/m3
		REAL*8 :: Radius(HowManyBins) !cm
		REAL*8 :: TermVel(HowManyBins) !m/s
		
		!! Internal Variables
		INTEGER :: I, J, NumAq, NumModes
		REAL*8 :: SaltConc(HowManyAqChems)
	
		TYPE(Particle), POINTER :: Cur
		
		! CMB add
        integer :: doAerosols != getAerosolProc()
		real*8 :: filler_value, ro2_t, rco3_t
        doAerosols = getAerosolProc()
		filler_value = 1.0E-34
		! end CMB add
		
		!Output Temperature and Pressure
		Temp = GetTemp()
		!if (Temp .lt. 220.) then
		!	Temp = 285.0
		!endif	

		Press = GetPress()
		
		! CMB (AER): Bug fix: Apply peroxy getters
		!ro2_t = GetTotalPeroxy()
		!rco3_t = GetTotalAcylPeroxy()
		!write(678,*) 'In OutputConcentrations(): ro2_t = ', ro2_t
		!write(678,*) 'In OutputConcentrationS(): rco3_t = ', rco3_t
		! end CMB add
	
		!Output Gas chemical concentrations
		DO I = 1, HowManyEvolveGasChems
			GasConc(I) = GridGasChem(I)

		END DO

        !if (doAerosols .gt. 0) then
            !Update aerosol concentrations
            !write(*,*) 'Calling SortAerosolAtGridPointForCoagulation() from OutputConcentrations'
            CALL SortAerosolAtGridPointForCoagulation () 
             NumModes = ubound(AerosolModes, 1)
             !WRITE(*,*) "NumModes: ", NumModes

            cur => particles%first
            I = 1
            NumAq = HowManyAqChems+HowManyAqCations+HowManyAqAnions
            DO WHILE(associated(cur))
                AeroNumConc(I) = cur%numberofparticles
                Radius(I) = cur%effectiveradius
				!write(*,*) "Radius (um): ", Radius(I)*1e4
                TermVel(I) = TerminalVelocity(cur)/100.0 !Convert from cm/s to m/s
                !write(*,*) "Terminal velocity (m/s)", TermVel(I)
                !Water
                AeroMassConc(I,1) = cur%AqChems(1)*(&
                                cur%Numberofparticles*AqMolecularMass(1)) &
                    /(1.0e-12)
                               
                !Force minimum water
				IF(AeroMassConc(I,1) .LT. 0.001) 	AeroMassConc(I,1) = 0.001						   
                
				!Do this step so that we only pass neutral electrolytes
                !to dynamics code - helps prevent odd charge balance problems
                ! CMB (uncomment write statements for now
                !write(*,*) 'AeroMassConc(i,1) = ', aeromassconc(i,1)
                !WRITE(*,*) "Before output HypotheticalElectrolyteConcentrations"
                SaltConc = HypotheticalElectrolyteConcentrations (cur)
		!WRITE(*,*) cur%AqChems(2), cur%AqChems(3), cur%AqChems(4), cur%AqChems(5), cur%AqChems(6) 
                !WRITE(*,*) findchem('NH3', 1)
                IF(cur%Numberofparticles .GT. filler_value) then !0.) THEN
                    !Be careful with NH3 to avoid double counting, MJA, 2016-06-14
                    DO J = 2, HowManyAqChems
                        if (J .EQ. findchem('NH3', 1)) then
						    AeroMassConc(I,J) = (Cur%AqChems(J))*(cur%Numberofparticles*AqMolecularMass(J)) &
                                            /(1.0e-12)						
						else 
						    AeroMassConc(I,J) = (SaltConc(J)+Cur%AqChems(J))*(cur%Numberofparticles*AqMolecularMass(J)) &
                                            /(1.0e-12)
					    endif
                    END DO
                    	! - looping over Anion and Cation Aq Chems - CRL
                    DO J = HowManyAqChems+1, NumAq
                        AeroMassConc(I,J) = 0.
                        !aeromassconc(i,j) = 1.0e-30
                    END DO

                    DO J = 1, HowManyOrgChems
                            AeroMassConc(I,J+NumAq) = Cur%OrgChems(J)*(cur%Numberofparticles*OrgMolecularMass(J)) &
                                            /(1.0e-12)
                    END DO
                    
                    DO J = 1, HowManyAqOrgChems
                            AeroMassConc(I,J+NumAq+HowManyOrgChems) = Cur%AqOrgChems(J)*(cur%Numberofparticles*AqOrgMolecularMass(J)) &
                                            /(1.0e-12)
                    END DO
                ELSE !No particles
                    AeroNumConc(I) = 0.
                    !aeronumconc(i) = 1.0e-30
                    DO J = 2, NumAq+HowManyOrgChems+HowManyAqOrgChems ! looping over all Aerosols (including Anions and Cations) -CRL
                        AeroMassConc(I,J) = 0.
                        !aeromassconc(i,j) = 1.0e-30
                    END DO				
                END IF
          
                I = I + 1
                cur => cur%next
                    
            END DO

            !Output Average Aerosol optical properties, 0 flag assumes core in shell mixing
            !print *, 'calling aerosoloptprop within outputconcentrations()'
            call flush(6)
            !CALL AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, SubExtCoeff, SubSSA, 0)
            !WRITE(*,*) ExtCoeff(1), SingScat(1), Assym(1)
            !PAUSE
       ! endif ! do aerosols
        
	END SUBROUTINE OutputConcentrations

	SUBROUTINE InitializeASP(doAerosols, Transmissivity, sza)

	!! Include modules
	!USE IFPORT,				 ONLY : RTC, ETIME
	USE ModelParameters
	USE InfrastructuralCode	
	USE GridPointFields
	USE Chemistry
	USE Aerosols		
	USE Condensation
	USE OutputRoutines
	USE Coagulation
	
	implicit none

    ! CMB add
    integer :: ii
    real*8, intent(in), optional :: sza, Transmissivity
    !logical, intent(in), optional :: doAerosols
    integer, intent(in), optional :: doAerosols
    !logical :: todoAerosols
    integer :: todoAerosols
    ! end CMB add
    
	!! Define the local variables
	integer :: i, j, k, l, seedarray(2),q, r, NumBins
    INTEGER :: time(8), sizerand
	REAL*4  :: etimearray(2)
	TYPE(Particle), POINTER :: Cur, Env
	CHARACTER (len = 8)     :: ErrorDate
	CHARACTER (len = 10)    :: ErrorTime

    ! CMB add
    todoAerosols = 1
    if (present(doAerosols)) todoAerosols = doAerosols
    call setAerosolProc(todoAerosols)
    ! end CMB add
    
	!! Initialize the random number generator 
    call DATE_AND_TIME(values=time)     
    
    ! Get the current time (for pgi compiler)
    seedarray(1) = time(4) * (360000*time(5) + 6000*time(6) + 100*time(7) + time(8))
    
    !seedarray(1) = rtc() !for ifort compiler
    seedarray(2) = FLOOR(SQRT(REAL(seedarray(1))))

    ! CMB (AER): Modify this to play well with gfortran
	!CALL Random_Seed(put=seedarray)
    sizerand = 2
    CALL Random_Seed(size=sizerand)!, put=seedarray)

	!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SET the DOMAIN SIZE !!
	CALL SetDomainSize ()
	
	!! Call subroutines to set module data
	CALL SetFileHandleCounter(12)
	CALL ReadMainInputDeck()
	
	!! Send a Welcome Message
	CALL Transcript ("**********************************************************")
	CALL Transcript ("** Welcome to the ASP Model developed by Matt Alvarado  **")
	CALL Transcript ("** (mjalvara@mit.edu).                                  **")
	CALL Transcript ("** ASP is an updated, expanded version of MELAM         **")
	CALL Transcript ("** by H. D. Steele (donnan@mit.edu).                    **")
	CALL Transcript ("**********************************************************")
	CALL Date_and_Time(ErrorDate, ErrorTime)
	CALL Transcript ("*********************************************")
	CALL Transcript ("** Run Started at                          **")
    
    ! CB: Modified the below call so that vim doesn't have a syntax
    ! highlighting heart attack
	CALL Transcript ("** "//ErrorTime(1:2)//":"//ErrorTime(3:4)//":"//&
                         ErrorTime(5:6)//" on "//ErrorDate(5:6)//&
			 "/"//ErrorDate(7:8)//"/"//ErrorDate(1:4)//" "//" "//"**")
	CALL Transcript ("******************************************")
				
	! INITIALIZE the CHEMISTRY
    !print *, 'CMB: Setting Chemistry Params' 
	CALL SetChemistryParams
	!print *, 'CMB: Done with SetChemistryParams'

	! CB: This is NEEDED and done only within
	! DevelopmentDriver.f90, which is the main().  We're
	! not calling that in STILT-ASP...
	!print *, 'CMB: Reading PHotolysis Rates'
	call ReadPhotolysisRates
	!print *, 'CMB: Recalculating Reaction Rates'
	call RecalculateReactionRates(ifgas, ifaq)
    ! end cb; we may need more...

    ! CMB: call recalculate photo rates sza after call of recalculatereactionrates
    !      with "FirstTime" set to TRUE
    if (present(sza)) call RecalculatePhotoRatesSZA(sza, Transmissivity, .TRUE.)
    ! end CMB add
    
	!! INITIALIZE the DISSOLUTION routine
    !print *, 'CMB: Setting all dissolution...'
	CALL SetAllDissolution

	!! INITIALIZE the PARTICLES
    !print *, 'CMB: initializing the particles'        
	CALL InitializeParticlesSectional
				
	!!WARNING Force initial water equilibrium between gas and aerosol
	!!(Uses open system for water eq.)
	!CALL EquilibrateInternallyAtGridPoint(1,1,1, EquilibrateWater = .TRUE., WaterOpenSystem = .TRUE.)

	!write(*,*) 'Before RegridAerosol in InitializeASP calling spillbeans...'
	!call spillbeans()
	!write(*,*) 'Regridding aerosol'
	CALL RegridAerosol ()
	!write(*,*) 'After RegridAerosol in InitializeASP calling spillbeans...'
	!call spillbeans()    
	END SUBROUTINE InitializeASP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ASPInterface(Timestep, Temp, Press, Dens, GasConc, &
				NumConc, MassConc, &
				!PhotoRates, &
				ExtCoeff, SingScat, Assym, Radius, TermVel, sza, Transmissivity, do_aerosols)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine connects ASP to an external 3D Eulerian model!! 
	!! MJA 062707							!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			
		!USE IFPORT,				ONLY : RTC, ETIME
		USE ModelParameters,	ONLY : HowManyBins
		USE InfrastructuralCode	
		USE GridPointFields
		USE Chemistry
		USE Aerosols		
		USE Condensation
		USE OutputRoutines
		USE Coagulation

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: Timestep !(in s)
		REAL*8 :: Temp  !(in K)
		REAL*8 :: Press !(in mbar)
		REAL*8 :: Dens  !(in kg/m3 dry air)
		REAL*8 :: RH_ASP 
		REAL*8 :: GasConc(HowManyEvolveGasChems)
		REAL*8 :: NumConc(HowManyBins)
		REAL*8 :: MassConc(HowManyBins, HowManyAqChems+HowManyAqCations+HowManyAqAnions+HowManyOrgChems+HowManyAqOrgChems)
		REAL*8 :: Radius(HowManyBins)
		REAL*8 :: TermVel(HowManyBins)
		!REAL*8 :: PhotoRates(12)
		
		! CMB add: 
		real*8, intent(in) :: sza, Transmissivity
		integer, intent(in), optional :: do_aerosols
		!logical :: want_aerosols
		! end CMB add

                ! CB: Ifort complains unless these dimensions match
		!REAL*8, DIMENSION(18) :: ExtCoeff, SingScat, Assym
                REAL*8, DIMENSION(451):: ExtCoeff, SingScat, Assym
		INTEGER :: I,mm
	    ! CMB add
	    !want_aerosols = .false.
		RH_ASP = GetRelativeHumidity ()
	    if (present(do_aerosols)) call setAerosolProc(do_aerosols)
	    ! end CMB add
		CALL InputConcentrations(Temp, Press, Dens, GasConc, &
									NumConc, MassConc, &
									!PhotoRates, 
									sza, Transmissivity)

        call flush(6)									
		CALL StepASPOnce (TimeStep)
		CALL OutputConcentrations(Temp, Press, Dens, GasConc, &
									NumConc, MassConc, &
									ExtCoeff, SingScat, Assym, Radius, TermVel)
	END SUBROUTINE ASPInterface

!---------------------------
      subroutine convert_gas_concs_units(tpar,ppar,cbsum_ppb,numtypes,&
                                         cbsum_molec_cm3, direction)
      implicit none

      real*8, intent(IN) :: ppar ! Pa
      real*8, intent(IN) :: tpar
      real*8, intent(INOUT) :: cbsum_ppb(:)
      real*8, intent(INOUT) :: cbsum_molec_cm3(:)
      integer, intent(IN) :: numtypes
      real*8 :: factor
      integer :: direction, i ! 1: ppb->molec/cm3, -1: reverse, others: no action
      real, parameter :: R__J_K_mol = 8.3144598  ! J mol-1 K-1 or other units
      integer, parameter :: M3_TO_CM3 = 100000
      real*8, parameter :: avogadrosnum = 6.022140857E23

      factor = (ppar*avogadrosnum/R__J_K_mol/tpar)/ &
                M3_TO_CM3/1.0E9 ! extra 1e6 is for ppm, 1e9 for ppb
!      print *, 'multiplicative factor = ', factor
!      print *, 'this value converts ppmv to molec/cm3'
     !print *, '    factor = ', factor 
     !print *, '    numtypes = ', numtypes

      do i = 2, numtypes
        !print *, 'i = ', i
        if (direction .eq. 1) then
            !print *, '    PPB = ', cbsum_ppb(i) 
            cbsum_molec_cm3(1) = cbsum_ppb(1)
            cbsum_molec_cm3(i) = cbsum_ppb(i)*factor 
            !print *, '    Molec/cm3 = ', cbsum_molec_cm3(i)
        else if (direction .eq. -1) then
            cbsum_ppb(1) = cbsum_molec_cm3(1)
            cbsum_ppb(i) = cbsum_molec_cm3(i)/factor
            !print *, '    PPM = ', cbsum_ppb(i)
        endif
      enddo

      return
      end subroutine




!---------------------------
      subroutine SAM_wrapper(numtyp, naddp, & !rdirt, numtyp, radirt, naddp, &
                                 solarzenithangle, &
                                 cbsum_ppb, &
                                 tempK, ppar, &   
                                 density, timestep, &
                                 casum_ppb, &
                                 MassConc, NumConc, TermVel, Transmissivity,iden,zden)
      !use GridPointFields
      !use domain, ONLY: ntracers
      use Chemistry    ! already imported above
      use Aerosols

      !use Condensation
      use OutputRoutines
      !use Coagulation
      !use StepASP
      use ModelParameters, ONLY : HowManyBins
      !use InfrastructuralCodes   
      implicit none
      !include 'DEFCONC.INC'
      !implicit none
      integer mm, ntr,iden,zden
      !! ------------- from ASP ----------------- External Variables
      REAL*8 :: timestep !(in s)
      !REAL*8 :: Temp  !(in K)
      REAL*8 :: Press !(in mbar)
      REAL*8 :: Dens  !(in kg/m3 dry air)
      REAL*8 :: GasConc(HowManyEvolveGasChems)
      REAL*8 :: NumConc(HowManyBins)
      REAL*8 :: MassConc(HowManyBins, HowManyAqChems+&
                         HowManyAqCations+HowManyAqAnions+& ! in # particles/cm3
                         HowManyOrgChems+HowManyAqOrgChems) ! in ug/m3
      REAL*8 :: AfterNum(HowManyBins)
      REAL*8 :: AfterMass(HowManyBins, HowManyAqChems+&
                         HowManyAqCations+HowManyAqAnions+& ! in # particles/cm3
                         HowManyOrgChems+HowManyAqOrgChems) ! in ug/m3
      !REAL*8 :: MassConc(HowManyBins,65)
      REAL*8 :: TermVel(HowManyBins)
      REAL*8 :: Radius(HowManyBins)
      !REAL*8 :: PhotoRates(12)
      ! CB: Ifort complains unless these dimensions match
      !REAL*8, DIMnENSION(18) :: ExtCoeff, SingScat, Assym
      REAL*8, DIMENSION(451):: ExtCoeff, SingScat, Assym
      ! ---------------------- end from ASP
      integer, intent(in) :: numtyp, naddp
      !logical, intent(inout) :: to_initialize_asp
      !type(pset), allocatable, intent(inout) :: dirt(:)
      real*8 :: factor
      integer :: direction, i ! 1: ppm->molec/cm3, -1: reverse, others: no action
      real, parameter :: R__J_K_mol = 8.3144598  ! or other units
      integer, parameter :: M3_TO_CM3 = 100000
      real*8, parameter :: avogadrosnum = 6.022140857E23
      real*8, intent(INOUT) :: cbsum_ppb(:)    ! before (ignore ppm, actually in ppb)
      real*8, intent(INOUT) :: casum_ppb(:)    ! after 
      real*8, intent(in) :: ppar, solarzenithangle, density
      real*8 :: press_mbar, tempK
      REAL*8 :: cbsum_molec_cm3(HowManyEvolveGasChems)
      REAL*8 :: casum_molec_cm3(HowManyEvolveGasChems)
      integer :: allo_status, iter
      real*8 :: derp
	  real*8, intent(in) :: Transmissivity !Estimated cloud cover 
	  
      ! bug fix: pressures already in mb
      press_mbar = ppar/100.0   ! convert to mbar from Pa

!
! Convert the input concentrations 
!
!     GAS (ppb to molec/cm3)
      call convert_gas_concs_units(tempK, ppar, cbsum_ppb, numtyp, &
                                   cbsum_molec_cm3, 1)

!     AEROSOL in ( to ug/m3)

      ! Reset photolysis parameters based on solar zenith angle
      ! UPDATE: Move this after the concentration mapper
      ! TODO: Do I need to include a recalc of ReactionRates?
	!write(*,*) 'In SAM Wrapper'
          !do ntr=1,HowManyEvolveGasChems
		!write(*,*) GasPhaseChemicalNames(ntr),cbsum_ppb(ntr)
	  !enddo
       if (iden.eq.11.and.zden.eq.32) then
	write(*,*) 'In SAM Wrapper Before ASPInterface  '
        do mm=1,HowManyOrgChems
           write(*,*) OrgPhaseChemicalNames(mm),MassConc(:,mm+HowManyAqChems+HowManyAqCations+HowManyAqAnions)
        enddo
       endif

	  call ASPInterface(timestep, tempK, press_mbar, density, cbsum_molec_cm3, &
                        NumConc, MassConc, &
                        !PhotoRates, &
                        ExtCoeff, SingScat, Assym, Radius, TermVel, &
                        solarzenithangle, Transmissivity,1) 

       if (iden.eq.11.and.zden.eq.32) then
	write(*,*) 'In SAM Wrapper After ASPInterface  '
        do mm=1,HowManyOrgChems
           write(*,*) OrgPhaseChemicalNames(mm), MassConc(:,mm+HowManyAqChems+HowManyAqCations+HowManyAqAnions)
        enddo
       endif
      casum_molec_cm3(1:HowManyEvolveGasChems) = cbsum_molec_cm3(1:HowManyEvolveGasChems)
	!write(*,*) 'In SAM Wrapper After ASPInterface printing waterconc ',casum_molec_cm3(1), cbsum_molec_cm3(1), cbsum_ppb(1) 
      
      ! Convert the returned concentrations back into ppb and deallocate
      ! the array containing CB4 species concentrations in molec/cm3
      !print *, 'Convert ASP from molec/cm3 back to ppb'

      call convert_gas_concs_units(tempK, ppar, casum_ppb, numtyp, &
                                   casum_molec_cm3, -1)

      return
      end subroutine
END MODULE StepASP
