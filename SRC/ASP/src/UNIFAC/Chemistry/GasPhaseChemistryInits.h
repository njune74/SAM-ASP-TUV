!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! GasPhaseChemistryInits.h
!! Reads in and stores data on the gas phase chemicals and the
!! gas phase kinetic reaction mechanism.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!     								             !!
!! Month  Year   Name            Description				     !!
!! 07     2006   Matt Alvarado   Began Update History			     !!
!! 07/24  2006   Matt Alvarado   Change "Pointer => NULL()" to	             !!
!!				   NULLIFY(POINTER) to fit pgf90	     !!
!! 09/07  2006   Matt Alvarado   Added in Deposition Velocity for gases      !!
!! 06/27  2007   Matt Alvarado   Removed Sensitivity Matrix Calculation      !!
!! 07/30  2007   Matt Alvarado   Fixed background water concentration in     !!
!!				 SetGasPhaseChemistry - removed eps from calc!!
!! 02/15  2012   Matt Alvarado     Removing Eulerian grids, making ASP       !!
!!                                 a one-box model or subroutine.            !!
!! 01/24  2013   Matt Alvarado   Added GasPhasePeroxy array to subroutine    !!
!!                                 SetGasPhaseChemistry to read and store    !!
!!                                 peroxy radical flag.                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	!!
!! 1. SUBROUTINE SetGasPhaseChemistry ()			!!
!! 2. SUBROUTINE SetGasPhaseODEandJacobian ()			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Gas Phase Chemical Information from 'GasPhaseChemistryInputDeck.in',!!
!! an input deck, defines the appropriate chemical arrays in this module,   !!
!! and calls the appropriate eulerian grid subroutines to establich the	    !!
!! chemical fields					      		    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetGasPhaseChemistry ()

	USE ModelParameters
	USE InfrastructuralCode
	USE GridPointFields, ONLY: SetHomogeneousChemFieldPPB,	&
				   SetHomogeneousChemFieldMGperM3,&
				   AllocateGasChemistryGrid,	&
				   SetHowManyGasChems,		&
				   GridGasChem,			&
				   GetM,			&
				   GetPress,			&
				   GetSatVapPress,		&
				   GetRelativeHumidity,		&
				   GetSatMixingRatio
	IMPLICIT NONE

	!! Local Infrastructural Variables
	integer				:: i, j, k,		&
					   allocation_error,    &
					   index1, index2, currindex
	logical				:: Scaffolding, Transcribe
	character (len=512)	:: CurrLine, Token, Name
	character (len=12)	:: Numbers = "0123456789+-"

	!! These Local Variables Arrive From the Input Deck
	integer	:: NumbChems
	REAL*8	:: Mgas, MixH2O, RH, SPress

	!! And there is a file handle to be used
	integer	:: FH

        ! MJA 2016-09-13: Variable to temporarily store initial concentration
        REAL*8	:: init_conc

        ! MJA 2016-09-13: Variable to set initial concentrations to 0 if a STILT run (hardcoded for now!)
        logical	:: STILT    

        ! MJA 2017-09-22: Variable to assume the emission factor space is actually initial concentrations,
        !                 like in ASP v2.1 and earlier
        logical :: OLDASP    

        ! MJA 2016-09-13: Variables to stroe CO index and CO initial concentration (hardcoded for now!)
        integer	:: CO_index
        REAL*8	:: CO_init

        ! CMB (AER, Inc): Insert this parameter as a "filler" value for
        !                 which to analyze 
        real*8 :: filler_conc
        filler_conc = 1.0e-12

	FH = GetFileHandle()
	HowManyGasChems=1
        ! Water is the first, which is why these arent 0
	HowManyEvolveGasChems=1

	IF (TranscriptFH > 0) THEN
		Transcribe = .TRUE.
	ELSE 
		Transcribe = .FALSE.
	END IF

	!! Use Scaffolding Code?
	!Scaffolding = .TRUE.
	Scaffolding = .FALSE.

	!! SCSCSCSC -- Announce Arrival
	IF (Scaffolding) CALL TRANSCRIPT(">>Entering SetGasPhaseChemistry()<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Most of the Parameters Arrive from an Input Deck !!
	!! That deck is "GasPhaseChemistryInputDeck.in" and !!
	!! is located in the main directory of the program  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !
        ! CMB (AER, Inc): Trying to determine reason for ordering
        !                 differences between mine and ASP Native.

	OPEN(UNIT=FH, FILE=TRIM(InputDeckSubDir)//"GasPhaseChems.in", STATUS='OLD')

	!! The first line is a flag for input concentration type which equals:
	!!  0 if inputs are in ppb
	!!  1 if inputs are in mg / m3
        !WRITE(*,*) "Open GasPhaseChems.in"
	CurrLine = GetLine(FH)
 

        !MJA 2016-09-13: First token of first uncommeted line is:
        ! 0 if initial concentration should be calculated from the emission factors
        ! 1 if they should be set to a negligible value, as in STILT mode, and 
        ! 2 if they should be read in from the emission factor slot (to have back-compatability with older input files)
        STILT = .FALSE.         
        CALL GetToken(CurrLine, ";", Token)
        k = STR2INT (Token)

        IF (k .EQ. 1) STILT = .TRUE.
        IF (k .EQ. 2) OLDASP = .TRUE.
        IF (k .GT. 2 .OR. k .LT. 0) &
	  CALL ERROR("In GasPhaseChems.in, the first non-commentary line should be a flag telling the model how you are "// &
	  "specifying the initial concentrations.  It must either be 0, 1, or 2, with no decimal points.  Your input "// &
	  "is out of range.")
 
        !MJA 2016-09-13: Second token of first uncommeted line is initial CO concentration in ppbv        
        CALL GetToken(CurrLine, ";", Token)
        CO_init = STR2REAL (Token)

	GasConcInputTypeFlag = 0 !Force ppb, mja 08/27/2010 

	!! Make a first pass through the file to determine the number of 
	!! chemicals specified therein.
	!! Need this number early to allocate all of the arrays.
10	CurrLine = GetLine(FH)
            
        !WRITE(*,*) CurrLine
	IF (CurrLine .NE. EOF) THEN
		HowManyGasChems = HowManyGasChems+1
		GOTO 10
	END IF
        
	!! Send this value to the grid point module
	CALL SetHowManyGasChems(HowManyGasChems)

	!! Back up for a second pass to parse the chemical traits
	REWIND(FH)

	!! Skip the input flag line on this pass
	CurrLine = GetLine(FH)

	!! SCSCSCSC
	IF (Scaffolding) CALL TRANSCRIPT(">>About to allocate<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Now ALLOCATE each of the arrays that are size HowManyGasChems !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!! -- up to 15 character CHEMICAL NAME which should be identical to that -- !!!
	!!! -- given in other input decks if it is to be matched across phases    -- !!!
	ALLOCATE (GasPhaseChemicalNames(HowManyGasChems), stat = allocation_error)
	if (allocation_error > 0) &
		CALL ERROR("Allocation of GasPhaseChemicalNames could not proceed in SetChemistryParams()")

	!!! -- ENTHALPY OF VAPORIZATION (Latent Heat) in J/g -- !!! For H2O only!
	ALLOCATE (EnthalpyOfVaporization(1), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of EnthalpyOfVaporization could not proceed in SetChemistryParams()")

	!!! -- MOLECULAR MASS in g / mole -- !!!
	ALLOCATE (GasMolecularMass(HowManyGasChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of GasMolecularMass could not proceed in SetChemistryParams()")

        !!! -- PEROXY RADICAL FLAG --!!
        !!! -- 2 if should be included in RO2_T and RCO3_T sums (acyl peroxy) -- !!!
        !!! -- 1 if should be included in RO2_T sum only, 0 if not -- !!!
	ALLOCATE (GasPhasePeroxy(HowManyGasChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of GasPhasePeroxy could not proceed in SetChemistryParams()")
		
	!!! -- GAS PHASE EMISSION FACTOR in g/kg dry matter burned -- !!!
	ALLOCATE (GasPhaseEmisFactor(HowManyGasChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of GasPhaseEmisFactor could not proceed in SetChemistryParams()")

	!!! -- GAS PHASE LATERAL BACKGROUND CONCENTRATION in ppb -- !!!
	ALLOCATE (GasPhaseBackground(HowManyGasChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of GasPhaseBackground could not proceed in SetChemistryParams()")

	!!! -- DEPOSITION VELOCITY in cm/s -- !!!
	ALLOCATE (GasPhaseDepVel(HowManyGasChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of GasPhaseDepVel could not proceed in SetChemistryParams()")

		
	!! Make the grid of gas phase chemicals
	CALL AllocateGasChemistryGrid(HowManyGasChems)

	!! SCSCSCSCSC
	IF (Scaffolding) CALL Transcript(">>Importing Gas Phase Species Info<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! -- IMPORT OTHER GAS PHASE SPECIES INFORMATION FROM THE DECK -- !!
	!!  Evolving Species Go First, Non-Evolving Last	          !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Hard Code the Information related to Water
	GasPhaseChemicalNames(1)    = "H2O"
	EnthalpyOfVaporization(1)   = 2501. * joules / grams ! (Iribarne and Godson)
	GasMolecularMass(1)         = 18.016 * grams
	GasPhaseDepVel(1)	    = 0.0
	
	! CMB (AER): Bug fix
	GasPhasePeroxy(1) = 0.0

	index1 = 2		! Evolving Species are Placed with this index
	index2 = HowManyGasChems  ! Non Evolving Species are Placed with this index

	!! Transcribe
	IF (Transcribe) THEN
		CALL Transcript("")
		CALL Transcript("__Initializing_Gas_Phase_Chemicals_("//TRIM(INT2STR(HowManyGasChems))//"_Species)_")
		CALL Transcript("_Species_Names_and_Initial_Mixing_Ratio_")
	END IF

	DO i = 2, HowManyGasChems

		CurrLine = GetLine(FH)

		!! Tokenize the input line !!
		!! The token order is:
		!! ---	1. Chemical Name (15 Characters Max)
		CALL GetToken(CurrLine, ";", Name)
                !print *, 'CMB: Name = ', name(1:len(trim(name)))
		! Hold until know where to place it

		!! Don't want the user to specify water
		IF (TRIM(Name) .EQ. "H2O" .OR.		&
			TRIM(Name) .EQ. "h2o" .OR.		&
			TRIM(Name) .EQ. "water" .OR.	&
			TRIM(Name) .EQ. "Water" .OR.	&
			TRIM(Name) .EQ. "WATER") &
			CALL ERROR ("You specified "//TRIM(Name)//" in GasPhaseChems.in.  This model hardcodes in the values for "// &
			"water.  Please don't specify anything about it and use 'H2O' when referring to it in chemical "// &
			"mechanisms.")

		!! ---	2. Evolve Chemical or Keep Constant?
		CALL GetToken(CurrLine,";",Token) 
		k = STR2INT (Token)

		IF (k == 1) THEN
                        !print *, '  CMB: K == 1'
			HowManyEvolveGasChems = HowManyEvolveGasChems + 1
			CurrIndex             = index1
			index1		      = index1 + 1
                        !print *, '  CMB: HowManyEvolveGasChems = ', &
                        !        howmanyevolvegaschems
                        !print *, '  CMB: CurrIndex = ', currindex
                        !print *, '  CMB: index1 = ', index1
		ELSE IF (k == 0) THEN
                        !print *, '  CMB: K = 0'
			CurrIndex             = index2
			index2		      = index2 - 1
                        !print *, '  CMB: CurrIndex = ', CurrIndex
                        !print *, '  CMB: Index2 = ', index2
		ELSE
		        CALL ERROR("Don't Understand "//Trim(StripToken(Token))//" as a value for the Constant or evolve flag"// &
		        "For Chemical "//Trim(GasPhaseChemicalNames(CurrIndex))//" in GasPhaseChems.in")
		END IF

		!! Now Place Name in Correct Slot

		GasPhaseChemicalNames(CurrIndex) = Trim(StripToken(Name))

                IF (GasPhaseChemicalNames(CurrIndex) .EQ. "CO") CO_index = CurrIndex !MJA 2016-09-13

		!! ---  3. Peroxy Radical Flag
		CALL GetToken(CurrLine, ";", Token)
		GasPhasePeroxy(CurrIndex) = STR2REAL (Token)

		!! ---  4. Molecular Mass of the Gas
		CALL GetToken(CurrLine, ";", Token)
		GasMolecularMass(CurrIndex) = STR2REAL (Token) * grams

		!! ---	5. Emission factor (g/kg dry matter burned) and Initial Concentration (mixing ratio)
		CALL GetToken(CurrLine,";",Token) 
                GasPhaseEmisFactor(CurrIndex) = STR2REAL(Token) !g/kg dry matter burned, 09-13-16 MJA
 		
		!! ---  6. Background Conc. of the Gas
		!
                ! 2/4/2016 CB: Hard-set this to 0 for now, since we are
                !              getting background concentrations from
                !              STILT-Chem inputs
				!
				! Update: Don't need to do this; it's not used in our application
                Mgas = GetM   ()
		CALL GetToken(CurrLine, ";", Token)
		GasPhaseBackground(CurrIndex) = STR2REAL (Token) * Mgas/1e9
                ! convert from ppb to molecules/cm3
                !GasPhaseBackground(CurrIndex) = 0.0
		
		!! ---  7. Deposition Velocity of the gas (cm/s)
                !
		CALL GetToken(CurrLine, ";", Token)
                IF (STILT) THEN
                    ! 2/4/2016 CB: Hard-set this to 0 for now; a specified
                    !              deposition velocity turns off a few of
                    !              the desired deposition parameters in 
                    !              STILT-Chem
                    GasPhaseDepVel(CurrIndex) = 0.0
                ELSE
		    GasPhaseDepVel(CurrIndex) = STR2REAL (Token)
                ENDIF
	END DO

        !2016-09-13 MJA: Set initial fire concentrations - for now, use a logical flag to decide between STILT-ASP (set to ~0)
        ! and regular method 
        ! NOTE: CO MUST BE THE FIRST SPECIES IN GasPhaseChems.in AND MUST EVOLVE!
	write(*,*) 'CO_init = ',CO_init

        DO i = 2, HowManyGasChems
            IF (STILT) THEN                
                ! 2/4/2016 CB: We are setting these concentrations to 0
                !              since we are getting the boundary conditions
                !              from STILT-Chem
		CALL SetHomogeneousChemFieldPPB (i, filler_conc)            
            ELSE IF (OLDASP) THEN 
                !Assume Emission factor space is initial concentrations in ppb, like ASP v2.1 and earlier
                init_conc = GasPhaseEmisFactor(i)
                CALL SetHomogeneousChemFieldPPB (i, init_conc)    
            ELSE !Calculate initial concentration from CO initial concentration and EFs
		init_conc = (GasPhaseEmisFactor(i)/GasPhaseEmisFactor(CO_index))*(GasMolecularMass(CO_index)/GasMolecularMass(i)) &
                               *(CO_init - (GasPhaseBackground(CO_index)*1e9/Mgas))+(GasPhaseBackground(i)*1e9/Mgas) !in ppbv

                IF (GasPhaseChemicalNames(i) .EQ. "O2") init_conc = 210000000 !CRL 2017-03-21
                IF (GasPhaseChemicalNames(i) .EQ. "N2") init_conc = 780000000 !CRL 2017-03-21
                CALL SetHomogeneousChemFieldPPB (i, init_conc)     
            END IF

	    !! Transcript
	    IF (Transcribe) CALL Transcript(GasPhaseChemicalNames(i)//" "//TRIM(REAL2STR(init_conc))//" (ppb)")

        END DO

	!Set background concentration of water
        !
        ! 2/4/2016 CB: For now let's remove this, though we will have
        !              to account for it in the StepASP steps I think?
	RH = Back_RH
	SPress = GetSatVapPress()
	GasPhaseBackground(1) = ChemScale * SPress * RH * GetM() / (GetPress())
        !MixH2O*RH*Mgas/1e9 !molecules/cm3
        !GasPhaseBackground(1) = 0

	!! SCSCSCSC -- Double Check ordering and values
	IF (Scaffolding) THEN
	CALL Transcript("")
	DO i = 1, HowManyGasChems
		CALL Transcript(TRIM(INT2STR(i))//" "//Trim(GasPhaseChemicalNames(i))//" "//TRIM(REAL2STR(EnthalpyOfVaporization(i))) &
		//" "//TRIM(REAL2STR(GasMolecularMass(i))))
	END DO
	END IF

	CLOSE (FH)
	CALL ReturnFileHandle(FH)
	
        !! SCSCSCSC -- Announce Departure
	IF (.TRUE.) CALL Transcript(">>Exiting SetGasPhaseChemistry()<<")

	RETURN
END SUBROUTINE SetGasPhaseChemistry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Jacobian and ODE set for the chemical solver is formulated
!! as an array of linked lists.  This function populates interprets
!! the input deck and populates that array.
SUBROUTINE SetGasPhaseODEandJacobian ()
	use ModelParameters
	USE InfrastructuralCode, ONLY : INT2STR, ERROR, Transcript
	implicit none

	!! Input Vars
	integer :: NumbChem, ListLength
	
	!! Internal Variables
	integer	:: i, j, allocation_error
	logical	:: Scaffolding, Transcribe

	!! Use Scaffolding Code?
	!Scaffolding = .TRUE.
	Scaffolding = .FALSE.

	IF (TranscriptFH > 0) THEN
		Transcribe = .TRUE.
	ELSE
		Transcribe = .FALSE.
	END IF

	!!! -- The ODE Set for Gas Phase Chemistry -- !!
	ALLOCATE (GasPhaseODESet(HowManyEvolveGasChems), stat = allocation_error)
	if (allocation_error > 0) &
	CALL ERROR("Allocation of the Gas Phase ODE Linked List Array could not proceed in SetGasPhaseODEandJacobian()")

	!!! -- The Jacobian for Gas Phase Chemistry -- !!
	ALLOCATE (GasPhaseJacobian(HowManyEvolveGasChems,HowManyEvolveGasChems), stat = allocation_error)
	if (allocation_error > 0) &
	CALL ERROR("Allocation of the Gas Phase Jacobian Array could not proceed in SetGasPhaseODEandJacobian()")

	!! Initialize the List Arrays
	DO i = 1, HowManyEvolveGasChems
		NULLIFY(GasPhaseODESet(i)%FirstTerm)
		!GasPhaseODESet(i)%FirstTerm => Null()
		DO j = 1, HowManyEvolveGasChems
			NULLIFY(GasPhaseJacobian(i,j)%FirstTerm)
			!GasPhaseJacobian(i,j)%FirstTerm => Null()
		END DO
	END DO

	!! SCSCSC
	IF (Scaffolding) CALL TRANSCRIPT("Trying to find water in the system: "//TRIM(INT2STR(FindChem("H2O", GasPhase))))

	!! Now Create the appropriate ODE Set
	CALL MakeODESet('GasChemicalMechanism.in',GasPhaseODESet,HowManyEvolveGasChems,HowManyGasChems,"+",GasPhase)
	!! And Print its Results
	CALL PrintODESet (GasPhaseODESet,HowManyEvolveGasChems,GasPhase,"GasPhaseODESet.txt")

	IF (Transcribe) THEN
		CALL Transcript ("")
		CALL Transcript ("Printing Gas Phase ODE Set to (GasPhaseODESet.txt)")
	END IF

	!! Now Use that ODE Set to Make a Jacobian
	CALL MakeJacobian(GasPhaseJacobian, GasPhaseODESet, HowManyEvolveGasChems, GasPhase)

	!! And Print its Results
	CALL PrintJacobian (GasPhaseJacobian,HowManyEvolveGasChems,GasPhase,"GasPhaseJacobian.txt")
	IF (Transcribe) CALL Transcript ("Printing Gas Phase Jacobian to (GasPhaseJacobian.txt)")

	RETURN
END SUBROUTINE SetGasPhaseODEandJacobian
