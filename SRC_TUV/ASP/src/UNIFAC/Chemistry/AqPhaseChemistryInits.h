!! ASP (c), 2004-2017, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! AqPhaseChemistryInits.h
!! Reads in and stores data on the aqueous phase inorganic chemicals, 
!! aqueous phase kinetic reactions, and aqueous phase 
!! inorganic equilibrium reactions.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 07/24  2006   Matt Alvarado     Change "Pointer => NULL()" to	     !!
!!					NULLIFY(POINTER) to fit pgf90	     !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!!                                 Also fixed type-conversion errors         !!
!! 02/27  2012   Matt Alvarado     Removed Refractive Index Data, now        !!
!!                                  hardcoded in SUBROUTINE                  !!
!!                                  ShellRefIndAndRad(InParticle) in file    !!
!!                                  ParticleAttributes.h                     !!
!! 02/29  2012   Matt Alvarado     Added flag to AqPhaseChem.in and          !!
!!                                  SetAqPhaseChemistry to indicate what     !!
!!                                  refractive index should be used for an   !!
!!                                  electrolyte.                             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:  !!
!! 1. SUBROUTINE SetAqPhaseChemistry ()				!!
!! 2. SUBROUTINE SetAqPhaseODEandJacobian ()			!!
!! 3. SUBROUTINE SetAqPhaseEquilibriumReactions ()		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Aq Phase Chemical Information from 'AqPhaseChemistryInputDeck.in',  !!
!! an input deck, defines the appropriate chemical arrays in this module,   !!
!! and calls the appropriate eulerian grid subroutines to establich the	    !!
!! chemical fields							    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetAqPhaseChemistry ()

	USE ModelParameters, ONLY : SetHowManyVariousAqChems, &
                                    SetProtonIndex, SetHydroxyIndex, cm, &
				    grams, moles, WaterMolecMass, & 
                                    InputDeckSubDir, TranscriptFH,   &
				    AqChemScale, Atmospheres
	USE InfrastructuralCode

	IMPLICIT NONE

		!! Local Infrastructural Variables
		INTEGER				:: i, j, k,			 &
							   allocation_error, &
							   index1, index2, currindex
		CHARACTER(len=512)	:: CurrLine, Token, Name
		CHARACTER (len=12)	:: Numbers = "0123456789+-"

		!! These Local Variables Arrive From the Input Deck
		INTEGER							 :: NumbChems
		REAL*8, ALLOCATABLE			     :: TempMassStorage(:)
		CHARACTER(len=512), ALLOCATABLE  :: TempNameStorage(:)

		!! Transcribe?  Use Scaffolding Code?
		LOGICAL			   :: Transcribe
		LOGICAL, PARAMETER :: Scaffolding =  .FALSE. 

		!! And there is a file handle to be used
		integer	:: FH,FH2
		FH  = GetFileHandle()
		FH2 = GetFileHandle()
		HowManyAqChems		  = 1	! These are set to 1 instead of 0 because water is hard-coded in
		HowManyEvolveAqChems  = 1
		HowManyAqAnions		  = 0
		HowManyAqCations	  = 0

		!! SCSCSCSC -- Announce Arrival
		IF (Scaffolding) CALL Transcript(">>Entering SetAqPhaseChemistry()<<")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! There are several Aqueous Phase Input Decks Required for Initialization: !!!!!!!!!!!!!!!!!!!!!!
		!!   1). AqPhaseChems.in (List of Undissociate Chemicals Present in Aq Phase.)			!!
		!!	 2). AqPhaseIons.in  (List of Ions Dissociated from Chemicals Present in Aq Phase.)	!!
		!!	 3). AqChemicalMechanism.in (Forward Reactions in Aqueous Phase.)			!!
		!!   4). AqEquilibriumReactions.in (Equilibrium Dissociation Reactions in Aqueous Phase.)	!!
		!!   5). AerosolThermoKusikMeissner.in (Thermodynamic data for equilib rxns)			!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! In this subroutine, we will create listings of all of the chemicals 
	    OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'AqPhaseChems.in', STATUS='OLD')	! the incoming chemicals
	    OPEN(UNIT=FH2, FILE=TRIM(InputDeckSubDir)//'AqPhaseIons.in', STATUS='OLD')	! the incoming ions

		!! Make a first pass through the file to determine the number of 
		!! chemicals specified therein.
		!! Need this number early to allocate all of the arrays.
10		CurrLine = GetLine(FH)
		IF (CurrLine .NE. EOF) THEN
			HowManyAqChems = HowManyAqChems+1
			GOTO 10
		END IF

		!! Back up for a second pass to parse the chemical traits
		REWIND(FH)

20		CurrLine = GetLine(FH2)
		IF (CurrLine .NE. EOF) THEN

			!! Page through name and grab the charge
			CALL GetToken(CurrLine,";",Token)
			CALL GetToken(CurrLine,";",Token)
			IF (STR2REAL(Token) .LE. 0) THEN
				HowManyAqAnions  = HowManyAqAnions+1	
			ELSE 
				IF (STR2REAL(Token) .EQ. 0) &
				CALL ERROR ("In AqPhaseIons.in one of the ions was given a zero charge.  This shouldn't happen.")
				HowManyAqCations = HowManyAqCations+1
			END IF
			GOTO 20
		END IF
		REWIND(FH2)

		!! SCSCSCSC
		IF (Scaffolding) CALL Transcript(">>About to allocate<<")

		!! Transcribe
		IF (TranscriptFH > 0) THEN
			Transcribe = .TRUE.
		ELSE
			Transcribe = .FALSE.
		END IF
		Transcribe = .false.	! cmb mod, don't want this right now
		
		IF (Transcribe) THEN
			CALL Transcript("")
			CALL Transcript("__Initializing_Aqueous_Phase_Chemicals_("//TRIM(INT2STR(HowManyAqChems))//"_Species)_")
		END IF

		!! Tell the ModelParameter Module How Many of Each Type of Chems There Are
		CALL SetHowManyVariousAqChems(HowManyAqChems, HowManyAqCations, HowManyAqAnions)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Now ALLOCATE each of the arrays that are size HowManyAqChems !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!! -- up to 15 character CHEMICAL NAME which should be identical to that -- !!!
		!!! -- given in other input decks if it is to be matched across phases    -- !!!
		ALLOCATE (AqPhaseChemicalNames(HowManyAqChems), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqPhaseChemicalNames could not proceed in SetAqPhaseChemistry()")

		!!! -- MOLECULAR MASS (grams / mole) -- !!!
		ALLOCATE (AqMolecularMass(HowManyAqChems), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqMolecularMass could not proceed in SetChemistryParams()")


		!!! -- SOLID SALT DENSITY (grams / cm3) -- !!!
		ALLOCATE (SolidSaltDensity(HowManyAqChems), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of SolidSaltDensity could not proceed in SetChemistryParams()")

		!!! -- REFRACTIVE INDEX FLAG -- !!!
		ALLOCATE (RefIndFlag(HowManyAqChems), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of RefIndFlag could not proceed in SetChemistryParams()")
		
		!! SCSCSCSCSC
		IF (Scaffolding) CALL Transcript (">>Importing Aq Phase Species Info<<")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! -- IMPORT OTHER GAS PHASE SPECIES INFORMATION FROM THE DECK -- !!
		!!  Evolving Species Go First, Non-Evolving Last				  !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		index1 = 2				 ! Evolving Species are Placed with this index
		index2 = HowManyAqChems  ! Non Evolving Species are Placed with this index

		!! Hard Code the Water 
		AqPhaseChemicalNames(1)		= "H2O"
		AqMolecularMass(1)			= WaterMolecMass
		
		DO i = 2, HowManyAqChems

			CurrLine = GetLine(FH)

			!! Tokenize the input line !!
			!! The token order is:
			!! ---	1. Chemical Name (15 Characters Max)
			CALL GetToken(CurrLine, ";", Name)				! Hold until know where to place it

			!! Don't want the user to specify water
			IF (TRIM(Name) .EQ. "H2O" .OR.		&
				TRIM(Name) .EQ. "h2o" .OR.		&
				TRIM(Name) .EQ. "water" .OR.	&
				TRIM(Name) .EQ. "Water" .OR.	&
				TRIM(Name) .EQ. "WATER")        &
				CALL ERROR ("You specified "//TRIM(Name)//" in AqPhaseChems.in.  This model hardcodes in the values for water."// &
				            "  Please don't specify anything about it and use 'H2O' when referring to it in chemical mechanisms.")

			!! ---	2. Evolve Chemical or Keep Constant?
			CALL GetToken(CurrLine,";",Token) 

			IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In AqPhaseChems.in, the Evolve or Not flag appears to be missing on "//  &
								TRIM(AqPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")

			k = STR2INT (Token)

			IF (k == 1) THEN
				HowManyEvolveAqChems = HowManyEvolveAqChems + 1
				CurrIndex             = index1
				index1				  = index1 + 1
			ELSE IF (k == 0) THEN
				CurrIndex             = index2
				index2				  = index2 - 1
			ELSE
				CALL ERROR("Don't Understand "//Trim(StripToken(Token))//" as a value for the Constant or evolve flag"//	&
						   "For Chemical "//Trim(AqPhaseChemicalNames(CurrIndex))//" in AqPhaseChems.in. ")
			END IF

			!! Now Place Name in Correct Slot
			AqPhaseChemicalNames(CurrIndex) = Trim(StripToken(Name))

			!! --- 3. Molecular Mass (g / mole)
			CALL GetToken(CurrLine,";",Token)
			IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In AqPhaseChems.in, the molecular mass appears to be missing on "//   &
							TRIM(AqPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
			AqMolecularMass(CurrIndex) = STR2REAL(Token) * grams / moles

			
			!! 4. SOLID SALT DENSITY (g/cm3)
			CALL GetToken(CurrLine,";",Token) 
			IF (LEN_TRIM(Token) .EQ. 0) &
			CALL ERROR("In AqPhaseChems.in, the solid salt density term appears to be missing on "//	&
					   TRIM(AqPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
			SolidSaltDensity(CurrIndex) = STR2REAL (Token)*grams/cm/cm/cm

	                !! 4. REFRACTIVE INDEX FLAG: 0 if the H2O refractive index is used
                        !!                           1 if the sulfate droplet refractive index is used
                        !!                           2 if the sea salt refractive index is used
                        !!                           3 if the mineral dust refractive index is used
			CALL GetToken(CurrLine,";",Token) 
			IF (LEN_TRIM(Token) .EQ. 0) &
			CALL ERROR("In AqPhaseChems.in, the refractive index flag appears to be missing on "//	&
					   TRIM(AqPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
			RefIndFlag(CurrIndex) = STR2INT(Token)
			
			!! Transcript
			IF (Transcribe) CALL Transcript(AqPhaseChemicalNames(CurrIndex))

		END DO

		!! SCSCSCSC -- Announce Departure
		IF (Scaffolding) CALL Transcript(">>Exiting SetAqPhaseChemistry()<<")
		
		CLOSE (FH)
		CALL ReturnFileHandle(FH)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Then ALLOCATE all of the ION related distributions !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!! -- MOLECULAR NAMES, revisited to include everything in one array -- !!!
		ALLOCATE (TempNameStorage(HowManyAqChems), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of TempNameStorage could not proceed in SetChemistryParams()")

		DO I = 1, HowManyAqChems
			TempNameStorage(I) = TRIM(AqPhaseChemicalNames(I)) 
		END DO
		DEALLOCATE (AqPhaseChemicalNames)

		ALLOCATE (AqPhaseChemicalNames(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqPhaseChemicalNames could not proceed in SetChemistryParams()")

		DO I = 1, HowManyAqChems
			AqPhaseChemicalNames(I) = TRIM(TempNameStorage(I))
		END DO
		DEALLOCATE (TempNameStorage)

		!!! -- CATION NAMES -- !!!
		ALLOCATE (AqCationNames(HowManyAqCations), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqCationNames could not proceed in SetAqPhaseChemistry()")

		!!! -- ANION NAMES -- !!!
		ALLOCATE (AqAnionNames(HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqAnionNames could not proceed in SetAqPhaseChemistry()")

		!!! -- CHARGE INDEX FOR ALL CHEMICALS -- !!!
		ALLOCATE (AqChemicalCharge(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqChemicalCharge could not proceed in SetAqPhaseChemistry()")

		DO I = 1, HowManyAqChems
			AqChemicalCharge(I) = 0
		END DO

		!!! -- CATION CHARGES -- !!!
		ALLOCATE (AqCationCharge(HowManyAqCations), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqCationCharge could not proceed in SetAqPhaseChemistry()")

		!!! -- ANION CHARGES -- !!!
		ALLOCATE (AqAnionCharge(HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqAnionCharge could not proceed in SetAqPhaseChemistry()")

		!!! -- CATION MASS -- !!!
		ALLOCATE (AqCationMass(HowManyAqCations), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqCationMass could not proceed in SetAqPhaseChemistry()")

		!!! -- ANION MASS -- !!!
		ALLOCATE (AqAnionMass(HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqAnionMass could not proceed in SetAqPhaseChemistry()")

		!!! -- MOLECULAR (grams) -- !!!
		!!! -- Revisit this array to include Cations and Anions in it!
		ALLOCATE (TempMassStorage(HowManyAqChems), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of TempMassStorage could not proceed in SetChemistryParams()")

		DO I = 1, HowManyAqChems
			TempMassStorage(I) = AqMolecularMass(I)
		END DO
		DEALLOCATE (AqMolecularMass)

		ALLOCATE (AqMolecularMass(HowManyAqChems+HowManyAqCations+HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqMolecularMass could not proceed in SetChemistryParams()")

		ALLOCATE (CationIonicRefraction(HowManyAqCations), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqMolecularMass could not proceed in SetChemistryParams()")

		ALLOCATE (AnionIonicRefraction(HowManyAqAnions), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqMolecularMass could not proceed in SetChemistryParams()")
		
		DO I = 1, HowManyAqChems
			AqMolecularMass(I) = TempMassStorage(I)
		END DO
		DEALLOCATE (TempMassStorage)

		!! SCSCSCSCSC
		IF (Scaffolding) CALL Transcript (">>            ---            <<")
		IF (Scaffolding) CALL Transcript (">>Importing Aq Phase Ion Info<<")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! -- Loop Over the Ion Inputs -- !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		!! Index1 is for Cations ; Index2 is for Anions
		Index1 = 1 ; Index2 = 1

		!! Tokenize the input line !!
		DO i = 1, HowManyAqCations + HowManyAqAnions 

			!! Parse from Second File
			CurrLine = GetLine(FH2)

			!! ---	1. Name (15 Characters Max)
			CALL GetToken(CurrLine, ";", Name)	! Hold until know where to place it

			!! ---	2. Charge
			CALL GetToken(CurrLine,";",Token) 
			k = STR2INT (Token)

			IF (k .GT. 0) THEN
				AqCationNames(Index1)  = Trim(StripToken(Name))
				AqCationCharge(Index1) = k
				CALL Transcript(TRIM(Name)//"		(+"//TRIM(INT2STR(k))//")")
			ELSE 
				AqAnionNames(Index2)  = Trim(StripToken(Name))
				AqAnionCharge(Index2) = k
				CALL Transcript(TRIM(Name)//"		("//TRIM(INT2STR(k))//")")
			END IF

			!! --- 3. Mass
			CALL GetToken(CurrLine,";",Token) 

			IF (k .GT. 0) THEN
				AqCationMass(Index1) = STR2REAL(Token) * grams / moles
			ELSE 
				AqAnionMass(Index2)  = STR2REAL(Token) * grams / moles
			END IF
			
			!! --- 4. Ionic Refraction
			CALL GetToken(CurrLine,";",Token)

			IF (k .GT. 0) THEN
				CationIonicRefraction(Index1) = STR2REAL(Token)
			ELSE 
				AnionIonicRefraction(Index2)  = STR2REAL(Token)
			END IF

			!! --- Iterate Counters
			IF (k .GT. 0) THEN
				index1 = index1 + 1
			ELSE 
				index2 = index2 + 1
			END IF

		END DO

		!! Now refill the Molecular Mass and Name Array to include both full molecules and ions
		DO I = 1, HowManyAqCations
			AqMolecularMass(I+HowManyAqChems)      = AqCationMass(I)
			AqPhaseChemicalNames(I+HowManyAqChems) = TRIM(AqCationNames(I))
			AqChemicalCharge(I+HowManyAqChems)	   = AqCationCharge(I)
		
		END DO
		DO I = 1, HowManyAqAnions
			AqMolecularMass(I+HowManyAqChems+HowManyAqCations)		= AqAnionMass(I)
			AqPhaseChemicalNames(I+HowManyAqChems+HowManyAqCations) = TRIM(AqAnionNames(I))
			AqChemicalCharge(I+HowManyAqChems+HowManyAqCations)	    = AqAnionCharge(I)
		
		END DO

		!! Tell ModelParameters which of the chemicals is a proton
		!! Will just give a general error if don't find H+, but every model scenario
		!! should include H+ as an option
		CALL SetProtonIndex(findchem("H+",AqPhase))
		CALL SetHydroxyIndex(findchem("OH-",AqPhase))

		CLOSE (FH2)
		CALL ReturnFileHandle(FH2)

		RETURN
	END SUBROUTINE SetAqPhaseChemistry 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The Jacobian and ODE set for the chemical solver is formulated	!!
	!! as an array of linked lists.  This function populates interprets !!
	!! the input deck and populates that array.							!!
	SUBROUTINE SetAqPhaseODEandJacobian ()

		USE InfrastructuralCode, ONLY : ERROR, Transcript
		USE ModelParameters

		IMPLICIT NONE

		!! Input Vars
		integer :: NumbChem, ListLength
		
		!! Internal Variables
		integer	:: i, j, allocation_error
		logical	:: Scaffolding, Transcribe

		!! Use Scaffolding Code?  Transcribe?
		Scaffolding = .FALSE. ! .true. ! cmb
		IF (TranscriptFH > 0) THEN
			Transcribe = .TRUE.
		ELSE
			Transcribe = .FALSE.
		END IF
		transcribe = .false.	! cmb
		
		!!! -- The ODE Set for Gas Phase Chemistry -- !!
		ALLOCATE (AqPhaseODESet(HowManyEvolveAqChems), stat = allocation_error)
		if (allocation_error > 0) &
		CALL ERROR("Allocation of the Aqueous Phase ODE Linked List Array could not proceed in SetAqPhaseODEandJacobian()")

		!!! -- The Jacobian for Gas Phase Chemistry -- !!
		ALLOCATE (AqPhaseJacobian(HowManyEvolveAqChems,HowManyEvolveAqChems), stat = allocation_error)
		if (allocation_error > 0) &
		CALL ERROR("Allocation of the Aqueous Phase Jacobian Array could not proceed in SetAqPhaseODEandJacobian()")

		!! Initialize the List Arrays
		DO i = 1, HowManyEvolveAqChems
			NULLIFY(AqPhaseODESet(i)%FirstTerm)
			!AqPhaseODESet(i)%FirstTerm => Null()
			DO j = 1, HowManyEvolveAqChems
				NULLIFY(AqPhaseJacobian(i,j)%FirstTerm)
				!AqPhaseJacobian(i,j)%FirstTerm => Null()
			END DO
		END DO

		!! Now Create the appropriate ODE Set
		IF (IfAq) THEN
			CALL MakeODESet('AqChemicalMechanism.in',AqPhaseODESet,HowManyEvolveAqChems,HowManyAqChems,"+",AqPhase)

			!! And Print its Results
			CALL PrintODESet (AqPhaseODESet,HowManyEvolveAqChems,AqPhase,"AqPhaseODESet.txt")
			IF (Transcribe) THEN
				CALL Transcript("")
				CALL Transcript("Printing Aqueous Phase ODE Set to (AqPhaseODESet.txt)")
			END IF

			!! Now Use that ODE Set to Make a Jacobian
			CALL MakeJacobian(AqPhaseJacobian, AqPhaseODESet, HowManyEvolveAqChems, AqPhase)

			!! And Print its Results
			CALL PrintJacobian (AqPhaseJacobian,HowManyEvolveAqChems,AqPhase,"AqPhaseJacobian.txt")
			IF (Transcribe) CALL Transcript("Printing Aqueous Phase Jacobian to (AqPhaseJacobian.txt)")
		END IF

		RETURN
	END SUBROUTINE SetAqPhaseODEandJacobian

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Read in the Equilibrium Reaction and Set up the Appropriate	!!
	!! Arrays for dealing with them.								!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetAqPhaseEquilibriumReactions ()

		USE InfrastructuralCode, ONLY : Warn,						&
										ERROR,						&
										Transcript,					&
										GetFileHandle,				&
										ReturnFileHandle,GetLine,	&
										GetToken,					&
										IsReal, STR2REAL,			&
										STR2INT, INT2STR,			&
										REAL2STR, StripToken, EOF

		USE ModelParameters,	 ONLY : TranscriptFH,				&
										InputDeckSubDir,			&
										grams, cm, m, moles,		&
										kilograms, dyncm,			&
										SetWaterEquilibriumIndex,	&
										WaterEquilibriumIndex,		&
										Metastable

		IMPLICIT NONE

		!! Internal Variables
		INTEGER				:: i, j, k, l, mm, n, O, P, Q, R, IonI(2), allocation_error, AnionIndex, CationIndex, electrolyteRxn
		LOGICAL				:: Scaffolding, Transcribe
		CHARACTER (LEN=256)	:: CurrLine, FullLine
		CHARACTER (LEN=64)	:: AssociatedForm, DissociatedForm, Token, Ion
		REAL*8			    :: II, JJ

		!! And there is a file handle to be used
		integer	:: FH
		FH  = GetFileHandle()

		!! Use Scaffolding Code?  Transcribe?
		Scaffolding = .FALSE. ! .true. ! cmb

		IF (TranscriptFH > 0) THEN
			Transcribe = .TRUE.
		ELSE
			Transcribe = .FALSE.
		END IF
		transcribe = .false.	! cmb
		
		!! SCSCSCSC -- Announce Arrival
		IF (Scaffolding) CALL Transcript("")
		IF (Scaffolding) CALL Transcript(">>Entering SetAqPhaseEquilibriumReactions()<<")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! There are several Aqueous Phase Input Decks Required for Initialization: !!!!!!!!!!!!!!!!
		!!   1). AqPhaseChems.in (List of Undissociate Chemicals Present in Aq Phase.)			  !!
		!!	 2). AqPhaseIons.in  (List of Ions Dissociated from Chemicals Present in Aq Phase.)	  !!
		!!	 3). AqChemicalMechanism.in (Forward Reactions in Aqueous Phase.)					  !!
		!!   4). AqEquilibriumReactions.in (Equilibrium Dissociation Reactions in Aqueous Phase.) !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'AqEquilibriumReactions.in', STATUS='OLD')

		!! SCSCSCSC
		IF (SCAFFOLDING) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("__Counting aqueous phase equilibrium reactions...")
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Count the Number of Reactions    !!
		!! (making this a two pass process) !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! Initialize the index.  A single electrolyte may be doubly counted in this list
		!! if it may dissociate into multiple sets of ions.  
		!! Equally a single set of ions may be doubly counted, meaning two electrolytes 
		!! dissociating into 
		!! So it is a # of Equilibrium Rxns, really.
		HowManyAqEqReactions = 0
		
5		CurrLine = GetLine(FH)
		IF (TRIM(CurrLine) .NE. EOF) THEN
			HowManyAqEqReactions = HowManyAqEqReactions + 1
			GOTO 5
		END IF

		REWIND FH

		!! SCSCSCSC
		IF (SCAFFOLDING) THEN
			CALL TRANSCRIPT("Counted "//TRIM(INT2STR(HowManyAqEqReactions))//" Aqueous Equilibrium Reactions")
		END IF


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Make a Tally List of the Equilibria So Don't Have to Check Every Mixture !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This list has three data values:											!!
		!!		1). The Aq Chemical Index of the Electrolyte						!!
		!!		2). The Aq Cation Index of the Cation								!!
		!!		3). The Aq Anion Index of the Anion									!!
		ALLOCATE (AqEquilibriaList(HowManyAqEqReactions,25), stat = allocation_error)
		if (allocation_error > 0) &
		CALL ERROR("Allocation of AqEquilibriaList could not proceed in SetAqPhaseEquilibriumReactions()")

		DO I = 1, HowManyAqEqReactions
			DO J = 1, 25
				AqEquilibriaList(I,J) = 0.
			END DO
		END DO

		!!! -- ION-ELECTROLYTE MAP and EQUILIBRIUM DATA -- !!!
		!!! --
		!!! -- DECODING KEY FOR AqEquilibriaList (a,b)
		!!! -- 
		!!! --   (i,1) : Electrolyte  : this is -1 if water is the only dissociator (e.g., H2O <=> H+ & OH-)
		!!! --   (i,2) : Cation 
		!!! --   (i,3) : Anion
		!!! --   (i,4) : Cation Stoicheometric Coefficient
		!!! --   (i,5) : Anion  Stoicheometric Coefficient
		!!! --   (i,6) : Flag for solid compound
		!!! --   (i,7) : Blank 
		!!! --   (i,8) : Number of H2O Molecules Needed on Left Hand Side of Dissociation
		!!! --   (i,9) : Equilibrium Coefficient at 298 K (either true or infinite, which would
		!!! --					imply that it condenses directly into a dissociated form)
		!!! --   (i,10): - Delta Ho0 / R T0   for Equilibrium Coefficient Temperature Dependence
		!!! --   (i,11): - Delta Cp0 / R      for Equilibrium Coefficient Temperature Dependence
		!!! --   (i,12): First Density Parameter for dissociated ion solution
		!!! --   (i,13): Second Density Parameter for dissociated ion solution
		!!! --   (i,14): Surface Tension Surface Excess at Saturation (mol / m2)
		!!! --   (i,15): Surface Tension Adsorption coefficient (dimentionless)
		!!!
		!!! ----------- Indexing of activities for various reaction types:
		!!! --- for simple reactions where everything fits KM, (I,17) = (I,18) = 0 and
		!!! ---  (I,16): Index into KMRef for thermodynamic reaction data
		!!!
		!!! --- for dissociating electrolytes (like bisulfate) (I,18) = 0 and
		!!! ---  (I,16): Index into AqEquilibriaList of numerator reaction
		!!! ---  (I,17): Index into KMRef of denominator reaction (which is this actual reaction)
		!!!
		!!! --- for complex reactions in which the activity is a composite of several other reactions
		!!! --   (i,16): Numerator Activity for Equilibrium Reaction (stored as reaction+exponent/10,000)
		!!! --   (i,17): Denominator Activity for Equilibrium Reaction (stored as reaction+exponent/10,000)
		!!! --   (i,18): Another Numerator Activity 
		!!! --- for levitocite (NH4)3H(SO4)3, follow Kim et al., 1993 (I,18) = 0 and
		!!! --   (i,16): Index into AqEquilibriaList for Ammonium Sulfate
		!!! --   (i,17): Index into AqEquilibriaList for Sulfuric Acid 
		!!!
		!!! --   (i,19): A flag: if this is a positive number then this reaction will be assigned the mixed coefficient of 
		!!!              another reaction (for bisulfate mixing and Levitocite)
		!!! --   (i,20): DRH at 298 K
		!!! --   (i,21): Temperature dependence of DRH
		!!! --   (i,22): Third dissociated species index (for levitocite)
		!!! --   (i,23): Third dissociated species stoichiometry
		!!! --   (i,24): Third dissociated species index (for hydrates, should be water)
		!!! --   (i,25): Third dissociated species stoichiometry

		!!! -- HYDRATE REACTIONS NOT YET IMPLEMENTED
		!!! -- will need a data structure and other issues

		! CMB: don't want these print statements right now
		!CALL TRANSCRIPT("")
		!CALL TRANSCRIPT("_About_to_Import_Aqueous_Phase_Equilibrium_Reactions_")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!! READ THE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Grab the next equilibrium reaction and process the inputs !!
		DO ElectrolyteRxn = 1, HowManyAqEqReactions

			CurrLine = GetLine(FH)

			IF (TRIM(CurrLine) .EQ. EOF) &
			CALL ERROR ("Somehow we miscounted the number of lines in AqEquilibriumReactions.in SetAqPhaseEquilibriumReactions "// &
						"().  This is kind of boffo.  Shouldn't happen.  Sorry.")

			FullLine = CurrLine

			!! Parse to separate Electrolyte, Ions, and Parameters
			CALL GetToken(CurrLine, "<=>", AssociatedForm)
			CALL GetToken(CurrLine, ";",   DissociatedForm)

			!! Record the Reaction
			! CMB: don't print this every time 
			!CALL TRANSCRIPT(TRIM(AssociatedForm)//" <=> "//TRIM(DissociatedForm))

			!! Find the Index Point in terms of Anions & Cations
			K = 0; L = 0; MM = 0; N = 0; O = 0; P = 0; Q = 0; R = 0;
			DO I = 1, 3

				J = 0

				CALL GetToken(DissociatedForm, "&", Ion)
				IF (len_trim(Ion) == 0) EXIT

				!! Check Token for Leading Factor
				CALL GetToken(Ion," ",Token)

				IF (IsReal(Token)) THEN
					J     = STR2REAL(Token)			! Store
					Token = StripToken(Ion)			! Push to Ion Name
				ELSE
					J     = 1
					Ion   = TRIM(TRIM(Ion) // " " // TRIM(Token))
				END IF

				IonI = FindIon(TRIM(Ion),SuppressError=.TRUE.)

				!! Cation
				IF (IonI(1) .GT. 0) THEN
					IF (K .GT. 0) CALL ERROR("Two Cations have been specified in the dissociation reaction of "//  &
											  TRIM(AssociatedForm)//" in AqEquilibriumReactions.in.  This is not allowed.")
					K = J
					L = IonI(2)
				END IF

				!! Anion
				IF (IonI(1) .LT. 0) THEN
					IF (MM .GT. 0) THEN
						CALL WARN("Two Anions have been specified in the dissociation reaction of "//TRIM(AssociatedForm)// &
							   " in AqEquilibriumReactions.in.  This should only be done for levitocite (NH4)3H(SO4)2.")
						O = J !Stoich
						P = IonI(2) !Second Anion
					
					ELSE	
						MM = J
						N = IonI(2)
					END IF
				END IF

				!! HYDRATION CHECK
				IF (IonI(1) .EQ. 0) THEN
				
					IonI(2) = FindChem(TRIM(Ion), AqPhase, SuppressError=.TRUE.)
					
					IF(AqPhaseChemicalNames(IonI(2)) .EQ. "H2O" .AND. I .EQ. 3) THEN
						!Hydrate Reaction
						Q = J		!Stoich
						R = IonI(2) !Index for water
					ELSE
						CALL ERROR("Couldn't locate the ion ("//trim(Ion)//") which you specified in "// &
												   "AqEquilibriumReactions.in.  Please make sure you've specified it in the "// &
												   "AqPhaseIons.in input deck.")
						
					END IF
				END IF
			
			END DO

			!! If there are more dissociated elements in the pipe then get confused
			IF (LEN_TRIM(StripToken(DissociatedForm)) .NE. 0) &
			CALL WARN ("  Warning!! The following Ion: ("//TRIM(StripToken(DissociatedForm))// &
					   ") Appeared to be Extra in Line... "//"  >>"//TRIM(FullLine)//"<<  Ignoring and Continuing...")

			!! Fill the Stoicheometric Coefficients and identifications
			AqEquilibriaList(electrolyteRxn,2) = REAL(L)
			AqEquilibriaList(electrolyteRxn,3) = REAL(N)
			AqEquilibriaList(electrolyteRxn,4) = REAL(K)
			AqEquilibriaList(electrolyteRxn,5) = REAL(MM)
			AqEquilibriaList(electrolyteRxn,22) = REAL(P)
			AqEquilibriaList(electrolyteRxn,23) = REAL(O)
			AqEquilibriaList(electrolyteRxn,24) = REAL(R)
			AqEquilibriaList(electrolyteRxn,25) = REAL(Q)
			IF(AqEquilibriaList(electrolyteRxn,22) .NE. 0) AqEquilibriaList(electrolyteRxn,19) = 1

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Now Parse the Electrolyte !!
			K = 0; MM = 0
			AqEquilibriaList(electrolyteRxn,1) = -1
			AqEquilibriaList(electrolyteRxn,8) = 0  ! Assume no water unless find some
			AqWaterDissociationReaction        = 0  ! And no water dissociation reaction

			DO I = 1, 2
				J = 0
				CALL GetToken(AssociatedForm, "&", Ion)

				IF (len_trim(Ion) == 0) EXIT

				!! Check Token for Leading Factor
30				CALL GetToken(Ion," ",Token)

				!! Should be water if Stoich > 1
				IF (IsReal(Token)) THEN

					!! error check
					IF (STR2REAL(Token) .EQ. 1) &
					CALL ERROR (" Don't specify a stoicheometric coefficient of one!  Simply leave it out.  Found such a "// &
								"gross thing in: >>"//TRIM(FullLine)//"<< ")

					AqEquilibriaList(electrolyteRxn,8) = STR2REAL(Token)
					Token = StripToken(Ion)				! Push to Chemical Name

					!! Only Water Can Have Multiple Instances
					IF (FindChem(TRIM(Token),AqPhase) .EQ. 1) THEN
						GOTO 30
					ELSE
						CALL ERROR ("Indicated a Stoicheometric Coefficient Above One for the Electrolyte in the following "// &
									"line.  This is not acceptable; we can only take one of the electrolyte.  Sorry.  The "//  &
									"offending line: "//TRIM(FullLine))
					END IF
				ELSE
					Ion = TRIM(Ion) // " " // TRIM(Token)
				END IF

				!! Check to see if it is water
				IF (FindChem(TRIM(Token),AqPhase) .EQ. 1) THEN
					AqEquilibriaList(electrolyteRxn,8) = 1
					CYCLE
				END IF
				
				!! Place Electrolyte Index in the Ion Grid
				AqEquilibriaList(electrolyteRxn,1) = FindChem(TRIM(Token),AqPhase)
			END DO

			IF (AqEquilibriaList(ElectrolyteRxn,1) .EQ. -1) CALL SetWaterEquilibriumIndex (ElectrolyteRxn)

			!! If there are more associated elements in the pipe then get confused
			IF (LEN_TRIM(StripToken(AssociatedForm)) .NE. 0) THEN
				CALL WARN ("  Warning!! The following Electrolyte: ("//TRIM(StripToken(AssociatedForm))//	&
					   ") Appeared to be Extra in Line... "//"  >>"//TRIM(FullLine)//&
					   "<<  Perhaps you "// &
					   "specified individual water molecules multiply?  If so, the program "//&
					   " doesn't understand it; you should use a stoicheometric coefficient.  Ignoring and Continuing...")

			END IF
			
			
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Now Parse Individual Data Elements !!

			!! 1. -- Flag for aqueous electrolyte (0) or solid electrolyte (1)
			CALL GetToken(CurrLine, ";", Token)
			IF (STR2INT(Token) .EQ. 0) THEN
				AqEquilibriaList(electrolyteRxn,6) = 0.0
			ELSE IF (STR2INT(Token) .EQ. 1) THEN
				AqEquilibriaList(electrolyteRxn,6) = 1.0
			ELSE
				CALL ERROR("Bad Input for the Aq/Solid flag in the following line "// &
						   "of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			END IF	
			!WRITE(*,*) "Electrolyte Reaction #", electrolyteRxn, " Flag: ", AqEquilibriaList(electrolyteRxn,6)

			!! 2. -- Empty
			AqEquilibriaList(ElectrolyteRxn,7) = -99

			!! 3. -- Third is the Reference Temperature Equilibrium Coefficient
			CALL GetToken(CurrLine, ";", Token)
			IF ((TRIM(Token) .EQ. "infinite" .OR. &
				TRIM(Token) .EQ. "INFINITE" .OR. &
				TRIM(Token) .EQ. "Infinite" .OR. &
				TRIM(Token) .EQ. "Inf" .OR. &
				TRIM(Token) .EQ. "INF" .OR. &
				TRIM(Token) .EQ. "inf") .OR. &
				(Metastable .AND. AqEquilibriaList(ElectrolyteRxn,6) .EQ. 1.0)) THEN
				AqEquilibriaList(ElectrolyteRxn,9) = -1.

			ELSE
				IF (.NOT. ISREAL(Token)) &
				CALL ERROR("Bad Input for the Equilibrium Coefficient at 298 K in the following line "// &
						   "of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
				AqEquilibriaList(ElectrolyteRxn,9) = STR2REAL(TRIM(Token))
			END IF

			!! 4. -- Fourth is the first thermodynamic composite quantity for temperature dependence of K298
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the First Coefficient of Temperature Dependence of the Equilibrium Coefficient "// &
					   "in the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,10) = STR2REAL(TRIM(Token))

			!! 5. -- Fifth is the second thermodynamic composite quantity for temperature dependence of K298
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Second Coefficient of Temperature Dependence of the Equilibrium Coefficient in "// &
					   "the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,11) = STR2REAL(TRIM(Token))

			!! 6. -- Sixth is the First Density Parameter (units of density, g / cm^3)
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Second Coefficient of Temperature Dependence of the Equilibrium Coefficient "// &
					   "in the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,12) = STR2REAL(TRIM(Token)) * grams / cm / cm / cm

			!! 7. -- Seventh is the Second Density Parameter (units of inverse molality, kg / moles)
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Second Coefficient of Temperature Dependence of the Equilibrium Coefficient in "// &
			           "the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,13) = STR2REAL(TRIM(Token)) !! / moles  -- probably molality should be counted 
																		!! in moles (duh) but right now it isn't reported as such, 
																		!! or at least this factor is lost somewhere.
			!! 8. -- Eighth is the Surface Tension Surface Excess at Saturation (mol / m2)
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Surface Excess of the Equilibrium Coefficient in the following line of "// &
			           "AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,14) = STR2REAL(TRIM(Token)) * moles / m / m 

			!! 9. -- Ninth is the Equilibrium Adsorption Coefficient 
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Adsorption Equilibrium Constantof Temperature Dependence of the Equilibrium "// &
					   "Coefficient in the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,15) = STR2REAL(TRIM(Token))

			!! 10. -- Tenth is the DRH at 298K
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Adsorption Equilibrium Constantof Temperature Dependence of the Equilibrium "// &
					   "Coefficient in the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,20) = STR2REAL(TRIM(Token))

			!! 11. -- Eleventh is the temperature dependence of DRH
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for the Adsorption Equilibrium Constantof Temperature Dependence of the Equilibrium "// &
					   "Coefficient in the following line of AqEquilibriumReacions.in: >>"//TRIM(FullLine)//"<<")
			AqEquilibriaList(ElectrolyteRxn,21) = STR2REAL(TRIM(Token))
		
		END DO

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! We don't want the user to double specify a reaction,	  !!
		!! so check to make sure that they're not and err if they !!
		!! are.	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO ElectrolyteRxn = 1, HowManyAqEqReactions !!!!!!!!!!!!!!!!

			IF (Transcribe .AND. AqEquilibriaList(ElectrolyteRxn,9) .EQ. -1.) THEN
				CALL Transcript ("") ; CALL Transcript(TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(ElectrolyteRxn,1))))//  &
				" is set to condense directly into its dissociated form, and will remain fully dissociated while in solution.")
			END IF

			DO I = ElectrolyteRxn+1, HowManyAqEqReactions
				IF (AqEquilibriaList(ElectrolyteRxn,1) .EQ. AqEquilibriaList(I,1) .AND.		&
					AqEquilibriaList(ElectrolyteRxn,2) .EQ. AqEquilibriaList(I,2) .AND.		&
					AqEquilibriaList(ElectrolyteRxn,3) .EQ. AqEquilibriaList(I,3))			&
					CALL ERROR ("It looks like the user double specied the following reaction: <<"//TRIM(AqPhaseChemicalNames(  &
					            INT(AqEquilibriaList(I,1))))//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4)))//" "//             &
								TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR(AqEquilibriaList(I,5)))//" "// &
								TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))//">> as specified in AqEquilibriumReactions.in.  "//  &
								"This is a problem worthy of stopping the program.")
			END DO
		END DO

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Warn if think we are underdefining the possible electrolytes    !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! It is more complicated than this, however, since this could not !!
		!! fail when it should.  This assumes that we will have only 1st   !!
		!! order electrolytes...										   !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		IF (HowManyAqEqReactions .LT. HowManyAqCations * HowManyAqAnions) &
		CALL WARN("MELAM believes that there are more combinations (forming electrolytes) of cation and anions possible than "// &
		          "you have defined equilibrium reactions for.  You defined: "//TRIM(int2str(HowManyAqCations))// &
				  " Cations and "//TRIM(int2str(HowManyAqAnions))//" Anions, leading me to believe that there were at least "// &
				  TRIM(int2str(HowManyAqCations*HowManyAqAnions))//" (C*A) electrolytes possible.  Yet you only defined "// &
				  trim(int2str(HowManyAqEqReactions))//" equilibria.")


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check Mass and Charge Balances !!
		DO I = 1, HowManyAqEqReactions !!!!!

			!! Mass Balance
			IF (AqEquilibriaList(I,1) .GT. 0) THEN
				II = AqMolecularMass(INT(AqEquilibriaList(I,1))) + AqMolecularMass(1) * AqEquilibriaList(I,8)
			ELSE  !! Here is a water only reaction
				II = AqMolecularMass(1) * AqEquilibriaList(I,8)

				!! Tag the water-only dissociation reaction
				IF (AqWaterDissociationReaction .GT. 0)			&
					CALL ERROR ("You specified two water-only aqueous dissociation reactions, which is a little weird.  "// &
								"MELAM can't handle this.  Sorry.")

				!! A water-only dissociation cannot be infinitely dissociating.  
				!! This would be a (big) problem.
				IF (AqEquilibriaList(I,9) .LT. 0)				&
					CALL ERROR("You asked that a water-only aqueous dissociation reaction infinitely dissociate.  If you do "// &
					           "this, all of the water would be ionic, which would really, really suck.  Don't do this.  Fix "// &
							   "it in AqEquilibriumReactions.in.")

				AqWaterDissociationReaction = I
			END IF
			
			IF (AqEquilibriaList(I,22) .NE. 0 .AND. AqEquilibriaList(I,24) .EQ. 0) THEN
				!!(Levitocite)
				JJ = AqEquilibriaList(I,4) * AqCationMass(INT(AqEquilibriaList(I,2))) + AqEquilibriaList(I,5) * &
					AqAnionMass(INT(AqEquilibriaList(I,3))) + AqEquilibriaList(I,23) * AqAnionMass(INT(AqEquilibriaList(I,22)))
			ELSE IF (AqEquilibriaList(I,22) .EQ. 0 .AND. AqEquilibriaList(I,24) .NE. 0) THEN
				!!Hydrates
					JJ = AqEquilibriaList(I,4) * AqCationMass(INT(AqEquilibriaList(I,2))) + AqEquilibriaList(I,5) * &
					AqAnionMass(INT(AqEquilibriaList(I,3))) + AqEquilibriaList(I,25) * AqMolecularMass(INT(AqEquilibriaList(I,24)))			
			ELSE
				!!Simple salts
				JJ = AqEquilibriaList(I,4) * AqCationMass(INT(AqEquilibriaList(I,2))) + AqEquilibriaList(I,5) * &
					AqAnionMass(INT(AqEquilibriaList(I,3)))
			END IF		 
			
			!! Err if More Egregious
			IF (II/JJ .GT. 1.001 .OR. II/JJ .LT. 0.999) THEN
				!write(*,*) II, JJ
				CALL ERROR("Mass More than 0.1% Out of Balance for the Reaction "//TRIM(AqPhaseChemicalNames(  &
							INT(AqEquilibriaList(I,1))))//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" "//   &
							TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" & "//TRIM(REAL2STR(AqEquilibriaList(I,5),1))// &
							" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))//" Aqueous Phase Equilibrium Reacion.  This is "// &
							"more than is acceptable; please rectify.")
			END IF

			!! Warn if Only Slightly Out of Balance
			IF (II/JJ .GT. 1.0001 .OR. II/JJ .LT. 0.9999) THEN
				CALL WARN("Mass More than 0.01% Out of Balance (but less than 0.5% Out of Balance) for the "//  &
						  TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))// &
						  " "//TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" & "//TRIM(REAL2STR(AqEquilibriaList(I,5),1))//" "// &
						  TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))//  &
						  " Aqueous Phase Equilibrium Reacion.  Make sure you are aware of this.")
			END IF

			!! Calculate the net charge of the reaction
			IF (AqEquilibriaList(I,22) .NE. 0) THEN
				!!Note that this formula assumes that only 1 cation and 2 anions are allowed! (Levitocite)
				II = AqEquilibriaList(I,4) * AqCationCharge(INT(AqEquilibriaList(I,2))) + AqEquilibriaList(I,5) * &
					 AqAnionCharge(INT(AqEquilibriaList(I,3))) + AqEquilibriaList(I,23) * AqAnionCharge(INT(AqEquilibriaList(I,22)))
			ELSE
				!Simple salts and their hydrates
				II = AqEquilibriaList(I,4) * AqCationCharge(INT(AqEquilibriaList(I,2))) + AqEquilibriaList(I,5) * &
					 AqAnionCharge(INT(AqEquilibriaList(I,3)))
			END IF
			IF (AqEquilibriaList(I,1) .GT. 0) II = II - AqChemicalCharge(INT(AqEquilibriaList(I,1)))

			!! Charge Balance
			IF (II .NE. 0) THEN
				CALL ERROR ("Charge out of balance for the reaction: <<"//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1)))) &
							//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4)))//" "//TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//  &
							" + "//TRIM(REAL2STR(AqEquilibriaList(I,5)))//" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))//       &
							">> as specified in AqEquilibriumReactions.in.  This is a problem worthy of stopping the program.")
			END IF
		END DO

		!! Report if all is copacetic
		! CMB: Comment this out right now
		!CALL TRANSCRIPT("")
		!CALL TRANSCRIPT("Mass and Charge Balance of Equilibrium Reactions OK")
		!CALL TRANSCRIPT("")

		CLOSE(FH)
		CALL ReturnFileHandle(FH)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Now Open and Deal with the Kusik Meissner Information
		!!
		!! Thermodynamic information need not be available for each
		!! reaction, as long as a suitable composite may be found.
		!! This is why the files are separated.
	    OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'AerosolThermoKusikMeissner.in', STATUS='OLD')

		!! SCSCSCSC
		IF (SCAFFOLDING) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("__About_to_Load_Kusik_Meissner_Data__")
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("__Counting aqueous phase equilibrium reactions...")
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Count the Number of Reactions    !!
		!! (making this a two pass process) !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!! Initialize the index.  A single electrolyte may be doubly counted in this list
		!! if it may dissociate into multiple sets of ions.  
		!! Equally a single set of ions may be doubly counted, meaning two electrolytes 
		!! dissociating into 
		!! So it is a # of Equilibrium Rxns, really.
		HowManyKMParams = 0
		
6		CurrLine = GetLine(FH)
		IF (TRIM(CurrLine) .NE. EOF) THEN
			HowManyKMParams = HowManyKMParams + 1
			GOTO 6
		END IF

		REWIND FH

		!! SCSCSCSC
		IF (SCAFFOLDING) THEN
			CALL TRANSCRIPT("Counted "//TRIM(INT2STR(HowManyKMParams))//" Kusik-Meissner Specifications")
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Make a Tally List of the Equilibria So Don't Have to Check Every Mixture !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This list has three data values:											!!
		!!		1). The Aq Chemical Index of the Electrolyte						!!
		!!		2). The Aq Cation Index of the Cation								!!
		!!		3). The Aq Anion Index of the Anion									!!
		ALLOCATE (KMRef(HowManyKMParams,8), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of KMRef could not proceed in SetAqPhaseEquilibriumReactions()")

		
		!!! -- ION-ELECTROLYTE MAP and EQUILIBRIUM DATA -- !!!
		!!! --
		!!! -- DECODING KEY FOR KMRef (a,b)
		!!! -- 
		!!! --   (i,1) : Electrolyte  : this is -1 if water is the only dissociator (e.g., H2O <=> H+ & OH-)
		!!! --   (i,2) : Cation 
		!!! --   (i,3) : Anion
		!!! --   (i,4) : Cation Stoicheometric Coefficient
		!!! --   (i,5) : Anion  Stoicheometric Coefficient
		!!! --   (i,6) : Kusik-Meissner Coefficient
		!!! --   (i,7) : Kusik-Meissner Temperature Dependence 
		!!! --   (i,8) : Number of H2O Molecules Needed on Left Hand Side of Dissociation

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!! READ THE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Grab the next equilibrium reaction and process the inputs !!
		DO ElectrolyteRxn = 1, HowManyKMParams

			CurrLine = GetLine(FH)

			IF (TRIM(CurrLine) .EQ. EOF) &
			CALL ERROR ("Somehow we miscounted the number of lines in AerosolThermoKusikMeissner.in " // &
					    "SetAqPhaseEquilibriumReactions ().  This is kind of boffo.  Shouldn't happen.  Sorry.")

			FullLine = CurrLine

			!! Parse to separate Electrolyte, Ions, and Parameters
			CALL GetToken(CurrLine, "<=>", AssociatedForm)
			CALL GetToken(CurrLine, ";",   DissociatedForm)

			!! Record the Reaction
			! CMB: comment this out right now
			!CALL TRANSCRIPT(TRIM(AssociatedForm)//" <=> "//TRIM(DissociatedForm))

			!! Find the Index Point in terms of Anions & Cations
			K = 0; L = 0; MM = 0; N = 0
			DO I = 1, 2

				J = 0

				CALL GetToken(DissociatedForm, "&", Ion)
				IF (len_trim(Ion) == 0) EXIT

				!! Check Token for Leading Factor
				CALL GetToken(Ion," ",Token)

				IF (IsReal(Token)) THEN
					J     = STR2REAL(Token)			! Store
					Token = StripToken(Ion)			! Push to Ion Name
				ELSE
					J     = 1
					Ion   = TRIM(Ion) // " " // TRIM(Token)
				END IF

				IonI = FindIon(TRIM(Ion))

				!! Cation
				IF (IonI(1) .GT. 0) THEN
					IF (K .GT. 0) CALL ERROR("Two Cations have been specified in the dissociation reaction of "//  &
											 TRIM(AssociatedForm)//" in AerosolThermoKusikMeissner.in.  This is not allowed.")
					K = J
					L = IonI(2)
				END IF

				!! Anion
				IF (IonI(1) .LT. 0) THEN
					IF (MM .GT. 0) CALL ERROR("Two Anions have been specified in the dissociation reaction of "//  &
											  TRIM(AssociatedForm)//" in AerosolThermoKusikMeissner.in.  This is not allowed.")
					MM = J
					N = IonI(2)
				END IF

				!! Error -- HYDRATION NOT YET IMPLEMENTED
				IF (IonI(1) .EQ. 0) &
				CALL ERROR("Hydration Reactions Not Yet Implemented, can't form the hydrate in AerosolThermoKusikMeissner.in.  "// &
						   "Please comment it out for the time being.")

			END DO

			!! If there are more dissociated elements in the pipe then get confused
			IF (LEN_TRIM(StripToken(DissociatedForm)) .NE. 0) &
			CALL WARN ("  Warning!! The following Ion: ("//TRIM(StripToken(DissociatedForm))//") Appeared to be Extra in Line... "// &
					   "  >>"//TRIM(FullLine)//"<<  Ignoring and Continuing...")

			!! Fill the Stoicheometric Coefficients and identifications
			KMRef(electrolyteRxn,2) = REAL(L)
			KMRef(electrolyteRxn,3) = REAL(N)
			KMRef(electrolyteRxn,4) = REAL(K)
			KMRef(electrolyteRxn,5) = REAL(MM)

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Now Parse the Electrolyte !!
			K = 0; MM = 0
			KMRef(electrolyteRxn,1) = -1
			KMRef(electrolyteRxn,8) = 0  ! Assume no water unless find some

			DO I = 1, 2
				J = 0
				CALL GetToken(AssociatedForm, "&", Ion)

				IF (len_trim(Ion) == 0) EXIT

				!! Check Token for Leading Factor
60				CALL GetToken(Ion," ",Token)

				!! Should be water if Stoich > 1
				IF (IsReal(Token)) THEN

					!! error check
					IF (STR2REAL(Token) .EQ. 1) &
					CALL ERROR (" Don't specify a stoicheometric coefficient of one!  Simply leave it out.  Found such a "// &
								"gross thing in: >>"//TRIM(FullLine)//"<< ")

					KMRef(electrolyteRxn,8) = STR2REAL(Token)
					Token = StripToken(Ion)				! Push to Chemical Name

					!! Only Water Can Have Multiple Instances
					IF (FindChem(TRIM(Token),AqPhase) .EQ. 1) THEN
						GOTO 60
					ELSE
						CALL ERROR ("Indicated a Stoicheometric Coefficient Above One for the Electrolyte in the following "// &
									"line.  This is not acceptable; we can only take one of the electrolyte.  Sorry.  The "//  &
									"offending line: "//TRIM(FullLine))
					END IF
				ELSE
					Ion = TRIM(Ion) // " " // TRIM(Token)
				END IF

				!! Check to see if it is water
				IF (FindChem(TRIM(Token),AqPhase) .EQ. 1) THEN
					KMRef(electrolyteRxn,8) = 1
					CYCLE
				END IF

				!! Place Electrolyte Index in the Ion Grid
				KMRef(electrolyteRxn,1) = FindChem(TRIM(Token),AqPhase)
			END DO

			IF (KMRef(ElectrolyteRxn,1) .EQ. -1) &
			CALL ERROR("I don't believe you that you think you have Kusik-Meissner parameters for the <<H2O <=> OH- & H+>> "// &
			           "reaction.  You say you do, but I frankly don't believe it.  You'll have to edit the code to get this "// &
					   "to work.  Start in SetAqPhaseEquilibriumReactions().")

			!! If there are more associated elements in the pipe then get confused
			IF (LEN_TRIM(StripToken(AssociatedForm)) .NE. 0) &
			CALL WARN ("  Warning!! The following Electrolyte: ("//TRIM(StripToken(AssociatedForm))// &
					   ") Appeared to be Extra in Line... "//&
					   "  >>"//TRIM(FullLine)//"<<  Perhaps you specified individual water molecules multiply?  If so, "//&
					   "the program doesn't understand it; you should use a stoicheometric coefficient.  "// &
					   "Ignoring and Continuing...")

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!! Now Parse Individual Data Elements !!

			!! 1. -- First is the basic Kusik-Meissner Coefficient
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) &
			CALL ERROR("Bad Input for Kusik-Meissner Coefficient in the following line of AerosolThermoKusikMeissner.in: >>"// &
					   TRIM(FullLine)//"<<")
			KMRef(electrolyteRxn,6) = STR2REAL(TRIM(Token))

			!! 2. -- Second is the Kusik-Meissner Temperature Coefficient
			CALL GetToken(CurrLine, ";", Token)
			IF (.NOT. ISREAL(Token)) CALL ERROR("Bad Input for Kusik-Meissner Temperature Dependence in the following line "// &
												"of AerosolThermoKusikMeissner.in: >>"//TRIM(FullLine)//"<<")
			KMRef(ElectrolyteRxn,7) = STR2REAL(TRIM(Token))

		END DO


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! We don't want the user to double specify a reaction,	  !!
		!! so check to make sure that they're not and err if they !!
		!! are.	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO ElectrolyteRxn = 1, HowManyKMParams !!!!!!!!!!!!!!!!!!!!!

			DO I = ElectrolyteRxn+1, HowManyKMParams
				IF (KMRef(ElectrolyteRxn,1) .EQ. KMRef(I,1) .AND.		&
					KMRef(ElectrolyteRxn,2) .EQ. KMRef(I,2) .AND.		&
					KMRef(ElectrolyteRxn,3) .EQ. KMRef(I,3))			&
					CALL ERROR ("It looks like the user double specied the following reaction: <<"//	&
								TRIM(AqPhaseChemicalNames(INT(KMRef(I,1))))//" <=> "//TRIM(REAL2STR(KMRef(I,4),1))// &
								" "//TRIM(AqCationNames(INT(KMRef(I,2))))//" & "//TRIM(REAL2STR(KMRef(I,5),1))//" "// &
								TRIM(AqAnionNames(INT(KMRef(I,3))))//">> as specified in AerosolThermoKusikMeissner.in.  "// &
								"This is a problem worthy of stopping the program.")
			END DO
		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check Mass and Charge Balances !!
		DO I = 1, HowManyKMParams !!!!!

			!! Mass Balance
			IF (KMRef(I,1) .GT. 0) THEN
				II = AqMolecularMass(INT(KMRef(I,1))) + AqMolecularMass(1) * KMRef(I,8)
			ELSE
				II = AqMolecularMass(1) * KMRef(I,8)
			END IF

			JJ = KMRef(I,4) * AqCationMass(INT(KMRef(I,2))) + KMRef(I,5) * AqAnionMass(INT(KMRef(I,3)))

			!! Err if More Egregious
			IF (II/JJ .GT. 1.001 .OR. II/JJ .LT. 0.999) THEN
				!WRITE(*,*) II, JJ
				CALL ERROR("Mass More than 0.1% Out of Balance for the Thermodynamic Specification for Reaction "// &
						   TRIM(AqPhaseChemicalNames(INT(KMRef(I,1))))//" <=> "//TRIM(REAL2STR(KMRef(I,4),1))//" "//     &
						   TRIM(AqCationNames(INT(KMRef(I,2))))//" & "//TRIM(REAL2STR(KMRef(I,5),1))//" "//              &
						   TRIM(AqAnionNames(INT(KMRef(I,3))))//" Aqueous Phase Equilibrium Reacion in "// &
						   "AerosolThermoKusikMeissner.in.  This is more than is acceptable; please rectify.")
			END IF

			!! Warn if Only Slightly Out of Balance
			IF (II/JJ .GT. 1.0001 .OR. II/JJ .LT. 0.9999) THEN
				!WRITE(*,*) II, JJ
				CALL WARN("Mass More than 0.01% Out of Balance for the Thermodynamic Specification for Reaction "// &
						  TRIM(AqPhaseChemicalNames(INT(KMRef(I,1))))//" <=> "//TRIM(REAL2STR(KMRef(I,4),1))//" "//      &
						  TRIM(AqCationNames(INT(KMRef(I,2))))//" & "//TRIM(REAL2STR(KMRef(I,5),1))//" "//               &
						  TRIM(AqAnionNames(INT(KMRef(I,3))))//" in AerosolThermoKusikMeissner.in.  This is more than "// &
						  "is acceptable; please rectify.")
			END IF

			!! Calculate the net charge of the reaction
			II = KMRef(I,4) * AqCationCharge(INT(KMRef(I,2))) + KMRef(I,5) * AqAnionCharge(INT(KMRef(I,3)))
			IF (KMRef(I,1) .GT. 0) II = II - AqChemicalCharge(INT(KMRef(I,1)))

			!! Charge Balance
			IF (II .NE. 0) THEN
				CALL ERROR ("Charge out of balance for the reaction: <<"//TRIM(AqPhaseChemicalNames(INT(KMRef(I,1))))// &
							" <=> "//TRIM(REAL2STR(KMRef(I,4)))//" "//TRIM(AqCationNames(INT(KMRef(I,2))))//" + "//          &
							TRIM(REAL2STR(KMRef(I,5)))//" "//TRIM(AqAnionNames(INT(KMRef(I,3))))//">> as specified in "//    &
							"AerosolThermoKusikMeissner.in.  This is a problem worthy of stopping the program.")
			END IF
		END DO

		!! Report if all is copacetic
		! cmb: comment this out right now
		!CALL TRANSCRIPT("")
		!CALL TRANSCRIPT("Mass and Charge Balance of Equilibrium Reactions OK")
		!CALL TRANSCRIPT("")


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Some of the Aqueous Electrolyte Dissociation Reactions !!
		!! will have direct matches in the Kusik Meissner chart.  !!
		!! Let us identify those and flag the others.			  !!
		DO I = 1, HowManyAqEqReactions	!!!!!!!!!!!!!!!!!!!!!!!!!!!!

			!! Initialize the indices of this thing...
			AqEquilibriaList(I,16) = -1
			AqEquilibriaList(I,17) = 0.
			AqEquilibriaList(I,18) = 0.

			DO J = 1, HowManyKMParams
				IF ((AqEquilibriaList(I,1) .EQ. KMRef(J,1)) .AND.    &
				    (AqEquilibriaList(I,2) .EQ. KMRef(J,2)) .AND.    &
				    (AqEquilibriaList(I,3) .EQ. KMRef(J,3)) .AND.    &
				    (AqEquilibriaList(I,4) .EQ. KMRef(J,4)) .AND.    &
				    (AqEquilibriaList(I,5) .EQ. KMRef(J,5)) .AND.    &
				    (AqEquilibriaList(I,8) .EQ. KMRef(J,8)))		 &
				  AqEquilibriaList(I,16) = J + (AqEquilibriaList(I,4) + AqEquilibriaList(I,5)) / 10000.
			END DO
		END DO


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! For bisulfate and other related systems of partial dissociation, MELAM allow    !!
		!! ions to dissociate into ions.  (e.g., HSO4- <=> H+ + SO4--).  You may either    !!
		!! specify the data directly, in which case it will be picked up in the last loop, !!
		!! or by specifying it indirectly (e.g., H2SO4 <=> H+ + HSO4-), which will then    !!
		!! be the thermodynamics used to equilibrate the first reaction.  This has the     !!
		!! advantage of allowing mixed activity coefficients to be calculated for it       !!
		!! (necessary for the inclusion of NaHSO4, NH4HSO4, etc).						   !!
		!!																				   !!
		!! It is negative to let KusikMeissner() know that it is to take the post-mixing   !!
		!! value																		   !!
		DO I = 1, HowManyAqEqReactions  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			IF (AqEquilibriaList(I,1)  .GT. HowManyAqChems .AND.	&   ! Ionic
				AqEquilibriaList(I,16) .LT. 0) THEN						! Unmatched

			DO J = 1, HowManyAqEqReactions

				!! Cationic
				IF (AqEquilibriaList(I,1) .LE. HowManyAqChems + HowManyAqCations) THEN

					IF ((AqEquilibriaList(I,1) .EQ. AqEquilibriaList(J,2)+HowManyAqChems) .AND.    &
						(AqEquilibriaList(I,3) .EQ. AqEquilibriaList(J,3)) .AND.    &
						(1					   .EQ. AqEquilibriaList(J,4)) .AND.    &
						(AqEquilibriaList(I,5) .EQ. AqEquilibriaList(J,5)) .AND.    &
						(AqEquilibriaList(I,8) .EQ. AqEquilibriaList(J,8)) .AND.	&
						(AqEquilibriaList(I,8) .EQ. 0).AND.	(I .NE. J)) THEN

					  AqEquilibriaList(I,16) = AqEquilibriaList(J,16)
					  AqEquilibriaList(I,17) = AqEquilibriaList(J,17)
					  AqEquilibriaList(I,18) = AqEquilibriaList(J,18)
					  AqEquilibriaList(I,19) = J
					  CALL TRANSCRIPT ("*********************************************************")
					  CALL TRANSCRIPT ("For the dissociation of "//TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))// &
									   ", MELAM is using the post-mixing thermodynamics for the dissociation of the fully "//  &
									   "associated species into that ion, following San Martini (2004), Kim et al (1993) "//   &
									   "and others.")
					  CALL TRANSCRIPT ("*********************************************************")
					END IF

				!! Anionic
				ELSE
					IF ((AqEquilibriaList(I,1) .EQ. AqEquilibriaList(J,3)+HowManyAqChems+HowManyAqCations) .AND.    &
						(AqEquilibriaList(I,2) .EQ. AqEquilibriaList(J,2)) .AND.    &
						(AqEquilibriaList(I,4) .EQ. AqEquilibriaList(J,4)) .AND.    &
						(1					   .EQ. AqEquilibriaList(J,5)) .AND.    &
						(AqEquilibriaList(I,8) .EQ. AqEquilibriaList(J,8)) .AND.	&
						(AqEquilibriaList(I,8) .EQ. 0) .AND. (I .NE. J))	THEN

					  AqEquilibriaList(I,16) = AqEquilibriaList(J,16)
					  AqEquilibriaList(I,17) = AqEquilibriaList(J,17)
					  AqEquilibriaList(I,18) = AqEquilibriaList(J,18)
					  AqEquilibriaList(I,19) = J
					  CALL TRANSCRIPT ("*********************************************************")
					  CALL TRANSCRIPT ("For the dissociation of "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,1))-  &
									   HowManyAqChems-HowManyAqCations))//", MELAM is using the post-mixing "// &
									   "thermodynamics for the dissociation of the fully associated species "// &
									   "into that ion, following San Martini (2004), Kim et al (1993) and others.")
					  CALL TRANSCRIPT ("*********************************************************")
					END IF
				END IF
			END DO

			END IF
		END DO


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Consider the mean activity for reactions for which we don't  !!
		!! have appropriate Kusik-Meissner parameters.				    !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Now consider the effective mean activity coefficients for    !!
		!! reactions for which no mean activity Kusik-Meissner params   !!
		!! have been specified.  In this case, MELAM will search for a  !!
		!! means way to form it from a composite of mean activities from!!
		!! other reactions.  A good example of what I mean by this is   !!
		!! the following reaction:										!!
		!!															    !!
		!!			H2O <=> OH- + H+								    !!
		!!																!!
		!! For which no KM data exists (or could exist, I think).  In   !!
		!! that case, an effective mean activity can be formed using    !!
		!! those of three other reactions:								!!
		!!																!!
		!!																!!
		!!  g(H+) g(OH-)     g(H+) g(Cl-)  g(K) g(OH-)					!!
		!!  ------------- = --------------------------------		    !!
		!!       1                  g(K) g(Cl-)							!!
		!!																!!
		!!    (g(HCl))^2 (g(KOH))^2										!!
		!! = --------------------------  in the parlance of this model	!!
		!!          (g(KCl))^2											!!
		!!																!!
		!! Mean activity coefficients for all three species are			!!
		!! available.  Something that allows this type of formulation   !!
		!! must exist in order for the water dissociation reaction to   !!
		!! be used.  An error will be reported if MELAM can't find an   !!
		!! appropriate thing to use.									!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO I = 1, HowManyAqEqReactions

			!! If this one is well identified, then go on
99			IF (AqEquilibriaList(I,16) .NE. -1) CYCLE

			!! Check for Levitocite
			IF (AqEquilibriaList(I,22) .NE. 0) THEN
				IF (AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))) .EQ. "(NH4)3H(SO4)2") THEN
					DO J = 1, HowManyAqEqReactions
						IF (AqEquilibriaList(J,1) .GT. 0) THEN !Exclude water as electrolyte
							
							!Find Ammonium Sulfate
							IF(AqPhaseChemicalNames(INT(AqEquilibriaList(J,1))) .EQ. "(NH4)2SO4") THEN
								AqEquilibriaList(I,16) = J
							END IF

							!Find Sulfuric Acid (H2SO4 <=> 2 H+ & SO4--)
							IF(AqPhaseChemicalNames(INT(AqEquilibriaList(J,1))) .EQ. "H2SO4") THEN
								AqEquilibriaList(I,17) = J
							END IF							
						END IF
					END DO
					CALL TRANSCRIPT ("*********************************************************")
					CALL TRANSCRIPT ("The levitocite (NH4)3H(SO4)2 reaction ")
					CALL TRANSCRIPT ("has no thermo defined for it.  But MELAM will use the formula from ")
					CALL TRANSCRIPT ("Kim et al, 1993:")
					CALL TRANSCRIPT ("  g((NH4)3H(SO4)2)^4 = ("//  &
												 "g("//"(NH4)2SO4 "//")^3"//" * "//   &
 												 "g("//"H2SO4"//")")
					CALL TRANSCRIPT ("*********************************************************")
					CALL TRANSCRIPT ("")
				
				END IF	
			END IF
			
			!! Try to match all of the pieces with other reactions
			!! for which we have thermodynamic data
			DO J = 1, HowManyKMParams
				
				!! Do the Cations Match?
				IF (AqEquilibriaList(I,2) .EQ. KMRef(J,2)) THEN

					!! Find one for which the Anions match as well
					DO K = 1, HowManyKMParams
						
						!! We should now have one for which the cations 
						!! and anions match.  Then look for something to 
						!! cancel the extra terms appropriately
						IF (AqEquilibriaList(I,3) .EQ. KMRef(K,3)) THEN
								
							DO L = 1, HowManyKMParams

								!! This should be the cancellation term if 
								!! its cations and anions cancel those of the
								!! other necessary species and all of the 
								!! stoicheometry matches
								IF ((KMRef(J,3)			   .EQ. KMRef(L,3)) .AND.  & !Anion Match
									(KMRef(K,2)			   .EQ. KMRef(L,2)) .AND.  & !Cation Match
								    (KMRef(J,5)			   .EQ. KMRef(L,5)) .AND.  & !Anion Stoich matches
									(KMRef(K,4)			   .EQ. KMRef(L,4)) .AND.  & !Cation Stioch matches
									(AqEquilibriaList(I,4) .EQ. KMRef(J,4)) .AND.  &
									(AqEquilibriaList(I,5) .EQ. KMRef(K,5))) THEN

									CALL TRANSCRIPT ("")

									!! Get the name of the electrolyte
									IF (AqEquilibriaList(I,1) .GT. 0) THEN
										IF (AqEquilibriaList(I,1) .LE. HowManyAqChems) THEN
											CurrLine = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))
										ELSE IF (AqEquilibriaList(I,1) .LE. HowManyAqChems+HowManyAqCations) THEN
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))
										ELSE
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations))
										END IF
									ELSE
										IF (AqEquilibriaList(I,8) .GT. 0) CurrLine = "H2O"
									END IF

									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("The electrolyte reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(CurrLine)//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" " &
												 //TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR  &
												 (AqEquilibriaList(I,5),1))//" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3)))))
									CALL TRANSCRIPT ("Has no thermo defined for it.  But MELAM found an ")
									CALL TRANSCRIPT ("appropriate derivative form:")
									CALL TRANSCRIPT ("  g("//TRIM(CurrLine)//")^"//TRIM(REAL2STR((AqEquilibriaList(I,4)+  &
								                 AqEquilibriaList(I,5)),1))//" = ("//  &
												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(J,1))))//")^"//TRIM(REAL2STR((KMRef(J,4)+ &
												 KMRef(J,5)),1))//" * "//   &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(K,1))))//")^"//TRIM(REAL2STR((KMRef(K,4)+ &
												 KMRef(K,5)),1))//") / "//  &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(L,1))))//")^"//TRIM(REAL2STR((KMRef(L,4)+ &
													 KMRef(L,5)),1)))
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("")

									!! Now tell AqEqList how to think of the thermodynamics going forward
									!! The format is: Index + exponent / 10000 
									AqEquilibriaList(I,16) = J + (KMRef(J,4)+KMRef(J,5)) / 10000. 
									AqEquilibriaList(I,17) = L + (KMRef(L,4)+KMRef(L,5)) / 10000. 
									AqEquilibriaList(I,18) = K + (KMRef(K,4)+KMRef(K,5)) / 10000. 


									GOTO 99  !! skip out of all of the loops
								END IF
								
								!For cases where Cation stoichiometry is off by a factor of 2
								IF ((KMRef(J,3)			   .EQ. KMRef(L,3)) .AND.  & !Anions Match
									(KMRef(K,2)			   .EQ. KMRef(L,2)) .AND. &  !Cations match
									(2*KMRef(K,4)			   .EQ. KMRef(L,4)) .AND.  & !Cation Stioch twice original
									(KMRef(J,5)			   .EQ. KMRef(L,5)) .AND. &
									(AqEquilibriaList(I,22) .EQ. 0)) THEN   !Not Levitocite
									
									CALL TRANSCRIPT ("")

									!! Get the name of the electrolyte
									IF (AqEquilibriaList(I,1) .GT. 0) THEN
										IF (AqEquilibriaList(I,1) .LE. HowManyAqChems) THEN
											CurrLine = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))
										ELSE IF (AqEquilibriaList(I,1) .LE. HowManyAqChems+HowManyAqCations) THEN
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))
										ELSE
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations))
										END IF
									ELSE
										IF (AqEquilibriaList(I,8) .GT. 0) CurrLine = "H2O"
									END IF

									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("The electrolyte reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(CurrLine)//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" " &
												 //TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR  &
												 (AqEquilibriaList(I,5),1))//" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3)))))
									CALL TRANSCRIPT ("Has no thermo defined for it.  But MELAM found an ")
									CALL TRANSCRIPT ("appropriate derivative form:")
									CALL TRANSCRIPT ("  g("//TRIM(CurrLine)//")^"//TRIM(REAL2STR((AqEquilibriaList(I,4)+  &
								                 AqEquilibriaList(I,5)),1))//" = ("//  &
												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(J,1))))//")^"//TRIM(REAL2STR((KMRef(J,4)+ &
												 KMRef(J,5)),1))//" * "//   &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(K,1))))//")^"//TRIM(REAL2STR((KMRef(K,4)+ &
												 KMRef(K,5)),1))//") / "//  &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(L,1))))//")^"//TRIM(REAL2STR((KMRef(L,4)+ &
													 KMRef(L,5)),1)))
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("")
									!WRITE(*,*) 
									
									!! Now tell AqEqList how to think of the thermodynamics going forward
									!! The format is: Index + exponent / 10000 
									AqEquilibriaList(I,16) = J + (KMRef(J,4)+KMRef(J,5)) / 10000. 
									AqEquilibriaList(I,17) = L + (KMRef(L,4)+KMRef(L,5)) / 10000. 
									AqEquilibriaList(I,18) = K + (KMRef(K,4)+KMRef(K,5)) / 10000. 


									GOTO 99  !! skip out of all of the loops
								
								END IF
							
								!For cases where Anion stoichiometry is off by a factor of 2, denominator high
								IF ((KMRef(J,3)			   .EQ. KMRef(L,3)) .AND.  & !Anions Match
									(KMRef(K,2)			   .EQ. KMRef(L,2)) .AND. &  !Cations match
									(KMRef(K,4)			   .EQ. KMRef(L,4)) .AND.  & !Cation Stioch 
									(2*KMRef(J,5)			   .EQ. KMRef(L,5)) .AND. & !Anion Stoich Twice Original
									(AqEquilibriaList(I,22) .EQ. 0)) THEN !Not Levitocite
									
									CALL TRANSCRIPT ("")

									!! Get the name of the electrolyte
									IF (AqEquilibriaList(I,1) .GT. 0) THEN
										IF (AqEquilibriaList(I,1) .LE. HowManyAqChems) THEN
											CurrLine = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))
										ELSE IF (AqEquilibriaList(I,1) .LE. HowManyAqChems+HowManyAqCations) THEN
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))
										ELSE
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations))
										END IF
									ELSE
										IF (AqEquilibriaList(I,8) .GT. 0) CurrLine = "H2O"
									END IF

									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("The electrolyte reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(CurrLine)//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" " &
												 //TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR  &
												 (AqEquilibriaList(I,5),1))//" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3)))))
									CALL TRANSCRIPT ("Has no thermo defined for it.  But MELAM found an ")
									CALL TRANSCRIPT ("appropriate derivative form:")
									CALL TRANSCRIPT ("  g("//TRIM(CurrLine)//")^"//TRIM(REAL2STR((AqEquilibriaList(I,4)+  &
								                 AqEquilibriaList(I,5)),1))//" = ("//  &
												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(J,1))))//")^"//TRIM(REAL2STR((KMRef(J,4)+ &
												 KMRef(J,5)),1))//" * "//   &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(K,1))))//")^"//TRIM(REAL2STR((KMRef(K,4)+ &
												 KMRef(K,5)),1))//") / "//  &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(L,1))))//")^"//TRIM(REAL2STR((KMRef(L,4)+ &
													 KMRef(L,5)),1)))
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("") 
									
									!! Now tell AqEqList how to think of the thermodynamics going forward
									!! The format is: Index + exponent / 10000 
									AqEquilibriaList(I,16) = J + (KMRef(J,4)+KMRef(J,5)) / 10000. 
									AqEquilibriaList(I,17) = L + (KMRef(L,4)+KMRef(L,5)) / 10000. 
									AqEquilibriaList(I,18) = K + (KMRef(K,4)+KMRef(K,5)) / 10000. 


									GOTO 99  !! skip out of all of the loops
								
								END IF
							
								!For cases where Anion and cation stoichiometry is off by a factor of 2, denominator low (i.e. CaCO3)
								IF ((KMRef(J,3)			   .EQ. KMRef(L,3)) .AND.  & !Anions Match
									(KMRef(K,2)			   .EQ. KMRef(L,2)) .AND. &  !Cations match
									(KMRef(K,4)			   .EQ. 2*KMRef(L,4)) .AND.  & !Cation Stioch 
									(KMRef(J,5)			   .EQ. 2*KMRef(L,5)) .AND. & !Anion Stoich Half Original
									(AqEquilibriaList(I,22) .EQ. 0) .AND. &    !Not Levitocite
										(AqEquilibriaList(I,24) .EQ. 0) .AND. & !Not CaSO4*2H2O 
										(AqCationNames(INT(AqEquilibriaList(I,2))) .NE. "Mg++")) THEN !Not a Mg salt
									
									!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), AqAnionNames(AqEquilibriaList(I,3)), "Hi!"
									IF (AqPhaseChemicalNames(INT(KMRef(K,1))) .EQ. "HCO3-") CYCLE
									CALL TRANSCRIPT ("")

									!! Get the name of the electrolyte
									IF (AqEquilibriaList(I,1) .GT. 0) THEN
										IF (AqEquilibriaList(I,1) .LE. HowManyAqChems) THEN
											CurrLine = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))
										ELSE IF (AqEquilibriaList(I,1) .LE. HowManyAqChems+HowManyAqCations) THEN
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))
										ELSE
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations))
										END IF
									ELSE
										IF (AqEquilibriaList(I,8) .GT. 0) CurrLine = "H2O"
									END IF

									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("The electrolyte reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(CurrLine)//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" " &
												 //TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR  &
												 (AqEquilibriaList(I,5),1))//" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3)))))
									CALL TRANSCRIPT ("Has no thermo defined for it.  But MELAM found an ")
									CALL TRANSCRIPT ("appropriate derivative form:")
									CALL TRANSCRIPT ("  g("//TRIM(CurrLine)//")^"//TRIM(REAL2STR((AqEquilibriaList(I,4)+  &
								                 AqEquilibriaList(I,5)),1))//" = ("//  &
												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(J,1))))//")^"//TRIM(REAL2STR((KMRef(J,4)+ &
												 KMRef(J,5)),1))//" * "//   &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(K,1))))//")^"//TRIM(REAL2STR((KMRef(K,4)+ &
												 KMRef(K,5)),1))//") / "//  &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(L,1))))//")^"//TRIM(REAL2STR(2*(KMRef(L,4)+ &
													 KMRef(L,5)),1)))
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("") 
									
									!! Now tell AqEqList how to think of the thermodynamics going forward
									!! The format is: Index + exponent / 10000 
									AqEquilibriaList(I,16) = J + (KMRef(J,4)+KMRef(J,5)) / 10000. 
									AqEquilibriaList(I,17) = L + 2*(KMRef(L,4)+KMRef(L,5)) / 10000. 
									AqEquilibriaList(I,18) = K + (KMRef(K,4)+KMRef(K,5)) / 10000. 
									
									GOTO 99  !! skip out of all of the loops
								
								END IF
								
									
								!For cases where Anion stoichiometry is off by a factor of 2, denominator low
								IF ((KMRef(J,3)			   .EQ. KMRef(L,3)) .AND.  & !Anions Match
									(KMRef(K,2)			   .EQ. KMRef(L,2)) .AND. &  !Cations match
									(KMRef(K,4)			   .EQ. KMRef(L,4)) .AND.  & !Cation Stioch 
									(KMRef(J,5)			   .EQ. 2*KMRef(L,5)) .AND. & !Anion Stoich Half Original
									(AqEquilibriaList(I,22) .EQ. 0) .AND. &    !Not Levitocite
										(AqEquilibriaList(I,24) .EQ. 0) .AND. & !Not CaSO4*2H2O 
										(AqCationNames(INT(AqEquilibriaList(I,2))) .NE. "Mg++")) THEN !Not a Mg salt
									
									!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), AqAnionNames(AqEquilibriaList(I,3)), "Hi!"
									IF (AqPhaseChemicalNames(INT(KMRef(K,1))) .EQ. "HCO3-") CYCLE
									CALL TRANSCRIPT ("")

									!! Get the name of the electrolyte
									IF (AqEquilibriaList(I,1) .GT. 0) THEN
										IF (AqEquilibriaList(I,1) .LE. HowManyAqChems) THEN
											CurrLine = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))
										ELSE IF (AqEquilibriaList(I,1) .LE. HowManyAqChems+HowManyAqCations) THEN
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))
										ELSE
											CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations))
										END IF
									ELSE
										IF (AqEquilibriaList(I,8) .GT. 0) CurrLine = "H2O"
									END IF

									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("The electrolyte reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(CurrLine)//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" " &
												 //TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR  &
												 (AqEquilibriaList(I,5),1))//" "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3)))))
									CALL TRANSCRIPT ("Has no thermo defined for it.  But MELAM found an ")
									CALL TRANSCRIPT ("appropriate derivative form:")
									CALL TRANSCRIPT ("  g("//TRIM(CurrLine)//")^"//TRIM(REAL2STR((AqEquilibriaList(I,4)+  &
								                 AqEquilibriaList(I,5)),1))//" = ("//  &
												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(J,1))))//")^"//TRIM(REAL2STR((KMRef(J,4)+ &
												 KMRef(J,5)),1))//" * "//   &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(K,1))))//")^"//TRIM(REAL2STR(2*(KMRef(K,4)+ &
												 KMRef(K,5)),1))//") / "//  &
 												 "g("//TRIM(AqPhaseChemicalNames(INT(KMRef(L,1))))//")^"//TRIM(REAL2STR(2*(KMRef(L,4)+ &
													 KMRef(L,5)),1)))
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("") 
									
									!! Now tell AqEqList how to think of the thermodynamics going forward
									!! The format is: Index + exponent / 10000 
									AqEquilibriaList(I,16) = J + (KMRef(J,4)+KMRef(J,5)) / 10000. 
									AqEquilibriaList(I,17) = L + 2*(KMRef(L,4)+KMRef(L,5)) / 10000. 
									AqEquilibriaList(I,18) = K + 2*(KMRef(K,4)+KMRef(K,5)) / 10000. 
									
									GOTO 99  !! skip out of all of the loops
								
								END IF
							
							END DO
						END IF
					END DO
				END IF
			END DO          
		END DO


		!! Check to see if any of the reactions are unaccounted for
		J = 0
		DO I = 1, HowManyAqEqReactions
			IF (AqEquilibriaList(I,16) .LE. 0) THEN

				!! Get the name of the electrolyte
				IF (AqEquilibriaList(I,1) .GT. 0) THEN
					IF (AqEquilibriaList(I,1) .LE. HowManyAqChems) THEN
						CurrLine = TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))
					ELSE IF (AqEquilibriaList(I,1) .LE. HowManyAqChems+HowManyAqCations) THEN
						CurrLine = TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))
					ELSE
						CurrLine = TRIM(AqAnionNames(INT(AqEquilibriaList(I,1))-HowManyAqChems-HowManyAqCations))
					END IF
				ELSE
					IF (AqEquilibriaList(I,8) .GT. 0) CurrLine = "H2O"
				END IF

				CALL WARN ("The program could not determine an appropriate way to form composite activity coefficients for "//  &
						   "the reaction <<"//TRIM(CurrLine)//" <=> "//TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" "//           &
						   TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR(AqEquilibriaList(I,5),1))//" "//    &
						   TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))//">> ")
				J = 1  ! flag for an error after listing all problems
			END IF

		END DO

		IF (J .EQ. 1) CALL ERROR("Failure to rectify cited thermodyanic issues means we cannot proceed.")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Consider the mean activity for partially dissociating rxns   !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! This is brought about by the sulfate example.  The correct	!!
		!! treatment, within the Kusik-Meissner framework, for			!!
		!! including bisulfate is to define each with respect to the	!!
		!! most fully dissociated state, e.g.:							!!
		!!																!!
		!!    H2SO4 <=> 2 H+ + SO4--									!!
		!!    HSO4- <=>   H+ + SO4--									!!
		!!																!!
		!! The equilibrium reactions depend on the ratio of activity	!!
		!! coeffs.  The denominator is always for a single electrolyte,	!!
		!! but, when that electrolyte is bisulfate it gets complicated. !!
		!! We want to rewrite this to be acceptable within the Kusik-   !!
		!! Meissner formulation, which only finds mean activity coeffs, !!
		!! not those for single ions.									!!
		!! For example, for the dissociation of bisulfate the ratio is:	!!
		!!																!!
		!!  gamma(H+) gamma(SO4--)    gamma(H+) gamma(H+) gamma(SO4--)  !!
		!!  ----------------------  = --------------------------------  !!
		!!         gamma(HSO4-)              gamma(H+) gamma(HSO4-)     !!
		!!																!!
		!!    (gamma(H2SO4))^3											!!
		!! = ------------------  in the parlance of this model			!!
		!!    (gamma(HSO4-))^2											!!
		!!																!!
		!! This requires several things:								!!
		!!																!!
		!! 1. that all reactions be written with respect to the most    !!
		!!    dissociated form.											!!
		!!																!!
		!! 2. that all mean activity coefficients for things charged    !!
		!!    electrolytes be defined to include implicit ions of the   !!
		!!    same type that dissociates from the given reaction.		!!
		!!	  (this is a limitation of the model, of course, as perhaps !!
		!!    other implicit cations would be appropriate.				!!
		!!																!!
		!! 3. that we identify and store the identity of the implicit   !!
		!!    and appropriate mean activity coefficients.  That we		!!
		!!    start with the first ratio and move to the third.			!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		DO I = 1, HowManyAqEqReactions

			!! If the index is higher than the number of 
			!! uncharged chemicals present, then it is an ion,
			!! which is what we're searching for.
			IF (AqEquilibriaList(I,1) .GT. HowManyAqChems) THEN

				!! Is it a Cation...?
				IF (AqEquilibriaList(I,1) .LE. HowManyAqChems + HowManyAqCations) THEN

					CALL ERROR ("Found an aqueous dissociation reaction in which the cation "//  &
							    TRIM(AqCationNames(INT(AqEquilibriaList(I,1))-HowManyAqChems))//      &
								" acts as an electrolyte.  This program is not prepared to handle "// &
								"this situation, and any charged electrolytes are taken to be anions.  "// &
								"If you cannot write this reaction w.r.t. the anion, then you are S.O.L..  "// &
								"You would have to rewrite everything to allow specification of implicit and "// &
								"explicit terms of the associated side of electrolytic dissociation reactions, "// &
								"and we only allow a single species.  Sorry.")

				!! ...or an Anion?
				ELSE

					CALL TRANSCRIPT ("*********************************************************")
					CALL TRANSCRIPT ("Identified that "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,1))-  &
								     HowManyAqChems-HowManyAqCations))//" acts as an electrolyte in solution... "// &
									 "and we are attempting to identify how to deal with the reaction.")
					CALL TRANSCRIPT ("")
					CALL TRANSCRIPT ("Trying to identify a matching reaction of the charged electrolyte...")
					CALL TRANSCRIPT ("")


					!! Find the matching reaction that will bring it to an undissociated form.
					!! The dissociated-side ANIONS will match...
					DO J = 1, HowManyAqEqReactions

						IF (I .EQ. J) CYCLE

						!! Does the reaction fit the criteria?
						IF (AqEquilibriaList(I,2) .EQ. AqEquilibriaList(J,2) .AND.  &
							AqEquilibriaList(I,3) .EQ. AqEquilibriaList(J,3) .AND.  &
							AqEquilibriaList(J,1) .LE. HowManyAqChems) THEN

							!! It cannot go through another partially dissociating reaction
							IF (AqEquilibriaList(J,17) .GT. 0) CYCLE


							CALL TRANSCRIPT ("We identified that Reaction (b) appears to allow a pathway")
							CALL TRANSCRIPT ("for reaction (a) to reach a fully associated form:")
							CALL TRANSCRIPT ("")
							CALL TRANSCRIPT ("  (a) "//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//" <=> "//  &
											 TRIM(REAL2STR(AqEquilibriaList(I,4),1))//" "//TRIM(AqCationNames(  &
											 INT(AqEquilibriaList(I,2))))//" + "//TRIM(REAL2STR(AqEquilibriaList(I,5),1))// &
											 " "//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3)))))
							CALL TRANSCRIPT ("  (b) "//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(J,1))))//" <=> "// &
											 TRIM(REAL2STR(AqEquilibriaList(J,4),1))//" "//TRIM(AqCationNames(  &
											 INT(AqEquilibriaList(J,2))))//" + "//TRIM(REAL2STR(AqEquilibriaList(J,5),1))//  &
											 " "//TRIM(AqAnionNames(INT(AqEquilibriaList(J,3)))))
							CALL TRANSCRIPT ("")
							CALL TRANSCRIPT ("Leading us to make the following assumptions regarding ")
							CALL TRANSCRIPT ("the use of mean activity coefficients (g) in equilibration:")
							CALL TRANSCRIPT ("")
							CALL TRANSCRIPT ("(g("//TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//")^"//TRIM(INT2STR( &
											  INT(AqEquilibriaList(I,4))))//    &
							                 " * g("//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))//")^"// &
											 TRIM(INT2STR(INT(AqEquilibriaList(I,5))))//  &
											 ") / g("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//") ")
							CALL TRANSCRIPT ("= (g("//TRIM(AqCationNames(INT(AqEquilibriaList(J,2))))//")^"// &
											 TRIM(INT2STR(INT(AqEquilibriaList(J,4))))//    &
							                 " * g("//TRIM(AqAnionNames(INT(AqEquilibriaList(J,3))))//")^"// &
											 TRIM(INT2STR(INT(AqEquilibriaList(J,5))))//  &
											 ") / (g("//TRIM(AqCationNames(INT(AqEquilibriaList(I,2))))//"))^"// &
											 TRIM(INT2STR(INT(AqEquilibriaList(J,4)-AqEquilibriaList(I,4))))//    &
											 " * g("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//")) ")
							CALL TRANSCRIPT ("= g("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(J,1))))//")^"// &
											 TRIM(INT2STR(INT(AqEquilibriaList(J,4)+AqEquilibriaList(J,5))))//   &
											 " / g("//TRIM(AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))))//")^"// &
										     TRIM(INT2STR(INT(AqEquilibriaList(I,4)+AqEquilibriaList(I,5)))))
							CALL TRANSCRIPT ("")
							CALL TRANSCRIPT ("Which will be used for reaction (a)")
							CALL TRANSCRIPT ("*********************************************************")

							!! Store the information about how the mean activity coefficients are 
							!! used in calculating equilibrium
							AqEquilibriaList(I,17) = AqEquilibriaList(I,16) !+ (AqEquilibriaList(I,4) + &
																			! AqEquilibriaList(I,5)) / 10000.
							AqEquilibriaList(I,16) = J + (AqEquilibriaList(J,4) + AqEquilibriaList(J,5)) / 10000.
						END IF
					END DO

					!! If you don't find an appropriate path up, err
					IF (AqEquilibriaList(I,16) .EQ. I .AND. AqEquilibriaList(I,17) .EQ. 0) &
						CALL ERROR("Did not find an appropriate path for this electrolyte ("//TRIM(AqPhaseChemicalNames( &
									INT(AqEquilibriaList(I,1))))//")to reach a fully associated state, which must be reached "// &
									"through a common dissociated form ("//TRIM(AqAnionNames(INT(AqEquilibriaList(I,3))))// &
									") using the same anion.  Please correct AqEquilibriumReactions.in")

				END IF
			END IF
		END DO

		!! SCSCSCSC -- Announce Arrival
		IF (Scaffolding) THEN
			CALL TRANSCRIPT("")
			CALL Transcript(">>Exiting SetAqPhaseEquilibriumReactions()<<")
		END IF

		CLOSE(FH)
		CALL ReturnFileHandle(FH)
		
		RETURN
	END SUBROUTINE SetAqPhaseEquilibriumReactions
