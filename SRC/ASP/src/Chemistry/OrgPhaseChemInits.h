!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! OrgPhaseChemInits.h
!! Reads in and stores data on the aqueous and hydrophobic phase 
!! organic chemicals.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 11/01  2006   Matt Alvarado     Added Refractive Index Data               !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/27  2012   Matt Alvarado     Removed Refractive Index Data, now        !!
!!                                  hardcoded in SUBROUTINE                  !!
!!                                  ShellRefIndAndRad(InParticle) in file    !!
!!                                  ParticleAttributes.h                     !!
!! 11/09  2012   Matt Alvarado     Changed BC density to 1.8 g/cm3 to match  !!
!!                                  density assumed for SP2 measurements     !!
!!                                  in ARCTAS                                !!
!! 01/24  2013   Matt Alvarado     Changed SetOrgPhaseChemistry to read      !!
!!                                  activity coefficient parameters.         !!
!! 01/30  2013   Matt Alvarado     Changed SetHydrophilicOrgPhaseChemistry   !!
!!                                  to read Kappa and                        !!
!!                                  activity coefficient parameters.         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:  !!
!! 1. SUBROUTINE SetOrgPhaseChemistry ()			!!
!! 2. SUBROUTINE SetHydrophilicOrgPhaseChemsitry ()		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

! **********
! Modified by CMB (AER): Make this play nice with gfortran (long lines)
! **********

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Org Phase Chemical Information from 'OrgPhaseChems.in',		    !!
!! an input deck, defines the appropriate chemical arrays in this module,   !!
!! and calls the appropriate eulerian grid subroutines to establich the	    !!
!! chemical fields							    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetOrgPhaseChemistry ()

	USE ModelParameters, ONLY : SetHowManyOrgChems, cm, &
				    grams, moles, InputDeckSubDir, &
                                    TranscriptFH,   &
				    AqChemScale, Atmospheres
	USE InfrastructuralCode

	IMPLICIT NONE

	!! Local Infrastructural Variables
	INTEGER	                :: i, j, k,			 &
		                   allocation_error,             &
		                   index1, index2, currindex
	CHARACTER(len=512)	:: CurrLine, Token, Name
	CHARACTER (len=12)	:: Numbers = "0123456789+-"

	!! These Local Variables Arrive From the Input Deck
	INTEGER			:: NumbChems
	REAL*8, ALLOCATABLE	:: TempMassStorage(:)
	CHARACTER(len=512), ALLOCATABLE  :: TempNameStorage(:)

	!! Transcribe?  Use Scaffolding Code?
	LOGICAL			:: Transcribe
	LOGICAL, PARAMETER      :: Scaffolding =  .FALSE. 

	!! And there is a file handle to be used
	integer	:: FH
	FH  = GetFileHandle()
	HowManyOrgChems		  = 1  
        !These are set to 1 because BC is hard coded
	HowManyEvolveOrgChems  = 1

	!! SCSCSCSC -- Announce Arrival
	IF (Scaffolding) CALL Transcript(">>Entering SetOrgPhaseChemistry()<<")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There are several Organic Phase Input Decks Required for Initialization: !!
!!   1). OrgPhaseChems.in (List of Chemicals Present in Organic Phase.)	    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! In this subroutine, we will create listings of all of the chemicals 
	OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'OrgPhaseChems.in', STATUS='OLD')	! the incoming chemicals
	    
	!! Make a first pass through the file to determine the number of 
	!! chemicals specified therein.
	!! Need this number early to allocate all of the arrays.
10	CurrLine = GetLine(FH)
	IF (CurrLine .NE. EOF) THEN
		HowManyOrgChems = HowManyOrgChems+1
		GOTO 10
	END IF

	!! Back up for a second pass to parse the chemical traits
	REWIND(FH)

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript(">>About to allocate<<")

	!! Transcribe
	IF (TranscriptFH > 0) THEN
		Transcribe = .TRUE.
	ELSE
		Transcribe = .FALSE.
	END IF
	Transcribe = .false.	! cmb add, don't want this anymore
	
	IF (Transcribe) THEN
		CALL Transcript("")
		CALL Transcript("__Initializing_Organic_Phase_Chemicals_("//TRIM(INT2STR(HowManyOrgChems))//"_Species)_")
	END IF

	!! Tell the ModelParameter Module How Many Org Chems There Are
	CALL SetHowManyOrgChems(HowManyOrgChems)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Now ALLOCATE each of the arrays that are size HowManyOrgChems!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!! -- up to 15 character CHEMICAL NAME which should be identical to that -- !!!
	!!! -- given in other input decks if it is to be matched across phases    -- !!!
	ALLOCATE (OrgPhaseChemicalNames(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgPhaseChemicalNames could not proceed in SetOrgPhaseChemistry()")

	!!! -- MOLECULAR MASS (grams / mole) -- !!!
	ALLOCATE (OrgMolecularMass(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgMolecularMass could not proceed in SetOrgChemistryParams()")

	!!! -- DENSITY (grams / cm3) -- !!!
	ALLOCATE (OrgDensity(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgDensity could not proceed in SetOrgChemistryParams()")

        !! ACTIVITY COEFF. FLAG
        !!         if 0. - number than follows is the index (starting at 1) 
        !!                 of the compound in unifacparams.h and the 
        !!                 subroutine UNIDRIV
        !!         if 1. - number that follows is a constant 
        !!                 organic activity coefficient (usually 1.0)
        !!         NOTE: It is safest to set these flags the same for all compounds
        !!         considered. Otherwise the UNIFAC compunds are treated as their 
        !!         own solution that doesn't inteact with the constant compounds.
	ALLOCATE (OrgActivityFlag(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgDensity could not proceed in SetOrgChemistryParams()")

	ALLOCATE (OrgActivityIndex(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgDensity could not proceed in SetOrgChemistryParams()")
			
	!! SCSCSCSCSC
	IF (Scaffolding) CALL Transcript (">>Importing Org Phase Species Info<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! -- IMPORT OTHER ORG PHASE SPECIES INFORMATION FROM THE DECK -- !!
	!!  Evolving Species Go First, Non-Evolving Last		  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	index1 = 2	! Evolving Species are Placed with this index
	index2 = HowManyOrgChems  
        ! Non Evolving Species are Placed with this index
		
	!Hard code values for BC
	OrgPhaseChemicalNames(1)	= "BC"
	OrgMolecularMass(1)		= 200.0* grams / moles
	OrgDensity(1)			= 1.8 * grams / cm / cm / cm 
        OrgActivityFlag(1)              = 1.0
        OrgActivityIndex(1)             = 1.0
        !ARCTAS SP2 observations used this density

	DO i = 2, HowManyOrgChems

		CurrLine = GetLine(FH)

		!! Tokenize the input line !!
		!! The token order is:
		!! ---	1. Chemical Name (15 Characters Max)
		CALL GetToken(CurrLine, ";", Name)
                !Hold until know where to place it

		!! Don't want the user to specify water
		IF (TRIM(Name) .EQ. "H2O" .OR.		&
			TRIM(Name) .EQ. "h2o" .OR.	&
			TRIM(Name) .EQ. "water" .OR.	&
			TRIM(Name) .EQ. "Water" .OR.	&
			TRIM(Name) .EQ. "WATER")        &
			CALL ERROR ("You specified "//TRIM(Name)//" in OrgPhaseChems.in.  This model hardcodes in the values for water."// &
			            "  Please don't specify anything about it and use 'H2O' when referring to it in chemical mechanisms.")

		!! ---	2. Evolve Chemical or Keep Constant?
		CALL GetToken(CurrLine,";",Token) 

		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In OrgPhaseChems.in, the Evolve or Not flag appears to be missing on "//  &
		TRIM(OrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")

		k = STR2INT (Token)

		IF (k == 1) THEN
			HowManyEvolveOrgChems = HowManyEvolveOrgChems + 1
			CurrIndex             = index1
			index1		      = index1 + 1
		ELSE IF (k == 0) THEN
			CurrIndex             = index2
			index2		      = index2 - 1
		ELSE
			CALL ERROR("Don't Understand "//Trim(StripToken(Token))//" as a value for the Constant or evolve flag"//	&
			"For Chemical "//Trim(OrgPhaseChemicalNames(CurrIndex))//" in OrgPhaseChems.in. ")
		END IF

		!! Now Place Name in Correct Slot
		OrgPhaseChemicalNames(CurrIndex) = Trim(StripToken(Name))

		!! --- 3. Molecular Mass (g / mole)
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In OrgPhaseChems.in, the molecular mass appears to be missing on "//   &
		                            TRIM(OrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		OrgMolecularMass(CurrIndex) = STR2REAL(Token) * grams / moles

		!! --- 4. Density (g / cm3)
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In OrgPhaseChems.in, the density appears to be missing on "//   &
						TRIM(OrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		OrgDensity(CurrIndex) = STR2REAL(Token) * grams / cm / cm / cm

		!! --- 5. Activity Coefficent Flag
                !!         if 0. - number than follows is the index (starting at 1) 
                !!                 of the compound in unifacparams.h and the 
                !!                 subroutine UNIDRIV
                !!         if 1. - number that follows is a constant 
                !!                 organic activity coefficient (usually 1.0)
                !!         NOTE: It is safest to set these flags the same for all compounds
                !!         considered. Otherwise the UNIFAC compunds are treated as their 
                !!         own solution that doesn't inteact with the constant compounds.
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In OrgPhaseChems.in, the acxtivity coeff. flag appears to be missing on "//   &
						TRIM(OrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		OrgActivityFlag(CurrIndex) = STR2REAL(Token)

		!! --- 6. Activity Coefficent Index
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In OrgPhaseChems.in, the acxtivity coeff. index appears to be missing on "//   &
						TRIM(OrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		OrgActivityIndex(CurrIndex) = STR2REAL(Token)			
			
		!! Transcript
		IF (Transcribe) CALL Transcript(OrgPhaseChemicalNames(CurrIndex))
	END DO

	!! SCSCSCSC -- Announce Departure
	IF (Scaffolding) CALL Transcript(">>Exiting SetOrgPhaseChemistry()<<")
		
	CLOSE (FH)
	CALL ReturnFileHandle(FH)

	RETURN
END SUBROUTINE SetOrgPhaseChemistry 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Hydrophilic Org Phase Chemical Information from 'AqOrgPhaseChems.in',!!
!! an input deck, defines the appropriate chemical arrays in this module,    !!
!! and calls the appropriate eulerian grid subroutines to establich the	     !!
!! chemical fields							     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetHydrophilicOrgPhaseChemistry ()

	USE ModelParameters, ONLY : SetHowManyAqOrgChems, cm,      &
				    grams, moles, InputDeckSubDir, &
                                    TranscriptFH,                  &
				    AqChemScale, Atmospheres
	USE InfrastructuralCode

	IMPLICIT NONE

	!! Local Infrastructural Variables
	INTEGER			:: i, j, k,			 &
			           allocation_error,             &
			           index1, index2, currindex
	CHARACTER(len=512)	:: CurrLine, Token, Name
	CHARACTER (len=12)	:: Numbers = "0123456789+-"

	!! These Local Variables Arrive From the Input Deck
	INTEGER			:: NumbChems
	REAL*8, ALLOCATABLE	:: TempMassStorage(:)
	CHARACTER(len=512), ALLOCATABLE  :: TempNameStorage(:)

	!! Transcribe?  Use Scaffolding Code?
	LOGICAL			   :: Transcribe
	LOGICAL, PARAMETER :: Scaffolding =  .FALSE. 

	!! And there is a file handle to be used
	integer	:: FH
	FH  = GetFileHandle()
	HowManyAqOrgChems	 = 0 
	HowManyEvolveAqOrgChems  = 0

	!! SCSCSCSC -- Announce Arrival
	IF (Scaffolding) CALL Transcript(">>Entering SetHyrdophilicOrgPhaseChemistry()<<")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There are several Organic Phase Input Decks Required for Initialization: !!
!!   1). AqOrgPhaseChems.in (List of Chemicals Present in Organic Phase.)   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! In this subroutine, we will create listings of all of the chemicals 
	OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'HydrophilicOrgChems.in', STATUS='OLD')	! the incoming chemicals
	    
	!! Make a first pass through the file to determine the number of 
	!! chemicals specified therein.
	!! Need this number early to allocate all of the arrays.
10	CurrLine = GetLine(FH)
	IF (CurrLine .NE. EOF) THEN
		HowManyAqOrgChems = HowManyAqOrgChems+1
		GOTO 10
	END IF

	!! Back up for a second pass to parse the chemical traits
	REWIND(FH)

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript(">>About to allocate<<")

	!! Transcribe
	IF (TranscriptFH > 0) THEN
		Transcribe = .TRUE.
	ELSE
		Transcribe = .FALSE.
	END IF
	Transcribe = .false.	! cmb
	
	IF (Transcribe) THEN
		CALL Transcript("")
!
! 29 September 2016 CMB: This is a bug.  This initialization statement is reporting the 
!                        incorrect # of hydrophilic organic species.
		!CALL Transcript("__Initializing_Hydrophilic_Organic_Chemicals_("//TRIM(INT2STR(HowManyOrgChems))//"_Species)_")
		CALL Transcript("__Initializing_Hydrophilic_Organic_Chemicals_("//TRIM(INT2STR(HowManyAqOrgChems))//"_Species)_")
	END IF

	!! Tell the ModelParameter Module How Many Aq Org Chems There Are
	CALL SetHowManyAqOrgChems(HowManyAqOrgChems)

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Now ALLOCATE each of the arrays that are size HowManyAqOrgChems!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! ** CMB (AER): Made a few of the conditionals into compound blocks; gfortran didn't like them
        !               before

	!!! -- up to 15 character CHEMICAL NAME which should be identical to that -- !!!
	!!! -- given in other input decks if it is to be matched across phases    -- !!!
	ALLOCATE (AqOrgPhaseChemicalNames(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) then
            CALL ERROR("Allocation of AqOrgPhaseChemicalNames could not proceed in SetHydrophilicOrgPhaseChemistry()")
        endif
	
	!!! -- MOLECULAR MASS (grams / mole) -- !!!
	ALLOCATE (AqOrgMolecularMass(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) then
            CALL ERROR("Allocation of AqOrgMolecularMass could not proceed in SetHydrophilicOrgChemistryParams()")
        endif

	!!! -- DENSITY (grams / cm3) -- !!!
	ALLOCATE (AqOrgDensity(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) then
            CALL ERROR("Allocation of AqOrgDensity could not proceed in SetHydrophilicOrgChemistryParams()")
        endif
		
	!!! -- Number of Carbons per molecule -- !!!
	ALLOCATE (AqOrgNumCarbon(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) then
            CALL ERROR("Allocation of AqOrgNumCarbon could not proceed in SetHydrophilicOrgChemistryParams()")
        endif

	!!! -- Kappa parameter to calculate associated water -- !!!
	ALLOCATE (Kappa(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) then
            CALL ERROR("Allocation of Kappa could not proceed in SetHydrophilicOrgChemistryParams()")
        endif

       !!! -- Raoult's law activity coefficient at infinte dilution -- !!!
       !!! -- Used to correct UNIFAC activity coefficient to Henry's law -- !!!
	ALLOCATE (GammaInf(HowManyAqOrgChems), stat = allocation_error)
	if (allocation_error > 0) then
            CALL ERROR("Allocation of Kappa could not proceed in SetHydrophilicOrgChemistryParams()")
        endif

        !! ACTIVITY COEFF. FLAG
        !!         if 0. - number than follows is the index (starting at 1) 
        !!                 of the compound in unifacparams.h and the 
        !!                 subroutine UNIDRIV
        !!         if 1. - number that follows is a constant 
        !!                 organic activity coefficient (usually 1.0)
        !!         NOTE: It is safest to set these flags the same for all compounds
        !!         considered. Otherwise the UNIFAC compunds are treated as their 
        !!         own solution that doesn't inteact with the constant compounds.
	ALLOCATE (AqOrgActivityFlag(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgDensity could not proceed in SetOrgChemistryParams()")

	ALLOCATE (AqOrgActivityIndex(HowManyOrgChems), stat = allocation_error)
	if (allocation_error > 0) CALL ERROR("Allocation of OrgDensity could not proceed in SetOrgChemistryParams()")
		
	!! SCSCSCSCSC
	IF (Scaffolding) CALL Transcript (">>Importing Hydrophilic Org Phase Species Info<<")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! -- IMPORT OTHER HYDROPHILIC ORG SPECIES INFORMATION FROM THE DECK !!
	!!  Evolving Species Go First, Non-Evolving Last		     !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	index1 = 1	! Evolving Species are Placed with this index
	index2 = HowManyAqOrgChems  ! Non Evolving Species are Placed with this index
		
	DO i = 1, HowManyAqOrgChems

		CurrLine = GetLine(FH)

		!! Tokenize the input line !!
		!! The token order is:
		!! ---	1. Chemical Name (15 Characters Max)
		CALL GetToken(CurrLine, ";", Name)			
                ! Hold until know where to place it

		!! Don't want the user to specify water
		IF (TRIM(Name) .EQ. "H2O" .OR.		&
			TRIM(Name) .EQ. "h2o" .OR.		&
			TRIM(Name) .EQ. "water" .OR.	&
			TRIM(Name) .EQ. "Water" .OR.	&
			TRIM(Name) .EQ. "WATER")        &
			CALL ERROR ("You specified "//TRIM(Name)//" in HydrophilicOrgChems.in.  This model hardcodes in the values for water."// &
			            "  Please don't specify anything about it and use 'H2O' when referring to it in chemical mechanisms.")

		!! ---	2. Evolve Chemical or Keep Constant?
		CALL GetToken(CurrLine,";",Token) 

		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In HydrophilicOrgChems.in, the Evolve or Not flag appears to be missing on "//  &
					    TRIM(AqOrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")

		k = STR2INT (Token)

		IF (k == 1) THEN
			HowManyEvolveAqOrgChems = HowManyEvolveAqOrgChems + 1
			CurrIndex               = index1
			index1			= index1 + 1
		ELSE IF (k == 0) THEN
			CurrIndex               = index2
			index2			= index2 - 1
		ELSE
			CALL ERROR("Don't Understand "//Trim(StripToken(Token))//" as a value for the Constant or evolve flag"//	&
				   "For Chemical "//Trim(AqOrgPhaseChemicalNames(CurrIndex))//" in HydrophilicOrgChems.in. ")
		END IF

		!! Now Place Name in Correct Slot
		AqOrgPhaseChemicalNames(CurrIndex) = Trim(StripToken(Name))

		!! --- 3. Molecular Mass (g / mole)
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In HydrophilicOrgChems.in, the molecular mass appears to be missing on "//   &
						TRIM(AqOrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		AqOrgMolecularMass(CurrIndex) = STR2REAL(Token) * grams / moles

		!! --- 4. Density (g / cm3)
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In HydrophilicOrgChems.in, the density appears to be missing on "//   &
						TRIM(AqOrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		AqOrgDensity(CurrIndex) = STR2REAL(Token) * grams / cm / cm / cm

		!! --- 5. Number of carbons per molecule
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) CALL ERROR("In HydrophilicOrgChems.in, the number of carbons appears to be missing on "//   &
						TRIM(AqOrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough items on that line?")
		AqOrgNumCarbon(CurrIndex) = STR2INT(Token)

		!! --- 6. Kappa (see Petters and Kreidenweis, ACP, 7, 1961, 2007.)
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) then
                    CALL ERROR("In HydrophilicOrgChems.in, the Kappa parameter appears to be missing on "//   &
			       TRIM(AqOrgPhaseChemicalNames(CurrIndex))//".  Perhaps you don't have enough "// &
                               "items on that line?")
                endif
		Kappa(CurrIndex) = STR2REAL(Token)

		!! --- 7. Raoult's Law Activity Coefficient at Infinite Dilution
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) then
                    CALL ERROR("In HydrophilicOrgChems.in, the activity coeff. at infinite dilution appears "//&
                               "to be missing on "//TRIM(AqOrgPhaseChemicalNames(CurrIndex))//&
                               ".  Perhaps you don't have enough items on that line?")
                endif
		GammaInf(CurrIndex) = STR2REAL(Token)
			
		!! --- 8. Activity Coefficent Flag
                !!         if 0. - number than follows is the index (starting at 1) 
                !!                 of the compound in unifacparams.h and the 
                !!                 subroutine UNIDRIV
                !!         if 1. - number that follows is a constant 
                !!                 organic activity coefficient FOR HENRY's LAW (usually 1.0)
                !!         NOTE: It is safest to set these flags the same for all compounds
                !!         considered. Otherwise the UNIFAC compunds are treated as their 
                !!         own solution that doesn't inteact with the constant compounds.
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) then
                    CALL ERROR("In HydrophilicOrgChems.in, the activity coeff. flag appears to be missing on "//   &
				TRIM(AqOrgPhaseChemicalNames(CurrIndex))//&
                                ".  Perhaps you don't have enough items on that line?")
                endif
		AqOrgActivityFlag(CurrIndex) = STR2REAL(Token)

		!! --- 9. Activity Coefficent Index
		CALL GetToken(CurrLine,";",Token)
		IF (LEN_TRIM(Token) .EQ. 0) then
                    CALL ERROR("In HydrophilicOrgChems.in, the activity coeff. index appears to be missing on "//   &
			       TRIM(AqOrgPhaseChemicalNames(CurrIndex))//&
                               ".  Perhaps you don't have enough items on that line?")
                endif
		AqOrgActivityIndex(CurrIndex) = STR2REAL(Token)

		!! Transcript
		IF (Transcribe) CALL Transcript(AqOrgPhaseChemicalNames(CurrIndex))
	END DO

	!! SCSCSCSC -- Announce Departure
	IF (Scaffolding) CALL Transcript(">>Exiting SetHydrophilicOrgPhaseChemistry()<<")
	
	CLOSE (FH)
	CALL ReturnFileHandle(FH)

	RETURN
END SUBROUTINE SetHydrophilicOrgPhaseChemistry 
