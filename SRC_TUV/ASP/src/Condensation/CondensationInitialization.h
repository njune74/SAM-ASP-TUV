!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! CondensationInitialization.h
!! Reads dissolution and gas-aerosol thermodynamic information from	
!! input decks, and stores it in the appropriate arrays.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!								             !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 06/27  2007   Matt Alvarado     Initially zero out DissolutionData        !!
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 01/23  2013   Matt Alvarado     Added flag to SetOrganicDissolution to    !!
!!                                  decide between two vaopr pressure methods!!
!! 01/24  2013   Matt Alvarado     Added flag to SetOrganicDissolution to    !!
!!                                  decide between two org act. coeff.methods!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	    !!
!! 1. SUBROUTINE SetAllDissolution()                                !!
!! 2. SUBROUTINE SetDissolution ()                                  !!
!! 3. SUBROUTINE SetOrganicDissolution()                            !!
!! 4. SUBROUTINE SetHydrophilicOrganicDissolution()		    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!! This allows all the intitialization routines to be called
!! with a single command
SUBROUTINE SetAllDissolution()

	IMPLICIT NONE

	CALL SetDissolution
	CALL SetOrganicDissolution
	CALL SetHydrophilicOrganicDissolution

END SUBROUTINE SetAllDissolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Aq Phase Chemical Information from				    !!
!! input decks, defines the appropriate chemical arrays in this module,	    !!
!! and calls the appropriate eulerian grid subroutines to establich the	    !!
!! chemical fields							    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetDissolution ()

		USE ModelParameters, ONLY : AqChemScale, Atmospheres,	&
					SetDissolutionEquilibriumOrFlux,&
					DissolutionEquilibriumOrFlux,	&
					ThermoBulkMode,                 &
					TranscriptFH, InputDeckSubDir,  &
                                        DoDissolutionFlag

		USE InfrastructuralCode

		USE Chemistry,   ONLY : GasPhase, AqPhase,		&
					FindChem, FindIon,		&
					GasPhaseChemicalNames,		&
					AqPhaseChemicalNames,		&
					AqCationNames,AqAnionNames,	&
					GasMolecularMass,		&
					AqMolecularMass,		&
					AqCationMass, AqAnionMass,	&
					AqEquilibriaList,		&
					HowManyAqChems,			&
					HowManyAqCations,		&
					HowManyAqEqReactions,		&
					GetAqCharge

		IMPLICIT NONE

		!! Local Infrastructural Variables
		INTEGER	:: I, J, K, allocation_error,IonIndex(2)
		REAL*8	:: MassTally, H2SO4GasIndex, SO4AqIndex, HAqIndex
		CHARACTER(len=512)	:: CurrLine, Token, GasChem1, &
                                           GasChem2, AqIon1, AqIon2
		CHARACTER (len=12)	:: Numbers = "0123456789+-"

		!! Transcribe?  Use Scaffolding Code?
		LOGICAL :: Transcribe, Match
                ! cmb: don't want this reported right now
		!LOGICAL, PARAMETER :: Scaffolding = .TRUE.
		LOGICAL, PARAMETER :: Scaffolding = .false.

		!! And there is a file handle to be used
		INTEGER :: FH
		FH  = GetFileHandle()

		!! SCSCSCSC -- Announce Arrival
		IF (Scaffolding) CALL Transcript(">>Entering SetDissolution()<<")


		!First is alwats water, second is always sulfate cond.
		HowManyDissolutionReactions = 2

		IF (DoDissolutionFlag) THEN
		!! In this subroutine, we will create listings of all of the chemicals 
			OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'Dissolution.in', STATUS='OLD')

			!! Make a first pass through the file to determine the number of 
			!! reactions specified therein.
			!! Need this number early to allocate all of the arrays.
			CurrLine = GetLine(FH)
			DO WHILE (CurrLine .NE. EOF)
				HowManyDissolutionReactions = HowManyDissolutionReactions+1
                                CurrLine = GetLine(FH)			
			END DO

			!! Back up for a second pass to parse the reactions
			REWIND(FH)

		END IF
			

		!! Transcribe
		IF (TranscriptFH > 0) THEN
			Transcribe = .TRUE.
		ELSE
			Transcribe = .FALSE.
		END IF
                Transcribe = .false.  ! cmb
		IF (Transcribe) CALL TRANSCRIPT("")
		IF (Transcribe) CALL TRANSCRIPT("_The_Following_Dissolution_Reactions_Are_Defined_")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Allocation the key DISSOLUTION Arrays				!!
!!									!!
!! The primary one is DissolutionData(i,j), in which i is the		!!
!! reaction and j is some data:						!!
!!							                !!
!! DissolutionData(i,1)  : Index in Gas Phase of Primay Chemical	!!
!! DissolutionData(i,2)  : Index in Gas Phase of Secondary Chemical	!!
!! DissolutionData(i,3)  : Index in Aqueous Phase of Primary Chemical	!!
!! DissolutionData(i,4)  : Index in Aqueous Phase of Secondary Chemical	!!
!! DissolutionData(i,5)  : Henrys Law Coefficient (H_298)		!!
!! DissolutionData(i,6)  : Henrys Law Coefficient (Delta H / R)		!!
!! DissolutionData(i,7)  : Henrys Law Coefficient (c_p / R)		!!
!! DissolutionData(i,8)  : For (i,2) is an ion, Numerator mean activity	!!
!!			   coefficient + Exponent / 10000		!!
!! DissolutionData(i,9)  : For (i,2) is an ion, Denominator mean activity!!
!!			   coefficient + Exponent / 10000		!!
!! DissolutionData(i,10) : Mass Accommodation Coefficient		!!
!! DissolutionData(i,11) : Reaction Type				!!
!! DissolutionData(i,12) : Index of AqEquilibriaList for related Aq reaction!!
!! DissolutionData(i,13) : Stoicheometry of r.h.s. slot 1		!!
!! DissolutionData(i,14) : Stoicheometry of r.h.s. slot 2 for type 2 or     !!
!!                                          l.h.s. slot 2 if  type 4	!!
!! DissolutionData(i,15) : Reaction I is coupled to through gas phase species !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! There are several reactions types:					     !!
!!									     !!
!! 0. Water Condensation.						     !!
!! 1. Direct dissolution of a single species not incorporating dissociation. !!
!! 2. Dissolution of a single species that then dissociates directly.	     !!
!! 4. Dissolution of a species that binds with an ion on the way in.	     !!
!! 6. Dissociation of two gas phase species to form 1 or 2 aqueous species.  !!
!! 7. Condensation of gas-phase H2SO4 to make 2 H+ and SO4--                 !!
!! Type 3 and 5 are reserved for coupled dissolution/dissociation            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		ALLOCATE (DissolutionData(HowManyDissolutionReactions,15), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of NumbReactions could not proceed in SetDissolution()")

		!Zero out DissolutionData to start
		DO I = 1, HowManyDissolutionReactions
			DO J = 1, 15
				DissolutionData(I,J) = 0.
			END DO
		END DO

		!! Hard Code the Water 
		DissolutionData(1,1)  = 1.
		DissolutionData(1,2)  = 0.
		DissolutionData(1,3)  = 1.
		DissolutionData(1,4)  = 0.
		DissolutionData(1,5)  = 0.
		DissolutionData(1,6)  = 0.
		DissolutionData(1,7)  = 0.
		DissolutionData(1,8)  = 0.
		DissolutionData(1,9)  = 0.
		DissolutionData(1,10) = 0.45
		DissolutionData(1,11) = 0.
		DissolutionData(1,12) = 0.
		DissolutionData(1,13) = 1.
		DissolutionData(1,14) = 0.
		DissolutionData(1,15) = 0.

		!! Rat out the hard coded water transfer
		IF (TRANSCRIBE) CALL TRANSCRIPT("H2O (v) <=> H2O (l)")

		!! Hard Code the Sulfate 
		H2SO4GasIndex = FindChem("H2SO4", 0)
		SO4AqIndex = FindChem("SO4--", 1)
		HAqIndex = FindChem("H+", 1)
		
		DissolutionData(2,1)  = H2SO4GasIndex
		DissolutionData(2,2)  = 0.
		DissolutionData(2,3)  = HAqIndex
		DissolutionData(2,4)  = SO4AqIndex
		DissolutionData(2,5)  = 0.
		DissolutionData(2,6)  = 0.
		DissolutionData(2,7)  = 0.
		DissolutionData(2,8)  = 0.
		DissolutionData(2,9)  = 0.
		DissolutionData(2,10) = 0.65
		DissolutionData(2,11) = 7.
		DissolutionData(2,12) = 0.
		DissolutionData(2,13) = 2.
		DissolutionData(2,14) = 1.
		DissolutionData(2,15) = 0.

		!! Rat out the hard coded sulfate transfer
		IF (TRANSCRIBE) CALL TRANSCRIPT("H2SO4 (v) <=> 2 H+ & SO4--")

		DO I = 3, HowManyDissolutionReactions-1
			DO J = 1, 15
				DissolutionData(i,j) = 0.
			END DO
		END DO

		DO I = 3, HowManyDissolutionReactions
			CurrLine = GetLine(FH)

			!! Tokenize the input line !!
			!! The token order is:
			!! ---	1. The reaction (which we must sub-divide)
			CALL GetToken(CurrLine, ";", Token)	

			CALL GetToken(Token, "<=>", GasChem2)	
			CALL GetToken(GasChem2, "&", GasChem1)

			AqIon2 = TRIM(Token)
			CALL GetToken(AqIon2, "&", AqIon1)

			!! Find the index for the first 
			IF (LEN_TRIM(GasChem1) .EQ. 0) &
				CALL ERROR("In Dissolution.in, either one of the reactions has no chemical specified or it is specificied "// &
						   "after an ampersand, leading me to think that the first of two chemicals is unspecified.  Both "// &
						   "of these are bad enough to ask you to fix it.")
			DissolutionData(i,1) = FindChem(TRIM(GasChem1),GasPhase,SuppressError=.TRUE.)	
			IF (DissolutionData(i,1) .EQ. 0) &
				CALL ERROR("I don't think you defined "//TRIM(GasChem1)// " in GasPhaseChems.in.  It must be there to "// &
						   "define a reaction using it in Dissolution.in, as you have done.")

			!! Presume reaction type 1 and change later if need be.
			DissolutionData(I,11) = 1. 

			!! Ask if there is a second species specified
			!! There are two possibilities for this chemical, if it exists:
			!!  a. It is a gas phase species
			!!  b. It is an aqueous phase ion 
			IF (LEN_TRIM(GasChem2) .GT. 0) THEN

				!! Use the FindChem and FindIon routines to identify it.
				DissolutionData(I,2)  = FindChem(TRIM(GasChem2),GasPhase,SuppressError=.TRUE.)
				DissolutionData(I,11) = 6.

				!! If it's an ion
				IF (DissolutionData(I,2) .EQ. 0) THEN

					IonIndex = FindIon(TRIM(GasChem2),SuppressError=.TRUE.)

					SELECT CASE(IonIndex(1))

					CASE(-1)
					DissolutionData(i,2) = IonIndex(2) + HowManyAqChems + HowManyAqCations

					CASE(1)
					DissolutionData(i,2) = IonIndex(2) + HowManyAqChems

					CASE(0)
					CALL ERROR("You asked for a second gas or ionic species ("//TRIM(GasChem2)//") on the left hand side of a "// &
							   "dissolution reaction that is not an aqueous phase ion or a gas phase species.  It must be one "// &
							   "of those to have a second species on the left hand side of the equation in Dissolution.in in "// &
							   "reaction ("//TRIM(GasChem1)//" & "//TRIM(GasChem2)//" <=> "//TRIM(Token)//")")

					END	SELECT

					!! An ion in this slot means it is reaction type 4.
					DissolutionData(i,11) = 4.
				END IF

			ELSE
				DissolutionData(i,2) = 0.
			END IF

			!! If there is a third species specified
			IF (LEN_TRIM(AqIon1) .EQ. 0) THEN
				IF (DissolutionData(i,2) .EQ. 0) THEN
					DissolutionData(i,3) = FindChem(TRIM(GasChem1),AqPhase,SuppressError=.TRUE.)
					IF (DissolutionData(i,3) .EQ. 0.) CALL ERROR("I think you defined "//TRIM(GasChem1)// &
						" in GasPhaseChems.in but not in AqPhaseChems.in and then tried to let it transfer between "// &
						"phases by writing a reaction for it in Dissolution.in, as you have done.  No dice.")
					DissolutionData(i,13) = 1.
				ELSE
					CALL ERROR("You can't define a two-body dissolution reaction without defining what it forms in the "// &
							   "aqueous phase, as you've done with "//trim(GasChem1)//" & "//trim(GasChem2)//" <=> "//		&
							   trim(AqIon1)//" & "//trim(AqIon2))
				END IF
			ELSE

				CALL GetToken(AqIon1, " ", Token)

				IF (LEN_TRIM(AqIon1) .EQ. 0) THEN
					AqIon1 = TRIM(Token)
					Token  = ""
				END IF

				IF (IsReal(Token)) THEN
					DissolutionData(I,13) = STR2REAL(Token)
				ELSE
					DissolutionData(I,13) = 1
				END IF

				DissolutionData(i,3) = FindChem(TRIM(StripToken(AqIon1)),AqPhase,SuppressError=.TRUE.)	
				IF (DissolutionData(i,3) .EQ. 0.) &
					CALL ERROR("I don't think you defined "//TRIM(Token)// " in AqPhaseChems.in.  It must be there to define "// &
							   "a reaction using it in Dissolution.in, as you have done.")
			END IF


			!! And if there is a fourth
			IF (LEN_TRIM(AqIon2) .NE. 0) THEN


				CALL GetToken(AqIon2, " ", Token)

				IF (LEN_TRIM(AqIon2) .EQ. 0) THEN
					AqIon2 = TRIM(Token)
					Token  = ""
				END IF

				IF (IsReal(Token)) THEN
					DissolutionData(I,14) = STR2REAL(Token)
				ELSE
					DissolutionData(I,14) = 1
				END IF

				DissolutionData(I,4) = FindChem(TRIM(StripToken(AqIon2)),AqPhase,SuppressError=.TRUE.)	
				IF (DissolutionData(I,4) .EQ. 0.) &
					CALL ERROR("I don't think you defined "//TRIM(Token)// " in AqPhaseChems.in.  It must be there to "// &
							   "define a reaction using it in Dissolution.in, as you have done.")

				!! i.d. the reaction type
				IF (DissolutionData(I,2) .EQ. 0.) THEN
					DissolutionData(I,11) = 2.
				ELSE
					DissolutionData(I,11) = 6.
				END IF
			END IF
			
			!! It is a problem if 
			IF (DissolutionData(I,3) .NE. FindChem(TRIM(GasChem1),AqPhase,SuppressError=.TRUE.) .AND. &
				DissolutionData(I,2) .EQ. 0 .AND. DissolutionData(I,4) .EQ. 0)						  &
				CALL WARN("In Dissolution.in, I am suspicious of the reaction"//trim(GasChem1)//" <=> "//trim(AqIon1)//	&
						  ", in which the chemical appears to change during the phase change.")

			!! Don't want the user to specify water
			IF ((TRIM(GasChem1) .EQ. "H2O"  .OR.	&
				TRIM(GasChem1) .EQ. "h2o"   .OR.	&
				TRIM(GasChem1) .EQ. "water" .OR.	&
				TRIM(GasChem1) .EQ. "Water" .OR.	&
				TRIM(GasChem1) .EQ. "WATER") .AND.  &
				DissolutionData(i,2) .EQ. 0.) &
					CALL ERROR ("You specified "//TRIM(GasChem1)//" in Dissolution.in.  This model hardcodes in the values for"// &
								" water.  Please don't specify anything about it and use 'H2O' when referring to it in "// &
								"chemical mechanisms.  If you want to change what is used for water, you'll have to go to "// &
								"SUBROUTINE SetDissolution() in CondensationInitialization.h.")

			!! ---	2. Mass Accommodation Coefficient
			CALL GetToken(CurrLine,";",Token) 
			DissolutionData(i,10) = STR2REAL (Token)
			IF (DissolutionData(i,10) .LT. 0. .OR. DissolutionData(i,10) .GT. 1.) THEN
				Token = TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) // " (g)"
				IF (DissolutionData(i,2) .GT. 0) &
					Token = TRIM(Token) // " + " // TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,2)))) // " (g)"
				Token = TRIM(Token) // " <=> " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) // " (aq)"
				CALL ERROR ("Specified an out of range mass accommodation coefficient in Dissolution.in for the reaction "// &
							TRIM(Token)//".  It should be between zero and unity.")
			END IF


			!! ---	3 & 4. HENRY'S LAW COEFFICIENT (three elements, H298 then Delta H/RT then Cp/R)
			!! ---  these are in appropriate units to negate molality and atmospheres
			CALL GetToken(CurrLine,";",Token) 
			DissolutionData(i,5) = STR2REAL (Token) / Atmospheres
			
			!For 2 Gas phase species, the units are in 1/atm2
			IF (DissolutionData(i,11) .EQ. 6.) DissolutionData(i,5) = DissolutionData(i,5) / Atmospheres

			CALL GetToken(CurrLine,";",Token) 
			DissolutionData(i,6) = STR2REAL (Token)

			CALL GetToken(CurrLine,";",Token) 
			DissolutionData(i,7) = STR2REAL (Token)


			!! Report the particular reaction
			Token = TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) // " (g)"

			IF (DissolutionData(i,2) .GT. 0 .AND. DissolutionData(i,11) .EQ. 4.) &
				Token = TRIM(Token) // " + " //TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,2)))) // " (aq)"
			IF (DissolutionData(i,2) .GT. 0 .AND. DissolutionData(i,11) .EQ. 6.) &
				Token = TRIM(Token) // " + " // TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,2)))) // " (g)"
			Token = TRIM(Token) // " <=> " 

			IF (DissolutionData(i,13) .GT. 1) Token = TRIM(Token) // " " //TRIM(REAL2STR(DissolutionData(i,13)))
			Token = TRIM(Token) // " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) // " (aq)"
			
			IF (INT(DissolutionData(i,4)) .GT. 0) THEN
				!WRITE(*,*) DissolutionData(i,4)
				Token = TRIM(Token) // " + " 
				IF (DissolutionData(i,14) .GT. 1) Token = TRIM(Token) // " " //TRIM(REAL2STR(DissolutionData(i,13)))
				Token = TRIM(Token) // " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,4)))) // " (aq)"
			END IF

			CALL TRANSCRIPT(TRIM(Token))			

		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check to make sure the same pair wasn't double specified !!
		IF (HowManyDissolutionReactions .GT. 2) THEN !!!!!!!!!!!!!!!!!

		DO I = 3, HowManyDissolutionReactions

		!! Check the reactions for errors of various types
		SELECT CASE (INT(DissolutionData(I,11)))
		CASE (1:4)


			DO J = 2,HowManyDissolutionReactions
				Match = .FALSE.
				IF (I .EQ. J) CYCLE
				IF (DissolutionData(I,1) .EQ. DissolutionData(J,1) .OR. &
					(DissolutionData(J,11) .EQ. 5. .AND.			    &
					 DissolutionData(I,1) .EQ. DissolutionData(J,2))) THEN

						IF (DissolutionData(J,11) .EQ. 5.) THEN
						   CALL WARN("It seems you specified one dissolution reaction in Dissolution.in that dissolves " &
						   // TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) //										 &
						   " (g) directly and one  that dissolves it in combination with another gas phase species.  "// &
						   "This may lead to an unequilatable system depending on contents and the relative equilibria " &
						   //"of these reactions.  You should revisit this if problems occur, and possibly otherwise as well.")

						ELSE IF (DissolutionData(I,11) .EQ. 4. .AND. DissolutionData(J,11) .EQ. 1.) THEN
							
						   Token = TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,2))))
						   CALL WARN("It seems you specified a dissolution reaction in Dissolution.in that dissolves "			&
						   //TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) // "(g)  in combination with " // TRIM(Token) // &
						   " (aq)" //	&
						   " and one  that dissolves  "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) // " directly. " //	&
						   "This is necessary (good job), but make sure the equilibria are consistent with eachother and an " // &
						   "aqueous equilibrium reaction that exchanges the forms.")

						   !! Identify that these reactions are coupled
						   !! make one negative (so we skip it in the cycling)
						   !! and one positive (which we key the integration off of)
						   DissolutionData(J,15) = -1 * I
						   DissolutionData(I,15) = J
						   MATCH = .TRUE.

						ELSE IF (DissolutionData(I,11) .EQ. 2. .AND. DissolutionData(J,11) .EQ. 1.) THEN

						   CALL WARN("It seems you specified a dissolution reaction in Dissolution.in that dissolves "			&
						   //TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) // "(g)  directly and one that dissociates it. "// &
						   "This is necessary (good job), but make sure the equilibria are consistent with eachother and an " // &
						   "aqueous equilibrium reaction that exchanges the forms.")

						   DissolutionData(J,15) = -1 * I
						   DissolutionData(I,15) = J
						   MATCH = .TRUE.
						END IF
						IF ((DissolutionData(I,11) .EQ. 1 .OR. DissolutionData(I,11) .EQ. 1) .AND.	&
							(DissolutionData(I,11) .NE. DissolutionData(J,11))) MATCH = .TRUE.
				END IF
			END DO
			IF (.NOT.MATCH) CALL WARN("In dissolution.in, you only specified one direct dissolution reaction for "//	&
									  TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1))))//".  You should specify reactions "// &
									  "that lead to both the associated and dissociated forms of the electrolyte unless: a). "// &
									  "it does not dissociate or b). only its dissociated form is allowed.")

		CASE (6)

			IF (DissolutionData(I,11) .EQ. 6) &
			CALL WARN ("The reaction "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1))))//" (g) + "//					&
									    TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,2))))//" (g) <=> "//				&
			TRIM(REAL2STR(DissolutionData(i,13)))//" "//TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3))))//" (s) "//	&
			"is a gas to solid reaction. No activity coefficients are necessary.")
		END SELECT

		END DO
		END IF


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check Charge and Mass Balance !!
		!! along with other matching and !!
		!! last pass type of things.     !!
		IF (HowManyDissolutionReactions .GT. 1) THEN
		DO I = 3, HowManyDissolutionReactions

			!! Figure Mass Balance
			MassTally = GasMolecularMass(INT(DissolutionData(I,1)))

			!! Add mass of 2nd l.h.s. chemical
			IF (DissolutionData(I,2) .GT. 0.) THEN
				IF (DissolutionData(i,11) .EQ. 4.) MassTally = MassTally + AqMolecularMass(INT(DissolutionData(I,2)))
				IF (DissolutionData(i,11) .EQ. 6.) MassTally = MassTally + GasMolecularMass(INT(DissolutionData(I,2)))
			END IF

			!! divide by the r.h.s. masses
			IF (INT(DissolutionData(I,4)) .EQ. 0.) THEN
				MassTally = MassTally / AqMolecularMass(INT(DissolutionData(I,3))) / INT(DissolutionData(I,13))
			ELSE 
				MassTally = MassTally / (AqMolecularMass(INT(DissolutionData(I,3))) * INT(DissolutionData(I,13)) + &
										 AqMolecularMass(INT(DissolutionData(I,4))) * INT(DissolutionData(I,14)))
			END IF


			!! err if it is not within 0.1% (which I arbitrarily set as a rounding error
			IF (MassTally .GT. 1.001 .OR. MassTally .LT. 0.999) THEN
				Token = TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) // " (g)"
				IF (DissolutionData(i,2) .GT. 0) THEN
					IF (DissolutionData(i,11) .EQ. 4.) &
					    Token = TRIM(Token) // " + " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,2)))) // " (aq)"
					IF (DissolutionData(i,11) .EQ. 6.) &
					    Token = TRIM(Token) // " + " // TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,2)))) // " (g)"
				END IF
				Token = TRIM(Token) // " <=> " // TRIM(REAL2STR(DissolutionData(i,13))) // " " // &
						TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) // " (aq)"
				IF (DissolutionData(i,4) .GT. 0) &
					Token = TRIM(Token) // " " // TRIM(REAL2STR(DissolutionData(i,14))) // " " // &
						    TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,4)))) // " (aq)"

				CALL ERROR("It seems that the mass balance is out of whack in dissolution reaction "//TRIM(Token)// &
				           " in Dissolution.in.  This is a problem.")
			END IF

			!! Figure Charge Balance and 
			!! Match Reactions
			SELECT CASE (INT(DissolutionData(i,11)))

			CASE (1)
				J = GetAqCharge(INT(DissolutionData(i,3)))
				IF (J .NE. 0) &
				CALL ERROR("While checking dissolution reactions, "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1))))// &
						   " (g) <=> "//TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) // " (aq) is not in charge balance.")

				!! Match Gas phase reaction to the appropriate aqueous phase reaction,
				!! if there is one (and the dissociation constant is high enough)
				DissolutionData(I,12) = 0.
				DO J = 1, HowManyAqEqReactions
					IF (DissolutionData(I,3) .EQ. AqEquilibriaList(J,1)) THEN
						DissolutionData(I,12) = J

						!! Error check this matchup.  Can't infinitely dissociate...
						IF (DissolutionData(I,11) .EQ. 1 .AND.	&
						    AqEquilibriaList(J,9) .EQ. -1)		&
							CALL ERROR ("Dissolution.in specifies that "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))// &
							            " dissolve directly and AqEquilibriumReactions.in states that it dissociate "//       &
										"immediately and infinitely once in solution.  This system is not equilabratable.  "// &
										"Please pick another setup.")
						EXIT
					END IF
				END DO

			CASE (2)
				J = DissolutionData(i,13)*GetAqCharge(INT(DissolutionData(i,3))) + &
				    DissolutionData(i,14)*GetAqCharge(INT(DissolutionData(i,4)))
				IF (J .NE. 0) &
				CALL ERROR("While checking dissolution reactions, "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1))))//   &
				           " (g) <=> "//TRIM(REAL2STR(DissolutionData(I,13))) // " " // &
						   TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) // " + "//TRIM(REAL2STR(DissolutionData(I,14))) &
						   // " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,4))))//" (aq) is not in charge balance.")

				!! Match Gas phase reaction to the appropriate aqueous phase reaction,
				!! if there is one (and the dissociation constant is high enough)
				DissolutionData(I,12) = 0.

				DO J = 1, HowManyAqEqReactions
					IF ((DissolutionData(I,3)-HowManyAqChems					.EQ. AqEquilibriaList(J,2) .AND. &
					     DissolutionData(I,4)-HowManyAqChems-HowManyAqCations   .EQ. AqEquilibriaList(J,3) .AND. &
						 DissolutionData(I,13)									.EQ. AqEquilibriaList(J,4) .AND. &
						 DissolutionData(I,14)									.EQ. AqEquilibriaList(J,5)).OR.  &
						(DissolutionData(I,3)-HowManyAqChems-HowManyAqCations   .EQ. AqEquilibriaList(J,3) .AND. &
					     DissolutionData(I,4)-HowManyAqChems					.EQ. AqEquilibriaList(J,2) .AND. &
						 DissolutionData(I,13)									.EQ. AqEquilibriaList(J,5) .AND. &
						 DissolutionData(I,14)									.EQ. AqEquilibriaList(J,4))) THEN
						DissolutionData(I,12) = J
						EXIT
					END IF
				END DO
				
				IF (DissolutionData(I,12) .EQ. 0.) &
				CALL ERROR("For a dissolution reaction (specified in Dissolution.in) that dissociates directly into ions, " // &
				           "there must also be an aqueous dissociation reaction with the same constituent ions forming an " // &
						   "aqueous electrolyte for which thermodynamic quantities are supplied.  You may make this an "    // &
						   "infinitely dissociating reaction, in which case there will be none of it present in aqueous "   // &
						   "associated form.  The problematic reaction is: "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1)))) &
						   // " (g) <=> " // TRIM(REAL2STR(DissolutionData(i,13))) // " " // &
						   TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) // " (aq)" // TRIM(REAL2STR(DissolutionData(i,14))) &
						   // " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,4)))) // " (aq)")

			CASE (4)
				J = DissolutionData(i,13)*GetAqCharge(INT(DissolutionData(i,3))) + DissolutionData(i,14)*  &
				    GetAqCharge(INT(DissolutionData(i,4))) - GetAqCharge(INT(DissolutionData(i,2)))

				IF (J .NE. 0) THEN
					Token = "While checking dissolution reactions, "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1))))// &
					        " (g) "//TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,2)))) // "(aq) <=> "//  &
							TRIM(REAL2STR(DissolutionData(I,13))) // " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) 
					IF(DissolutionData(i,4) .GT. 0) Token = TRIM(Token) // " + "//TRIM(REAL2STR(DissolutionData(I,14))) // &
					                                        " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,4))))
					Token = TRIM(Token) //" (aq) is not in charge balance."
					CALL ERROR(Token)
				END IF

			CASE (6)
				J = DissolutionData(i,13)*GetAqCharge(INT(DissolutionData(i,3))) + &
				    DissolutionData(i,14)*GetAqCharge(INT(DissolutionData(i,4)))

				IF (J .NE. 0) THEN
					Token = "While checking dissolution reactions, "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,1))))// &
					        " (g) "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(i,2)))) // "(g) <=> "//                    &
							TRIM(REAL2STR(DissolutionData(I,13))) // " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,3)))) 
					IF(DissolutionData(i,4) .GT. 0) Token = TRIM(Token) // " + "//TRIM(REAL2STR(DissolutionData(I,14))) // &
					                                        " " // TRIM(AqPhaseChemicalNames(INT(DissolutionData(i,4))))
					Token = TRIM(Token) //" (aq) is not in charge balance."
					CALL ERROR(Token)
				END IF
			END SELECT

			!! Work through the thermodynamics if there is a loose ion attaching in the 
			!! aqueous phase.  Match a ratio to the problem at hand.
			!!
			!! e.g.,  for NH3 (g) + H+ (aq) <=> NH4+ (aq)
			!! 
			!! you'll need Gamma^2(NH4/OH) / Gamma^2(H+/NO3-) or the like...
			!! Does case 4 only.
			IF (DissolutionData(I,11) .EQ. 4.) THEN
				
				!! The ions are both cations...
				IF (DissolutionData(I,2) .LE. HowManyAqChems+HowManyAqCations) THEN

					DO J = 1, HowManyAqEqReactions

						!! Match the first cation (the denominator)
						IF (DissolutionData(I,2) .EQ. AqEquilibriaList(J,2)+HowManyAqChems .AND. &
							AqEquilibriaList(J,4) .EQ. 1) THEN

							DO K = 1, HowManyAqEqReactions

								!! Then the second cation (the denominator)
								IF (DissolutionData(I,3) .EQ. AqEquilibriaList(K,2)+HowManyAqChems .AND. &
								    AqEquilibriaList(J,3) .EQ. AqEquilibriaList(K,3) .AND. &
									AqEquilibriaList(J,5) .EQ. 1					 .AND. &
									AqEquilibriaList(K,5) .EQ. 1 ) THEN
							
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("Found a thermodynamic treatment for the dissolution reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + " &
									                 //TRIM(AqCationNames(INT(DissolutionData(I,2))-HowManyAqChems))//" <=> " &
													 //TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems)))
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("  g("//TRIM(AqCationNames(INT(DissolutionData(I,2))-HowManyAqChems))//") / g("// &
									                 TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems))//") =")
									CALL TRANSCRIPT ("  g("//TRIM(AqCationNames(INT(AqEquilibriaList(J,2))))//"/"//  &
									                 TRIM(AqAnionNames(INT(AqEquilibriaList(J,3))))//")^2 / g("//    &
													 TRIM(AqCationNames(INT(AqEquilibriaList(K,2))))//"/"//          &
													 TRIM(AqAnionNames(INT(AqEquilibriaList(K,3))))//")^2")
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("")

									DissolutionData(I,8) = K + 2. / 10000.
									DissolutionData(I,9) = J + 2. / 10000.

									GOTO 20
								END IF
							END DO
						END IF
					END DO

20					IF (DissolutionData(I,8) .EQ. 0) &
							CALL ERROR("Couldn't find a thermodynamic treatment for the dissolution reaction "//	&
							TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + "//TRIM(AqCationNames			&
							(INT(DissolutionData(I,2))-HowManyAqChems))//" <=> "//TRIM(AqCationNames(INT(DissolutionData(I,3))- &
							HowManyAqChems))//" for which we wanted to match some mean activities to: g("//			&
							TRIM(AqCationNames(INT(DissolutionData(I,2))-HowManyAqChems))//") / g("//					&
							TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems))//") = ???.  You should have reactions "// &
							"in AqEquilibriaReactions.in that will provide correct fodder.")
			
				!! If the Ions are both Anions
				ELSE

					DO J = 1, HowManyAqEqReactions

						!! Match the first cation (the denominator)
						IF (DissolutionData(I,2) .EQ. AqEquilibriaList(J,3)+HowManyAqChems+HowManyAqCations .AND. &
							AqEquilibriaList(J,5) .EQ. 1) THEN

							DO K = 1, HowManyAqEqReactions

								!! Then the second cation (the denominator)
								IF (DissolutionData(I,3) .EQ. AqEquilibriaList(K,3)+HowManyAqChems+HowManyAqCations .AND. &
								    AqEquilibriaList(J,2) .EQ. AqEquilibriaList(K,2) .AND. &
									AqEquilibriaList(J,4) .EQ. 1					 .AND. &
									AqEquilibriaList(K,4) .EQ. 1 ) THEN

									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("Found a thermodynamic treatment for the dissolution reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + "// &
													 TRIM(AqCationNames(INT(DissolutionData(I,2))-HowManyAqChems))//" <=> "//  &
													 TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems)))
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("  g("//TRIM(AqAnionNames(INT(DissolutionData(I,2))-HowManyAqChems-	&
													 HowManyAqCations))//") / g("//TRIM(AqAnionNames(INT(DissolutionData(I,3))- &
													 HowManyAqChems-HowManyAqCations))//") =")
									CALL TRANSCRIPT ("  g("//TRIM(AqCationNames(INT(AqEquilibriaList(K,2))))//"/"//TRIM(  &
													 AqAnionNames(INT(AqEquilibriaList(K,3))))//")^2 / g("//TRIM(         &
													 AqCationNames(INT(AqEquilibriaList(J,2))))//"/"//TRIM(               &
													 AqAnionNames(INT(AqEquilibriaList(J,3))))//")^2")
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("This may be buggy -- the anion matching part of this.")
									CALL TRANSCRIPT ("If it doesn't work well check CondensationInitialization.h")
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("")

									DissolutionData(I,8) = K + 2. / 10000.
									DissolutionData(I,9) = J + 2. / 10000.

									GOTO 25
								END IF
							END DO
						END IF
					END DO
				END IF

25				IF (DissolutionData(I,8) .EQ. 0) &
					CALL ERROR("Couldn't find a thermodynamic treatment for the dissolution reaction "//	&
							 TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + "//                   &
							 TRIM(AqAnionNames(INT(DissolutionData(I,2))-HowManyAqChems-HowManyAqCations))//   &
							 " <=> "//TRIM(AqAnionNames(INT(DissolutionData(I,3))-HowManyAqChems-HowManyAqCations))// &
							 " for which we wanted to match some mean activities to: g("//                &
							 TRIM(AqAnionNames(INT(DissolutionData(I,2))-HowManyAqChems-HowManyAqCations))// &
							 ") / g("//TRIM(AqAnionNames(INT(DissolutionData(I,3))-HowManyAqChems-HowManyAqCations))// &
							 ") = ???.  You should have reactions in AqEquilibriaReactions.in that will provide correct fodder.")

			END IF

			!!Try to find matching dissociation reaction for a reaction
			!!where 2 gas-phase species become 2 aq phase ions, so that we can 
			!!calculate the mixed mean activity coefficient.
			
			IF (1 .EQ. 0) THEN
			!IF (DissolutionData(I,11) .EQ. 6.) THEN
				
				!! The first ion is a cation
				
				!WRITE(*,*) AqPhaseChemicalNames(INT(DissolutionData(I,3)))
				!WRITE(*,*) AqPhaseChemicalNames(INT(DissolutionData(I,4)))
				IF (DissolutionData(I,3) .LE. HowManyAqChems+HowManyAqCations) THEN
	
					DO J = 1, HowManyAqEqReactions

						!! Match the first cation (the denominator)
						IF (DissolutionData(I,3) .EQ. AqEquilibriaList(J,2)+HowManyAqChems .AND. & !Cation matches
							AqEquilibriaList(J,4) .EQ. 1) THEN !Cation stoichiometry = 1

							!WRITE(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(J,2)+HowManyAqChems))

							DO K = 1, HowManyAqEqReactions

								!! Then the anion
								IF (DissolutionData(I,4) .EQ. AqEquilibriaList(K,3)+HowManyAqChems+HowManyAqCations .AND. & !Anion matches
								    AqEquilibriaList(J,2) .EQ. AqEquilibriaList(K,2) .AND. & !Cation Matches
									AqEquilibriaList(J,5) .EQ. 1					 .AND. & !Anion stoichiometry = 1 in both reactions
									AqEquilibriaList(K,5) .EQ. 1 ) THEN
							
									
									!WRITE(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(K,3)+HowManyAqChems+HowManyAqCations))
									!WRITE(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(K,2)+HowManyAqChems))
									
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("Found a thermodynamic treatment for the dissolution reaction: ")
									CALL TRANSCRIPT ("      "//TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + " &
									                 //TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,2))))//" <=> " &
													 //TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems))//" + " &
													 //TRIM(AqAnionNames(INT(DissolutionData(I,4))-HowManyAqChems-HowManyAqCations)))
									CALL TRANSCRIPT ("")
									CALL TRANSCRIPT ("Will use Dissolution Reaction #"//TRIM(INT2STR(K))//" to calculate mean activity.")
									CALL TRANSCRIPT ("*********************************************************")
									CALL TRANSCRIPT ("")
									

									DissolutionData(I,12) = K !Index of associated dissolution reaction

									GOTO 30
								END IF
							END DO
						END IF
					END DO

30					IF (DissolutionData(I,12) .EQ. 0) &
							CALL ERROR("Couldn't find a thermodynamic treatment for the dissolution reaction "//	&
										TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + " &
									    //TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,2))))//" <=> " &
										//TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems))// " + " &
										//TRIM(AqAnionNames(INT(DissolutionData(I,4))-HowManyAqChems-HowManyAqCations)))
								
				ELSE	!Anion Listed First
					CALL ERROR("Couldn't find a thermodynamic treatment for the dissolution reaction "//	&
										TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,1))))//" + " &
									    //TRIM(GasPhaseChemicalNames(INT(DissolutionData(I,2))))//" <=> " &
										//TRIM(AqCationNames(INT(DissolutionData(I,3))-HowManyAqChems))// " + " &
										//TRIM(AqAnionNames(INT(DissolutionData(I,4))-HowManyAqChems-HowManyAqCations))//". "// &
										"Please list the cation first on the right hand side.")
							
				
				END IF
			END IF

		END DO
		END IF

		!! SCSCSCSC -- Announce Departure
		IF (Scaffolding) CALL Transcript(">>Exiting SetDissolution()<<")
		
		CLOSE (FH)
		CALL ReturnFileHandle(FH)

		RETURN
	END SUBROUTINE SetDissolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Hydrophobic Phase Chemical Information from			    !!
!! input decks, defines the appropriate chemical arrays in this module,	    !!
!! and calls the appropriate subroutines to establish the	            !!
!! chemical fields							    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetOrganicDissolution ()

		USE ModelParameters, ONLY : AqChemScale, Atmospheres,	&
					SetDissolutionEquilibriumOrFlux, &
					DissolutionEquilibriumOrFlux,	 &
					ThermoBulkMode, &
					TranscriptFH, InputDeckSubDir, &
					DoHydrophobicOrgDissolutionFlag

		USE InfrastructuralCode

		USE Chemistry,       ONLY : GasPhase, AqPhase,		&
					FindChem, FindIon,		&
					GasPhaseChemicalNames,		&
					OrgPhaseChemicalNames,		&
					OrgMolecularMass, &
					GasMolecularMass 


		IMPLICIT NONE

		!! Local Infrastructural Variables
		INTEGER	:: I, J, K, allocation_error,IonIndex(2)
		REAL*8	:: MassTally, H2SO4GasIndex, HSO4AqIndex, HAqIndex
		CHARACTER(len=512)	:: CurrLine, Token, GasChem1, &
                                           GasChem2, AqIon1, AqIon2
		CHARACTER (len=12)	:: Numbers = "0123456789+-"

		!! Transcribe?  Use Scaffolding Code?
		LOGICAL :: Transcribe, Match
		LOGICAL, PARAMETER :: Scaffolding = .TRUE.

		!! And there is a file handle to be used
		INTEGER :: FH
		FH  = GetFileHandle()

		!! SCSCSCSC -- Announce Arrival
		IF (Scaffolding) CALL Transcript(">>Entering SetOrganicDissolution()<<")


		HowManyOrganicDissolutionReactions = 0

		IF (DoHydrophobicOrgDissolutionFlag) THEN
		!! In this subroutine, we will create listings of all of the chemicals 
			OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'OrganicDissolution.in', STATUS='OLD')

			!! Make a first pass through the file to determine the number of 
			!! reactions specified therein.
			!! Need this number early to allocate all of the arrays.
			CurrLine = GetLine(FH)
			DO WHILE (CurrLine .NE. EOF)
				HowManyOrganicDissolutionReactions = HowManyOrganicDissolutionReactions+1
                                CurrLine = GetLine(FH)
			END DO

			!! Back up for a second pass to parse the reactions
			REWIND(FH)

		END IF
			

		!! Transcribe
		IF (TranscriptFH > 0) THEN
			Transcribe = .TRUE.
		ELSE
			Transcribe = .FALSE.
		END IF
                Transcribe = .false.    ! cmb
		IF (Transcribe) CALL TRANSCRIPT("")
		IF (Transcribe) CALL TRANSCRIPT("_The_Following_Organic_Dissolution_Reactions_Are_Defined_")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Allocation the key ORGANIC DISSOLUTION Arrays			    !!
		!!									    !!
		!! The primary one is OrganicDissolutionData(i,j), in which i is the	    !!
		!! reaction and j is some data:						    !!
		!!									    !!
		!! OrganicDissolutionData(i,1)  : Index in Gas Phase of Primary Chemical    !!
		!! OrganicDissolutionData(i,2)  : Flag1: 0 for Myrdal and Yalkowski (1997)  !!
                !!                                method; 1 for Psat(300K) and Hvap method  !!
		!! OrganicDissolutionData(i,3)  : Index in Organic Phase of Primary Chemical!!
		!! OrganicDissolutionData(i,4)  : Blank                                     !!
		!! OrganicDissolutionData(i,5)  : Boiling Point (K) (for flag1 = 0)	    !!
		!! OrganicDissolutionData(i,6)  : # of Effective Torsional Bonds            !!
                !!                                (for flag1 = 0)	                    !!
		!! OrganicDissolutionData(i,7)  : Hydrogen Bond Number (for flag1 = 0)      !!
		!! OrganicDissolutionData(i,8)  : Vapor Pressure Correction Factor	    !!
                !!                                (for flag1 = 0)	                    !!
		!! OrganicDissolutionData(i,9)  : Blank                 		    !!
		!! OrganicDissolutionData(i,10) : Mass Accommodation Coefficient	    !!
		!! OrganicDissolutionData(i,11) : Reaction Type	(For StepCondensation)	    !!
		!! OrganicDissolutionData(i,12) : Blank                                     !!
		!! OrganicDissolutionData(i,13) : Stoicheometry of r.h.s. slot 1	    !!
		!! OrganicDissolutionData(i,14) : Psat at 300 K (mbar) (for flag1 = 1)      !!
		!! OrganicDissolutionData(i,15) : DeltaHvap (kJ/mol) (for flag1 = 1)        !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! There is only one reaction type:					    !!
		!!									    !!
		!! 8. Direct dissolution of a single species not incorporating dissociation !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		ALLOCATE (OrganicDissolutionData(HowManyOrganicDissolutionReactions,15), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of OrganicDissolutionData could not proceed in SetOrganicDissolution()")

		DO I = 1, HowManyOrganicDissolutionReactions-1
			DO J = 1, 15
				OrganicDissolutionData(i,j) = 0.
			END DO
		END DO

		DO I = 1, HowManyOrganicDissolutionReactions
			CurrLine = GetLine(FH)

			!! Tokenize the input line !!
			!! The token order is:
			!! ---	1. The reaction (which we must sub-divide)
			CALL GetToken(CurrLine, ";", Token)	

			CALL GetToken(Token, "<=>", GasChem2)	
			CALL GetToken(GasChem2, "&", GasChem1)

			AqIon2 = TRIM(Token)
			CALL GetToken(AqIon2, "&", AqIon1)

			!! Find the index for the first 
			IF (LEN_TRIM(GasChem1) .EQ. 0) &
				CALL ERROR("In OrganicDissolution.in, either one of the reactions has no chemical specified or it is specificied "// &
						   "after an ampersand, leading me to think that the first of two chemicals is unspecified.  Both "// &
						   "of these are bad enough to ask you to fix it.")
			
			OrganicDissolutionData(i,1) = FindChem(TRIM(GasChem1),GasPhase,SuppressError=.TRUE.)	
			
			IF (OrganicDissolutionData(i,1) .EQ. 0) &
				CALL ERROR("I don't think you defined "//TRIM(GasChem1)// " in GasPhaseChems.in.  It must be there to "// &
						   "define a reaction using it in OrganicDissolution.in, as you have done.")

			!! Presume reaction type 8
			OrganicDissolutionData(I,11) = 8.
			OrganicDissolutionData(i,3) = FindChem(TRIM(StripToken(AqIon1)),2,SuppressError=.TRUE.)	
			OrganicDissolutionData(i,13) = 1.0
				IF (OrganicDissolutionData(i,3) .EQ. 0.) &
					CALL ERROR("I don't think you defined "//TRIM(Token)// " in OrgPhaseChems.in.  It must be there to define "// &
							   "a reaction using it in OrganicDissolution.in, as you have done.")

	
						
			!! It is a problem if 
			IF (OrganicDissolutionData(I,3) .NE. FindChem(TRIM(GasChem1),2,SuppressError=.TRUE.) .AND. &
				OrganicDissolutionData(I,2) .EQ. 0 .AND. OrganicDissolutionData(I,4) .EQ. 0)						  &
				CALL WARN("In OrganicDissolution.in, I am suspicious of the reaction"//trim(GasChem1)//" <=> "//trim(AqIon1)//	&
						  ", in which the chemical appears to change during the phase change.")

			!! Don't want the user to specify water
			IF ((TRIM(GasChem1) .EQ. "H2O"  .OR.	&
				TRIM(GasChem1) .EQ. "h2o"   .OR.	&
				TRIM(GasChem1) .EQ. "water" .OR.	&
				TRIM(GasChem1) .EQ. "Water" .OR.	&
				TRIM(GasChem1) .EQ. "WATER") .AND.  &
				OrganicDissolutionData(i,2) .EQ. 0.) &
					CALL ERROR ("You specified "//TRIM(GasChem1)//" in OrganicDissolution.in.  This model hardcodes in the values for"// &
								" water.  Please don't specify anything about it and use 'H2O' when referring to it in "// &
								"chemical mechanisms.  If you want to change what is used for water, you'll have to go to "// &
								"SUBROUTINE SetOrganicDissolution() in CondensationInitialization.h.")

			!! ---	2. Mass Accommodation Coefficient
			CALL GetToken(CurrLine,";",Token) 
			OrganicDissolutionData(i,10) = STR2REAL (Token)
			IF (OrganicDissolutionData(i,10) .LT. 0. .OR. OrganicDissolutionData(i,10) .GT. 1.) THEN
				Token = TRIM(GasPhaseChemicalNames(INT(OrganicDissolutionData(i,1)))) // " (g)"
				IF (OrganicDissolutionData(i,2) .GT. 0) &
					Token = TRIM(Token) // " + " // TRIM(GasPhaseChemicalNames(INT(OrganicDissolutionData(i,2)))) // " (g)"
				Token = TRIM(Token) // " <=> " // TRIM(OrgPhaseChemicalNames(INT(OrganicDissolutionData(i,3)))) // " (aq)"
				CALL ERROR ("Specified an out of range mass accommodation coefficient in OrganicDissolution.in for the reaction "// &
							TRIM(Token)//".  It should be between zero and unity.")
			END IF


                        !! --- 3. Flag:  0 for Myrdal and Yalkowski (1997) method; 1 for Psat(300K) and Hvap method  !!
                        CALL GetToken(CurrLine,";",Token) 
			OrganicDissolutionData(i,2) = STR2REAL (Token)

                        IF (OrganicDissolutionData(i,2) .LT. 0.005) THEN

			    !! ---	4, 5, 6, 7,  VAPOR PRESSURE PARAMETERS, following Myrdal and Yalkowski
			
			    !4. BOILING POINT in K
			    CALL GetToken(CurrLine,";",Token) 
			    OrganicDissolutionData(i,5) = STR2REAL (Token)
			
			    !5. TORSIONAL BOND NUMBER
			    CALL GetToken(CurrLine,";",Token) 
			    OrganicDissolutionData(i,6) = STR2REAL (Token)

			    !6. HYDROGEN BOND NUMBER
			    CALL GetToken(CurrLine,";",Token) 
			    OrganicDissolutionData(i,7) = STR2REAL (Token)

			    !7. Vapor Pressure Correction Factor
			    CALL GetToken(CurrLine,";",Token) 
			    OrganicDissolutionData(i,8) = STR2REAL (Token)
                                 
                        ELSE

			    !! ---	4, 5   Psat(300K) and Hvap (KJ/mol)
			
			    !4. Psat at 300 K in mbar
			    CALL GetToken(CurrLine,";",Token) 
			    OrganicDissolutionData(i,14) = STR2REAL (Token)
			
			    !5. Delta Hvap in KJ/mol
			    CALL GetToken(CurrLine,";",Token) 
			    OrganicDissolutionData(i,15) = STR2REAL (Token)
                            
                        ENDIF

			!! Report the particular reaction
			Token = TRIM(GasPhaseChemicalNames(INT(OrganicDissolutionData(i,1)))) // " (g)"

			Token = TRIM(Token) // " <=> " 

			Token = TRIM(Token) // " " // TRIM(OrgPhaseChemicalNames(INT(OrganicDissolutionData(i,3)))) // " (org)"

			CALL TRANSCRIPT(TRIM(Token))			

		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check to make sure the same pair wasn't double specified !!
		IF (HowManyOrganicDissolutionReactions .GE. 1) THEN !!!!!!!!!!!!!!!!!

		DO I = 1, HowManyOrganicDissolutionReactions

			DO J = 1,HowManyOrganicDissolutionReactions
				Match = .FALSE.
				IF (I .EQ. J) CYCLE
				IF (OrganicDissolutionData(I,1) .EQ. OrganicDissolutionData(J,1) ) THEN

					CALL ERROR("You double-specified a OrganicDissolution reaction in OrganicDissolution.in. Please fix this.")

				END IF
			END DO
		END DO
		END IF


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check Mass Balance !!
		!! along with other matching and !!
		!! last pass type of things.     !!
		IF (HowManyOrganicDissolutionReactions .GE. 1) THEN
		DO I = 1, HowManyOrganicDissolutionReactions

			!! Figure Mass Balance
			MassTally = GasMolecularMass(INT(OrganicDissolutionData(I,1)))

			!! divide by the r.h.s. masses
			MassTally = MassTally / OrgMolecularMass(INT(OrganicDissolutionData(I,3))) / INT(OrganicDissolutionData(I,13))

			!! err if it is not within 0.1% (which I arbitrarily set as a rounding error
			IF (MassTally .GT. 1.001 .OR. MassTally .LT. 0.999) THEN
				Token = TRIM(GasPhaseChemicalNames(INT(OrganicDissolutionData(i,1)))) // " (g)"
				Token = TRIM(Token) // " <=> " // TRIM(REAL2STR(OrganicDissolutionData(i,13))) // " " // &
						TRIM(OrgPhaseChemicalNames(INT(OrganicDissolutionData(i,3)))) // " (org)"
				CALL ERROR("It seems that the mass balance is out of whack in OrganicDissolution reaction "//TRIM(Token)// &
				           " in OrganicDissolution.in.  This is a problem.")
			END IF

		END DO
		END IF

		!! SCSCSCSC -- Announce Departure
		IF (Scaffolding) CALL Transcript(">>Exiting SetOrganicDissolution()<<")
		
		CLOSE (FH)
		CALL ReturnFileHandle(FH)

		RETURN
	END SUBROUTINE SetOrganicDissolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read Hydrophilic Phase Chemical Information from			    !!
!! input decks, defines the appropriate chemical arrays in this module,	    !!
!! and calls the appropriate eulerian grid subroutines to establich the	    !!
!! chemical fields							    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetHydrophilicOrganicDissolution ()

		USE ModelParameters, ONLY : AqChemScale, Atmospheres,	&
					SetDissolutionEquilibriumOrFlux, &
					DissolutionEquilibriumOrFlux,	 &
					ThermoBulkMode, &
					TranscriptFH, InputDeckSubDir, &
					DoHydrophilicOrgDissolutionFlag

		USE InfrastructuralCode

		USE Chemistry,   ONLY : GasPhase, AqPhase,		&
				        FindChem, FindIon,		&
					GasPhaseChemicalNames,		&
					AqOrgPhaseChemicalNames,	&
					AqOrgMolecularMass, &
					GasMolecularMass 


		IMPLICIT NONE

		!! Local Infrastructural Variables
		INTEGER	:: I, J, K, allocation_error,IonIndex(2)
		REAL*8	:: MassTally, H2SO4GasIndex, HSO4AqIndex, HAqIndex
		CHARACTER(len=512)	:: CurrLine, Token, GasChem1, &
                                           GasChem2, AqIon1, AqIon2
		CHARACTER (len=12)	:: Numbers = "0123456789+-"

		!! Transcribe?  Use Scaffolding Code?
		LOGICAL :: Transcribe, Match
		LOGICAL, PARAMETER :: Scaffolding = .TRUE.

		!! And there is a file handle to be used
		INTEGER :: FH
		FH  = GetFileHandle()

		!! SCSCSCSC -- Announce Arrival
		IF (Scaffolding) CALL Transcript(">>Entering SetHydrophilicOrganicDissolution()<<")


		HowManyAqOrganicDissolutionReactions = 0

		IF (DoHydrophilicOrgDissolutionFlag) THEN
		!! In this subroutine, we will create listings of all of the chemicals 
			OPEN(UNIT=FH,  FILE=TRIM(InputDeckSubDir)//'HydrophilicOrganicDissolution.in', STATUS='OLD')


			!! Make a first pass through the file to determine the number of 
			!! reactions specified therein.
			!! Need this number early to allocate all of the arrays.
			CurrLine = GetLine(FH)
			DO WHILE (CurrLine .NE. EOF)
				HowManyAqOrganicDissolutionReactions = HowManyAqOrganicDissolutionReactions+1
				CurrLine = GetLine(FH)
			END DO

			!! Back up for a second pass to parse the reactions
			REWIND(FH)

		END IF
			

		!! Transcribe
		IF (TranscriptFH > 0) THEN
			Transcribe = .TRUE.
		ELSE
			Transcribe = .FALSE.
		END IF
                Transcribe = .false.    ! cmb
		IF (Transcribe) CALL TRANSCRIPT("")
		IF (Transcribe) CALL TRANSCRIPT("_The_Following_Hydrophilic_Organic_Dissolution_Reactions_Are_Defined_")

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Allocation the key HYDROPHILIC ORGANIC DISSOLUTION Arrays									!!
		!!																			!!
		!! The primary one is AqOrganicDissolutionData(i,j), in which i is the				!!
		!! reaction and j is some data:												!!
		!!																			!!
		!! AqOrganicDissolutionData(i,1)  : Index in Gas Phase of Primay Chemical			!!
		!! AqOrganicDissolutionData(i,2)  : Empty											!!
		!! AqOrganicDissolutionData(i,3)  : Index in Aqueous Phase of Primary Chemical		!!
		!! AqOrganicDissolutionData(i,4)  : Empty											!!
		!! AqOrganicDissolutionData(i,5)  : Henry's Law Constant								!!
		!! AqOrganicDissolutionData(i,6)  : DelH of Henry's Law Const					!!
		!! AqOrganicDissolutionData(i,7)  : Acid K1								!!
		!! AqOrganicDissolutionData(i,8)  : DelH of Acid K1											!!
		!! AqOrganicDissolutionData(i,9)  : Acid K2											!!
		!! AqOrganicDissolutionData(i,10) : Mass Accommodation Coefficient					!!
		!! AqOrganicDissolutionData(i,11) : Reaction Type									!!
		!! AqOrganicDissolutionData(i,12) : DelH of Acid K2										!!
		!! AqOrganicDissolutionData(i,13) : Stoicheometry of r.h.s. slot 1					!!
		!! AqOrganicDissolutionData(i,14) : Blank											!!!
		!! AqOrganicDissolutionData(i,15) : Blank											 !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! There is only one reaction type:											!!
		!!																				!!														!!
		!! 8. Direct dissolution of a single species not incorporating dissociation.	!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
		ALLOCATE (AqOrganicDissolutionData(HowManyAqOrganicDissolutionReactions,15), stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("Allocation of AqOrganicDissolutionData could not proceed in SetHydrophilicOrganicDissolution()")

		DO I = 1, HowManyAqOrganicDissolutionReactions-1
			DO J = 1, 15
				AqOrganicDissolutionData(i,j) = 0.
			END DO
		END DO

		DO I = 1, HowManyAqOrganicDissolutionReactions
			CurrLine = GetLine(FH)

			!! Tokenize the input line !!
			!! The token order is:
			!! ---	1. The reaction (which we must sub-divide)
			CALL GetToken(CurrLine, ";", Token)	

			CALL GetToken(Token, "<=>", GasChem2)	
			CALL GetToken(GasChem2, "&", GasChem1)

			AqIon2 = TRIM(Token)
			CALL GetToken(AqIon2, "&", AqIon1)

			!! Find the index for the first 
			IF (LEN_TRIM(GasChem1) .EQ. 0) &
				CALL ERROR("In HydrophilicOrganicDissolution.in, either one of the reactions has no chemical specified or it is specificied "// &
						   "after an ampersand, leading me to think that the first of two chemicals is unspecified.  Both "// &
						   "of these are bad enough to ask you to fix it.")
			
			AqOrganicDissolutionData(i,1) = FindChem(TRIM(GasChem1),GasPhase,SuppressError=.TRUE.)	
			
			IF (AqOrganicDissolutionData(i,1) .EQ. 0) &
				CALL ERROR("I don't think you defined "//TRIM(GasChem1)// " in GasPhaseChems.in.  It must be there to "// &
						   "define a reaction using it in HydrophilicOrganicDissolution.in, as you have done.")

			!! Presume reaction type 9 
			AqOrganicDissolutionData(I,11) = 9. 
			AqOrganicDissolutionData(I,15) = 0.

			AqOrganicDissolutionData(i,3) = FindChem(TRIM(StripToken(AqIon1)),3,SuppressError=.TRUE.)	
			AqOrganicDissolutionData(i,13) = 1.0
				IF (AqOrganicDissolutionData(i,3) .EQ. 0.) &
					CALL ERROR("I don't think you defined "//TRIM(Token)// " in HydrophilicOrgChems.in.  It must be there to define "// &
							   "a reaction using it in HydrophilicOrganicDissolution.in, as you have done.")

	
						
			!! It is a problem if 
			IF (AqOrganicDissolutionData(I,3) .NE. FindChem(TRIM(GasChem1),3,SuppressError=.TRUE.) .AND. &
				AqOrganicDissolutionData(I,2) .EQ. 0 .AND. AqOrganicDissolutionData(I,4) .EQ. 0)						  &
				CALL WARN("In HydrophilicOrganicDissolution.in, I am suspicious of the reaction"//trim(GasChem1)//" <=> "//trim(AqIon1)//	&
						  ", in which the chemical appears to change during the phase change.")

			!! Don't want the user to specify water
			IF ((TRIM(GasChem1) .EQ. "H2O"  .OR.	&
				TRIM(GasChem1) .EQ. "h2o"   .OR.	&
				TRIM(GasChem1) .EQ. "water" .OR.	&
				TRIM(GasChem1) .EQ. "Water" .OR.	&
				TRIM(GasChem1) .EQ. "WATER") .AND.  &
				AqOrganicDissolutionData(i,2) .EQ. 0.) &
					CALL ERROR ("You specified "//TRIM(GasChem1)//" in HydrophilicOrganicDissolution.in.  This model hardcodes in the values for"// &
								" water.  Please don't specify anything about it and use 'H2O' when referring to it in "// &
								"chemical mechanisms.  If you want to change what is used for water, you'll have to go to "// &
								"SUBROUTINE SetDissolution() in CondensationInitialization.h.")

			!! ---	2. Mass Accommodation Coefficient
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,10) = STR2REAL (Token)
			IF (AqOrganicDissolutionData(i,10) .LT. 0. .OR. AqOrganicDissolutionData(i,10) .GT. 1.) THEN
				Token = TRIM(GasPhaseChemicalNames(INT(AqOrganicDissolutionData(i,1)))) // " (g)"
				IF (AqOrganicDissolutionData(i,2) .GT. 0) &
					Token = TRIM(Token) // " + " // TRIM(GasPhaseChemicalNames(INT(AqOrganicDissolutionData(i,2)))) // " (g)"
				Token = TRIM(Token) // " <=> " // TRIM(AqOrgPhaseChemicalNames(INT(AqOrganicDissolutionData(i,3)))) // " (aq)"
				CALL ERROR ("Specified an out of range mass accommodation coefficient in HydrophilicOrganicDissolution.in for the reaction "// &
							TRIM(Token)//".  It should be between zero and unity.")
			END IF


			!! ---	Henry's Law Constants and temp dependence

			!Henry's Law Constant
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,5) = STR2REAL (Token)
			
			!Del H of Henry's Law (kcal/mol)
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,6) = STR2REAL (Token)

			!Acid K1 (mol/kg)
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,7) = STR2REAL (Token)

			!DelH of Acid K1 (kcal/mol)
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,8) = STR2REAL (Token)

			!Acid K2 (mol/kg)
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,9) = STR2REAL (Token)

			!DelH of Acid K2 (kcal/mol)
			CALL GetToken(CurrLine,";",Token) 
			AqOrganicDissolutionData(i,12) = STR2REAL (Token)
			
			!! Report the particular reaction
			Token = TRIM(GasPhaseChemicalNames(INT(AqOrganicDissolutionData(i,1)))) // " (g)"

			Token = TRIM(Token) // " <=> " 

			Token = TRIM(Token) // " " // TRIM(AqOrgPhaseChemicalNames(INT(AqOrganicDissolutionData(i,3)))) // " (aq)"

			CALL TRANSCRIPT(TRIM(Token))			

		END DO

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check to make sure the same pair wasn't double specified !!
		IF (HowManyAqOrganicDissolutionReactions .GE. 1) THEN !!!!!!!!!!!!!!!!!

		DO I = 1, HowManyAqOrganicDissolutionReactions

			DO J = 1,HowManyAqOrganicDissolutionReactions
				Match = .FALSE.
				IF (I .EQ. J) CYCLE
				IF (AqOrganicDissolutionData(I,1) .EQ. AqOrganicDissolutionData(J,1) ) THEN

					CALL ERROR("You double-specified a AqOrganicDissolution reaction in HydrophilicOrganicDissolution.in. Please fix this.")

				END IF
			END DO
		END DO
		END IF


		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Check Mass Balance !!
		!! along with other matching and !!
		!! last pass type of things.     !!
		IF (HowManyAqOrganicDissolutionReactions .GE. 1) THEN
		DO I = 1, HowManyAqOrganicDissolutionReactions

			!! Figure Mass Balance
			MassTally = GasMolecularMass(INT(AqOrganicDissolutionData(I,1)))

			!! divide by the r.h.s. masses
			MassTally = MassTally / AqOrgMolecularMass(INT(AqOrganicDissolutionData(I,3))) / INT(AqOrganicDissolutionData(I,13))

			!! err if it is not within 0.1% (which I arbitrarily set as a rounding error
			IF (MassTally .GT. 1.001 .OR. MassTally .LT. 0.999) THEN
				Token = TRIM(GasPhaseChemicalNames(INT(AqOrganicDissolutionData(i,1)))) // " (g)"
				Token = TRIM(Token) // " <=> " // TRIM(REAL2STR(AqOrganicDissolutionData(i,13))) // " " // &
						TRIM(AqOrgPhaseChemicalNames(INT(AqOrganicDissolutionData(i,3)))) // " (org)"
				CALL ERROR("It seems that the mass balance is out of whack in AqOrganicDissolution reaction "//TRIM(Token)// &
				           " in HydrophilicOrganicDissolution.in.  This is a problem.")
			END IF

		END DO
		END IF

		!! SCSCSCSC -- Announce Departure
		IF (Scaffolding) CALL Transcript(">>Exiting SetHydrophilicOrganicDissolution()<<")
		
		CLOSE (FH)
		CALL ReturnFileHandle(FH)

		RETURN
	END SUBROUTINE SetHydrophilicOrganicDissolution
