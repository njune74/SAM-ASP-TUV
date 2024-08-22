!! ASP (c), 2004-2017, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ChemistryParametersInitRoutines.h
!! This contains the routines for making and printing the ODE,
!! Jacobian, and sensitivity arrays for the gas chemical mechanism.
!! Also includes a subroutine for reading time-varying photolysis 
!! rates from a file.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY			                                     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 07/24  2006   Matt Alvarado     Change "Pointer => NULL()" to	     !!
!!				     NULLIFY(POINTER) to fit pgf90	     !!
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 01/03  2013   Matt Alvarado     Fixed second order "+ M" reactions from   !!
!!                                   MCM v3.2 to act as pseudo-1st order     !!
!! 01/04  2013   Matt Alvarado     Fixed PrintODESet to allow                !!
!!                                   reaction numbers greater than 1000      !!
!! 01/04  2013   Matt Alvarado     Made MakeODESet print reaction and        !!
!!                                   rxn number to Transcript.TXT            !!
!! 01/18  2013   Matt Alvarado     Removed hardwired CO + OH reaction and    !!
!!                                  replaced with a K1 + K2 reaction         !!
!! 01/18  2013   Matt Alvarado    Addeda thermal 3rd order rate constant     !!
!! 01/24  2013   Matt Alvarado    Expanded Product arrays in MakeODESet to 40!!
!!                                  and increased size of string variables   !!
!! 01/25  2013   Matt Alvarado    Converted                                  !!
!!                                 + RO2_T reactions to special first order  !!
!! 01/25  2013   Matt Alvarado    Added reactions like HO2 + HO2 + H2O       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		     !!
!! 1. SUBROUTINE MakeODESet(InputDeck,ODEs,Size,ChemSize,Delimitor,Phase)     
!!		NOTE:	Any changes made to rate constant functions here 
!!				must be reflected in subroutine 
!!				ResetGasPhaseReactionRates 
!!				in EulerianEnvironment/GasPhaseChemistry.h
!! 2. SUBROUTINE PrintODESet (ODEs,Size,Phase,OutputFile)
!! 3. SUBROUTINE MakeJacobian(Jac, ODEs, Size, Phase)
!! 4. FUNCTION DerivativeOfODEElement(InEl, dChem) RESULT (OutEl)
!! 5. SUBROUTINE PrintJacobian (Jac,Size,Phase,OutputFile)
!! 6. FUNCTION CopyODEElement(InEl) RESULT (OutEl)
!! 7. SUBROUTINE PrintODEElement(PrintEqTerm)
!! 8. SUBROUTINE MakeSensitivityMatrix
!! 9. FUNCTION ParamDerivativeOfODEElement
!! 10. SUBROUTINE PrintSensitivityMatrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Reactions to be Parsed are stored in this input deck.		!!
!! Chemicals used therein must be detailed in the GasPhaseChems.in	!!
!! deck as well.							!!
!! Size is the number of elements in the ODE deck			!!
!! ChemSize is the number of total chemicals in this phase	        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MakeODESet(InputDeck,ODEs,Size,ChemSize,Delimitor,Phase)

	USE InfrastructuralCode
	USE ModelParameters
	USE GridPointFields, ONLY : AllocateGasReactionGrid
	
	! CMB add
	!use StringIO, ONLY : fortranStr2Float
	! end CMB
	
	IMPLICIT NONE

	!! Input Variables
	CHARACTER*(*)	  :: InputDeck, Delimitor
	TYPE (DiffEqList) :: ODEs(:)
	INTEGER	          :: ChemSize,Size,Phase

	!! A Linked List Type Used Only in this Subroutine
	TYPE ReactRate
	    integer	  :: NumbBodies	! Number of bodies in the reaction
	    integer	  :: NumbRateFacs
	    real*8	  :: Rate(20)	  ! Elements needed for rate
	    type (ReactRate), pointer :: Next
	END TYPE ReactRate

	!! Local Variables
	INTEGER		  :: FH, allocation_error, I, J
	INTEGER		  :: Reacts(10),Prods(40),ReactionType,NumbRateFacs
	INTEGER		  :: NumbProducts,NumbReactants,NumbReactions, & 
                             MaxNumbRateFacs
	REAL*8		  :: ReactFactors(10),ProdFactors(40), II, RateFacs(20)
	character (len=1024)      :: FullLine,CurrLine
	character (len=1024)	  :: Token,Reactants,Reactant
	character (len=1024)	  :: Products,Product
	character (len=1024)	  :: Rates,Rate
	TYPE(DiffEqTerm), POINTER :: Current, New
	TYPE(ChemTerm),   POINTER :: CurrentChem, NewChem
	TYPE(ReactRate),  POINTER :: FirstRate,CurrentRate,NewRate

	integer :: iter ! CMB add, debugging
	
	logical, parameter :: scaffolding = .false.	! cmb in
	!LOGICAL, PARAMETER		  :: Scaffolding = .TRUE.	! cmb out

	!! SCSCSCSCSC
	IF (Scaffolding) THEN
		CALL TRANSCRIPT("")
		CALL TRANSCRIPT("**** Entering MakeODESet() ****")
		CALL TRANSCRIPT("")
	END IF

	!! Initialize Reaction-Array relevent stuff
	NumbReactions   = 0
	MaxNumbRateFacs = 0
	NULLIFY(FirstRate)
	!FirstRate	   => Null()

	FH = GetFileHandle()
        OPEN(UNIT=FH,FILE=TRIM(InputDeckSubDir)//TRIM(InputDeck),STATUS='OLD')

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! GET as many REACTIONS from the input deck as it has !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DO 

	   !! Get a new data line (non-comment) from the input deck
10	   CurrLine = GetLine(FH)

	   !! Exit the DO loop if the end of the file has been reached
	   IF (CurrLine .EQ. EOF) EXIT

	   NumbReactions = NumbReactions + 1
	   FullLine = CurrLine

	   !! Break the input line into three parts
	   CALL GetToken(CurrLine,"=>", Reactants)
	   CALL GetToken(CurrLine,";",  Products)
	   Rates = Trim(CurrLine)
	   	   !print *, '****** REACTANTS (trimmed) = "', reactants(1:len(trim(reactants))), '".'
		   !print *, '****** PRODUCTS (trimmed) = "', products(1:len(trim(products))), '".'
	   	   !print *, '****** RATES (trimmed) = "', rates(1:len(trim(rates))), '".'
	   !! Check to make sure the line parsed plausibly
	   ! Modified by CMB (AER, Inc): LEN_TRIM can be broken
	   if (len(trim(striptoken(reactants))) == 0 .or. &
	           len(trim(striptoken(products))) == 0 .or. &
			   len(trim(striptoken(rates))) == 0) then
			   
	   !IF (Len_Trim(StripToken(Reactants)) == 0 .OR.		&
	!	Len_Trim(StripToken(Products))  == 0 .OR.		&
	!	Len_Trim(StripToken(Rates))		== 0)		&
		CALL ERROR("The following input Line to MakeODESet() from "//TRIM(InputDeck)//" Parsed incorrectly:"// &
		"    "//TRIM(FullLine))
	    endif
		
        !!Print Reaction to transcript with reaction number
		! cmb: comment out
        !CALL TRANSCRIPT("Rxn #"//TRIM(INT2STR(NumbReactions))//"    "//TRIM(FullLine))
	   
        !!!!!!!!!!!!!!!
	   !! REACTANTS !!
	   !!!!!!!!!!!!!!!
	
	   ReactionType  = 0		! the number of bodies in the reaction
	   NumbReactants = 0

	   !! Tokenize the Reactants (up to three)
	   DO I = 1, 3

		ii = 1.
		!print *, 'CMB: Delimitor = "', delimitor, '".'
		CALL GetToken(Reactants,Delimitor,Reactant)
		!IF (len_trim(Reactant) == 0) EXIT
		if (len(trim(reactant)) == 0) exit
		ReactionType = ReactionType + 1		
                ! Count the number of bodies in the reaction

		!! Check Token for Leading Factor
		CALL GetToken(Reactant," ",Token)

		IF (IsReal(Token)) THEN
			ii    = STR2REAL(Token)			! Store
			Token = StripToken(Reactant)	! Push to Chemical
		END IF

		!! If it not the generic M, then go ahead
		IF (Trim(StripToken(Token)) .NE. 'M' .AND.		&
			Trim(StripToken(Token)) .NE. 'm') THEN
			NumbReactants = NumbReactants + 1			
		IF (NumbReactants > 2) &
		    CALL ERROR("Too many reactants specified in the following reaction:"//&
	            "  "//TRIM(FullLine))

		    Reacts(NumbReactants) = FindChem(Trim(StripToken(Token)), Phase)
		    ReactFactors(NumbReactants) = ii
		END IF

		! CMB (AER, Inc): More replacements of len_trim
		!IF (len_trim(StripToken(Reactants)) == 0) EXIT
		if (len(trim(striptoken(reactants))) == 0) exit
	   END DO

	   !! If there is more reactant text in the pipe, then get confused
	   if (len(trim(striptoken(reactants))) .ne. 0) &
	   !IF (len_trim(StripToken(Reactants)) .NE. 0) &
		CALL WARN ("  Warning!! The following Reactant: ("//TRIM(StripToken(Reactants))// &
		     ") Appeared to be Extra in Line... "//"  >>"//TRIM(FullLine)//"<<  Ignoring and Continuing...")

	   !!!!!!!!!!!!!!
	   !! PRODUCTS !!
	   !!!!!!!!!!!!!!
	   NumbProducts = 0
	   
	   !print *, 'CMB: Tokenized reactants'
	   !print *, '------------------------'
	   !do iter = 1, NumbReactants
	   !	   print *, '    Reactant[', iter, '] = "', &
	!	   	   	   reacts(iter), '".  Factor = ', ReactFactors(iter), '.'
	 !  enddo
	   
	   
	   !! Tokenize the Products (up to three)
	   DO
		ii = 1.
		CALL GetToken(Products,Delimitor,Product)
		!IF (len_trim(Product) == 0) EXIT
		IF (len(trim(Product)) == 0) EXIT
		
		!! Check Token for Leading Factor
		CALL GetToken(Product," ",Token)
		!print *, '--- > Token = ', token
		IF (IsReal(Token)) THEN
			ii    = STR2REAL(Token)			! Store
			Token = StripToken(Product)	    ! Push to Chemical
		END IF

		!! If it not the generic M, then go ahead
		IF (Trim(StripToken(Token)) .NE. 'M' .AND.		&
			Trim(StripToken(Token)) .NE. 'm') THEN
			NumbProducts   = NumbProducts + 1
			Prods(NumbProducts) = FindChem(Trim(StripToken(Token)), Phase)
			ProdFactors(NumbProducts) = ii
		END IF

		!IF (len_trim(StripToken(Products)) == 0) EXIT
		IF (len(trim(StripToken(Products))) == 0) EXIT
	   END DO

	   !! If there is more reactant text in the pipe, then get confused
	   if (len(trim(striptoken(products))) .ne. 0) &
	   !IF (len_trim(StripToken(Products)) .NE. 0) &
		CALL WARN("Warning!! The following Product: ("//TRIM(StripToken(Products))//") Appeared to be Extra in Line... "//  &
		"  >>"//TRIM(FullLine)//"  Ignoring and Continuing...")

				 !print *, 'CMB: Tokenized products'
				!	   print *, '------------------------'
				!	   do iter = 1, NumbProducts
				!	   	   print *, '    Products[', iter, '] = "', Prods(iter), '".  Factor = ', ProdFactors(iter), '.'
				!	   enddo		
	   !!!!!!!!!!!
	   !! RATES !!
	   !!!!!!!!!!!

	   !! Count the number of rate terms in Rates
	   J = 1
	   !print *, 'CMB: Trying to count the # of rate terms...'
	   !print *, 'CMB: prior to that if statement: '
	   !print *, 'CMB: # products = ', NumbProducts
	   !	   print *, 'CMB: len(trim(Rates)) = "', rates(1:len(trim(rates))), '".'
	   IF (Rates(Len_Trim(StripToken(Rates))+1:Len_Trim(StripToken(Rates))+1) == ';') &
	    Rates(Len_Trim(StripToken(Rates))+1:Len_Trim(StripToken(Rates))+1) = ' '
	   ! IF (Rates(Len(Trim(StripToken(Rates))) + 1:Len(Trim(StripToken(Rates)))+1) == ';') &
	   ! 	Rates(Len(Trim(StripToken(Rates))) + 1:Len(Trim(StripToken(Rates)))+1) = ' '
	   !print *, 'CMB: after that if statement'
	   !print *, 'CMB: # products = ', NumbProducts
	   !print *, 'CMB: len(trim(Rates)) = "', rates(1:len(trim(rates))), '".'
	   DO I = 1, len(Rates)
	   	   !print *, '    CMB: i = ', i
		IF (Rates(i:i) .EQ. ';') then
			!print *, '        CMB: found semicolon at i = ', i
			J = J+1 
			endif
           END DO
		   
		   ! CMB :try removing this?
           !J = J-1 !mja, 09/16/2010
		   ! end CMB mods
	!print *, 'CMB: j = ', j	   
	   DO I = 1, 20
		RateFacs(I) = 0.
	   END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!! N.B. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!! Any Change to these reaction  !!!
	!!! parsings must be mirrored in  !!!
	!!! both the input file and in    !!!
	!!! the RecalculateReactionRates()!!!
	!!! subroutine and its sub-       !!!
	!!! subroutines.  DO NOT FORGET!  !!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Gas Phase Reactions Have Several Types: Photolytic, Two-Body, !!
!! Three-Body, each of which have different functional forms for !!
!! their respective reaction rates.  These are different for	 !!
!! each phase.							 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           IF (Phase .EQ. GasPhase) THEN
		
	   SELECT CASE (ReactionType)

	   !! An error in parsing the reaction type could have occured.  
           !! Check and exit.
	     CASE (:0)
		CALL ERROR("The following line from  "//Trim(InputDeck)//" produced an"// &
		" impossible reaction type (too few):   >>"//TRIM(FullLine))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  For ONE-BODY REACTIONS, assume that the reaction rate is proportional !!
!!  to concentration. Additional dependencies may later be added, such as !!
!!  interactive radiation (11/01).  So:					  !!
!!                1. Reaction Rate (sec^-1)				  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     CASE (1)

	        !Check for PAN, Degradation, or Heterogeneous flag
		CALL GetToken(Rates,";",Rate)

		!print *, '*** CMB: Rate = "', rate, '".'
		!print *, '*** CMB: Trimmed rate = "', rate(1:len(trim(rate))), '".'
		!print *, '*** CMB: StripToken(rate) = "', stripToken(rate), '".'
		!print *, '*** CMB: STR2REAL(stripToken(Rate)) = ', str2real(striptoken(rate)), '.'
		!print *, '*** CMB: INT of this = ', int(str2real(striptoken(rate))), '.'
		!print *, '*** CMB: J = (number of rate terms) = ', j
		
		!Photolysis Case (Flag = 0)
		IF (INT(STR2REAL(StripToken(Rate))) .EQ. 0) THEN
		        !! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 2	 
			!WRITE(*,*) 'J =', J
			IF (J .NE. NumbRateFacs) &
			CALL ERROR("The following line from "//Trim(InputDeck)//		 &
			" had the wrong number of rate terms.  Expected "//	 &
			 TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
			 //"<<")
                        !WRITE(*,*) Rate
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 1 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
				
			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)
                                !WRITE(*,*) Rate
				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 2 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO


		!NOT USED: Old PAN Degradation Case (Flag = 1)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 1) THEN	
	                !! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 6

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 3 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
				      CALL ERROR("The following line from "//Trim(InputDeck)//	&
						 " 4 apparently has a bad rate expression ("//		&
						 Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!Isomerization Case (Flag = 2)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 2) THEN	
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 5 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 6 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
			
		!Heterogeneous case (Flag = 3)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 3) THEN
		        !! This number should match the number of rate inputs 
                        !! for this reaction type
		        NumbRateFacs = 4

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 7 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 8 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO

		!Equilibrium degradation case (Flag = 4)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 4) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 8

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 9 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
				      CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 10 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO

		!Fixed Rate Case (Flag = 5)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 5) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 2

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 11 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
			
			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
			            RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
				    CALL ERROR("The following line from "//Trim(InputDeck)//	&
					 " 12 apparently has a bad rate expression ("//		&
					Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!HNO4 case - 
                !this is hard-wired in ResetGasPhaseReactionRates (Flag = 6)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 6) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 1

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 13 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
			
		!N2O5 case - 
                ! this is hard-wired in ResetGasPhaseReactionRates (Flag = 7)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 7) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 1

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 14 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
			
		!O1D case - 
                ! this is hard-wired in ResetGasPhaseReactionRates (Flag = 8)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 8) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 1

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 15 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
			
		!O + O2 + M case 
                ! - this is hard-wired in ResetGasPhaseReactionRates (Flag = 9)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 9) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 1

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 16 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
			
		!RAD* + O2 case - 
                ! this is hard-wired in ResetGasPhaseReactionRates (Flag = 10)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 10) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 1

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 17 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
			
			!MCM v3.2 PAN Degradation Case (Flag = 11)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 11) THEN	
	                !! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 6

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 18 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
				      CALL ERROR("The following line from "//Trim(InputDeck)//	&
						 " 100 apparently has a bad rate expression ("//		&
						 Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO	
		!RO2_T pseudo-first order case (Flag = 12) (MJA, 01-25-2013)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 12) THEN	
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 19 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 20 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
		END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For TWO-BODY REACTIONS (thermal or pseudo-second-order pressure dependent)!!
!!  assume the form:							     !!
!!			k = A * T^B exp(C/T)^D			 	     !!
!!  All units are appropriate combinations of moles, cm, seconds, and K      !!
!!  So, specify:  1. "A": Reaction Rate at 298				     !!
!!		  2. "B": Direct Temperature Dependence			     !!
!!		  3. "C": Arrhenius Temperature Dependence		     !!
!!		  4. "D": Arrhenius Off-Factor for Temperature Depedence     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     CASE (2)

		!Check for Nitrate flags (first term 1 or 2)
		CALL GetToken(Rates,";",Rate)			
			
		!Normal 2-Body Case
		IF (INT(STR2REAL(StripToken(Rate))) .EQ. 0) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 4

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
				   " 21 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN	
			               RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 22 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!Alkyl-Nirate Case (Flag = 1)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 1) THEN
		        !! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 9

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 23 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				     RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
				     CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 24 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!Corresponding to nitrate case (Flag = 2)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 2) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 9

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")
				
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 25 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF


			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
					RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 26 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!K = K1 + K2
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 3) THEN
			!! This number should match the number of rate inputs 
                        ! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 27 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN	
			               RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 28 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
		
		!k = k1 + k2*M case (Flag = 4)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 4) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 29 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN	
			               RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 30 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO

		!k = k1 + k3*M(1+k3*M/k2) case (Flag = 5)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 5) THEN
		!! This number should match the number of rate inputs for 
                !! this reaction type
			NumbRateFacs = 7

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 31 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
					RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 32 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!Heterogeneous case (Flag = 6)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 6) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 2

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 33 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 34 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO

		!k = A*exp(B/T)*(1-(1/(1+C*exp(D/T)))) case (Flag = 7)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 7) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 35 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   "36 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO

                !k = A*exp(B/T)*(1/(1+C*exp(D/T))) case (Flag = 8)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 8) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 37 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 38 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
	
		!k = A*exp(B/T)*(1-C*exp(D/T)) case (Flag = 9)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 9) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 39 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 40 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
		!Flag 10, k = C_H2O*(A1*exp(C1/T)+Mgas*A2*exp(C2/T)	
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 10) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 41 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 42 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO
			
		!Flag 11, k = TMP3*TMP4/100
                !!     TMP3 = A*exp(B/T)
                !!     TMP4 = [C/T + D*PRES+F] 
                !!     PRES = ATMOSPHERIC PRESSURE IN PASCAL
                ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 11) THEN
			!! This number should match the number of rate inputs 
                        !! for this reaction type
			NumbRateFacs = 6

			IF (J .NE. NumbRateFacs) &
				CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 43 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF

			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				       RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 44 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
				END IF
			END DO				
		END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!For THREE-BODY REACTIONS (kinetically pressure-dependent), assume the form:!!
!!			k = (k_0 * k_inf) / (k_0 + k_inf) * F		     !!
!!  With appropriate definitions of F.  So, specify with appropriate units:  !!
!!	      1. k_0							     !!
!!		  2. Exponent for k_0 Temp Dependence			     !!
!!		  3. k_inf						     !!
!!		  4. Exponent for k_inf Temp Dependence			     !!
!!		  5. F_c						     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     CASE (3)

		!Check for 3rd order rate constant flag (first term  = 1.0)
		CALL GetToken(Rates,";",Rate)		
			
		!If the first term is below 1, 
                ! assume it is a pseudo-2nd order reaction (JPL Style)
		IF(INT(STR2REAL(StripToken(Rate))) .EQ. 0) THEN
			
			!! This number should match the number of rate inputs 
                        !!for this reaction type
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR(" The following line from "//Trim(InputDeck)//		&
				" had the wrong number of rate terms.  Expected "//	&
				TRIM(INT2STR(NumbRateFacs))//"  >>"//TRIM(FullLine)//"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
				" 45 apparently has a bad rate expression ("//		&
				Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF
				
			DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
				   RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						" 46 apparently has a bad rate expression ("//Trim(StripToken(Rate))  &
						//"):  "//TRIM(FullLine))
				END IF
			END DO
			
		!Use third-order rate constant (as in O + O2 + M)
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 1) THEN
			NumbRateFacs = 3

			IF (J .NE. NumbRateFacs) &
				CALL ERROR(" The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//"  >>"//TRIM(FullLine)//"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
		        ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   "47 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
	        	END IF
				
		        DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
  				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 48 apparently has a bad rate expression ("//Trim(StripToken(Rate))  &
						   //"):  "//TRIM(FullLine))
				END IF
			END DO

		!Use third-order thermal rate constant 
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 2) THEN
			NumbRateFacs = 5

			IF (J .NE. NumbRateFacs) &
				CALL ERROR(" The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//"  >>"//TRIM(FullLine)//"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
		        ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 49 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
	        	END IF
				
		        DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
  				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 50 apparently has a bad rate expression ("//Trim(StripToken(Rate))  &
						   //"):  "//TRIM(FullLine))
				END IF
			END DO
	
	        !THIRD ORDER IUPAC Style
		ELSE IF (INT(STR2REAL(StripToken(Rate))) .EQ. 3) THEN
			NumbRateFacs = 6

			IF (J .NE. NumbRateFacs) &
				CALL ERROR(" The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//"  >>"//TRIM(FullLine)//"<<")

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(1) = STR2REAL(StripToken(Rate))
		        ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
					   " 51 apparently has a bad rate expression ("//		&
					   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
	        	END IF
				
		        DO I = 2, NumbRateFacs

				CALL GetToken(Rates,";",Rate)

				IF (IsReal(StripToken(Rate))) THEN
  				      RateFacs(I) = STR2REAL(StripToken(Rate))
				ELSE
					CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 52 apparently has a bad rate expression ("//Trim(StripToken(Rate))  &
						   //"):  "//TRIM(FullLine))
				END IF
			END DO

		ELSE
			CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " is a third-order reaction that doesn't fit the current code:  "//TRIM(FullLine))
		END IF

 	     CASE (4:)
		CALL ERROR("The following line from "//Trim(InputDeck)//			&
				   " produced an impossible reaction type (too many): >>"// &
				   TRIM(FullLine)//"<<")
           END SELECT
           END IF

	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   !! -- Now Consider Aqueous Phase -- !!
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	   IF (Phase .EQ. AqPhase) THEN

	   SELECT CASE (ReactionType)

	   !! An error in parsing the reaction type could have occured. 
           !! Check and exit.
	     CASE (:0)
		CALL ERROR("The following line from  "//Trim(InputDeck)//" produced an"// &
			   " impossible reaction type (too few):   >>"//TRIM(FullLine))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  For ONE-BODY REACTIONS, assume that the reaction rate is proportional !!
!!  to concentration. Additional dependencies may later be added, such as !!
!!  interactive radiation (11/01).  So:					  !!
!!                1. Reaction Rate (sec^-1)				  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     CASE (1)

		!! This number should match the number of rate inputs 
                !! for this reaction type
		NumbRateFacs = 1

		IF (J .NE. NumbRateFacs) &
			CALL ERROR("The following line from "//Trim(InputDeck)//		 &
					   " had the wrong number of rate terms.  Expected "//	 &
					   TRIM(INT2STR(NumbRateFacs))//":    >>"//TRIM(FullLine) &
					   //"<<")

		DO I = 1, NumbRateFacs

			CALL GetToken(Rates,";",Rate)
	
			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(I) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("  The following line from "//Trim(InputDeck)//	&
				          " 54 apparently has a bad rate expression ("//		&
						  Trim(StripToken(Rate))//"):"//TRIM(FullLine))
			END IF
		END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  For TWO-BODY OR MORE BODY REACTIONS,				     !!
!!  assume the form:							     !!
!!			k = A * exp(B*((To/T)-1))		             !!
!!  All units are appropriate combinations of moles, cm, seconds, and K      !!
!!  So, specify:  1. "A": Reaction Rate at 298				     !!
!!	          2. "B": Temperature Dependence	                     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	     CASE (2:)

		!! This number should match the number of rate inputs 
                !! for this reaction type
		NumbRateFacs = 2

		IF (J .NE. NumbRateFacs) &
			CALL ERROR("  The following line from "//Trim(InputDeck)//		&
					   " had the wrong number of rate terms.  Expected "//	&
					   TRIM(INT2STR(NumbRateFacs))//": >>"//TRIM(FullLine)//"<<")

		DO I = 1, NumbRateFacs

			CALL GetToken(Rates,";",Rate)

			IF (IsReal(StripToken(Rate))) THEN
				RateFacs(I) = STR2REAL(StripToken(Rate))
			ELSE
				CALL ERROR("The following line from "//Trim(InputDeck)//	&
						   " 55 apparently has a bad rate expression ("//		&
						   Trim(StripToken(Rate))//"):"//"    "//TRIM(FullLine))
			END IF
		END DO

           END SELECT
           END IF

	   !!!!!!!!!!!!!!!!!!!!!!!!!
	   !! ALLOCATE RATE TERMS !!!!!!!!!!!!!!!!!!!!!!
	   !! Allocate and load the new rate elements !!
	   ALLOCATE (NewRate, stat = allocation_error)
	   if (allocation_error > 0) CALL ERROR("Allocation of Rate Element NewRate Failed in MakeODESet()")

	   NULLIFY(NewRate%Next)
	   !NewRate%Next		=> Null()
	   NewRate%NumbBodies	 = ReactionType
	   NewRate%NumbRateFacs = NumbRateFacs
	   Do I = 1,NumbRateFacs
	 	NewRate%Rate(I) = RateFacs(I)
	   End Do

	   IF (Associated(FirstRate)) THEN
	 	CurrentRate%Next => NewRate
		CurrentRate      => CurrentRate%Next
	   ELSE
		FirstRate		 => NewRate
		CurrentRate		 => FirstRate
	   END IF

	   IF (NumbRateFacs > MaxNumbRateFacs) MaxNumbRateFacs = NumbRateFacs
	   Nullify(NewRate)

	   !!!!!!!!!!!!!!!!!!!!!!!!
	   !! ALLOCATE ODE TERMS !!
           !!!!!!!!!!!!!!!!!!!!!!!!
	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Each term of the chemical ODE set is held as an element in a linked list !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	   !! Allocate and load the new element
	   ALLOCATE (New, stat = allocation_error)
	   if (allocation_error > 0) CALL ERROR("  Allocation of Linked List Element New Failed in MakeODESet()")
	   NULLIFY(New%Next)
	   !New%Next	  => Null()
	   New%ReactionNumber = NumbReactions

	   !! Allocate, load, and dangle the chemicals involved
	   ALLOCATE (NewChem, stat = allocation_error)
	   if (allocation_error > 0) CALL ERROR("Allocation of Chemical Term fo a Linked List Element New Failed in MakeODESet()")
		
	   NewChem%WhichChem = Reacts(1)
	   NewChem%Stoich    = ReactFactors(1)
	   NULLIFY(NewChem%Next)
	   !NewChem%Next     => Null()
	   New%DerefChem    => NewChem
	   Nullify(NewChem,CurrentChem)

	   !!If there is more than one chemical, then load and dangle the extra
	   IF (ReactionType > 1) THEN
		DO I = 2, NumbReactants	

			!! Make and load the chemical
			ALLOCATE (NewChem, stat = allocation_error)
			if (allocation_error > 0) &
			CALL ERROR("Allocation of Chemical Term fo a Linked List Element New Failed in MakeODESet()")

			NewChem%WhichChem = Reacts(I)
			NewChem%Stoich    = ReactFactors(I)
			Nullify(NewChem%Next)
			CurrentChem		 => New%DerefChem

			!! Place the chemical in the chemical list (ordered by chem number)
			DO 
				IF (Associated(CurrentChem%Next)) THEN
					IF (CurrentChem%Next%WhichChem >= NewChem%WhichChem) THEN

			!! Simply increase the stoicheometric constant by one
			!! if the chemical is repeated.  Insert the chemical
			!! where appropriate if it is not repeated.
						IF (CurrentChem%Next%WhichChem .EQ. NewChem%WhichChem) THEN
							CurrentChem%Next%Stoich = CurrentChem%Stoich + 1
							Deallocate(NewChem)
							Nullify(NewChem)
						ELSE
							NewChem%Next     => CurrentChem%Next
							CurrentChem%Next => NewChem
						END IF
						EXIT
					ELSE
						CurrentChem => CurrentChem%Next
						CYCLE
					END IF
				ELSE
					!! Check for Duplicate Chemicals
					IF (CurrentChem%WhichChem .EQ. NewChem%WhichChem) THEN
						CurrentChem%Stoich = CurrentChem%Stoich + 1
						Deallocate(NewChem)
						Nullify(NewChem)
					ELSE
						CurrentChem%Next => NewChem
					END IF
					EXIT
				END IF							
			END DO
		END DO
	   END IF
	   Nullify(Current,NewChem,CurrentChem)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Add the term to the chemical array at the correct places,                 !!
!! citing one instance                                                       !!
!! for each of the reactants and products				     !!

!! First: REACTANTS
		DO I = 1, NumbReactants

			IF (Reacts(I) > ChemSize) CALL ERROR(" Reactant identified is out of bounds! While parsing "//  &
						trim(InputDeck)//" in MakeODESet()")

			IF (Reacts(I) > -1 .AND. Reacts(I) <= Size) THEN	!! Use -1 to signal a repeated reactant
				!! Walk to end of appropriate chain and attatch element
				IF(Associated(ODEs(Reacts(I))%FirstTerm)) THEN
					Current => ODEs(Reacts(I))%FirstTerm

22					IF (Associated(Current%Next)) THEN
							Current => Current%Next
							GOTO 22
					END IF
					Current%Next => CopyODEElement(New)		!! add a copy to the end of a list
					Current      => Current%Next
				ELSE
					ODEs(Reacts(I))%FirstTerm => CopyODEElement(New)
					Current => ODEs(Reacts(I))%FirstTerm
				END IF

				!! Calculate the correct leading factor
				Current%LeadingFactor = 1.
				IF (NumbReactants > I) THEN
					DO J = I+1,NumbReactants
						IF (Reacts(I) == Reacts(J)) THEN
							Reacts(J) = -1
							Current%LeadingFactor = Current%LeadingFactor + 1.
						END IF
					END DO
				END IF

                                ! CMB (AER): gfortran requires () around the -1
				Current%LeadingFactor = Current%LeadingFactor * (-1.)
				Nullify(Current)
			END IF
		END DO

		!!!!!!!!!!!!!!!!!!!!!!
		!! Second: PRODUCTS !!
		DO I = 1, NumbProducts

			IF (Prods(I) > ChemSize) CALL ERROR (" Reactant identified is out bounds!"//&
					"  While parsing "//trim(InputDeck)//" in MakeODESet()")

			IF (Prods(I) > -1 .AND. Prods(I) <= Size) THEN	!! Use -1 to signal a repeated reactant
				!! Walk to end of appropriate chain and attatch element
				IF(Associated(ODEs(Prods(I))%FirstTerm)) THEN
					Current => ODEs(Prods(I))%FirstTerm

25					IF (Associated(Current%Next)) THEN
							Current => Current%Next
							GOTO 25
					END IF
					Current%Next => CopyODEElement(New)		!! add a copy to the end of a list
					Current     => Current%Next
				ELSE
					ODEs(Prods(I))%FirstTerm => CopyODEElement(New)
					Current => ODEs(Prods(I))%FirstTerm
				END IF
				
				!! Calculate the correct leading factor
				Current%LeadingFactor = ProdFactors(I)
				IF (NumbProducts > I) THEN
					DO J = I+1,NumbProducts
						IF (Prods(I) == Prods(J)) THEN
							Prods(J) = -1
							Current%LeadingFactor = Current%LeadingFactor + 1.
						END IF
					END DO
				END IF
				Nullify(Current)
			END IF
		END DO
		!! Free the memory used by New when done with this
		Deallocate(New)
		Nullify(New)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Go On to Next Reaction  !!
		!! Or Close the Input Deck !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		END DO
		CLOSE (FH)
		CALL ReturnFileHandle(FH)

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!! N.B. !!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!! Any Change to these reaction  !!!
		!!! allocations must be mirrored  !!!
		!!! in both the input file and in !!!
		!!! the RecalculateReactionRates()!!!
		!!! subroutine and its sub-		  !!!
		!!! subroutines.  DO NOT FORGET!  !!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! DEFINE THE REACTION RATE ARRAYS !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SELECT CASE (Phase)

		CASE (GasPhase)
		ALLOCATE(GasPhaseReactionRates(NumbReactions,MaxNumbRateFacs+2),stat=allocation_error)
		if (allocation_error > 0) CALL ERROR("  Allocation of GasPhaseReactionRates failed in MakeODESet()")
		GasPhaseReactionRateArraySize(1) = NumbReactions
		GasPhaseReactionRateArraySize(2) = MaxNumbRateFacs

		CASE (AqPhase)
		ALLOCATE(AqPhaseReactionRateInfo(NumbReactions,MaxNumbRateFacs+2),stat=allocation_error)
		if (allocation_error > 0) CALL ERROR(" Allocation of AqPhaseReactionRateInfo failed in MakeODESet()")
		ALLOCATE(AqPhaseReactionRates(NumbReactions,NumbAqPhaseRateTemps),stat=allocation_error)
		if (allocation_error > 0) CALL ERROR(" Allocation of AqPhaseReactionRates failed in MakeODESet()")
		AqPhaseReactionRateArraySize(1) = NumbReactions
		AqPhaseReactionRateArraySize(2) = MaxNumbRateFacs

		END SELECT
 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Fill the RATE ARRAYS												!!
		!!																	!!
		!! The arrays have the following format:							!!
		!!   Term (1): ReactionType (1 2 or 3 bodied reation)				!!
		!!   Term (2): Number of Rate Factors								!!
		!!   Term (3)...(Term (n)): Rate Elements to be used in calculation !!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		CurrentRate => FirstRate
		I = 1

35		SELECT CASE (Phase)
		CASE (GasPhase)
			GasPhaseReactionRates(I,1) = CurrentRate%NumbBodies 
			GasPhaseReactionRates(I,2) = CurrentRate%NumbRateFacs
		CASE (AqPhase)
			AqPhaseReactionRates(I,1) = CurrentRate%NumbBodies 
			AqPhaseReactionRates(I,2) = CurrentRate%NumbRateFacs
		END SELECT

		DO J = 3,2+CurrentRate%NumbRateFacs
			SELECT CASE (Phase)
			CASE (GasPhase)
				GasPhaseReactionRates(I,J)   = CurrentRate%Rate(J-2)
			CASE (AqPhase)
				AqPhaseReactionRates(I,J)    = CurrentRate%Rate(J-2)
			END SELECT
		END DO

		IF (Associated(CurrentRate%Next)) THEN
			CurrentRate => CurrentRate%Next
			I = I+1
			GOTO 35
		END IF			

		!! Allocate the Gridded Array that will hold the functioning reaction rates over the domain
		IF (Phase .EQ. GasPhase) CALL AllocateGasReactionGrid(NumbReactions)

		!! SCSCSCSCSC
		IF (Scaffolding) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("**** Exiting MakeODESet() ****")
			CALL TRANSCRIPT("")
		END IF

	RETURN
END SUBROUTINE MakeODESet

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
	!! SCSC   Print the ODE set for a given phase  SCSC !!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE PrintODESet (ODEs,Size,Phase,OutputFile)

		USE ModelParameters, ONLY : OutputDeckSubDir
		USE InfrastructuralCode
		IMPLICIT NONE

		!! Input Variables
		TYPE (DiffEqList)   :: ODEs(:)
		INTEGER			    :: Size,Phase
		CHARACTER*(*)	    :: OutputFile
		
		TYPE (DiffEqTerm), POINTER :: Current
		TYPE (ChemTerm),   POINTER :: CurrentChem

		CHARACTER (len = 50) :: ChemName
		INTEGER				 :: FH,I

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Output the calculated ODE set to a given Output File !!
		FH = GetFileHandle()
	    OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFile))

		DO I = 1, Size

			SELECT CASE (Phase)
			CASE (GasPhase)
				ChemName = '['//TRIM(GasPhaseChemicalNames(I))//']'
			CASE (AqPhase)
				ChemName = '['//TRIM(AqPhaseChemicalNames(I))//']'
			END SELECT

			WRITE (FH,'(3a)',advance='no') "d ",Trim(ChemName)," / dt	="

			Current => ODEs(I)%FirstTerm

			!! Loop through elements
			IF (.NOT.Associated(Current)) THEN
				WRITE (FH,'(a)') " 0"
			ELSE
10				IF (Current%LeadingFactor > 0) WRITE (FH,'(a)',advance='no') " +"
				IF (Current%LeadingFactor == -1) WRITE (FH,'(a)',advance='no') " -"
				IF (Current%LeadingFactor == -2) WRITE (FH,'(a)',advance='no') " -2.0"
				IF (Current%LeadingFactor .NE. -2 .AND. Current%LeadingFactor .NE. -1 .AND. Current%LeadingFactor .NE. 1) &
				WRITE (FH,'(F5.3)',advance='no') Current%LeadingFactor 

				WRITE (FH,'(a,i4,a)',advance='no') "k(",Current%ReactionNumber,")"
				
				CurrentChem => Current%DerefChem
20				IF(Associated(CurrentChem)) THEN

					ChemName = ""
					SELECT CASE (Phase)
					CASE (GasPhase)
						ChemName = '['//TRIM(GasPhaseChemicalNames(CurrentChem%WhichChem))//']'
					CASE (AqPhase)
						ChemName = '['//TRIM(AqPhaseChemicalNames(CurrentChem%WhichChem))//']'
					END SELECT

					IF(CurrentChem%Stoich > 1) THEN
						WRITE (FH,'(a,a,f3.1)',advance='no') TRIM(ChemName),'^',CurrentChem%Stoich
					ELSE
						WRITE (FH,'(a)',advance='no') TRIM(ChemName)
					ENDIF

				    IF(Associated(CurrentChem%Next)) THEN
						CurrentChem => CurrentChem%Next
						Goto 20
					ENDIF
				ENDIF 
			
			    IF(Associated(Current%Next)) THEN
					Current => Current%Next
					Goto 10
				ENDIF
				WRITE (FH,'(a)') ""
			END IF

		END DO
		CLOSE (FH)

	END SUBROUTINE PrintODESet


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Once an ODE Set has been created, it is used to make an appropriate !!
	!! Jacobian matrix.  The appropriate Jacobian matrix has the following !! 
	!! form:															   !!
	!!			y1				y2				y3			...			   !!
	!!																	   !!
	!! x1	| d2[x1]/dy1dt	d2[x1]/dy2dt	d2[x1]/dy3dt ...	|		   !!
	!! x2	| d2[x2]/dy1dt	d2[x2]/dy2dt	d2[x2]/dy3dt ...	|		   !!
	!! x3	| d2[x3]/dy1dt	d2[x3]/dy2dt	d2[x3]/dy3dt ...	|		   !!
	!! ...	|	...				...				...				|		   !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE MakeJacobian(Jac, ODEs, Size, Phase)

		USE InfrastructuralCode
		USE ModelParameters
		IMPLICIT NONE

		!! Input Variables
		TYPE (DiffEqList) :: ODEs(:), Jac(:,:)
		INTEGER			  :: Size,Phase

		!! Local Variables
		INTEGER					  :: allocation_error, I, J, K, NonZeroElements
		TYPE(DiffEqTerm), POINTER :: Current, CurrentJac, DerivativeTerm

		LOGICAL, PARAMETER		  :: Scaffolding = .TRUE.
		TYPE(DiffEqTerm), POINTER :: TempTerm

		!! SCSCSCSCSC
		IF (Scaffolding) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("**** Entering MakeJacobian() ****")
			CALL TRANSCRIPT("")
		END IF

		!! Loop over all of the elements, tabulating the Jacobian
		DO I = 1, SIZE
			DO J = 1, SIZE

				Current => ODEs(I)%FirstTerm

				IF (Associated(Current)) THEN

					!! Start a loop over all of the terms in the ODE string
					!! Take the derivative of each term and attach it appropriately
10					DerivativeTerm => DerivativeOfODEElement(Current,J)

					!! Attach the derivative, if not NULL, to the Jacobian Structure
					IF (Associated(DerivativeTerm)) THEN
						CurrentJac => Jac(I,J)%FirstTerm
						IF (.NOT.Associated(CurrentJac)) THEN
							Jac(I,J)%FirstTerm => DerivativeTerm
						ELSE
20							IF(Associated(CurrentJac%Next)) THEN
								CurrentJac => CurrentJac%Next
								GOTO 20
							END IF
							CurrentJac%Next => DerivativeTerm
						END IF
					END IF
					IF (Associated(Current%Next)) THEN
						Current => Current%Next
						GOTO 10
					END IF
				ELSE
				NULLIFY(Jac(I,J)%FirstTerm)
				!Jac(I,J)%FirstTerm => Null()
				END IF
			END DO
		END DO

		!! LSODES wants to know the number of non-zero jacobian array
		!! elements.  Count them here.
		NonZeroElements = 0
		DO I = 1, SIZE
			DO J = 1, SIZE
				IF (Associated(Jac(I,J)%FirstTerm)) NonZeroElements = NonZeroElements + 1
			END DO
		END DO

		SELECT CASE (Phase)
		CASE (GasPhase)
			HowManyNonZeroGasJacTerms	= NonZeroElements

			!!! -- LSODE-related arrays for Gas Phase Chemistry -- !!
			!!! --      these specify the sparcity structure    -- !!
			ALLOCATE (GasIa(HowManyEvolveGasChems+1), stat = allocation_error)
			if (allocation_error > 0) CALL ERROR("Allocation of the GasIa could not proceed in MakeJacobian()")
			ALLOCATE (GasJa(NonZeroElements), stat = allocation_error)
			if (allocation_error > 0) CALL ERROR("Allocation of the GasIa could not proceed in MakeJacobian()")

			K = 1
			DO J = 1, SIZE
			    GasIa(J) = K
				DO I = 1, SIZE
					IF (Associated(Jac(I,J)%FirstTerm)) THEN
						GasJa(K) = I
						K		 = K + 1
					END IF
				END DO
			END DO
			GasIa(HowManyEvolveGasChems+1) = NonZeroElements+1

		CASE (AqPhase)
			HowManyNonZeroAqJacTerms	= NonZeroElements

			!!! -- LSODE-related arrays for Aqueous Phase Chemistry -- !!
			!!! --        these specify the sparcity structure      -- !!
			ALLOCATE (AqIa(HowManyEvolveAqChems+1), stat = allocation_error)
			if (allocation_error > 0) CALL ERROR("Allocation of the AqIa could not proceed in MakeJacobian()")
			ALLOCATE (AqJa(NonZeroElements), stat = allocation_error)
			if (allocation_error > 0) CALL ERROR("Allocation of the AqJa could not proceed in MakeJacobian()")

			K = 1
			DO J = 1, SIZE
			    AqIa(J) = K
				DO I = 1, SIZE
					IF (Associated(Jac(I,J)%FirstTerm)) THEN
						AqJa(K) = I
						K		 = K + 1
					END IF
				END DO
			END DO
			AqIa(HowManyEvolveAqChems+1) = NonZeroElements+1
		END SELECT

		!! SCSCSCSCSC
		IF (Scaffolding) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("**** Exiting MakeJacobian() ****")
			CALL TRANSCRIPT("")
		END IF

	END SUBROUTINE MakeJacobian

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! For the Jacobian calculation, the derivative of each		!!
	!! of the ODE chemical elements is taken with respect		!!
	!! to each of the evolving chemicals.  This function does	!!
	!! that.  IF the derivative is Zero, return NULL			!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION DerivativeOfODEElement(InEl, dChem) RESULT (OutEl)

		USE InfrastructuralCode, ONLY : Transcript

		IMPLICIT NONE

		TYPE (DiffEqTerm), POINTER  :: OutEl, InEl
		INTEGER						:: dChem, Count
		TYPE (ChemTerm), POINTER	:: LeadingChem, TrailingChem

		LOGICAL :: Scaffolding = .FALSE. !.TRUE.

		IF (Scaffolding) CALL TRANSCRIPT("*** Entering DerivativeOfEDEElement() ***")

		IF(.NOT.ASSOCIATED(InEl)) THEN
			NULLIFY(OutEl)
			!OutEl => Null()
			RETURN
		END IF

		OutEl		 => CopyODEElement(InEl)
		LeadingChem  => OutEl%DerefChem
		TrailingChem => LeadingChem

		!! Tally the number of instances of dChem exist

		!! Check to see if dChem is in the first position
		IF (LeadingChem%WhichChem == dChem) THEN
			!! IF Exponent is 1 then remove chemical
			IF (LeadingChem%Stoich == 1) THEN
				OutEl%DerefChem => LeadingChem%Next
				Deallocate(LeadingChem)
			ELSE
				OutEl%LeadingFactor = OutEl%LeadingFactor * LeadingChem%Stoich
				LeadingChem%Stoich = LeadingChem%Stoich - 1.
			END IF
			IF (Scaffolding) CALL TRANSCRIPT("*** Exiting DerivativeOfEDEElement() ***")
			Return
		END IF

		!! Could be that there is only one chem and it is not equal to dChem
		IF (.NOT.Associated(LeadingChem%Next)) THEN
			Deallocate(OutEl)
			Nullify(OutEl)
			IF (Scaffolding) CALL TRANSCRIPT("*** Exiting DerivativeOfEDEElement() ***")
			RETURN
		END IF

		!! Stagger Leading and Trailing Chem
		LeadingChem => LeadingChem%Next

		!! Find the place to take the derivative
10		IF (LeadingChem%WhichChem == dChem) THEN
			!! IF Exponent is 1 then remove chemical
			IF (LeadingChem%Stoich == 1) THEN
				TrailingChem%Next => LeadingChem%Next
				Deallocate(LeadingChem)
			ELSE
				OutEl%LeadingFactor = OutEl%LeadingFactor * LeadingChem%Stoich
				LeadingChem%Stoich = LeadingChem%Stoich - 1.
			END IF
			IF (Scaffolding) CALL TRANSCRIPT("*** Exiting DerivativeOfEDEElement() ***")
			Return
		END IF

		!! Loop again or return a null element
		IF (Associated(LeadingChem%Next)) THEN
			LeadingChem  => LeadingChem%Next
			TrailingChem => TrailingChem%Next
			Goto 10
		ELSE 

			Deallocate(OutEl)
			Nullify(OutEl)
			IF (Scaffolding) CALL TRANSCRIPT("*** Exiting DerivativeOfEDEElement() ***")
			RETURN
		END IF		

	END FUNCTION DerivativeOfODEElement


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS !!
	!! SCSC   Print the Jacobian for a given phase  SCSC !!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE PrintJacobian (Jac,Size,Phase,OutputFile)

		USE ModelParameters, ONLY : OutputDeckSubDir
		USE InfrastructuralCode
		IMPLICIT NONE

		!! Input Variables
		TYPE (DiffEqList)   :: Jac(:,:)
		INTEGER			    :: Size,Phase
		CHARACTER*(*)	    :: OutputFile
		
		TYPE (DiffEqTerm), POINTER :: Current
		TYPE (ChemTerm),   POINTER :: CurrentChem

		CHARACTER (len = 50) :: ChemName,dChem
		INTEGER				 :: FH,I,J

		FH = GetFileHandle()
	    OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFile))

		!! Loop over the Jacobian elements, row by row, and print them to the file
		DO I = 1, Size
						
			WRITE (FH,*) " "
			WRITE (FH,'(a,I5,a)') "*****Begin Row ",I," of Jacobian*****"
			DO J = 1, Size
				Current => Jac(I,J)%FirstTerm

				!! Only print non-zero elements
				IF (Associated(Current)) THEN
				
					SELECT CASE (Phase)
					CASE (GasPhase)
						ChemName = '['//TRIM(GasPhaseChemicalNames(I))//']'
						dChem	 = '['//TRIM(GasPhaseChemicalNames(J))//']'
					CASE (AqPhase)
						ChemName = '['//TRIM(AqPhaseChemicalNames(I))//']'
						dChem    = '['//TRIM(AqPhaseChemicalNames(J))//']'
					END SELECT

					WRITE (FH,'(5a)',advance='no') "d^2 ",Trim(ChemName)," / dt d",Trim(dChem),"	="

				
10					IF (Current%LeadingFactor > 0) WRITE (FH,'(a)',advance='no') " +"
					IF (Current%LeadingFactor == -1) WRITE (FH,'(a)',advance='no') " -"
					IF (Current%LeadingFactor .NE. -1 .AND. Current%LeadingFactor .NE. 1) &
					WRITE (FH,'(F7.3)',advance='no') Current%LeadingFactor 

					WRITE (FH,'(a,i5,a)',advance='no') "k(",Current%ReactionNumber,")"
					
					CurrentChem => Current%DerefChem
20					IF(Associated(CurrentChem)) THEN

						ChemName = ""
						SELECT CASE (Phase)
						CASE (GasPhase)
							ChemName = '['//TRIM(GasPhaseChemicalNames(CurrentChem%WhichChem))//']'
						CASE (AqPhase)
							ChemName = '['//TRIM(AqPhaseChemicalNames(CurrentChem%WhichChem))//']'
						END SELECT

						IF(CurrentChem%Stoich > 1) THEN
							WRITE (FH,'(a,a,f5.3)',advance='no') TRIM(ChemName),'^',CurrentChem%Stoich
						ELSE
							WRITE (FH,'(a)',advance='no') TRIM(ChemName)
						ENDIF

						IF(Associated(CurrentChem%Next)) THEN
							CurrentChem => CurrentChem%Next
							Goto 20
						ENDIF
					ENDIF 
				
					IF(Associated(Current%Next)) THEN
						Current => Current%Next
						Goto 10
					ENDIF
					WRITE (FH,'(a)') ""
				END IF
			END DO
		END DO
		CLOSE (FH)

	END SUBROUTINE PrintJacobian

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This function copies an ODE List Element, Nullifying !!
	!! the %Next pointer.									!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION CopyODEElement(InEl) RESULT (OutEl)

		USE InfrastructuralCode, ONLY : ERROR
		implicit none

		TYPE (DiffEqTerm), POINTER  :: OutEl, InEl
		TYPE (ChemTerm), POINTER	:: NewChem, CurrentInChem, CurrentOutChem
		Integer						:: allocation_error

		IF (.NOT.Associated(InEl)) CALL ERROR(" CopyODEElement() passed an unassociated pointer.  This shouldn't have happened.")

		!! Allocate and load the new element
		ALLOCATE (OutEl, stat = allocation_error)
		if (allocation_error > 0) CALL ERROR("  Allocation of Linked List Element New Failed in CopyODEElement()")

		OutEl%LeadingFactor  = InEl%LeadingFactor
		OutEl%ReactionNumber = InEl%ReactionNumber
		NULLIFY(OutEl%Next)
		!OutEl%Next		    => Null()

		NULLIFY(OutEl%DerefChem)
		!OutEl%DerefChem	=> Null()
		CurrentOutChem  => OutEl%DerefChem
		CurrentInChem   => InEl%DerefChem

10		IF (Associated(CurrentInChem)) THEN

			ALLOCATE (NewChem, stat = allocation_error)
			if (allocation_error > 0) &
			CALL ERROR(" Allocation of Chemical Term fo a Linked List Element New Failed in MakeODESet()")
			NewChem%WhichChem = CurrentInChem%WhichChem
			NewChem%Stoich    = CurrentInChem%Stoich
			NULLIFY(NewChem%Next)
			!NewChem%Next	 => Null()
			
			!! Attach and move on
			IF (Associated(OutEl%DerefChem)) THEN
				CurrentOutChem%Next => NewChem
				CurrentOutChem      => CurrentOutChem%Next
			ELSE
				OutEl%DerefChem => NewChem
				CurrentOutChem  => OutEl%DerefChem
			END IF

			!! Effect a loop
			CurrentInChem => CurrentInChem%Next
			Goto 10
		END IF

		NULLIFY(NewChem,CurrentInChem,CurrentOutChem)

		RETURN
	END FUNCTION CopyODEElement

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
	!! SCSC   This function dumps an ODE Element to the  SCSC !!
	!! SCSC   Standard Out, primarily for debugging.	 SCSC !!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSC !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE PrintODEElement(PrintEqTerm)

		USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, TRANSCRIPT
		implicit none

		Type (DiffEqTerm), POINTER :: PrintEqTerm
		Type (ChemTerm),   POINTER :: CurrChem
		Integer					   :: I

		CALL TRANSCRIPT("")
		CALL TRANSCRIPT("------------------------")
		CALL TRANSCRIPT("  Printing ODE Element")
		CALL TRANSCRIPT("------------------------")


		CALL TRANSCRIPT("LeadingFactor"//TRIM(REAL2STR(PrintEqTerm%LeadingFactor)))
		CALL TRANSCRIPT("ReactionNumber"//TRIM(INT2STR(PrintEqTerm%ReactionNumber)))

		CurrChem => PrintEqTerm%DerefChem
		I = 1

10		IF (Associated(CurrChem)) THEN

		!CALL TRANSCRIPT("WhichChem("//TRIM(INT2STR(I))//")"//TRIM(INT2STR(CurrChem%WhichChem)//  &
		!       " "//TRIM(GasPhaseChemicalNames(CurrChem%WhichChem))))
		!CALL TRANSCRIPT("Stoich("//TRIM(INT2STR(I)//")"//TRIM(REAL2STR(CurrChem%Stoich))))
		I		  = I +1

		IF (Associated(CurrChem%Next)) THEN
			CurrChem => CurrChem%Next
			Goto 10
		END IF
		END IF

		CALL TRANSCRIPT("------------------------")
		CALL TRANSCRIPT("")

		RETURN
	END SUBROUTINE PrintODEElement

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Once an ODE Set has been created, it is used to make an appropriate !!
	!! Sensitivity matrix.  The appropriate Jacobian matrix has the following !! 
	!! form:															   !!
	!!			p1				p2				p3			...			   !!
	!!																	   !!
	!! x1	| d2[x1]/dp1dt	d2[x1]/dp2dt	d2[x1]/dp3dt ...	|		   !!
	!! x2	| d2[x2]/dp1dt	d2[x2]/dp2dt	d2[x2]/dp3dt ...	|		   !!
	!! x3	| d2[x3]/dp1dt	d2[x3]/dp2dt	d2[x3]/dp3dt ...	|		   !!
	!! ...	|	...				...				...				|		   !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE MakeSensitivityMatrix(Sens, ODEs, NumChems, NumRxns, Phase)

		USE InfrastructuralCode
		USE ModelParameters
		IMPLICIT NONE

		!! Input Variables
		INTEGER			  :: NumChems, NumRxns, Phase
		TYPE (DiffEqList) :: ODEs(:), Sens(:,:)
		
		!! Local Variables
		INTEGER					  :: allocation_error, I, J, K, NonZeroElements
		TYPE(DiffEqTerm), POINTER :: Current, CurrentSens, DerivativeTerm

		LOGICAL, PARAMETER		  :: Scaffolding = .TRUE.
		TYPE(DiffEqTerm), POINTER :: TempTerm

		!! SCSCSCSCSC
		IF (Scaffolding) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("**** Entering MakeSensitivityMatrix() ****")
			CALL TRANSCRIPT("")
		END IF

		!! Loop over all of the elements, tabulating the Sensitivity Matrix
		!! I is the chemical index, J is the reaction index
		
		DO I = 1, NumChems
			DO J = 1, NumRxns

				Current => ODEs(I)%FirstTerm

				IF (Associated(Current)) THEN
					
					!! Start a loop over all of the terms in the ODE string
					!! Take the derivative of each term and attach it appropriately
10					DerivativeTerm => ParamDerivativeOfODEElement(Current,J)
					
					!! Attach the derivative, if not NULL, to the Sens Structure
					IF (Associated(DerivativeTerm)) THEN
						CurrentSens => Sens(I,J)%FirstTerm
						IF (.NOT.Associated(CurrentSens)) THEN
							Sens(I,J)%FirstTerm => DerivativeTerm
						ELSE
20							IF(Associated(CurrentSens%Next)) THEN
								CurrentSens => CurrentSens%Next
								GOTO 20
							END IF
							CurrentSens%Next => DerivativeTerm
						END IF
					END IF
					IF (Associated(Current%Next)) THEN
						Current => Current%Next
						GOTO 10
					END IF
				ELSE
				NULLIFY(Sens(I,J)%FirstTerm)
				!Sens(I,J)%FirstTerm => Null()
				END IF
			END DO
		END DO

		!! SCSCSCSCSC
		IF (Scaffolding) THEN
			CALL TRANSCRIPT("")
			CALL TRANSCRIPT("**** Exiting MakeSensitivityMatrix() ****")
			CALL TRANSCRIPT("")
		END IF

	END SUBROUTINE MakeSensitivityMatrix

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! For the Sensitivity matrix calculation, the derivative of each		!!
	!! of the ODE chemical elements is taken with respect		!!
	!! to each of the rate constants.  This function does	!!
	!! that.  IF the derivative is Zero, return NULL.
	!! If not, just return the ODE element. The evaluator function will then 
	!!  ignore the rate constant (set it to 1).			!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION ParamDerivativeOfODEElement(InEl, dRxn) RESULT (OutEl)

		USE InfrastructuralCode, ONLY : Transcript

		IMPLICIT NONE

		TYPE (DiffEqTerm), POINTER  :: OutEl, InEl
		INTEGER						:: dRxn, Count
		TYPE (ChemTerm), POINTER	:: LeadingChem, TrailingChem

		LOGICAL :: Scaffolding = .FALSE. !.TRUE.

		IF (Scaffolding) CALL TRANSCRIPT("*** Entering ParamDerivativeOfEDEElement() ***")

		IF(.NOT.ASSOCIATED(InEl)) THEN
			NULLIFY(OutEl)
			!OutEl => Null()
			RETURN
		END IF

		OutEl => CopyODEElement(InEl)
		IF (OutEl%ReactionNumber .EQ. dRxn) THEN 
			!This means that the derivative is non-zero, so just return element
			IF (Scaffolding) CALL TRANSCRIPT("*** Exiting ParamDerivativeOfEDEElement() ***")
			RETURN
		ELSE
			!In this case, return NULL 
			Deallocate(OutEl)
			Nullify(OutEl)
			IF (Scaffolding) CALL TRANSCRIPT("*** Exiting ParamDerivativeOfEDEElement() ***")
			RETURN
		END IF	

	END FUNCTION ParamDerivativeOfODEElement

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS !!
	!! SCSC   Print the Sensitivity Matrix for a given phase  SCSC !!
	!! SCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCSCS !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE PrintSensitivityMatrix (Sens,NumChems, NumRxns,Phase,OutputFile)

		USE ModelParameters, ONLY : OutputDeckSubDir
		USE InfrastructuralCode
		IMPLICIT NONE

		!! Input Variables
		TYPE (DiffEqList)   :: Sens(:,:)
		INTEGER			    :: NumChems, NumRxns, Phase
		CHARACTER*(*)	    :: OutputFile
		
		TYPE (DiffEqTerm), POINTER :: Current
		TYPE (ChemTerm),   POINTER :: CurrentChem

		CHARACTER (len = 50) :: ChemName, dChem,dRxn
		INTEGER				 :: FH,I,J

		FH = GetFileHandle()
	    OPEN(UNIT=FH, FILE=TRIM(OutputDeckSubDir)//TRIM(OutputFile))

		!! Loop over the Jacobian elements, row by row, and print them to the file
		DO I = 1, NumChems
			WRITE (FH,*) " "
			WRITE (FH,'(a,I5,a)') "*****Begin Row ",I," of Sensitivity Matrix*****"
			DO J = 1, NumRxns
				
				Current => Sens(I,J)%FirstTerm

				!! Print only non-zero elements
				IF (Associated(Current)) THEN
				
					SELECT CASE (Phase)
					CASE (GasPhase)
						dRxn     = '['//TRIM(INT2STR(J))//']'
						dChem	 = '['//TRIM(GasPhaseChemicalNames(I))//']'
					CASE (AqPhase)
						dRxn     = '['//TRIM(INT2STR(J))//']'
						dChem	 = '['//TRIM(GasPhaseChemicalNames(I))//']'
					END SELECT

					WRITE (FH,'(5a)',advance='no') "d^2 ",Trim(dChem)," / dt dk(",Trim(dRxn),")	="

				
10					IF (Current%LeadingFactor > 0) WRITE (FH,'(a)',advance='no') " +"
					IF (Current%LeadingFactor == -1) WRITE (FH,'(a)',advance='no') " -"
					IF (Current%LeadingFactor .NE. -1 .AND. Current%LeadingFactor .NE. 1) &
					WRITE (FH,'(F5.3)',advance='no') Current%LeadingFactor 
					
					CurrentChem => Current%DerefChem
20					IF(Associated(CurrentChem)) THEN

						ChemName = ""
						SELECT CASE (Phase)
						CASE (GasPhase)
							ChemName = '['//TRIM(GasPhaseChemicalNames(CurrentChem%WhichChem))//']'
						CASE (AqPhase)
							ChemName = '['//TRIM(AqPhaseChemicalNames(CurrentChem%WhichChem))//']'
						END SELECT

						IF(CurrentChem%Stoich > 1) THEN
							WRITE (FH,'(a,a,f5.3)',advance='no') TRIM(ChemName),'^',CurrentChem%Stoich
						ELSE
							WRITE (FH,'(a)',advance='no') TRIM(ChemName)
						ENDIF

						IF(Associated(CurrentChem%Next)) THEN
							CurrentChem => CurrentChem%Next
							Goto 20
						ENDIF
					ENDIF 
				
					IF(Associated(Current%Next)) THEN
						Current => Current%Next
						Goto 10
					ENDIF
					WRITE (FH,'(a)') ""
				END IF
			END DO
		END DO
		CLOSE (FH)

	END SUBROUTINE PrintSensitivityMatrix
