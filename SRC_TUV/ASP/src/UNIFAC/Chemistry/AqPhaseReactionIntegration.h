!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! AqPhaseReactionIntegration.h
!! WARNING!: This is where the integration of the aqueous-phase
!! kinetic reaction mechanism would go, but it hasn't been coded yet!
!! Right now, it can only store and reproduce reaction rates.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY								!!
!!										!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	    !!
!! 1. SUBROUTINE SetAqPhaseReactionRates (ReactionRates,ArraySize)  !!
!! 2. FUNCTION GetAqPhaseReactionRates (Temp) RESULT (RateVec)      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Reaction Rates for Aqueous Phase Reactions are Stored !!
	!! in an intepolatable array (against temperature).      !!
	!! The idea is to simply grab the vector in the least	 !!
	!! computationally expensive manner inside the LSODES	 !!
	!! call for each droplet.				 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetAqPhaseReactionRates (ReactionRates,ArraySize)

		Use ModelParameters
		USE GridPointFields,      ONLY : GetTemp
                IMPLICIT NONE

		REAL*8,  INTENT(in) :: ReactionRates(:,:)
		INTEGER, INTENT(in) :: ArraySize(2)

		INTEGER :: I,J
		REAL*8	:: Temp

		!! Loop over all gridcells and reactions
		DO I = 1, ArraySize(1)
			DO J = 1, NumbAqPhaseRateTemps

				!! Use Select to Differentiate Reaction Types
				SELECT CASE (INT(ReactionRates(I,1)))	! Number of Bodies in the reaction

				!! One body reactions (units are per-second)
				CASE(1)
				AqPhaseReactionRates(I,J) = ReactionRates(I,3)
						
				!! Two or More (use Arrheneous form, with Concentration units)
				CASE(2:)
				Temp = GetTemp()
				AqPhaseReactionRates(I,J) = ReactionRates(I,3) * EXP(ReactionRates(I,4)*((298.15/Temp)-1))

				END SELECT
			END DO
		END DO

	END SUBROUTINE SetAqPhaseReactionRates

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! "The idea is to simply grab the vector in the least	 !!
	!! computationally expensive manner inside the LSODES	 !!
	!! call for each droplet." (from last note)		 !!
	!!							 !!
	!! This is the function that grabs this vector (using    !!
	!! linear interpolation).				 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION GetAqPhaseReactionRates (Temp) RESULT (RateVec)

		IMPLICIT NONE

		!! External Variables
		REAL*8 :: RateVec(AqPhaseReactionRateArraySize(1))
		REAL*8 :: Temp

		!! Internal Variables
		INTEGER :: I, LowIndex
		REAL*8  :: Weighting, BinWidth

		BinWidth  = (MaxAqPhaseRateTemp-MinAqPhaseRateTemp)/(NumbAqPhaseRateTemps-1)
		LowIndex  = FLOOR((Temp - MinAqPhaseRateTemp)/BinWidth)
		Weighting = (Temp - MinAqPhaseRateTemp - LowIndex*BinWidth) / BinWidth


		!! Loop Over All Reactions and Fill the Output Vector
		DO I = 1, AqPhaseReactionRateArraySize(1)
			RateVec(I) = AqPhaseReactionRates(I,LowIndex) + Weighting *							&
						 (AqPhaseReactionRates(I,LowIndex+1)-AqPhaseReactionRates(I,LowIndex))
		END DO
		
	END FUNCTION GetAqPhaseReactionRates
