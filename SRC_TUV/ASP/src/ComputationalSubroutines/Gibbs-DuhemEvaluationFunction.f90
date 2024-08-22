!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Gibbs-DuhemEvaluationFunction.f90
!! Provides the needed integrad function for binary water activity
!! calculations.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES						        !!
!! 1. NONE									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. FUNCTION EvalIdLogGamma (IonicStr, Q, Stoich, Charge)         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! A portion of the integrated form of the Gibbs-Duhem Equation is			!!
!! in integral from 0 to I of I d log_10 Gamma, which is integrated using	!!
!! DQXGS.  DQXGS requires a function that evaluates the integrand at a		!!
!! given value.  This is it.  This uses the same form of the equations as   !!
!! Resch, 1995.																!!
!! ----------------------------------------------------------------------   !!
!!    dlogG     q*B*I*(1+0.1 I)^(q-1)      0.5107(0.5*Sqrt[I]+0.069(C-1)I^4)!!
!!  I ----- = ------------------------  -  -----------------------------	!!
!!     dI     10ln10[1+B(1+0.1I)^q -B]            (1+C*Sqrt[I])^2			!!
!! ----------------------------------------------------------------------   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This can't be housed within the Module Aerosols and be made EXTERNAL as  !!
!! well.  So it's is floating on its own.									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION EvalIdLogGamma (IonicStr, Q, Stoich, Charge)

	IMPLICIT NONE

	!! External Variables
	REAL*8  :: IonicStr, Q(3), Stoich(4), Charge(4)		! Total Ionic Strength and the Kusik-Meissner Q Value for the I,J pair

	!! Internal Variables
	REAL*8  :: B,C,D,SQRTI
	INTEGER :: I

	SQRTI = SQRT(IonicStr)

!print *, IonicStr,"Q :",Q(1),Q(2),Q(3)
!print *, IonicStr,"S :",Stoich(1),Stoich(2),Stoich(3)
!print *, IonicStr,"C :",Charge(1),Charge(2),Charge(3)

	IF ((Q(2) .EQ. 0.) .AND. (Q(3) .EQ. 0.) .AND. (Stoich(1) .EQ. 0.)) THEN

		B = 0.75-0.065*Q(1)
		C = 1. + 0.055*Q(1)*EXP(-0.023*IonicStr**3.)
		D = (1.+0.1*IonicStr)**(Q(1)-1.)

		EvalIdLogGamma = Q(1)*B*IonicStr*D / (10.*LOG(10.)*(1.0+B*D*(1.0+0.1*IonicStr)-B)) &
						 - 0.5107*(0.5*SQRTI+0.069*IonicStr**4.*(C-1.))					   &
						 / (1.+C*SQRTI) / (1.+C*SQRTI)
	ELSE

!print *, "--Entering EvalIdLogGamma--"

		I = 1
		B = 0.75-0.065*Q(I)
		C = 1. + 0.055*Q(I)*EXP(-0.023*IonicStr**3.)
		D = (1.+0.1*IonicStr)**(Q(I)-1.)

		IF (1.+B*(1.+0.1*IonicStr)**Q(I)-B .LE. 0.0) THEN
			WRITE(*,*) "IonicStr", IonicStr
			WRITE(*,*) "B: ", B, "Q(1): ", Q(I), "Q(2): ", Q(2), "Q(3): ", Q(3)
			write(*,*) 1.+B*(1.+0.1*IonicStr)**Q(I)-B
		END IF

		EvalIdLogGamma = Stoich(I) * Charge(I) *			        & !!-- Leading Factors
						  (LOG10(1.+B*(1.+0.1*IonicStr)**Q(I)-B) -	& 
						  0.5107*SqrtI/(1.0+SqrtI*C))

		I = 2
		B = 0.75-0.065*Q(I)
		C = 1. + 0.055*Q(I)*EXP(-0.023*IonicStr**3.)
		D = (1.+0.1*IonicStr)**(Q(I)-1.)

		EvalIdLogGamma = EvalIdLogGamma + Stoich(I) * Charge(I) *		&				
						 (LOG10(1.+B*(1.+0.1*IonicStr)**Q(I)-B) -		&
						 0.5107*SqrtI/(1.0+SqrtI*C))
!print *, "Eval2 ",EvalIdLogGamma

		I = 3
		B = 0.75-0.065*Q(I)
		C = 1. + 0.055*Q(I)*EXP(-0.023*IonicStr**3.)
		D = (1.+0.1*IonicStr)**(Q(I)-1.)

		EvalIdLogGamma = EvalIdLogGamma - Stoich(I) * Charge(I) *		& 			
						 (LOG10(1.+B*(1.+0.1*IonicStr)**Q(I)-B) -		& 
						 0.5107*SqrtI/(1.0+SqrtI*C))
!print *, "Eval3 ",EvalIdLogGamma

		EvalIdLogGamma = EvalIdLogGamma / Charge(4) / Stoich(4)

!print *, "4Eval: ",IonicStr, EvalIdLogGamma

!print *, "--Exiting EvalIdLogGamma--"

	END IF

!	WRITE(*,*) EvalIdLogGamma


	RETURN
END FUNCTION EvalIdLogGamma
