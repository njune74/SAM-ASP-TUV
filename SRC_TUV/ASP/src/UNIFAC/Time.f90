!     path:      $HeadURL: https://svn.aer.com/svn/aer/project/RD/ASP/branches/remove-eulerian-coords/src/Time.f90 $
!     author:    $Author: malvarad $
!     revision:  $Revision: 10800 $
!     created:   $Date: 2011-10-31 09:20:57 -0400 (Mon, 31 Oct 2011) $
!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Time.f90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES						        !!
!! 1. ModelParameters						!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!	
!! 1. SUBROUTINE ResetTime ()										!!
!! 2. SUBROUTINE ResetTimeStepSize (inTimeStepSize)					!!
!! 3. SUBROUTINE StepTime ()										!!
!! 4. REAL*8 FUNCTION GetTime ()									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! CENTRALIZE MODEL TIME !!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Current Time, Time Steps, Etc. !!
	!! are all controlled from here.  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	MODULE Time

		USE ModelParameters, ONLY : Seconds, Minutes, Hours, Days

		IMPLICIT NONE

		PRIVATE

		PUBLIC ::	TimeStep, BeginningTime,		&
					CurrentTime, TimeStepSize,		&
					ResetTime, StepTime,			&
					GetTime, ResetTimeStepSize

		REAL*8, PARAMETER :: BeginningTime = 0.
		REAL*8            :: TimeStepSize

		INTEGER :: TimeStep 
		REAL*8  :: CurrentTime 

	CONTAINS

	SUBROUTINE ResetTime ()

		IMPLICIT NONE

		TimeStep    = 0
		CurrentTime = BeginningTime

	END SUBROUTINE ResetTime


	SUBROUTINE ResetTimeStepSize (inTimeStepSize)

		IMPLICIT NONE

		REAL*8 :: inTimeStepSize

		TimeStepSize = inTimeStepSize

		RETURN
	END SUBROUTINE ResetTimeStepSize


	SUBROUTINE StepTime ()

		IMPLICIT NONE

		TimeStep = TimeStep + 1
		CurrentTime = CurrentTime + TimeStepSize

	END SUBROUTINE StepTime



	REAL*8 FUNCTION GetTime ()
		IMPLICIT NONE

		GetTime = CurrentTime

		RETURN
	END FUNCTION

	END MODULE Time
