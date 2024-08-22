!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Pressure.h
!! Contains subroutines for inputting and outputting pressure data,
!! changing temperature and pressure, calcualing M (total gas) 
!! concentration, and expanding grid point volumes.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 02/15  2012   Matt Alvarado     Removing Eulerian grids, making ASP       !!
!!                                 a one-box model or subroutine.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. SUBROUTINE SetPressField (, Press)
!! 2. FUNCTION   GetPress () RESULT (Press)
!! 3. FUNCTION   GetPressFromGrid () RESULT (Press)
!! 4. FUNCTION   GetM () RESULT (M)
!! 5. SUBROUTINE ChangeTempPress (NewTemperature, NewPressure)
!! 6. SUBROUTINE ExpandGPVolume (VolumeScale)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set the initial PRESSURE field.  This will be ammended or	  !!
!! forced eventually to reflect an actual environment.	          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetPressField (Press)

	USE ModelParameters, ONLY : PressureScale

	implicit none

	!! Inputs Variables:
	real*8  :: Press

	GridPress = Press * PressureScale

	RETURN
END SUBROUTINE SetPressField

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retrieve the PRESSURE    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetPress () RESULT (Press)

	implicit none

	Press = GridPress

	return
END FUNCTION GetPress

REAL*8 FUNCTION GetPressFromGrid () RESULT (Press)

	IMPLICIT NONE

	Press = GridPress
	return
END FUNCTION GetPressFromGrid


!!!!!!!!!!!!!!!!
!! Retrieve M !!
!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetM () RESULT (M)
		
	USE ModelParameters, ONLY : PressureScale

	implicit none

	!! Use the Ideal Gas Law to scale a reference M to local values
!	M =	2.688e19 * (GetPress()/1013./PressureScale)/(GetTemp ()/273.)
	M =	2.68758e19 * (GetPress()/1013./PressureScale)/(GetTemp ()/273.)
	return
END FUNCTION GetM

