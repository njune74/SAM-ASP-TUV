!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Initialize.h
!! This file sets the temperature and pressure grids
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 02/15  2012   Matt Alvarado     Removing most of these grids, making ASP  !!
!!                                 a one-box model or subroutine.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE InitializeGridPointFields (T, Press)			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE InitializeGridPointFields (T, Press)

	  USE InfrastructuralCode, ONLY : Error
	  implicit none

	  !! Define the input variables
	  real*8  :: T, Press

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! Initialize Primary Variables First!!
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  CALL SetGrids()	! This sets volume to 1 cm^3 

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! Set the TEMPERATURE           !!
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          CALL SetTempField (T)

	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! Populate the PRESSURE grid !!
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  CALL SetPressField (Press)
          AllSamePressure = .TRUE.

          RETURN
	END SUBROUTINE InitializeGridPointFields
