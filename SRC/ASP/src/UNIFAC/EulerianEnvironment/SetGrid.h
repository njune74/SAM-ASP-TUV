!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! SetGrid.h
!! This file sets grids and their gridp point volumes.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!							      	             !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History		             !!
!! 02/15  2012   Matt Alvarado     Removing most of these grids, making ASP  !!
!!                                 a one-box model or subroutine.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	 !!
!! 1. SUBROUTINE SetGrids ()					 !!
!! 2. SUBROUTINE SetGridPointVolumes ()			         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This subroutine sets the number of gridcells for each    !!
	!! of the grids in the Eulerian Section of the Array.  It   !!
	!! should eventually be updated to read from a deck.	    !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetGrids ()

	  IMPLICIT NONE

	  CALL SetGridPointVolumes

	  AllSameGrid = .TRUE.
	
	END SUBROUTINE SetGrids

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Set Grid Point Volumes !!
        !! 02-15-2012 MJA Always 1 cm^3
	SUBROUTINE SetGridPointVolumes ()

	  USE ModelParameters, ONLY : DomainX, DomainY, DomainZ
	
	  IMPLICIT NONE

	  GridPointVolume = DomainX*DomainY*DomainZ

	  RETURN
	END SUBROUTINE SetGridPointVolumes
