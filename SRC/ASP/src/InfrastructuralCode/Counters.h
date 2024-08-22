!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Counters.h
!! This module provides counters for file handles   
!! and so-forth										
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY								!!
!!										!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:   !!
!! 1. SUBROUTINE SetFileHandleCounter (InitialFileHandle)        !!
!! 2. FUNCTION GetFileHandle() RESULT (CFH)                      !!
!! 3. SUBROUTINE  ReturnFileHandle(CFH)                          !!
!! 4. SUBROUTINE SetParticleID (inAID)                           !!
!! 5. FUNCTION GetParticleID ()	                                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This set of functions is used to provide unique filehandles !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetFileHandleCounter (InitialFileHandle)

		INTEGER :: InitialFileHandle

		FH = InitialFileHandle

	END SUBROUTINE SetFileHandleCounter

	!! -- Get the current file handle -- !!
	INTEGER FUNCTION GetFileHandle() RESULT (CFH)

		CFH = FH
		FH  = FH + 1

	END FUNCTION

	!! -- If the file handle was the last one to be handed out, then	-- !!
	!! -- it may be returned, thus preventing endless increases in FH.	-- !!
	SUBROUTINE  ReturnFileHandle(CFH)

		INTEGER :: CFH

		IF (CFH == FH - 1) FH = FH-1

	END SUBROUTINE ReturnFileHandle

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This set of functions is used to provide unique Aerosol IDs !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SetParticleID (inAID)

		INTEGER :: inAID

		AID = inAID

	END SUBROUTINE SetParticleID

	INTEGER FUNCTION GetParticleID ()

		GetParticleID = AID
		AID		      = AID + 1

	END FUNCTION GetParticleID
