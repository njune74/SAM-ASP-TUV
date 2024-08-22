!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! InfrastructuralCode.f90
!! This is the main source file for functions useful
!! in programming, like counters, input and output commands,
!! and commands to convert numbers to strings and vice versa
!! This is the source file linking many header files together.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						    !!
!!								    !!
!! Month  Year   Name              Description			    !!
!! 07     2006   Matt Alvarado     Began Update History		    !!
!! 07/19  2013   Matt Alvarado     Made Module private by default   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES				    !!
!! 1. ModelParameters			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES						    !!
!! 1. InputOutputCommands.h					    !!
!! 2. StringCommands.h						    !!
!! 3. Counters.h					            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE InfrastructuralCode

!    use StringIO

PRIVATE
PUBLIC ::  GetLine,              &
           Error,                &
           Warn,                 &
           Transcript,           &
           REAL2STR,             &
           INT2STR,              &
           STR2REAL,             &
           STR2INT,              &
           GetToken,             &
           StripToken,           &
           IsReal,               &
           IsInteger,            &
           SetFileHandleCounter, &
           GetFileHandle,        &
           ReturnFileHandle,     &
           GetParticleID,        &
           SetParticleID,        &
           EOF,                  &
           IsNaN

!! This is the current file handle
INTEGER :: FH

!! This is the current aerosol id
INTEGER :: AID

!! This is the generic end-of-file statement
CHARACTER (len=33), PARAMETER :: EOF = "END-OF-FILE-012343210-ELIF-FO-DNE"

CONTAINS

INCLUDE "InputOutputCommands.h"
INCLUDE "StringCommands.h"
INCLUDE "Counters.h"

! CMB: Modify to account for floating-point issues?
! Nah, keep this for now.
logical function IsNaN(a)
   IMPLICIT NONE
   real ::  a
!   if (a + 1.0 .eq. a) then
   if (a.ne.a) then
       isnan = .true.
   else
       isnan = .false.
   end if
   return
end function

END MODULE InfrastructuralCode
