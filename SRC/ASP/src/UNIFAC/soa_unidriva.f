!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! soa_unidriva.f
!! Calls unifac.f for aqueous phase organics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																	!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES												!!
!! 1. unifacparamMattAq.h													!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	  !!
!! 1. SUBROUTINE UNIDRIVA(XPASS,GAMMA,NUM,T)						  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE UNIDRIVA(XPASS,GAMMA,NUM,T)
C
C     XPASS - passed mole fraction array, indexed starting at 0
C     GAMMA - activity coefficient, passed back as act...
C     NUM - number of possible species present - primary + secondary
C       - XPASS of any more than those present are 0
C
      REAL XPASS(0:NUM), TINY, X(0:NUM), SUMX, Z
      REAL GAMMA(0:NUM),T
      INTEGER N, DIMFUN, DIMMOL, I, J, NMOL, NFUNC
C
      INCLUDE 'unifacparamMattAq.h'  
C            uses DIMMOL and DIMFUN
C
C   Notes:
C      1. this routine uses indices of 0 .. n-1 and is the only one that 
C          does not follow the NR tradition of indexing arrays using 1 .. n
C
C      2. only one set of concentrations (X) are inputted 
C          at each call.  RMOL is a one-dimension array.
C       Therefore, Nsolutions = No. of UNIFAC calls = 1.
C
C      3. input X(i) = passed mole fraction of each component, 
C          n = NSP(6), output GAMMA = act coeff.  File input parameters!
C       No file output.
C
C       4. xpass is used to receive x inputs, x is local and modified 
C          (normalized) in unidriver.  while xpass, pointing to x, in typea
C          is not changed because it continues to be used by typea.
C        
C     put mole frac into temp array
C
c
      DO 20 I = 0,NUM-1
         X(I) = XPASS(I)
 20   CONTINUE    
C
C  check that x sums to 1 or renormalize
C
      SUMX = 0.0
      DO 30 J = 0,NMOL-1
         SUMX=SUMX+X(J)
 30   CONTINUE   
C
C
      IF ((SUMX - 1.0) .GT. TINY .OR. (SUMX - 1.0) .LT.
     +        (-1.0*TINY)) THEN
         DO 40 I = 0,NMOL-1
            X(I) = X(I)/SUMX
 40      CONTINUE   
      ENDIF
C
C
C
C	WRITE(*,*) "In Temp: ", T, NUM
	CALL UNIFAC(NMOL,NFUNC,NU,X,A,RG,QG,Z,T,GAMMA)
C
C
C
      RETURN
      END







