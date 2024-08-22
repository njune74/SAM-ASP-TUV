!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! unifac.f
!! Calculates activity coefficienct for organics in aqueous and hydrophobic phases
!! Developed by Pradeep Saxena and Betty Pun, provided by Rob Griffin (UNH)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																	!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	  !!
!! 1. SUBROUTINE UUNIFAC(NMOL,NFUNC,NU,X,A,RG,QG,Z,TEMP,GAMA)		  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE UNIFAC(NMOL,NFUNC,NU,X,A,RG,QG,Z,TEMP,GAMA)
C
C********************************************************
C Function: unifac
C           SUBROUTINE TO ESTIMATE ACTIVITY COEFFICIENTS &
C           AW OF MULTICOMPONENT MIXTURES FROM UNIFAC
C
C Preconditions required: called by unidriver (b in this case)
C
C Return values:  GAMA(I): ACTIVITY COEFFICIENT OF COMPONENT I ON
C                 MOLE FRACTION SCALE
C
C
c Notes: 1. FOR DIMENSIONAL PURPOSES, WE HAVE ASSUMED THAT
C           MAXIMUM NUMBER OF FUNCTIONAL GROUPS =50 AND
C           MAXIMUM NUMBER OF MOLECULAR ENTITIES =50.
C     
C        2. INPUT SUMMARY: 1) GROUP VOLUME (RG), GROUP SURFACE
C           AREA (QG) OF EACH FUNCTIONAL GROUP IN EACH
C           COMPONENT; 2) MOLE FRACTION OF ALL COMPONENTS
C           3) INTER-FUNCTIONAL GROUP INTERACTION PARAMETERS
C           A(J1,J2) FOR ALL FUNCTIONAL GROUP COMBOS.
C           "COMPONENT" MEANS A MOLECULAR LEVEL ENTITY
c
c        3. AVOID SENDING ZERO MOLE FRACTIONS IN ARRAY X(I).  THE
C           PROGRAM CAN HANDLE X(I) = 0 (IT SETS ACTIVITY COEFF TO 1)
C           BUT BETTER TO AVOID IT. (note from P.Saxena.)
c!!!           (BKP NOTE: activity coefficient is not set to 1 at x = 0)
C
C Revision History: Developed by Pradeep Saxena, EPRI, 95
C                   Revised by Betty Pun, AER, Nov 99 Under EPRI/CARB
C                   funding to comply with models-3 standards
c
c  Revised 1/29/2016 by CMB (AER): Make compatible with gfortran; changed the 
c                                  of conditional statements (.eq. .true. is 
c                                  not valid in gfortran)
c
C ***************************************************************** 
C........Arguments and their description
C
      INTEGER NMOL         ! TOTAL NO. OF MOLECULAR ENTITIES IN THE MIXTURE
      INTEGER NFUNC        ! TOTAL NO. OF FUNCTIONAL GROUPS IN THE MIXTURE
C
      INTEGER NU(NMOL,NFUNC)    ! NU(I,J): VECTOR OF NO. OF A PARTICULAR 
                             ! FUNCTIONAL GROUP J IN MOLECULE I 
                             ! (FROM STOICHIOMETRY)
      REAL X(NMOL)           ! MOLE FRACTION (or amount) OF MOLECULE I
      REAL A(NFUNC,NFUNC)        ! INTERACTION PARAMETER FOR GROUPS J1 & J2
      REAL RG(NFUNC)          ! VAN DER VAAL VOLUME FOR GROUP J
      REAL QG(NFUNC)          ! VAN DER VAAL SURFACE AREA FOR GROUP J
      REAL GAMA(NMOL)        ! calculated activity coefficients
      REAL Z               ! COORDINATION NUMBER FOR THE SOLVENT = 10
      REAL TEMP            ! TEMPERATURE IN DEGREES KELVIN
C 
C ....  Parameters and their descriptions 
C
      REAL R(50)           ! total volume (R) for each molecule
      REAL Q(50)           ! total surface area (Q) for each molecule 
      REAL RL(50)          ! L for each molecule
      REAL RX(50)          ! R(I) * X(I) 
      REAL QX(50)          ! Q(I) * X(I)
      REAL XL(50)          ! RL(I) * X(I)
      REAL PHI(50)
      REAL THETA(50)
      REAL XGM(50)
      REAL XGP(50,50)
      REAL THTAGP(50,50)
      REAL THTAGM(50)
      REAL TTSIM(50)
      REAL TTSIP(50)
      REAL SI(50,50)		!exp(-A/Temp)
      REAL GAMMLN(50)
      REAL GAMPLN(50,50)
      REAL YCLN(50)
      REAL YRLN(50)

      REAL ZHALF           ! 0.5 * Z
      REAL ZOT1
      REAL ZOT2
      REAL ZOT4
      REAL ZOT5
      REAL ZOT6
      REAL ZOT7
      REAL ZOT8
      REAL ZOT9
      REAL ZOT10
      REAL ZOT11
      REAL SUMXGP
      REAL SUMXGM
      REAL SUMXL          ! SUM(RL(I) * X(I))
      REAL SUMJX
      REAL CTEMP
      REAL SUMTGM
      REAL SUMTGP
      REAL SUMQX          ! SUM (Q(I) * X(I))
      REAL SUMRX          ! SUM (R(I) * X(I))
 
      INTEGER I, J, J1, J2  ! loop counters
	LOGICAL SCAF
C
C
C    Subroutine writes to error file, need to open
      SCAF = .FALSE.
      if (scaf) write(*,*) "UNIFAC Check 1"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 1"
	  OPEN(7,FILE='unifac.err')
C
      ZHALF = 0.5 * Z
C
C     **************************************************************
C     COMBINATORIAL PART--THIS IS FOR ORIGINAL UNIFAC
C     NEED TO REVISE IT FOR MODIFIED UNIFAC
C     **************************************************************
C     COMPUTE R AND Q FOR EACH MOLECULE
C
      DO 20 I = 1,NMOL
      R(I) = 0.0
      Q(I) = 0.0
C
C
      DO 10 J = 1, NFUNC
      R(I) = R(I) + FLOAT(NU(I,J)) * RG(J)
      Q(I) = Q(I) + FLOAT(NU(I,J)) * QG(J)
10    CONTINUE
C
20    CONTINUE
C     WRITE(7,9900)Q,R
9900  FORMAT('Q AND R in b',4F12.5)
C
      DO 30 I = 1,NMOL
      RL(I) = ZHALF * (R(I) - Q(I)) - R(I) + 1.
      RX(I) = R(I) * X(I)
      QX(I) = Q(I) * X(I)
      XL(I) = X(I) * RL(I)
30    CONTINUE
C
      SUMRX = 0.
      SUMQX = 0.
      SUMXL = 0.
      DO 50 I = 1,NMOL
      SUMRX = SUMRX + RX(I)
      SUMQX = SUMQX + QX(I)
      SUMXL = SUMXL + XL(I)
50    CONTINUE
C
      ZOT1 = 0.0
      IF (SUMRX .NE. 0.) THEN
      ZOT1 = 1./SUMRX
      ELSE
      WRITE (7,5000)
5000  FORMAT('ZERO DIVISION IN COMBINATORIAL--SUMRX=0. in b')
      STOP
      END IF
      ZOT2 = 0.0
      IF (SUMQX .NE. 0.) THEN
      ZOT2 = 1./ SUMQX
      ELSE
      WRITE(7,5010)
5010  FORMAT('ZERO DIVISION IN COMBINATORIAL-SUMQX=0. in b')
      STOP
      END IF
C
      DO 60 I = 1,NMOL
      PHI(I) = RX(I) * ZOT1
      THETA(I) = QX(I) * ZOT2
60    CONTINUE
C
C     SOLVE FOR THE COMBINATORIAL LN(GAMAC)
      DO 70 I = 1,NMOL
      YCLN(I) = 0.0
      IF (X(I) .EQ. 0.0) THEN
      GO TO 70
      ELSE
      END IF
      IF((PHI(I)) .NE. 0. .AND. (THETA(I)) .NE. 0.) THEN
      ZOT4 = PHI(I)/X(I)
      ZOT5 = THETA(I)/PHI(I)
      YCLN(I) = ALOG(ZOT4) + ZHALF * Q(I) * ALOG(ZOT5) + RL(I) -
     $ ZOT4 * SUMXL
      ELSE
      WRITE(7,5015)
5015  FORMAT('THETA(I) ,OR. PHI(I) .EQ. 0.0 in b')
      STOP
      END IF
70    CONTINUE
      WRITE(7,9991)YCLN(1),YCLN(2)
9991  FORMAT('YCLN in b',2E12.5)

      if (scaf) write(*,*) "UNIFAC Check 2"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 2"
C     *****************************************************************
C     NOW COMPUTE THE RESIDUAL PART--THIS IS FOR ORIGINAL
C     UNIFAC.  NEED TO REVISE IT FOR MODIFIED UNIFAC
C     *****************************************************************
C
C       COMPUTE GROUP MOLE FRACTIONS
C
C       COMPUTE TOTAL GROUP MOLES AND GROUP MOLE FRACTIONS
C       IN PURE COMPOUNDS AND MIXTURE
C
      SUMXGM = 0.0
      DO 150 I = 1, NMOL
C
      SUMXGP = 0.0
      DO 110 J = 1, NFUNC
      SUMXGM = SUMXGM + FLOAT(NU(I,J)) * X(I)
      SUMXGP =SUMXGP + FLOAT(NU(I,J))
110   CONTINUE
C
C     NOW LOOP THRU TO GET MOLE FRACTIONS FOR PURE C COMPOUNDS
      IF (SUMXGP .NE. 0.) THEN
      ZOT6 = 1./ SUMXGP
      ELSE
      WRITE(7,5020)
5020  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMXGP=0. in b')
      STOP
      END IF
      DO 120 J = 1, NFUNC
      XGP(I,J) = FLOAT(NU(I,J)) * ZOT6
C     WRITE(6,7770)I,J,XGP(I,J)
7770  FORMAT('XGP in b',2I10,2E12.5)
120   CONTINUE
C
C
150   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 3"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 3"
C
C
C     NOW WE HAVE TOTAL MOLES OF ALL GROUPS IN THE MIX
C     LOOP ONCE THRU TO COMPUTE GROUP MOLE FRACTION
C     IN THE MIXTURE
      IF (SUMXGM .LE. 0.) THEN
      WRITE(7,5030)
5030  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMXGM=0. in b')
      STOP
      ELSE
      ZOT7 = 1./SUMXGM
      END IF
      DO 180 J = 1,NFUNC
      SUMJX = 0.0
      DO 170 I=1,NMOL
      SUMJX = FLOAT(NU(I,J)) * X(I) + SUMJX
170   CONTINUE
      XGM(J) = SUMJX * ZOT7
C     WRITE(6,7771)J,XGM(J)
7771  FORMAT('XGM in b',I10,E12.5)
180   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 4"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 4"
C
C     COMPUTE TOTAL GROUP AREA & GROUP AREA FRACTIONS
C     IN PURE COMPOUNDS
C
      DO 250 I = 1, NMOL
      SUMTGP = 0.0
C
      DO 210 J = 1, NFUNC
      SUMTGP =SUMTGP + XGP(I,J) * QG(J)
210   CONTINUE
C
C     NOW LOOP THRU TO GET AREA FRACTIONS FOR PURE
C     COMPOUNDS
      IF (SUMTGP .NE. 0.) THEN
      ZOT8 = 1./ SUMTGP
      ELSE
      WRITE(7,5040)
5040  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMTGP=0. in b')
      STOP
      END IF
C
      DO 220 J = 1, NFUNC
      THTAGP(I,J) = QG(J) * XGP(I,J)* ZOT8
C     WRITE(6,7772)I,J,THTAGP(I,J)
7772  FORMAT('THTAGP in b',2I10,E12.5)
220   CONTINUE
C
C
250   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 5"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 5"
C
C     COMPUTE TOTAL GROUP AREA IN MIXTURE
      SUMTGM = 0.0
      DO 260 J = 1, NFUNC
      SUMTGM = SUMTGM + XGM(J) * QG(J)
260   CONTINUE
C
C     NOW WE HAVE TOTAL AREA OF ALL GROUPS IN THE MIX
C     LOOP ONCE THRU TO COMPUTE GROUP AREA FRACTION
C     IN THE MIXTURE
      IF (SUMTGM .LE. 0.) THEN
      WRITE(7,5050)
5050  FORMAT('ZERO DIVISION IN RESIDUAL PART-SUMTGM=0. in b')
      STOP
      ELSE
      ZOT9 = 1./SUMTGM
      END IF
      DO 280 J = 1,NFUNC
      THTAGM(J) = QG(J) * XGM(J) * ZOT9
C     WRITE(6,7773)J,THTAGM(J)
7773  FORMAT('THTAGM in b',I10,E12.5)
280   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 6"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 6"
C
C
C
C     COMPUTE SI VALUES FROM A(J1,J2) AND TEMPERATURE
      IF (TEMP .NE. 0.)THEN
      CTEMP = 1./TEMP
C	WRITE(*,*) "Temp: ", TEMP
      ELSE
      WRITE(7,5060)
5060  FORMAT('ZERO DIVIDE IN RESIDUAL-TEMP=0.0 in b')
      STOP
      END IF
      DO 300 J1 = 1,NFUNC
      DO 300 J2 = 1,NFUNC
      SI(J1,J2) = EXP (-A(J1,J2)*CTEMP)
C     WRITE(6,7774)J1,J2,SI(J1,J2)
7774  FORMAT('SI in b',2I10,E12.5)
300   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 7"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 7"
C
C
C
C     NOW COMPUTE CAPITAL GAMMAS FOR PURE COMPOUNDS
C     AND MIXTURES
C
C     CAPITAL GAMAS FOR FUNCTIONAL GROUPS IN MIXTURES
C     ARE INDEPENDENT OF MOLECULAR ENTITIES
C
      DO 500 J1 = 1, NFUNC
      TTSIM(J1) = 0.0
      DO 400 J2 = 1,NFUNC
      TTSIM(J1) = TTSIM(J1) + THTAGM(J2) * SI(J2,J1)
400   CONTINUE
500   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 8"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 8"
C
C
      DO 550 J1 = 1,NFUNC
      ZOT10 = 0.0
      DO 520 J2 = 1,NFUNC
      IF (TTSIM(J2) .NE. 0.0) THEN
      ZOT10 = ZOT10 + THTAGM(J2) * SI(J1,J2)/TTSIM(J2)
      ELSE
      WRITE(7,5070)
5070  FORMAT('ZERO DIVIDE IN RESIDUAL TTSIM = 0.0 in b')
      STOP
      END IF
520   CONTINUE
      IF (TTSIM(J1) .NE. 0.0) THEN
      GAMMLN(J1) = QG(J1) * (1. - ALOG(TTSIM(J1)) - ZOT10)
C     WRITE(6,7775)J1,GAMMLN(J1)
7775  FORMAT('GAMMLN in b',I10,E12.5)
      ELSE
      WRITE(7,5080)
5080  FORMAT('ZERO ARGUMENT OF LOG IN RES TTSIM in b')
      STOP
      END IF
550   CONTINUE
      if (scaf) write(*,*) "UNIFAC Check 9"
c	IF (SCAF .EQ. .TRUE.) WRITE(*,*) "UNIFAC Check 9"
C
C
C     CAPITAL GAMAS FOR FUNCTIONAL GROUP IN PURE
C     COMPOUNDS ARE DEPENDENT ON MOLECULE--SO
C     WE LOOP THRU EACH MOLECULE-FUNCTIONAL GROUP
C     PAIR
C
      DO 900 I = 1, NMOL
C
C
      DO 700 J1 = 1, NFUNC
      TTSIP(J1) = 0.0
      DO 600 J2 = 1,NFUNC
      TTSIP(J1) = TTSIP(J1) + THTAGP(I,J2) * SI(J2,J1)
600   CONTINUE
700   CONTINUE
C
C
      DO 750 J1 = 1,NFUNC
      ZOT11 = 0.0
      DO 720 J2 = 1,NFUNC
      IF (TTSIP(J2) .NE. 0.0) THEN
      ZOT11 = ZOT11 + THTAGP(I,J2) * SI(J1,J2)/TTSIP(J2)
      ELSE
      WRITE(7,5090)
5090  FORMAT('ZERO DIVIDE IN RESIDUAL TTSIP = 0.0 in b')
      STOP
      END IF
720   CONTINUE
      IF (TTSIP(J1) .NE. 0.0) THEN
      GAMPLN(I,J1) = QG(J1) * (1. - ALOG(TTSIP(J1)) - ZOT11)
C     WRITE(6,7776)I,J1,GAMPLN(I,J1)
7776  FORMAT('GAMPLN in b',2I10,E12.5)
      ELSE
      WRITE(7,6000)
6000  FORMAT('ZERO ARGUMENT OF LOG IN RES TTSIP in b')
      STOP
      END IF
750   CONTINUE
C
c
900   CONTINUE
C
C     COMPUTE RESIDUAL PART AND ADD TO
C     COMBINATORIAL PART
C
      DO 1100 I = 1, NMOL
C
      YRLN(I) = 0.0
      DO 1000 J = 1,NFUNC
      YRLN(I) = YRLN(I) + FLOAT(NU(I,J)) * (GAMMLN(J) - GAMPLN(I,J))
1000  CONTINUE
C
C     WRITE(*,*) I, YCLN(I), YRLN(I)
	GAMA(I) = EXP ( YCLN(I) + YRLN(I))
C
1100  CONTINUE
C
      RETURN
      END
