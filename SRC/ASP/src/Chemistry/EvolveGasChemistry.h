!! ASP (c), 2004-2017, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! EvolveGasChemistry.h
!! This Module Deals with the Gas Phase Chemical	
!! integration required in the system.  It has a	
!! controlling driver function and definitions of	
!! the LSODES-related subroutines.	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY						            !!
!!									    !!
!! Month  Year   Name              Description				    !!
!! 07     2006   Matt Alvarado     Began Update History			    !!
!! 10/26  2006   Matt Alvarado     Forced LRW to be above 500000            !!
!! 02/15  2012   Matt Alvarado     Removing Eulerian grids, making ASP      !!
!!                                 a one-box model or subroutine.           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	        !!
!! 1. SUBROUTINE SetEvolveGasChemConstants ()				!!
!! 2. SUBROUTINE StepGasChemistry(TimeStepSize,NumTimeSteps)	        !!
!! 3. SUBROUTINE GasODEEvaluator (neq, t, y, ydot)			!!
!! 4. SUBROUTINE GasJacobianEvaluator (neq, t, y, I, ia, ja, pdj)       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
SUBROUTINE SetEvolveGasChemConstants ()

	USE GridPointFields, ONLY : SetYLEN,SetHowManyEvolveGasChems

	IMPLICIT NONE

	GasLIW	 = 31+HowManyEvolveGasChems+HowManyNonZeroGasJacTerms
	YLEN = HowManyGasChems+3
	LRW  = 40 + 34*HowManyEvolveGasChems + 5*HowManyNonZeroGasJacTerms
	IF (LRW .LT. 5000000) THEN
		LRW = 5000000
	END IF
	
	!! Replicate this value in the eulerian module
	CALL SetYLEN(YLEN)
	CALL SetHowManyEvolveGasChems(HowManyEvolveGasChems)

END SUBROUTINE SetEvolveGasChemConstants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This Subroutine Steps forward the Gas Phase Chemical Concentrations	     !!
!! by a certain timestep.  It is intended to be called from the controlling  !!
!! driver function.  It controls access to LSODES, and all of the tolerences !!
!! and related parameters are set or read in here.			     !!
!!								             !!
!! XP, YP, and ZP are the gridpoint numbers in gas grid space.		     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE StepGasChemistry(TimeStepSize,NumTimeSteps)

	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR, Warn, &
                                        Transcript, GetFileHandle
        USE ModelParameters,	 ONLY : ppbscale
	
	USE GridPointFields,	 ONLY : GetGasChemVector,		&
					getM,				&
					UpdateChemicalConcentrations
		
	IMPLICIT NONE

	!! External Variables
	REAL*8  :: TimeStepSize		!! Should be in Seconds
	INTEGER :: NumTimeSteps

	! 8 October 2016 CMB: AS per the migration from lsodes to dlsodes,
	!					  we need to declare and store these common blocks
	!					  on each routine that calls dlsodes.
	real*8 :: rowns, ccmax, el0, h, hmin, hmxi, hu, rc, tn, uround
	integer :: init, mxstep, mxhnil, nhnil, NSLAST, NYH, IOWNS, ICF, IERPJ, &
			   IERSL, JCUR, JSTART, KFLAG, L, LYH, LEWT, LACOR, LSAVF, &
			   LWM, LIWM, METH, MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, &
			   NQ, NST, NFE, NJE, NQU
	COMMON /DLS001/ ROWNS(209), CCMAX, EL0, H, HMIN, HMXI, HU, RC, &
					TN, UROUND, INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, &
					NYH, IOWNS(6), ICF, IERPJ, IERSL, JCUR, JSTART, &
					KFLAG, L, LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, &
					MITER, MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, &
					NFE, NJE, NQU
	real*8 :: con0, conmin, ccmxj, psmall, rbig, seth
	integer :: IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, &
			   IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,&
			   LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,&
			   NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU			
	COMMON /DLSS01/ CON0, CONMIN, CCMXJ, PSMALL, RBIG, SETH, &
					IPLOST, IESP, ISTATC, IYS, IBA, IBIAN, IBJAN, IBJGP, &
					IPIAN, IPJAN, IPJGP, IPIGP, IPR, IPC, IPIC, IPISP, IPRSP, IPA,&
					LENYH, LENYHM, LENWK, LREQ, LRAT, LREST, LWMIN, MOSS, MSBJ,&
					NSLJ, NGP, NLU, NNZ, NSP, NZL, NZU				
	! end CMB
	
	!! Internal Variables
        ! CB: debugging
	LOGICAL, PARAMETER :: Scaffolding = .FALSE.
	!LOGICAL, PARAMETER :: Scaffolding = .TRUE.
        REAL*8  :: BeginTime,EndTime
	INTEGER	::	q, OutFH

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!! LSODES Related Internal Variables !!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!Length of RWORK is >= 20 + (2 + 1./lenrat)*nnz + (11 + 9./lenrat)*neq
	!! where..							    !!
	!! nnz    = the number of nonzero elements in the sparse	    !!
	!!          jacobian (if this is unknown, use an estimate), and      !!
	!! lenrat = the real to integer wordlength ratio (usually 1 in       !!
        !!          single precision and 2 in double precision).	     !!
	!!								     !!
	!!  Length of Y is the Total Number of Gas Chemicals                 !!
        !!      (Evolving + Non-Evolving) Plus 3                             !!
        !!      (those last three elements are the gridcell location)	     !!
	!!								     !!
	!!  Length of IWork is >= 31 + neq + nnz			     !!
	!!								     !!
	!!  Value of LRW is the declared length of RWORK		     !!
	!!								     !!
	!!  Value of LIW is the declared length of IWORK		     !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8, ALLOCATABLE	:: y(:), rwork(:)
	REAL*8			:: t,outtime,rtol,atol
	INTEGER, ALLOCATABLE	:: iwork(:)
	INTEGER			:: itol,itask,istate,iopt,mf,neq
		
	!! Non-LSODES Local Variables
	INTEGER :: I,J,allocation_error
	REAL*8  :: M

	!! SCSCSCSC
	IF (Scaffolding) CALL Transcript(">>>Entering StepGasChemistry<<<")

	IF (.NOT. IfGas) &
	  CALL ERROR ("You are trying to evolve the gas phase chemistry but have set in MELAMInputDeck.in to not "//  &
	  "run gas phase chemistry.  So we did not read the necessary input files.  This won't work.  "// &
	  "(in 	SUBROUTINE StepGasChemistry)")

	!! SCSCSCSCSC
	IF (Scaffolding) THEN
	  CALL Transcript("Calculated array lengths:"//trim(int2str(ylen)))
	END IF

	!! Allocate the LSODES working arrays
	ALLOCATE (Y(ylen), stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of Y could not proceed in StepGasChemistry()")

	ALLOCATE (rwork(LRW), &
			  stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR ("Allocation of RWORK could not proceed in StepGasChemistry()")

	ALLOCATE (iwork(GasLIW), &
			  stat = allocation_error)
	IF (allocation_error > 0) CALL ERROR("Allocation of IWORK could not proceed in StepGasChemistry()")

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set LSODES-Related Variable (Notes Excerpted !!
!! from LSODES file)			        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! itask  = an index specifying the task to be performed.		!!
!!          input only.  itask has the following values and meanings.	!!
!!          1  means normal computation of output values of y(t) at	!!
!!             t = tout (by overshooting and interpolating).		!!
!!          2  means take one step only and return.			!!
!!          3  means stop at the first internal mesh point at or	!!
!!             beyond t = tout and return.				!!
!!          4  means normal computation of output values of y(t) at	!!
!!             t = tout but without overshooting t = tcrit.		!!
!!             tcrit must be input as rwork(1).  tcrit may be equal to	!!
!!             or beyond tout, but not behind it in the direction of	!!
!!             integration.  this option is useful if the problem	!!
!!             has a singularity at or beyond t = tcrit.		!!
!!          5  means take one step, without passing tcrit, and return.	!!
!!             tcrit must be input as rwork(1).				!!
!!								     	!!
!!          note..  if itask = 4 or 5 and the solver reaches tcrit	!!
!!          (within roundoff), it will return t = tcrit (exactly) to	!!
!!          indicate this (unless itask = 4 and tout comes before tcrit,!!
!!          in which case answers at t = tout are returned first).	!!
itask  = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! istate = an index used for input and output to specify the
!!         the state of the calculation.
!!
!!         on input, the values of istate are as follows.
!!          1  means this is the first call for the problem
!!             (initializations will be done).  see note below.
!!          2  means this is not the first call, and the calculation
!!             is to continue normally, with no change in any input
!!             parameters except possibly tout and itask.
!!             (if itol, rtol, and/or atol are changed between calls
!!             with istate = 2, the new values will be used but not
!!             tested for legality.)
!!          3  means this is not the first call, and the
!!             calculation is to continue normally, but with
!!             a change in input parameters other than
!!             tout and itask.  changes are allowed in
!!             neq, itol, rtol, atol, iopt, lrw, liw, mf,
!!             the conditional inputs ia and ja,
!!             and any of the optional inputs except h0.
!!             in particular, if miter = 1 or 2, a call with istate = 3
!!             will cause the sparsity structure of the problem to be
!!             recomputed (or reread from ia and ja if moss = 0).
!!          note..  a preliminary call with tout = t is not counted
!!          as a first call here, as no initialization or checking of
!!          input is done.  (such a call is sometimes useful for the
!!          purpose of outputting the initial conditions.)
!!          thus the first call for which tout .ne. t requires
!!          istate = 1 on input.
!!
!!          on output, istate has the following values and meanings.
!!           1  means nothing was done, as tout was equal to t with
!!              istate = 1 on input.  (however, an internal counter was
!!              set to detect and prevent repeated calls of this type.)
!!           2  means the integration was performed successfully.
!!          -1  means an excessive amount of work (more than mxstep
!!              steps) was done on this call, before completing the
!!              requested task, but the integration was otherwise
!!              successful as far as t.  (mxstep is an optional input
!!              and is normally 500.)  to continue, the user may
!!              simply reset istate to a value .gt. 1 and call again
!!              (the excess work step counter will be reset to 0).
!!              in addition, the user may increase mxstep to avoid
!!              this error return (see below on optional inputs).
!!          -2  means too much accuracy was requested for the precision
!!              of the machine being used.  this was detected before
!!              completing the requested task, but the integration
!!              was successful as far as t.  to continue, the tolerance
!!              parameters must be reset, and istate must be set
!!              to 3.  the optional output tolsf may be used for this
!!              purpose.  (note.. if this condition is detected before
!!              taking any steps, then an illegal input return
!!              (istate = -3) occurs instead.)
!!          -3  means illegal input was detected, before taking any
!!              integration steps.  see written message for details.
!!              note..  if the solver detects an infinite loop of calls
!!              to the solver with illegal input, it will cause
!!              the run to stop.
!!          -4  means there were repeated error test failures on
!!              one attempted step, before completing the requested
!!              task, but the integration was successful as far as t.
!!              the problem may have a singularity, or the input
!!              may be inappropriate.
!!          -5  means there were repeated convergence test failures on
!!              one attempted step, before completing the requested
!!              task, but the integration was successful as far as t.
!!              this may be caused by an inaccurate jacobian matrix,
!!              if one is being used.
!!          -6  means ewt(i) became zero for some i during the
!!              integration.  pure relative error control (atol(i)=0.0)
!!              was requested on a variable which has now vanished.
!!              the integration was successful as far as t.
!!          -7  means a fatal error return flag came from the sparse
!!              solver cdrv by way of prjs or slss (numerical
!!              factorization or backsolve).  this should never happen.
!!              the integration was successful as far as t.
!!
!!          note.. an error return with istate = -1, -4, or -5 and with
!!          miter = 1 or 2 may mean that the sparsity structure of the
!!          problem has changed significantly since it was last
!!          determined (or input).  in that case, one can attempt to
!!          complete the integration by setting istate = 3 on the next
!!          call, so that a new structure determination is done.
!!
!!          note..  since the normal output value of istate is 2,
!!          it does not need to be reset for normal continuation.
!!          also, since a negative input value of istate will be
!!          regarded as illegal, a negative output value requires the
!!          user to change it, and possibly other inputs, before
!!          calling the solver again.
istate = 1	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! iopt   = an integer flag to specify whether or not any optional      !!
!!          inputs are being used on this call.  input only.		!!
!!          the optional inputs are listed separately below.		!!
!!          iopt = 0 means no optional inputs are being used.		!!
!!                   default values will be used in all cases.		!!
!!          iopt = 1 means one or more optional inputs are being used.	!!
iopt   = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! mf     = the method flag.  used only for input.		    !!
!!          the standard choices for mf are..			    !!
!!            mf = 10  for a nonstiff problem,			    !!
!!            mf = 21 or 22 for a stiff problem with ia/ja supplied !!
!!                     (21 if jac is supplied, 22 if not),	    !!
!!            mf = 121 for a stiff problem with jac supplied,	    !!
!!                     but not ia/ja,				    !!
!!            mf = 222 for a stiff problem with neither ia/ja nor   !!
!!                     jac supplied.				    !!
!!          the sparseness structure can be changed during the	    !!
!!          problem by making a call to lsodes with istate = 3.	    !!
mf     =  21 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This is the maximum number of internal steps allowed. !!
	!! The default is 500 and that's not enough.		 !!
	iwork(6) = 2000 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!!!!!!!!!!!!!!!!!!!!
	!! Specify Ia and Ja !!
	DO I = 1,HowManyEvolveGasChems+1
		iwork(30+I) = GasIa(I)
	END DO
	DO I = 1,HowManyNonZeroGasJacTerms
		iwork(31+HowManyEvolveGasChems+I) = GasJa(I)
	END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! atol   = an absolute error tolerance parameter, either a scalar or	!!
!!          an array of length neq.  input only.			!!
!!							       		!!
!!             the input parameters itol, rtol, and atol determine	!!
!!          the error control performed by the solver.  the solver will	!!
!!          control the vector e = (e(i)) of estimated local errors	!!
!!          in y, according to an inequality of the form		!!
!!                      rms-norm of ( e(i)/ewt(i) )   .le.   1,		!!
!!          where       ewt(i) = rtol(i)*abs(y(i)) + atol(i),		!!
!!          and the rms-norm (root-mean-square norm) here is		!!
!!          rms-norm(v) = sqrt(sum v(i)**2 / neq).  here ewt = (ewt(i))	!!
!!          is a vector of weights which must always be positive, and	!!
!!          the values of rtol and atol should all be non-negative.	!!
!!          the following table gives the types (scalar/array) of	!!
!!          rtol and atol, and the corresponding form of ewt(i).	!!
!!									!!
!!             itol    rtol       atol          ewt(i)			!!
!!              1     scalar     scalar     rtol*abs(y(i)) + atol	!!
!!              2     scalar     array      rtol*abs(y(i)) + atol(i)	!!
!!              3     array      scalar     rtol(i)*abs(y(i)) + atol	!!
!!              4     array      array      rtol(i)*abs(y(i)) + atol(i)	!!
!!									!!
!!          when either of these parameters is a scalar, it need not	!!
!!          be dimensioned in the user-s calling program.		!!
!!									!!
!!          if none of the above choices (with itol, rtol, and atol	!!
!!          fixed throughout the problem) is suitable, more general	!!
!!          error controls can be obtained by substituting		!!
!!          user-supplied routines for the setting of ewt and/or for	!!
!!          the norm calculation.  see part iv below.			!!
!!									
!!          if global errors are to be estimated by making a repeated	!!
!!         run on the same problem with smaller tolerances, then all	!!
!!          components of rtol and atol (i.e. of ewt) should be scaled	!!
!!          down uniformly.					        !!
itol   = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 12 April 2016 CMB(AER): Request to run with rtol of 1e-3; let's see what happens
! UPDATE: Had minimal effect on the run time
rtol   = 1.0e-6
!rtol = 1.0e-3

	M = getM()
	atol   = 1.e-6*ppbscale*M

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialize the Concentrations from the local chemical concentrations: !! 
y = GetGasChemVector()

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Tell the program which gridpoint it's operating at: !!
	y(HowManyGasChems+1) = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	y(HowManyGasChems+2) = 1
	y(HowManyGasChems+3) = 1

	! Reset Time
	t         = 0.
	outtime   = 0.

	IF (Scaffolding) CALL TRANSCRIPT("About to loop calls to LSODES")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Integrate Forward the Specified Number of TimeSteps !!
	DO J = 0,NumTimeSteps !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  IF (Scaffolding) then
            CALL TRANSCRIPT("Before LSODES; calling it on GasODEEvaluator and GasJacobianEvaluator")
            print *, 'istate = ', istate
            print *, 'howmanyevolvegaschems = ', howmanyevolvegaschems
          endif
		  
	  ! October 7 CMB (AER, Inc): Replacing LSODES with DLSODES
10	  CALL dLSODES (GasODEEvaluator,HowManyEvolveGasChems,y,t,outtime, &
	                itol,rtol,atol,itask,&
	  	            istate,iopt,rwork,lrw,iwork,GasLIW,GasJacobianEvaluator,mf)	  
!10	  CALL LSODES (GasODEEvaluator,HowManyEvolveGasChems,y,t,outtime, &
!               itol,rtol,atol,itask,&
!	       istate,iopt,rwork,lrw,iwork,GasLIW,GasJacobianEvaluator,mf)
			
	  IF (Scaffolding) CALL TRANSCRIPT("After LSODES : "//INT2STR(ISTATE))
	  !! If the maximum number of timesteps occurred, then warn the user 
          !!(in Scaffolding mode) and push the integration back into LSODES 
          !! until the appropriate end time is reached.
	  IF (ISTATE .EQ. -1 .AND. T < OUTTIME) THEN
		IF (Scaffolding) &
		CALL WARN("Warning! Nominal maximum number of steps in LSODES was exceeded in Gas Phase Integration... Recycling")
		ISTATE = 2
		GOTO 10
	  END IF

	  !! An Error will occur if the array lengths were too small 
          !! (check to make sure IWORK(17) is reasonable!)
	  IF (IWORK(17) > LRW .AND. IWORK(17) < 1000000)  &
		CALL ERROR("The Given Length for LSODES' Working Array RWORK (accessing in StepGasChemistry())"//		&
		" Should have been of size "//TRIM(INT2STR(iwork(17)))//" but was only "//					&
		TRIM(INT2STR(LRW))//" instead.  Fix the assignment in EvolveGasChemistry's SetEvolveGasChemConstants().")
			

	  !! iterate timestep forward
	  outtime = (J+1) * TimeStepSize
	END DO

	IF (Scaffolding) CALL TRANSCRIPT("Before Update Gas Chems")
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Replace the new concentrations into the correct grid cell !!
	CALL UpdateChemicalConcentrations(Y(1:HowManyEvolveGasChems))
	!! SCSCSCSC
	IF (Scaffolding) CALL TRANSCRIPT (">>>Exiting StepGasChemistry<<<")
        DEALLOCATE (Y, RWORK, IWORK, STAT = allocation_error)

RETURN
END SUBROUTINE StepGasChemistry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MELAM -- 2001 -- Donnan Steele			!!
!!							!!
!! The ODE and Jacobian subroutines to be used by	!!
!! LSODES.  They calculate the rate constants and	!!
!! dereference the ODE and Jacobian linked lists	!!
!! contained in the ChemistryParameters module.		!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This program provides the ODE ydot vector to LSODES. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! These take the linked list form of the ODE sets and  !!
!! Jacobian and create the numerical form of the LSODES-!!
!! related functions.					!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GasODEEvaluator (neq, t, y, ydot)

	USE GridPointFields,	 ONLY :	GasPhaseChemicalRates
	USE InfrastructuralCode, ONLY : Transcript

	IMPLICIT NONE

	!! External Variable (prescribed by LSODES)
	REAL*8		:: t, y, ydot
	INTEGER		:: neq
	DIMENSION	:: y(1), ydot(1)

	!! Internal Variables
	REAL*8	:: Term, T_mix, m, b
	INTEGER :: I,GP(3)

	TYPE (DiffEqTerm), POINTER :: Current
	TYPE (ChemTerm),   POINTER :: CurrentChem

        ! CB: Make true for debugging
	LOGICAL, PARAMETER :: Scaffolding =  .FALSE. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Structure of the Y(1) input vecor is taken to be:		     !!
!!    (1 - neq) : evolving chemicals				             !!
!! (neq+1  - HowManyGasChems)			  : non-evolving chems	     !!
!!(HowManyGasChems+1 - HowManyGasChems+3) : Grid Point Specification (x,y,z) !!
!! 02-15-2012 MJA Grid point specification is now not needed,                !!
!! but left in to not screw things up. Always 1, 1, 1                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! SCSCSCSC
	IF (Scaffolding) CALL TRANSCRIPT (">>>Entering GasODEEvaluator()<<<")

	!! The Gas Cell Reference for the Gas Phase Chemistry
	!! is Passed as the first three elements after the GridCell
	GP(1) = INT(y(HowManyGasChems+1))
	GP(2) = INT(y(HowManyGasChems+2))
	GP(3) = INT(y(HowManyGasChems+3))

	!! Initialize variables (set to 0)
	DO I = 1, neq
		ydot(I) = 0.
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! First evaluate Gas Phase Chemical Contributions, if desired !!
	IF (IFGAS) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !print *, 'Evaluating Gas Phase Chemical Contributions'

	   !! Evaluate contributions to each term of Y
	   DO I = 1, neq

	        Current => GasPhaseODESet(I)%FirstTerm
		!! Loop through elements
		IF (Associated(Current)) THEN
			!! Accumulate current derivative in term

10			Term = Current%LeadingFactor *			&
			       GasPhaseChemicalRates(Current%ReactionNumber)
				
			!! Multiply By Chemicals
			CurrentChem => Current%DerefChem
20			IF(Associated(CurrentChem)) THEN
				IF(CurrentChem%Stoich .NE. 1) THEN
					Term = Term * (Y(CurrentChem%WhichChem)**CurrentChem%Stoich)
				ELSE
					Term = Term * Y(CurrentChem%WhichChem)
				ENDIF

				IF(Associated(CurrentChem%Next)) THEN
					CurrentChem => CurrentChem%Next
					Goto 20
				ENDIF
			ENDIF 
			
			!! Add term to overall derivative
			ydot(I) = ydot(I) + Term

			IF(Associated(Current%Next)) THEN
				Current => Current%Next
				Goto 10
			ENDIF
		END IF
	   END DO
	END IF

	!! SCSCSCSC
	IF (Scaffolding) CALL TRANSCRIPT(">>>Exiting GasODEEvaluator()<<<")

	RETURN
END SUBROUTINE GasODEEvaluator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the JACOBIAN from the Linked List		  !!
!! Array Contained in the ChemistryParameters module      !!		
!! The Subroutine, as proscribed by LSODES, returns a     !!
!! vector of the Jacobian, set as the Ith vertical vector !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GasJacobianEvaluator (neq, t, y, I, ia, ja, pdj)

	USE GridPointFields, ONLY: GasPhaseChemicalRates

        !
        ! Modified by CMB (AER): Needs to Use the logical variables from Chemistry
        !                        module 
        !
        !use chemistry, only: ifgas, ifaq, ifsolid

	IMPLICIT NONE

	!! External Variables
	REAL*8		:: t, y, pdj
	INTEGER		:: I, neq, ia, ja
	DIMENSION	:: y(1), ia(1), ja(1), pdj(1)

	!! Internal Variables
	REAL*8	:: GP(3), Term, RunningTotal, T_mix, m, b
	INTEGER :: J

	TYPE (DiffEqTerm), POINTER :: Current
	TYPE (ChemTerm),   POINTER :: CurrentChem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The Structure of the Y(1) input vecor is taken to be:		     !!
!!    (1 - neq) : evolving chemicals				             !!
!! (neq+1  - HowManyGasChems)			  : non-evolving chems	     !!
!!(HowManyGasChems+1 - HowManyGasChems+3) : Grid Point Specification (x,y,z) !!
!! 02-15-2012 MJA Grid point specification is now not needed,                !!
!! but left in to not screw things up. Always 1, 1, 1                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	!! The Gas Cell Reference for the Gas Phase Chemistry	    !!
	!! is Passed as the first three elements after the GridCell !!
	GP(1) = y(HowManyGasChems+1)
	GP(2) = y(HowManyGasChems+2)
	GP(3) = y(HowManyGasChems+3)

	!! Initialize variables
	DO J = 1, neq
		pdj(J) = 0.
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! If Gas Phase Chemistry Occurs !!
	IF (IFGAS) THEN !!!!!!!!!!!!!!!!!!!
            !print *, 'ifgas is true inside of GasJacobianEvaluator'
	   DO J = 1, neq

		RunningTotal = 0
		Current => GasPhaseJacobian(J,I)%FirstTerm

		!! Loop through elements
		IF (Associated(Current)) THEN
			!! Accumulate current derivative in term
10			Term = Current%LeadingFactor *			&
			       GasPhaseChemicalRates(Current%ReactionNumber)
				
			!! Multiply By Chemicals
			CurrentChem => Current%DerefChem
20			IF(Associated(CurrentChem)) THEN
				IF(CurrentChem%Stoich .NE. 1) THEN
					Term = Term * (Y(CurrentChem%WhichChem)**CurrentChem%Stoich)
				ELSE
					Term = Term * Y(CurrentChem%WhichChem)
				ENDIF

				IF(Associated(CurrentChem%Next)) THEN
					CurrentChem => CurrentChem%Next
					Goto 20
				ENDIF
			ENDIF 
			
			!! Add term to overall derivative
			RunningTotal = RunningTotal + Term

			IF(Associated(Current%Next)) THEN
				Current => Current%Next
				Goto 10
			ENDIF

			pdj(J) = RunningTotal

		END IF
	   END DO
	END IF

! CMB: The IFSOLID branch causes an error.
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	!! If Aqueous Phase Chemistry Occurs !!
!	IF (IFAQ) THEN !!!!!!!!!!!!!!!!!!!!!!!!
!            print *, 'ifaq = true in GasJacobianEvaluator; TODO CB'
!	END IF
!
!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	!! If Solid Phase Chemistry Occurs !!
!	IF (IFSOLID) THEN !!!!!!!!!!!!!!!!!!!
!            print *, 'ifsolid = true in GasJacobianEvaluator; TODO CB'
!	END IF					  

RETURN
END SUBROUTINE GasJacobianEvaluator
