!! ASP (c), 2004-2006, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ChemistryDerivativeStructures.h
!! This defines the elements of the ODE and Jacobian linked lists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY																!!
!!																				!!
!! Month  Year   Name              Description									!!
!! 07     2006   Matt Alvarado     Began Update History							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file defines the following data types:	!!
!! 1. TYPE ChemTerm								!!
!! 2. TYPE DiffEqTerm							!!
!! 3. TYPE DiffEqList							!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Each sub-term of an ODE or of the Jacobian is an element  !!
	!! of a linked list.  The external subroutines provided to   !!
	!! LSODE or its cousins will be calculated from these lists. !!
	!!															 !!
	!! Each term of a chemical ODE (or Jacobian) element has the !!
	!! following form:											 !!
	!!															 !!
	!!		(+/-) Rate * Factor * PI_n (Y(i)^n)					 !!
	!!															 !!
	!! Which is represented in this linked list by one element of!!
	!! DiffEqTerm, with:										 !!
	!!		 +1 or -1 in PlusMinus representing (+/-)			 !!!!!!!!!!!!!!!!!!!!!!!
	!!		 Rate holding the rate constant to be used in appropriate units			  !!
	!!		 Factor holding the leading factor (integer) for the term				  !!
	!!		 A linked list of ChemTerms representing each multiple in the PI notation !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	TYPE ChemTerm
		integer					    :: WhichChem
		real*8						:: Stoich
		type (ChemTerm), pointer    :: Next
	END TYPE ChemTerm

	TYPE DiffEqTerm
		real*8						:: LeadingFactor ! Contains + or -, as appropriate, and stoichiometric coeff
		integer						:: ReactionNumber
		type (ChemTerm), pointer	:: DerefChem
		type (DiffEqTerm), pointer  :: Next
	END TYPE DiffEqTerm

	TYPE DiffEqList
		type (DiffEqTerm), pointer	:: FirstTerm
	END TYPE DiffEqList
