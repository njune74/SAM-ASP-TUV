!! ASP (c), 2004-2012, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Condensation.f90
!! This is the main source file for all gas to aerosol transfers,
!! including the equilibrium and flux-limited kinetic routines.
!! It mainly links all the header files together.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY								!!
!!										!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP           !!
!!                                 a one-box model or subroutine.               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES		                    !!
!! 1. ModelParameters			    !!
!! 2. InfrastructuralCode		    !!
!! 3. GridPointFields			    !!
!! 4. Chemistry				    !!
!! 5. Aerosols				    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES				        !!
!! 1. CondensationInitialization.h		        !!
!! 2. CondensationIntegrator.h			        !!
!! 3. CondensationRelatedFunctions.h			!!
!! 4. HydrophobicCondensationFunctions.h		!!
!! 5. binsolu.h (Organic water content data)		!!
!! 6. OrgCondensationIntegrator.h			!!
!! 7. AqOrgCondensationIntegrator.h			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!	
!! 1. NONE								!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Condensation

USE Aerosols, ONLY : Particle

PRIVATE
PUBLIC :: StepCondensation,		   &
          StepCondensationAll,             &
          SetAllDissolution,	           &
          EquilibrateGridPoint,            &
          EqAllWater,                      &
          EquilibrateOrganicParticle,      &
          EquilibrateInternallyatGridPoint,&
          OrganicDissolutionData

        !!These should be replaced by Get functions, so
        !!outside routines can't change them MJA, 20130719
	INTEGER :: HowManyDissolutionReactions, &
                   HowManyOrganicDissolutionReactions, &
	           HowManyAqOrganicDissolutionReactions

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The primary array is DissolutionData(i,j), defined in 
	!! CondensationInitialization.h, in which i is the
	!! reaction and j is some data:
	!!
	!! DissolutionData(i,1) : Index in Gas Phase of Primay Chemical
	!! DissolutionData(i,2) : Index in Gas Phase of Secondary Chemical
	!! DissolutionData(i,3) : Index in Aqueous Phase
	!! DissolutionData(i,4) : Mass Accommodation Coefficient
	!! DissolutionData(i,5) : Henrys Law Coefficient (H_298)
	!! DissolutionData(i,6) : Henrys Law Coefficient (Delta H / R)
	!! DissolutionData(i,7) : Index of AqEquilibriaList for related Aq reaction
	!! These are loaded in from the input files Dissolution.in, 

        !!These should be replaced by Get functions, so
        !!outside routines can't change them MJA, 20130719
	REAL*8,	ALLOCATABLE :: DissolutionData(:,:)
	REAL*8,	ALLOCATABLE :: OrganicDissolutionData(:,:)
	REAL*8,	ALLOCATABLE :: AqOrganicDissolutionData(:,:)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This loads the input decks and the appropriate arrays !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE "CondensationInitialization.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The routines that step Condensation forward !!
!! and interface with LSODES                   !!
INCLUDE "CondensationIntegrator.h" !!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Random other routines, including saturation !!
!! pressure calculations.                      !!
INCLUDE "CondensationRelatedFunctions.h" !!!!!!!!

INCLUDE "HydrophobicCondensationFunctions.h"

INCLUDE "binsolu.h"

INCLUDE "OrgCondensationIntegrator.h"

INCLUDE "AqOrgCondensationIntegrator.h"

END MODULE Condensation

