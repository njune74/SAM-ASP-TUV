!! ASP (c), 2004-2013, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! binsolu.h
!! This file calculates the amount of water associated with the aqueous organic
!! species. The data for each compound comes from a UNIFAC calculation of the 
!! water activity in a binary prganic-water solution. Values are linearly 
!! interpolated.
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY								!!
!!										!!
!! Month  Year   Name              Description					!!
!! 07     2006   Matt Alvarado     Began Update History				!!
!! 01/30  2013   Matt Alvarado     Changed HydrophilicOrganicWaterContent to    !!
!!                                   use Kappa                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. FUNCTION HydrophilicOrganicWaterContent(InParticle, Temp, RH)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!This function calculates the water content (in mol/particle) caused by the
!hydrophilic oprganic using Kappa param
REAL*8 FUNCTION HydrophilicOrganicWaterContent(InParticle, Temp, RH)
      
        USE Chemistry,  ONLY : Kappa, AqOrgDensity, AqOrgMolecularMass,&
                               AqMolecularMass, HowManyAqOrgChems 
	
        IMPLICIT NONE

        !External Variables
        REAL*8 :: Temp, RH
        TYPE(Particle),POINTER :: InParticle

        !Internal Variables
        INTEGER :: I
        REAL*8 :: Sum 
        REAL*8, ALLOCATABLE :: MOLALRH(:), ORGCONC(:)

        ALLOCATE(MOLALRH(1:HowManyAqOrgChems),ORGCONC(1:HowManyAqOrgChems))

	Sum = 0.0
	DO I = 1, HowManyAqOrgChems
                IF (Kappa(I) .GT. 1.0e-5) THEN !For very small kappa, no water is associated with the organic
                   !Assume water density = 1.0 g/cm3
		   MOLALRH(I) = (AqOrgDensity(I)/1.0)*((1.0-RH)/RH)/Kappa(I)/AqOrgMolecularMass(I) 
                   !umol/ug = mol org/g H2O
	
 		   ORGCONC(I) = InParticle%AqOrgChems(I)*InParticle%NumberofParticles !Mol org/cm3

                   !ZSR formula
		   Sum = Sum + (ORGCONC(I)/MOLALRH(I)) !g H2O /cm3 air
		ENDIF	
	END DO

	
	Sum = Sum/AqMolecularMass(1) !Water content in mol/cm3 air
	!WRITE(*,*) "Water mass: ", AqMolecularMass(1)
	HydrophilicOrganicWaterContent = Sum/InParticle%NumberofParticles !water content in mol/particle
	!WRITE(*,*) "OrgWater: ", HydrophilicOrganicWaterContent*InParticle%NumberofParticles, "RH: ", RH

END FUNCTION HydrophilicOrganicWaterContent



