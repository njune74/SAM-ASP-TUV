!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Temperature.h
!! This file contains the calculation and retrieval     
!! functions for various thermodynamics fields of the gridded
!! environment.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							    !!
!!									    !!
!! Month  Year   Name              Description				    !!
!! 07     2006   Matt Alvarado     Began Update History			    !!
!! 02/15  2012   Matt Alvarado     Removing Eulerian grids, making ASP      !!
!!                                 a one-box model or subroutine.           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. SUBROUTINE SetTempField (Temp)
!! 2. FUNCTION GetThermalVelocityOfAir ()
!! 3. FUNCTION MolecularThermalDiffusivity () 
!! 4. FUNCTION MolecularThermalConductivity ()
!! 5. FUNCTION GetDynamicViscosityOfAir () RESULT (DynVisc)
!! 6. FUNCTION GetTemp () RESULT (Temp)
!! 7. FUNCTION GetTempFromGrid () RESULT (Temp)
!! 8. FUNCTION MeanFreePathOfAir ()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set the initial TEMPERATURE field.  This will be ammended or	  !!
!! forced eventually to reflect an actual environment.	          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetTempField (Temp)

	implicit none

	!! Input Variables
	real*8  :: Temp

	GridTemp = Temp

	RETURN
END SUBROUTINE SetTempField

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The THERMAL VELOCITY of AIR is simply proportional to the sq root	!!
!! of the temperature.  Derivation is from gas kinetic theory, setting	!!
!! average kinetic energy of air molec to 4/pi kB T.			!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetThermalVelocityOfAir ()

	use ModelParameters, ONLY : ThmVelFac, cm, AirMolecMass

	implicit none

	GetThermalVelocityOfAir = ThmVelFac * SQRT(GetTemp() / AirMolecMass )

	RETURN
END FUNCTION GetThermalVelocityOfAir 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The MOLECULAR THERMAL DIFFUSIVITY, calculated from List (1984)'s !!
!! Smithsonian Meteorological Tables fit.			    !!
!!								    !!
!! Dh = Kappa_d / (rho cpd)				            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION MolecularThermalDiffusivity () 

	USE ModelParameters, ONLY : m, joules, cpd

	IMPLICIT NONE


	!! The molecular thermal diffusivity (J cm-1 s-1 K-1)
	MolecularThermalDiffusivity = MolecularThermalConductivity () &
                                      / GetAirDensity() / cpd

END FUNCTION MolecularThermalDiffusivity 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The MOLECULAR THERMAL CONDUCTIVITY, calculated from List (1984)'s!!
!! Smithsonian Meteorological Tables fit.			    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION MolecularThermalConductivity () 

  USE ModelParameters, ONLY : joules, cm

  IMPLICIT NONE

  !! List 1984's empirical fit to Kd (J cm-1 s-1 K-1):
  MolecularThermalConductivity = &
       ( 2.3807d-4 + 7.1128d-7 * (GetTemp()-273.16) )* joules / cm

END FUNCTION MolecularThermalConductivity 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The DYNAMIC VISCOSITY of AIR is calculated using Sutherland's Equation,!! 
!! taken from List, 1984 (Smithsonian Met Tables) via Jacobson's Modeling !!
!! book.  The quantity itself is the ratio of shearing stress to shear.	  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetDynamicViscosityOfAir () RESULT (DynVisc)

	use ModelParameters
	implicit none

	real*8			   :: T			 ! internal

	T = GetTemp ()

	DynVisc = 1.4962864d-5 * T**1.5 / (T + 120.) * grams / cm

END FUNCTION GetDynamicViscosityOfAir 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Retrieve the TEMPERATURE			                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION GetTemp () RESULT (Temp)

	implicit none

	Temp =  GridTemp

        return
END FUNCTION GetTemp
 
!! Same as GetTemp 
 REAL*8 FUNCTION GetTempFromGrid () RESULT (Temp)
		
	USE InfrastructuralCode, ONLY : ERROR

	IMPLICIT NONE

	Temp = GridTemp

	RETURN
END FUNCTION GetTempFromGrid



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Find the Mean Free Path at a Particular Location !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION MeanFreePathOfAir ()

        USE ModelParameters, ONLY : grams, cm, lengthscale, massscale

	IMPLICIT NONE

	!! Internal
	REAL*8 :: T 

	T = GetTemp()

	!! The mean free path is equal to 2 * the dynamic viscosity of air
	!! divided by (density of air * the thermal velocity of air)
	MeanFreePathOfAir = 2*GetDynamicViscosityOfAir()/ &
					(GetThermalVelocityOfAir() *	&
					GetAirDensity()*lengthscale)

	!write(*,*) 'CRL: In Temperature.h: MeanFreePathOfAir = ',MeanFreePathOfAir
	!write(*,*) 'CRL: GetDynamicViscosityOfAir() = ', GetDynamicViscosityOfAir()
	!write(*,*) 'CRL: GetThermalVelocityOfAir() = ',GetThermalVelocityOfAir() 
	!write(*,*) 'CRL: GetAirDensity() = ',GetAirDensity()
	!write(*,*) 'CRL: lengthscale = ',lengthscale 
	!write(*,*) ''
	!! This result does include an implicit factor of lengthscale!
	RETURN
END FUNCTION MeanFreePathOfAir

