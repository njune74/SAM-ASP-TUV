!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! Water.h
!! This file contains the calculation and retrieval     
!! functions for contains the calculation and retrieval 
!! functions for the Saturation Vapor Pressure of Water in Air

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name            Description				     !!
!! 07     2006   Matt Alvarado   Began Update History			     !!
!! 07/13  2007   Matt Alvarado   Fixed Relative Humidity functions           !!
!! 07/30  2007   Matt Alvarado   Fixed more relative humidity functions      !!
!! 08/29  2007   Matt Alvarado   Fixed mixing ratio and air density functions!!
!! 02/15  2012   Matt Alvarado   Removing Eulerian grids, making ASP         !!
!!                                 a one-box model or subroutine.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. FUNCTION GetSatVapPress () RESULT (SatVapPress)
!! 2. FUNCTION GetSatVapBurden()
!! 3. FUNCTION GetSatVapConcentration ()
!! 4. FUNCTION GetSatMixingRatio () RESULT (SatMixRatio)
!! 5. FUNCTION GetMixingRatio (RHin) RESULT (MixRatio)
!! 6. FUNCTION GetAirDensity () RESULT (Density)
!! 7. FUNCTION GetRelativeHumidity ()
!! 8. FUNCTION GetRelativeHumidityFromGrid ()
!! 9. FUNCTION GetRelativeHumidityFromBurden (WaterBurden)
!!10. SUBROUTINE SetAllRelativeHumidities (RH)
!!11. SUBROUTINE SetRelativeHumidity (RH)
!!12. FUNCTION SurfaceTensionOfWater ()
!!13. FUNCTION DiffusionCoefficientOfWater ()
!!14. FUNCTION GetCpm ()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ****************************************************************
! Modified by CMB (AER): gfortran requires some explicit casting
!                        to double precision (dexp)
! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Saturation Vapor Pressure Values.  At some point,
	!! this may be updated to take interpolated values.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!! This function checks out (07/13/07, MJA)
!#ifdef __GFORTRAN__
        REAL*8 function GetSatVapPress()
!#else
!	REAL*8 FUNCTION GetSatVapPress()  RESULT (SatVapPress)
!#endif
	  USE ModelParameters, ONLY : mbar

	  implicit none
 
	  real*8  :: Temp

          ! CB: Modify to make it look like a function call so ifort doesn't complain 
	  !Temp  = GetTemp 
          Temp = GetTemp()
!#ifdef __GFORTRAN__
          GetSatVapPress = 6.112 * dexp(17.67*(Temp-273.15)/(Temp-29.65)) * mbar
!#else
!	  SatVapPress = 6.112 * dexp(17.67*(Temp-273.15) /(Temp-29.65)) * mbar
!#endif

	RETURN
	END FUNCTION GetSatVapPress

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Saturation Vapor Concentration Values.  At some point,
	!! this may be updated to take interpolated values.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetSatVapBurden()

	  IMPLICIT NONE
          
          ! CB: make look like a func call...
	  !GetSatVapBurden = GetSatVapConcentration 
	  GetSatVapBurden = GetSatVapConcentration() 

	RETURN
	END FUNCTION GetSatVapBurden

	REAL*8 FUNCTION GetSatVapConcentration()  

	  USE ModelParameters, ONLY : eps, ChemScale

	  IMPLICIT NONE

	  real*8 :: SatVapPress

	  SatVapPress			   = GetSatVapPress ()

	  GetSatVapConcentration = ChemScale * SatVapPress/	&
				   (GetPress()) * GetM()

	RETURN
	END FUNCTION GetSatVapConcentration

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Saturation Mass Mixing Ratio of water wv_sat
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        REAL*8 FUNCTION GetSatMixingRatio () RESULT (SatMixRatio)
        REAL*8 FUNCTION GetSatMixingRatio()

		use ModelParameters, ONLY : eps, ppb
		implicit none

		real*8			   :: SatVapPress	! internal

                ! CB: add parentheses at end so that ifort doesn't complain
		SatVapPress = GetSatVapPress()

                ! CB: add parentheses after getpress
!		SatMixRatio = SatVapPress*eps/	& ! ppb; See Emanuel, 1994
!			      (GetPress() - SatVapPress) * 1.0d9 * ppb
		GetSatMixingRatio = SatVapPress*eps/	& ! ppb; See Emanuel, 1994
			      (GetPress() - SatVapPress) * 1.0d9 * ppb
                
	RETURN
	END FUNCTION GetSatMixingRatio

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Mass Mixing Ratio  of Water Vapor, wv
	!! See 2.31 of Jacobson, 2005 (2nd ed.)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	REAL*8 FUNCTION GetMixingRatio (RHin) RESULT (MixRatio)
	REAL*8 FUNCTION GetMixingRatio(RHin)

		use ModelParameters, ONLY : eps, ppb
		implicit none

		real*8, optional   :: RHin
		real*8	           :: SatVapPress, RH

		IF (Present(RHin)) THEN
			RH = RHin
		ELSE
                        ! CB: Add parenthesis so that ifort doesn't complain
			RH = GetRelativeHumidity()
		END IF

                ! CB: add parentheses so that ifort doesnt complain
		SatVapPress = GetSatVapPress()
!		MixRatio = RH*SatVapPress*eps /(GetPress()-RH*SatVapPress) 
		GetMixingRatio = RH*SatVapPress*eps /(GetPress()-RH*SatVapPress) 

                ! Mixing ratio (-)
		!! removed a factor of 1e9 from this, 11/12/2002

	RETURN
	END FUNCTION GetMixingRatio

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Retrieve the Air Density				!!
	!!							!!
	!! Inputs include location and Relative Humidity (RH)	!!
	!! Output is Air Density in (g / lengthscale^3)		!!
	!! Tested against Example 2.3 in Jacobson	        !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	REAL*8 FUNCTION GetAirDensity () RESULT (Density)
	REAL*8 FUNCTION GetAirDensity()

		use ModelParameters
		implicit none

		real*8			   :: MixRatio, RH

		MixRatio = GetMixingRatio ()
!		Density = GetPress () * (1. + MixRatio) /		&
!				  (GetTemp () * Rdry * (1. + MixRatio/eps))
		GetAirDensity = GetPress () * (1. + MixRatio) /		&
				  (GetTemp () * Rdry * (1. + MixRatio/eps))

	RETURN
	END FUNCTION GetAirDensity


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Invert the Concentration of Water to get the                     !!
        !! local relative humidity                                          !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetRelativeHumidity()

		USE ModelParameters, ONLY : eps, ChemScale

		IMPLICIT NONE

		!! Local Variables
		REAL*8 :: WaterConc

		WaterConc = GetGasChem(1)

		GetRelativeHumidity = GetPress() * WaterConc / &
					(GetSatVapPress() * (GetM()))



		!write(*,*) 'CRL: GetRelativeHumidity = ',GetRelativeHumidity 

		!write(*,*) 'CRL: GetPress() = ',GetPress()
		!write(*,*) 'CRL: WaterConc = ',WaterConc 
		!write(*,*) 'CRL: GetSatVapPress()= ',GetSatVapPress()
		!write(*,*) 'CRL: GetM()= ',GetM()
		!write(*,*) ''

		RETURN
	END FUNCTION GetRelativeHumidity

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Same as GetRelativeHumidity, but at a particular grid point instead!!
      !! of location.							    !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetRelativeHumidityFromGrid ()

		USE InfrastructuralCode, ONLY : ERROR
		USE ModelParameters, ONLY : eps, ChemScale

		IMPLICIT NONE 

		!! Local Variables
		REAL*8 :: T, WaterConc

		WaterConc = GetGasChem(1)

		GetRelativeHumidityFromGrid = GetPress() * WaterConc / &
				(GetSatVapPress() * GetM())


		RETURN
	END FUNCTION GetRelativeHumidityFromGrid

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! There are times in the condensation routine that we are holding a   !!
      !! water burden for a particular gridcell in a vector that is running  !!
      !! through LSODES.  We want to be able to get a relative humidity value!!
      !! from this without actually altering the Eulerian Array. So this     !!
      !! takes that burden, diagnoses the other environmental parameters, and!!
      !! then delivers a relative humidity value.			     !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetRelativeHumidityFromBurden (WaterBurden)

		USE ModelParameters, ONLY : eps, ChemScale

		!! Local Variables
		REAL*8 :: T, WaterBurden
		
		GetRelativeHumidityFromBurden =	GetPress() * WaterBurden  / &
				(GetSatVapPress() * (GetM()))

		RETURN
	END FUNCTION GetRelativeHumidityFromBurden

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Reset all of the relative humidities at once !!
!! primarily for debugging but also for simple  !!
!! experiments.									!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetAllRelativeHumidities (RH)

	IMPLICIT NONE

	REAL*8  :: RH

        CALL SetRelativeHumidity (RH)


END SUBROUTINE SetAllRelativeHumidities 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Reset the relative humidity at a given  !!
!! gridcell location to a particular value !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SetRelativeHumidity (RH)

	USE ModelParameters, ONLY : eps, ChemScale, ppbscale

	IMPLICIT NONE

	REAL*8  :: RH

	!! Local Variables
	REAL*8 :: T, Press, WaterMR, SPress

	SPress = GetSatVapPress()

	!! Maybe put a MM on the top?
	GridGasChem(1) = ChemScale * SPress * RH * GetM() / GetPress()

	RETURN

END SUBROUTINE SetRelativeHumidity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This approximation for SURFACE TENSION of Water (dyn / cm) is    !!
!! from Jacobson (15.2)  from Prupp & Klett fit of a polynomian	    !!
!! expression to Dorsch and Hacker (1951) see books for ref.	    !!
!! It is only really good for [-40,40] but use it everywhere here.  !!
!! Further dependence of surface tension on electrolytic properties !!
!! is treated by SurfaceTension() in ParticleAttributes.h           !! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION SurfaceTensionOfWater ()

	  USE InfrastructuralCode, ONLY : Warn, REAL2STR
	  USE ModelParameters,     ONLY : dyncm

	  IMPLICIT NONE

	  !! Local Variables
	  REAL*8 :: Tc

	  Tc = GetTemp() - 273.15

	  !! Supercooled or not Supercooled
	  IF (Tc < 0) THEN
            SurfaceTensionOfWater = 75.93+0.115*Tc+0.06818*(Tc**2) &
                                    +6.511d-3*(Tc**3)+		   &
			            2.933d-4*(Tc**4)+6.283d-6*(Tc**5) &
                                    +5.285d-8*(Tc**6)
	  ELSE
	     SurfaceTensionOfWater = 76.1 - 0.155*Tc
	  END IF

	  SurfaceTensionOfWater = SurfaceTensionOfWater * dyncm

	  RETURN
	END FUNCTION SurfaceTensionOfWater

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIFFUSIVITY OF WATER VAPOR as a function of temperature and pressure. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION DiffusionCoefficientOfWater ()

		USE ModelParameters, ONLY : Atmospheres, cm

		IMPLICIT NONE



		!! Calculate the Diffusion Coefficient of Water (cm^2/s)
		!! See Seinfeld and Pandis p.801 or Prupp & Klett
		DiffusionCoefficientOfWater = (0.211*Atmospheres/GetPress())*(GetTemp()/273.)**1.94 * cm * cm

		RETURN
	END FUNCTION DiffusionCoefficientOfWater

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! SPECIFIC HEAT OF MOIST AIR AT CONSTANT PRESSURE !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION GetCpm ()

		USE Modelparameters, ONLY : Cpd, CpV

		IMPLICIT NONE

		!! Internal Variables
		REAL*8 :: MR

		MR = GetMixingRatio()

		GetCpm = (Cpd + MR * CpV) / (1. + MR)

		RETURN
	END FUNCTION GetCpm
