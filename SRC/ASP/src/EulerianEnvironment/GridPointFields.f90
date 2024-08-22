!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! GridPointFields.f90
!! This Module Provides chemical, thermodynamic, dynamical and 
!! other particle-independent quantities to subroutines.  Each 
!! are updated after each time step.
!! This is the source file linking many header files together.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY					                    !!
!!									    !!
!! Month  Year   Name              Description				    !!
!! 07     2006   Matt Alvarado     Began Update History			    !!
!! 02/15  2012   Matt Alvarado     Removing most of these grids, making ASP !!
!!                                 a one-box model or subroutine.           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEPENDENCIES                             !!
!! 1. ModelParameters			    !!
!! 2. InfrastructuralCode		    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! LINKED HEADER FILES						    !!
!! 1. SetGrid.h							    !!
!! 2. Initialize.h						    !!
!! 3. Temperature.h						    !!
!! 4. Pressure.h						    !!
!! 5. Water.h							    !!
!! 6. GasPhaseChemistry.h				            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:      !!
!! 1. SUBROUTINE SetHowManyGasChems(inHowManyGasChems)              !!
!! 2. SUBROUTINE SetYLEN(inYLEN)                                    !!
!! 3. SUBROUTINE SetHowManyEvolveGasChems(inHowManyEvolveGasChems)  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	module GridPointFields

	implicit none

	!! Everything is private unless declared otherwise
	PRIVATE

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Now specify which subroutines may be accessed by the outside     !!
	!! most will be retrieval (GET---)                                  !!
        !! or initialization (SET--- or INIT---)                            !!
	!! functions.							    !! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	PUBLIC	::	AllocateGasChemistryGrid,		& 
			AllocateGasReactionGrid,		&
			DiffusionCoefficientOfWater,	        &
			GasPhaseChemicalRates,			&
			GetAirDensity,				&
			GetCpm,					&
			GetDynamicViscosityOfAir,		& 
			GetGasChem,				&
			GetGasBurden,				&
			GetGasChems,				&
			GetGasChemVector,			&
			GetGasChemFromBurden,			&
			GetM,					&
			GetMixingRatio,				&
			GetPartialPressure,			&
			GetPress, GetPressFromGrid,		& 
			GetRelativeHumidity,			&
			GetRelativeHumidityFromGrid,	        &
			GetRelativeHumidityFromBurden,	        &
			GetSatVapPress,				& 
			GetSatMixingRatio,			&
			GetSatVapConcentration,			&
			GetSatVapBurden,			&
			GetTemp, GetTempFromGrid,		& 
			GetThermalVelocityOfAir,		&  
			GridPointVolume,			&
			InitializeGridPointFields,		& 
			MeanFreePathOfAir,			&
			MolecularThermalConductivity,	        &
			MolecularThermalDiffusivity,	        &
			ResetGasPhaseReactionRates,		&
			SetGridPointVolumes,			&
			SetHomogeneousChemFieldPPB,		&
			SetHomogeneousChemFieldMGperM3,	        &
			SetHowManyGasChems,			&
			SetWaterVaporField,			&
			SetPressField,				&
			SetTempField,				&
			SetYLEN,SetHowManyEvolveGasChems,	&
			SurfaceTensionOfWater,			&
			SetRelativeHumidity,			&
			SetAllRelativeHumidities,		&
			UpdateChemicalBurdens,			&
			UpdateChemicalConcentrations,	        &
			GetGridCellChemBurden,			&
			AddToGridCellChemBurden,		&
			ReplaceGridCellChemBurden,		&
			AllSameGrid,				&
                        !gridgaschem    ! cb out
                        gridgaschem, linear2dInterp ! cb in

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Each of the EULERIAN GRIDS have their own GRID SPACING !!
	!! We define the several sets of those here.	          !!
	!!							  !!
	!!   1 is TEMPERATURE					  !!
	!!   2 is PRESSURE					  !!
	!!   3 is WATER SATURATION VAPOR PRESSURE		  !!
	!!   4 is GAS PHASE CHEMISTRY and Particles		  !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	integer, parameter :: HowManyGrids	         = 3
!        integer, parameter :: TempGridIndex		 = 1
!	integer, parameter :: PressGridIndex		 = 2
!	integer, parameter :: GasPhaseChemGridIndex      = 3  
        !! Also the grid for particles

	!! These arrays hold the number of grid points in each cartesian 
	!! coordinate for each grid mesh.
!	integer :: XGridPoints(HowManyGrids)
!	integer :: YGridPoints(HowManyGrids)
!	integer :: ZGridPoints(HowManyGrids)
	real*8  :: GridPointVolume

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The folling arrays are PRIMARY !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8              :: GridTemp 	! Air Temperature (K)
	real*8              :: GridPress	! Air Pressure (mb)
	real*8, allocatable :: GridGasChem(:)   
        ! Gas Phase Chemistry (mol per cc for each molecule)
	real*8              :: TemporaryWaterArray
        ! Relative Humidity (between 0 and 1) is stored temporarily

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! The folling arrays hold SECONDARY THERMODYNAMIC variables, 
	!! and are diagnosed from the primary variables above.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real*8, allocatable :: GasPhaseChemicalRates(:)

	!! If all the grids are the same scale, then this is true
	LOGICAL :: AllSameGrid
	LOGICAL :: AllSamePressure !! Same if common pressure

	!! replicated from chemistry module
	INTEGER :: HowManyGasChems, ylen, HowManyEvolveGasChems

	contains

	SUBROUTINE SetHowManyGasChems(inHowManyGasChems)
		INTEGER :: inHowManyGasChems
		HowManyGasChems = inHowManyGasChems
	END SUBROUTINE

	SUBROUTINE SetYLEN(inYLEN)
		INTEGER :: inYLEN
		YLEN = inYLEN
	END SUBROUTINE

	SUBROUTINE SetHowManyEvolveGasChems(inHowManyEvolveGasChems)
		INTEGER :: inHowManyEvolveGasChems
		HowManyEvolveGasChems = inHowManyEvolveGasChems
	END SUBROUTINE

        ! CB: Need to access this header file...
        !     Nope, let's just include that source file here...
        !include "InterpolateGrid.h"
        ! CB: Add simple linear interpolation routine.  Yeah, it 
        !     uses floating-point equality...I know that's bad, but
        !     it will suit us for the purposes of this example.
        !
        ! UPDATE: The sqrt(-1.0) does not work for a NaN set under
        ! gfortran.  Let's try tricking the compiler...
        subroutine linear2dInterp(x1, y1, x2, y2, x3, y3)
        implicit none

        real*8, intent(in) :: x1, y1, x2, y2, x3
        real*8, intent(out) :: y3
        real*8 :: slope, intercept
        
        ! CMB: Hack here
        real*4 :: nan
        nan = 0.

        if (x1 .eq. x2) then
            if (y1 .ne. y2) then
                y3 = nan/nan !sqrt(-1.0)     ! force NaN
            else
                y3 = y2
            endif
        else
            slope = (y2 - y1)/(x2 - x1)
            intercept = y2 - slope*x2
            y3 = slope*x3 + intercept
        endif

        end subroutine linear2dInterp
        ! end cb add

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! This series of includes specify target files that contain       !!
        !! grid point field subrountines to be included in this module     !!
        !! (and thus not compiled on their own but rather drawn into here. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	include "SetGrid.h"
	include "Initialize.h"
	include "Temperature.h"
	include "Pressure.h"
	include "Water.h"
	include "GasPhaseChemistry.h"

        end module GridPointFields
