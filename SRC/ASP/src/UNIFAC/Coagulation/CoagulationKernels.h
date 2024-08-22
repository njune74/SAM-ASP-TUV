!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! CoagulationKernels.h
!! This file contains the functions for calculating Brownian diffusion
!! and gravitational collection coagulation kernels.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 08/29  2007   Matt Alvarado   Fixed mean free path calculation in         !!
!!				    CoagKernBrownianDiff                     !!
!!				 Fixed calculation of Schmidt number in      !!
!!				    CoagKernBrownianDiff                     !!
!! 02/17  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:	!!
!! 1. FUNCTION CoagKernBrownianDiff (Particle1, Particle2)	!!
!! 2. FUNCTION CoagKernGravCollection (Particle1, Particle2)	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Brownian Diffusion Coagulation Kernel for a		    !!
!! two given sizes of particles.  Ther kernel has the form:		    !!
!!									    !!
!!		KBD(r1,r2) = 4*pi (ri + rj) (Di + Dj) Beta		    !!
!!									    !!
!! In which ri is the radius of i, Di is the diffusion coefficient,	    !!
!! and Beta is a transition-region correction factor based on the Knudsen   !!
!! number.  The inputs rib and rit signify the bottom and top of a	    !!
!! bin's size range.  The bins have flat distributions.  The Mass values    !!
!! are the full masses associated with a size bin (we track both mass and   !!
!! number).								    !!
!!									    !!
!! We also calculate the Coagulation Kernel of Brownian Convective          !!
!! Enhancement, which is a separate kernel addressing the pull of larger    !!
!! particles' wakes on smaller particles.				    !!
!!									    !!
!! The routine returns the SUM of the two kernels, appropriate given the    !!
!! way these will be solved.						    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CoagKernBrownianDiff (Particle1, Particle2)

	USE GridPointFields, ONLY : GetDynamicViscosityOfAir,		&
				    GetTemp,				&
				    GetAirDensity

	USE Chemistry, ONLY        : AqMolecularMass

	USE Aerosols,		 ONLY : Particle,			&
					ParticleKnudsen,	&
					ParticleMass,		&
					ReynoldsNumber	

	USE ModelParameters, ONLY : Pi, kB, ThmVelFac

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle), POINTER :: Particle1, Particle2

	!! Local Variables
	REAL*8 :: DynVisc, Temp, AirDensity	!! Environmental
	REAL*8 :: D1, D2, Knudsen1, Knudsen2, Cs1, Cs2,	&	!! Particle
			  Mass1, Mass2, MFP1, MFP2, Del1, Del2, CE, BetaInv, &
			  ReNumb, ScNumb

!print *, "Rad1: ",Particle1%EffectiveRadius,Particle2%EffectiveRadius

	!! CE is a stand-in for collection efficiency, set here to 1
	CE = 1.

	!! Use the location of one of those as a stand-in.
	!! First Calculate the Di's.  These depend on kB, T, the Cunningham Correction
	!! Factor, the Dynamic Viscosity, and the radius.
	DynVisc    = GetDynamicViscosityOfAir ()
	Temp       = GetTemp ()
	AirDensity = GetAirDensity ()

	Knudsen1   = ParticleKnudsen (Particle1)
	Knudsen2   = ParticleKnudsen (Particle2)

!print *, "Knud 1 2 ",Knudsen1,Knudsen2

	!! The Cunniongham Slip Correction Factor (Cs) has three possible forms.
	!! The first comes from Miliken via Jacobson's Atm. Modeling Book:
	IF (.FALSE.) THEN
	!! From Kasten 1968 (Jacobson)
		Cs1 = 1. + Knudsen1 * (1.249 + 0.42 * DEXP(-0.87 / Knudsen1))
		Cs2 = 1. + Knudsen2 * (1.249 + 0.42 * DEXP(-0.87 / Knudsen2))
	ELSEIF (.FALSE.) THEN
	!! From Allen and Raabe 1982
		Cs1 = 1. + Knudsen1 * (1.257 + 0.4 * DEXP(-1.1 / Knudsen1))
		Cs2 = 1. + Knudsen2 * (1.257 + 0.4 * DEXP(-1.1 / Knudsen2))
	ELSEIF (.TRUE.) THEN
	!! From Fuchs 1964 (Seinfeld and Pandis, p. 661)
		Cs1 = (5. + 4.*Knudsen1 + 6.*Knudsen1*Knudsen1 + 18.*Knudsen1*Knudsen1*Knudsen1)/&
			  (5. - Knudsen1    + (8.+pi)*Knudsen1*Knudsen1)
		Cs2 = (5. + 4.*Knudsen2 + 6.*Knudsen2*Knudsen2 + 18.*Knudsen2*Knudsen2*Knudsen2)/&
			  (5. - Knudsen2    + (8.+pi)*Knudsen2*Knudsen2)
	END IF

	!! Which provides enough to calculate Di (cm^2 / s)
	D1 = kB*Temp*Cs1 / (6.*pi*DynVisc*Particle1%EffectiveRadius)
	D2 = kB*Temp*Cs2 / (6.*pi*DynVisc*Particle2%EffectiveRadius)

	!! Get particle mass
	Mass1 = ParticleMass(Particle1)
	Mass2 = ParticleMass(Particle2)
!print *, "Mass 1 2 ",Mass1, Mass2 

	!! Calculate the mean free paths
	MFP1 = 8.*D1*Sqrt(Mass1/Temp) / Pi / ThmVelFac
	MFP2 = 8.*D2*Sqrt(Mass2/Temp) / Pi / ThmVelFac

	!! We next define the inverse of beta, the transition factor
	!! This comes from Fuchs' 1964 book.  This depends on delta, which
	!! we calculate first:
	Del1 = ((2.*Particle1%EffectiveRadius + MFP1)**3. -										&
			(4.*Particle1%EffectiveRadius * Particle1%EffectiveRadius + MFP1*MFP1)**1.5)	&
		   /(6.*Particle1%EffectiveRadius * MFP1) - 2.*Particle1%EffectiveRadius

	Del2 = ((2.*Particle2%EffectiveRadius + MFP2)**3. -										&
			(4.*Particle2%EffectiveRadius * Particle2%EffectiveRadius + MFP2*MFP2)**1.5)	&
		   /(6.*Particle2%EffectiveRadius * MFP2) - 2.*Particle2%EffectiveRadius

	!! The Fuchs interpolation form for the transition factor (here, the inverse):
	BetaInv = (Particle1%EffectiveRadius + Particle2%EffectiveRadius) / (Particle1%EffectiveRadius +						 &
			  Particle2%EffectiveRadius + Sqrt(Del1*Del1+Del2*Del2)) +														 &
			  4.*CE*(D1+D2)/(Particle1%EffectiveRadius + Particle2%EffectiveRadius)/ThmVelFac/Sqrt(Temp/Mass1 + Temp/Mass2)

	!! Coagulation Kernel
	CoagKernBrownianDiff = 4. * Pi * (Particle1%EffectiveRadius + Particle2%EffectiveRadius) * (D1 + D2) / BetaInv


!print *, "CoagKernBrownianDiff ",CoagKernBrownianDiff 

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Now add the terms to do with the Convective Enhancement  !!
	!! of Brownian Diffusion Kernel.  This is the effect of the !!
	!! larger particle's  pulling in of smaller ones with its   !!
	!! turbulent wake.  cf. Jacobson p. 446 and erratta.        !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Get the Reynolds Number of the Larger Particle and the
	!! Schmidt Number of the smaller.
	IF (Particle1%EffectiveRadius .GE. Particle2%EffectiveRadius) THEN
		ReNumb = ReynoldsNumber (Particle1)
		ScNumb = DynVisc / AirDensity / D2
	ELSE
		ReNumb = ReynoldsNumber (Particle2)
		ScNumb = DynVisc / AirDensity / D1
	END IF

	!! Make the correction
	IF (ReNumb .LE. 1) THEN
		CoagKernBrownianDiff = CoagKernBrownianDiff * (1. + 0.45 * (ReNumb * ScNumb)**(1./3.))
	ELSE
		CoagKernBrownianDiff = CoagKernBrownianDiff * (1. + 0.45 * ReNumb**0.5 * ScNumb**(1./3.))
	END IF

	RETURN
END FUNCTION CoagKernBrownianDiff


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Gravitational Collection Coagulation Kernel for       !!
!! two given sizes of particles.  Ther kernel has the form:	       !!
!!								       !!
!!		KGC(r1,r2) = Ec pi (r1 + r2)^2 (Vf_t,1 - Vf_t,2)       !!
!!								       !!
!! The difficulty with this form is properly writing the collection    !!
!! efficiency, Ec.   						       !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CoagKernGravCollection (Particle1, Particle2)

	USE GridPointFields, ONLY : GetDynamicViscosityOfAir,	&
					GetAirDensity

	USE Chemistry,       ONLY : AqMolecularMass

	USE Aerosols,	     ONLY : Particle,		&
				    TerminalVelocity

	USE ModelParameters, ONLY : Pi

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle), POINTER :: Particle1, Particle2

	!! Local Variables
	REAL*8 :: dVel, Ec
	REAL*8 :: DynVisc, AirDensity  !! Environmental   


	!! Both particles will be associated with the same gridpoint, so simply use 
	!! the location of one of those as a stand-in.
	!! First Calculate the Di's.  These depend on kB, T, the Cunningham Correction
	!! Factor, the Dynamic Viscosity, and the radius.
	DynVisc    = GetDynamicViscosityOfAir ()
	AirDensity = GetAirDensity ()

	!!From Jacobson, "Fundamentals of Atmospheric Modeling", p. 510
	Ec = CollisionEfficiency(Particle1, Particle2)

	!! Both sectional means the velocity differential is of the fall velocities
	IF (Particle1%Sectional .AND. Particle2%Sectional) THEN

		dVel = DABS(TerminalVelocity(Particle1)-TerminalVelocity(Particle2))

	!! One or the other is sectional
	ELSE IF (Particle1%Sectional) THEN


		dVel = SQRT((Particle2%Velocity(3)-TerminalVelocity(Particle1))**2 +	&
					Particle2%Velocity(1)*Particle2%Velocity(1) +				&
					Particle2%Velocity(2)*Particle2%Velocity(2))
	
	ELSE IF (Particle2%Sectional) THEN

		dVel = SQRT((Particle1%Velocity(3)-TerminalVelocity(Particle2))**2 +	&
					Particle1%Velocity(1)*Particle1%Velocity(1) +				&
					Particle1%Velocity(2)*Particle1%Velocity(2))

	!! Neither are sectional
	ELSE 

		dVel = SQRT((Particle1%Velocity(1)-Particle2%Velocity(1))**2 +	&
					(Particle1%Velocity(2)-Particle2%Velocity(2))**2 +	&
					(Particle1%Velocity(3)-Particle2%Velocity(3))**2)

	END IF

	!! The kernel is then defined as a simple product
	CoagKernGravCollection = Ec * CE * Pi * dVel *		&
							 (Particle1%EffectiveRadius+Particle2%EffectiveRadius)**2
	
	RETURN
END FUNCTION CoagKernGravCollection

