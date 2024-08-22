!! ASP (c), 2004-2013, Matt Alvarado (malvarad@aer.com)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! ParticleAttributes.h
!! Contains functions for calculating particle properties 
!! (density, radius, mass
!! curvature correction, Knudsen number, Reynolds number, etc.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 07/24  2006   Matt Alvarado     Removed "ISNAN" to fit pgf90		     !!
!! 09/07  2006   Matt Alvarado   Changed RecalculateRadius to give lower edge!!
!!				      when number conc. = 0.		     !!
!! 11/01  2006   Matt Alvarado     Created ShellRefIndAndRad
!! 05/18  2007   Matt Alvarado     Revised ShellRefIndAndRad
!!				   Created AerosolOptProp
!! 07/13  2007   Matt Alvarado     Added line to  ShellRefIndAndRad
!!	 			      to keep core rad below shell rad
!! 07/17  2007   Matt Alvarado     Added line to SurfaceTension
!!				      to keep value above that of pure Ethanol
!! 07/18  2007   Matt Alvarado     Removed calls to ParticleDensity from
!!				CurvatureCorrection and OrgCurvatureCorrection
!!			 Just use values in particle data structure instead
!!				   In ParticleDensity, if ion loop fails, force
!!				      density of ionic solution to 1.0 g/cm3
!! 07/26  2007   Matt Alvarado     Fixed AbsCoreRad calculation in 
!!                                 ShellRefIndAndRad
!! 09/07  2007   Matt Alvarado     Fixed shell refractive index calcultion
!! 09/17  2007   Matt Alvarado     Removed write statement from ParticleDensity
!! 10/01  2007   Matt Alvarado     Set ParticleDensity to skip bins 
!!                                 with 0 particles
!! 10/03  2007   Matt Alvarado     Set Denomsum limit in SurfaceTension to 
!!                                 1.0e40
!! 10/04  2007   Matt Alvarado     Protected ShellRefIndAndRad from near zero 
!!                                  number concentrations
!! 10/11  2007   Matt Alvarado     Set HypotheticalElectrolyteConcentrations 
!!                                 to warn of charge
!!                                 imbalance, but not stop.
!! 10/15  2007   Matt Alvarado     Fixed Assymetry parameter in 
!!                                 ShellRefIndAndRad
!! 08/30  2010   Matt Alvarado     Changed ShellRefIndandRad, 
!!                                 AerosolOpticalDepth, and AerosolOptProp 
!!                                 to allow 62 photolysis wavelengths
!! 02/16  2012   Matt Alvarado     Removed Eulerian grids, making ASP        !!
!!                                 a one-box model or subroutine.            !!
!! 02/27  2012   Matt Alvarado     Changed ShellRefIndandRad, 
!!                                  AerosolOpticalDepth, and AerosolOptProp 
!!                                  to allow 451 bins between 250 nm and 700 nm
!!                                 ShellRefIndandRad now used OPAC wavelength
!!                                  dependent complex refractive indices
!! 02/29  2012   Matt Alvarado     Changed ShellRefIndandRad to use different
!!                                  wavelength dependent refractive indices
!!                                  for different electrolytes based on flags
!!                                  in AqPhaseChems.in. Also added refractive
!!                                  index for dust.
!! 05/03  2012   Matt Alvarado     Changed AerosolOptProp to calculate
!!                                  backscattering coefficient.
!! 05/08  2012   Matt Alvarado     Added special routine for dry BC only particles
!!                                  to particle density.
!! 05/03  2012   Matt Alvarado     Fixed calculation of hemispherical
!!                                  backscattering coefficient in AerosolOptProp
!  08/17  2012   Matt Alvarado     Added capability for external mixtures to ShellRefIndAndRad
!! 11/08  2012   Matt Alvarado     Added Maxwell-Garnett mixing rule to ShellRefIndAndRad
!! 01/24  2013   Matt Alvarado     Updated OrgCurvatureCorrection to give 1.0
!!                                  if particle number = 0.0
!! 07/10  2013   Matt Alvarado     Made it so only RecalculateRadius resets 
!!                                  the particle%surfacetension by calling 
!!                                  SurfaceTension, and CurvatureCorrection
!!                                  uses the particle%surfacetension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:		!!
!! 1. SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel)
!! 2. FUNCTION ParticleKnudsen (InParticle)
!! 3. FUNCTION ParticleKnudsenForEnergy (InParticle)
!! 4. FUNCTION CunninghamSlipCorrectionFactor (InParticle)
!! 5. FUNCTION AerosolPH (InParticle)
!! 6. FUNCTION ParticleDensity (InParticle) !Sets radius and effective radius
!! 7. FUNCTION ParticleMass (InParticle)
!! 8. FUNCTION GetPatwardhandAndKumarWeightings (InParticle)
!! 9. FUNCTION HypotheticalElectrolyteConcentrations (InParticle)
!!10. FUNCTION TerminalVelocity (InParticle)
!!11. SUBROUTINE RecalculateRadius (InParticle) !Calls particle density
!!12. FUNCTION CurvatureCorrection (InParticle)
!!13. FUNCTION SurfaceTension (InParticle)
!!14. FUNCTION ReynoldsNumber (InParticle)
!!15. REAL*8 FUNCTION OrgCurvatureCorrection (InParticle)
!!16. REAL*8 FUNCTION MolalityAqCarbon(InParticle)
!!17. REAL*8 FUNCTION AqueousSolutionMass
!!18. SUBROUTINE ShellRefIndAndRad(InParticle)
!!19. REAL*8 FUNCTION AerosolOpticalDepth()
!!20. SUBROUTINE AerosolOptProp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
! Modified by CMB (AER): Gfortran complains when "Particle" is used both
!                        as a structure type declarator and a local variable
!                        name; I have replaced the variable names with
!                        "p_Particle"

SUBROUTINE ResetAllOriginPressureLevels (OriginPressureLevel)

	IMPLICIT NONE

	Type(Particle),POINTER :: p_Particle
	REAL    :: OriginPressureLevel

	p_Particle => Particles%First

	IF (ASSOCIATED(p_Particle)) p_Particle%OriginPressureLevel = OriginPressureLevel
		
	DO WHILE (ASSOCIATED(p_Particle%Next))
		p_Particle => p_Particle%Next
		IF (ASSOCIATED(p_Particle)) p_Particle%OriginPressureLevel = OriginPressureLevel
	END DO

	RETURN
END SUBROUTINE ResetAllOriginPressureLevels

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number of the Particle	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsen (InParticle)

	USE GridPointFields, ONLY : MeanFreePathOfAir

	IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsen = MeanFreePathOfAir ()/ InParticle%EffectiveRadius
	
	!write(*,*) 'CRL: ParticleKnudsen = ',ParticleKnudsen
	!write(*,*) 'CRL: MeanFreePathOfAir = ',MeanFreePathOfAir ()
	!write(*,*) 'CRL: InParticle%EffectiveRadius= ',InParticle%EffectiveRadius
 	!write(*,*) ''

END FUNCTION ParticleKnudsen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the Knudsen Number Of Energy of the Particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleKnudsenForEnergy (InParticle)

	USE GridPointFields, ONLY : MolecularThermalDiffusivity, & 
                                    GetThermalVelocityOfAir

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle

	ParticleKnudsenForEnergy = 3 *MolecularThermalDiffusivity() &
			/ GetThermalVelocityOfAir()   &
			/ InParticle%EffectiveRadius

	RETURN
END FUNCTION ParticleKnudsenForEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the CUNNINGHAM SLIP CORRECTION FACTOR, which !!
!! has a number of possible forms.  This function selects !!
!! the appropriate one and returns the value.			  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CunninghamSlipCorrectionFactor (InParticle)

    IMPLICIT NONE

	Type(Particle),POINTER :: InParticle
	REAL*8		   :: Knudsen

	Knudsen = ParticleKnudsen (InParticle)

	!! From Kasten 1968 (via Jacobson, 1997)
	CunninghamSlipCorrectionFactor = 1 + Knudsen * (1.249 + 0.42 * dexp(-0.87 / Knudsen))

	!write(*,*) 'CRL: CunninghamSlipCorrectionFactor = ',CunninghamSlipCorrectionFactor
	!write(*,*) 'CRL: Knudsen = ',Knudsen
	
	RETURN
END FUNCTION CunninghamSlipCorrectionFactor



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns the pH of an aerosol. !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AerosolPH (InParticle)

	USE ModelParameters, ONLY : ProtonIndex

	IMPLICIT NONE

	!! External Variables
	TYPE(PARTICLE),POINTER :: InParticle

	!!Internal Variables
	REAL*8 :: ProtonMolarity


	IF (InParticle%NumberOfParticles .EQ. 0. .OR. InParticle%Dry) THEN
		AerosolPH = 14.

		RETURN
	END IF

	IF (InParticle%AqChems(ProtonIndex) .GT. 0) THEN
		!Molarity in mol/L solution
		ProtonMolarity = InParticle%AqChems(ProtonIndex)*InParticle%SolutionDensity*1000/AqueousSolutionMass(InParticle)
		
		AerosolPH = -1.*LOG10(ProtonMolarity)

	ELSE
		AerosolPH = 14.

	END IF


END FUNCTION AerosolPH



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE DENSITY for a lagrangian particle !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleDensity (InParticle)
!MJA 05-008-2012 Added special routine for dry BC only particles

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, &
                                        moles, pi, avogadro,    &
                                        SmallestAerosolPossible, & 
									micron ! cmb add

	USE Chemistry,		ONLY : IFSOLID, AqMolecularMass,       &
                                       HowManyAqChems,HowManyAqCations,&
				       HowManyAqAnions,                &
                                       HowManyAqEqReactions,           & 
                                       AqEquilibriaList,	       &
				       AqAnionCharge, AqCationCharge,  &
				       HowManyOrgChems,                &
                                       OrgMolecularMass, OrgDensity,   &
                                       SolidSaltDensity,               &
                                       HowManyAqOrgChems,              &
                                       AqOrgMolecularMass,             &
                                       AqOrgDensity

	USE InfrastructuralCode, ONLY : ERROR,REAL2STR, WARN, IsNaN, int2str

	! cmb add
	!use outputroutines, only : spillbeans
	! 
	
	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I,K, SkipFlag, C
	REAL*8  ::Volume,Mass,ionMass,PSIij, TempElectrolytes(HowManyAqChems),&
	BinaryDensity, NumerSum, DenomSum, TempMolality, &
  	PatKumWeights(HowManyAqEqReactions,3), XX, YY, Residual, dR, &
	SolutionDensity, InsolMass, InsolVolume, InsolDensity, InorgSolnDensity

	REAL*8, PARAMETER :: ResidualThreshold = 0.0025 !**3.

        LOGICAL :: ForceSolnDens

	real*8 :: sum_electrolytes

        !Use this flag to force the aqueous solution density to be 1.0 g/cm3 and skip a lot of code
        ForceSolnDens = .TRUE.
        
        ! CB: Modify this equality check to have a designated error 
        !     tolerance        
        !IF (InParticle%NumberOfParticles .LE. 0. ) THEN
        if (inparticle%numberofparticles .le. 1.0E-40) then
			ParticleDensity = 1.0
			InParticle%EmbryoRadius    = 0.
			InParticle%EffectiveRadius = 0.
			InParticle%InsolubleRadius = 0.
			InParticle%ParticleDensity = 1.
			InParticle%SolutionDensity = 1.
			InParticle%InsolubleDensity = 1.
			InParticle%InorgSolnDensity = 1.
			RETURN
		Else
			! cb add
			!call system("/bin/rm -f fort.2222")
			!call spillbeans(tout_unit=2222)
			!write(*,*) 'Nonzero # of particles (', inparticle%numberofparticles, ')'
			!write(*,*) ''
			!write(*, *)"Number ",trim(INT2STR(I)),",  ",inparticle%numberofparticles,inparticle%sectional
			!write(*, *)"Effective Radius: ",trim(real2str(inparticle%EffectiveRadius/micron)),",  (",trim(real2str(inparticle%Edges(1)/micron))," to ",trim(real2str(inparticle%Edges(2)/micron)),")"
			!write(*, *)"Embryo Radius: ",trim(real2str(inparticle%EmbryoRadius/micron))
			!write(*, *)"Insoluble Radius: ",trim(real2str(inparticle%InsolubleRadius/micron))
			!write(*, *)"Absorbing Core Radius: ",trim(real2str(inparticle%AbsCoreRad/micron))
			!write(*, *)"Shell Real Ref. Ind. 550 nm",trim(real2str(inparticle%ShellRealRefracInd(1)))
			!write(*, *)"Shell Imag. Ref. Ind. 550 nm",trim(real2str(inparticle%ShellImagRefracInd(1)))
			!write(*, *)"Extinction Efficiency: 550 nm",trim(real2str(inparticle%ExtEff(1)))
			!write(*, *)"Scattering Efficiency: 550 nm",trim(real2str(inparticle%ScaEff(1)))
			!write(*, *)"BackScattering Efficiency: 550 nm",trim(real2str(inparticle%BackScaEff(1)))
			!write(*, *)"Assymetry Parameter: 550 nm",trim(real2str(inparticle%AssymParam(1)))
		
			!write(*, *)"Temperature ",trim(real2str(inparticle%Temperature))
			!write(*, *)"Ionic Strength ",trim(real2str(inparticle%IonicStr))
			!write(*, *)"Water Activity ",trim(real2str(inparticle%WaterActivity))
			!write(*, *)"SolutionDensity ",trim(real2str(inparticle%SolutionDensity))
			!!write(*, *)"InsolubleDensity ",trim(real2str(inparticle%InsolubleDensity))
			!write(*, *)"ParticleDensity ",trim(real2str(inparticle%ParticleDensity))
			!write(*, *)"Ph ", trim(real2str(AerosolPh(inparticle)))
			!write(*, *)"Curvature Correction ", trim(real2str(CurvatureCorrection(inparticle)))
			!write(*, *)"OrgCurvature Correction ", trim(real2str(OrgCurvatureCorrection(inparticle)))	
			!write(*, *)"Solution Surface Tension", trim(real2str(inparticle%SurfaceTension))
			
			!write(*, *)"Temperature ", trim(real2str(inparticle%temperature))
			
			if (inparticle%numberofparticles .GT. 0.0) then
				!write(*,*) '--> calculating HypotheticalElectrolyteConcentrations <--'
				!		SaltConc = HypotheticalElectrolyteConcentrations (inparticle)
			endif   
			!write(*,*) ''
			!write(*,*) 'AqChems:'
			!write(*,*) '--------'
			!do i = 1, size(inparticle%aqchems)
		!		write(*,*) 'i = ', i, ', inparticle%aqchems(i) = ', inparticle%aqchems(i)
		!	enddo
		!	write(*,*) ''
			!write(*,*) 'OrgChems:'
			!write(*,*) '---------'
			!do i = 1, size(inparticle%orgchems)
		!		write(*,*) 'i = ', i, ', inparticle%orgchems(i) = ', inparticle%orgchems(i)
		!	enddo
		!	write(*,*) ''
		!	write(*,*) 'AqOrgChems:'
		!	write(*,*) '-----------'
		!	do i = 1, size(inparticle%aqorgchems)
		!		write(*,*) 'i = ', i, ', inparticle%aqorgchems(i) = ', inparticle%aqorgchems(i)
		!	enddo
			! end cb add
		endif

	!MJA 05-08-2012 Special routine for dry (BC only) particles
        IF (InParticle%Dry) THEN
		ParticleDensity = OrgDensity(1) !BC density
		InsolVolume = InParticle%OrgChems(1)*OrgMolecularMass(1)/OrgDensity(1) 
                InParticle%EmbryoRadius    = 0.0
		InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333
                InParticle%EffectiveRadius = InParticle%InsolubleRadius

		InParticle%ParticleDensity = OrgDensity(1)
		InParticle%SolutionDensity = OrgDensity(1)
		InParticle%InsolubleDensity = OrgDensity(1)
		InParticle%InorgSolnDensity = OrgDensity(1)
		RETURN
	END IF

        Mass = 0 ; Volume = 0; ionMass = 0; NumerSum = 0; DenomSum = 0
        IF (ForceSolnDens) THEN
            InorgSolnDensity = 1.0
            DO I = 1, HowManyAqChems+HowManyAqAnions+HowManyAqCations
                Mass = Mass    + AqMolecularMass(I)*InParticle%AqChems(I)
            ENDDO
            Volume = Mass*InorgSolnDensity
        ELSE        


           !! Get electrolyte equivalents of solution containing dissociated ions,
           !! which assumes a two-parameter fit proscribed in Resch, 1995 
           !! (assuming that the density of water does not change with Temperature).
           !! This is a SOLUTION density, which includes the incorporated water.
           !! So we scale up by the mass of the water as well!
           TempElectrolytes = HypotheticalElectrolyteConcentrations (InParticle)
           ! cmb in: 
           !write(*,*) 'CMB: Output from HypotheticalElectrolyteConcentrations...'
           call flush(6)
           ! end cmb
           do i = 1, size(TempElectrolytes)
              !write(*,*) 'i = ', i, ', TempElectrolytes = ', TempElectrolytes(i)
              sum_electrolytes = sum_electrolytes + TempElectrolytes(i)
           enddo
           ! end cmb in
	
           !For some reason, This routine returns NaN if there isn't a write
           !statement here. MJA, 070507
           WRITE(1234,*) TempElectrolytes(1)
	   ! cb add
           close(1234)

           !! Get the weightings for the Patwardhan and Kumar mixing rule
           !WRITE(*,*)'Just Before PatKumWeights spilling beans'
           !call spillbeans()
           PatKumWeights = GetPatwardhandAndKumarWeightings (InParticle)

           !! Loop over each reaction and consider the appropriate equation
           DO I = 1, HowManyAqEqReactions
              !write(*,*) 'i (of HowManyAqEqReactions) = ', i

              !! Some quasi-electrolytes (such as HSO4-) are not quantifyable in this
              !! density parameterization and should have been pushed to upper level
              !! electrolytes anyways (such as H2SO4).  So we ignore those electrolytes!
              ! CMB (AER, Inc): floating-point inequality check again?
              IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .ge. 1.0e-40) THEN ! GT 0 disallows water
                 !IF (AqEquilibriaList(I,1) .LE. HowManyAqChems .AND. AqEquilibriaList(I,1) .GT. 0.) THEN ! GT 0 disallows water

                 !! This is the equivalent molality of the associated electrolyte
                 !write(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(I,1))),I
                 !Write(*,*) "Temp Electrolytes: ", TempElectrolytes(AqEquilibriaList(I,1))
                 !Write(*,*) "InParticle%AqChems(1) ", InParticle%AqChems(1)
                 !Write(*,*) "moles ", moles
                 !Write(*,*) "grams ", grams
                 !Write(*,*) "AqMolecularMass(1) ", AqMolecularMass(1)			! cmb: changed int to nint?  not sure if necessary...
                 TempMolality = TempElectrolytes(nINT(AqEquilibriaList(I,1))) / InParticle%AqChems(1) / &
                      (moles / grams * AqMolecularMass(1)) * 1000.

                 !Mass    = Mass    + AqMolecularMass(INT(AqEquilibriaList(I,1)))*TempElectrolytes(INT(AqEquilibriaList(I,1)))
                 !ionMass = ionMass + AqMolecularMass(INT(AqEquilibriaList(I,1)))*TempElectrolytes(INT(AqEquilibriaList(I,1)))
                 Mass    = Mass    + AqMolecularMass(nINT(AqEquilibriaList(I,1)))*TempElectrolytes(nINT(AqEquilibriaList(I,1)))
                 ionMass = ionMass + AqMolecularMass(nINT(AqEquilibriaList(I,1)))*TempElectrolytes(nINT(AqEquilibriaList(I,1)))

                 IF(IsNaN (TempMolality)) THEN
                    WRITE(*,*) "Problem in Temp Molality"
                    WRITE(*,*) InParticle%NumberOfParticles, I, AqEquilibriaList(I,13)
                    WRITE(*,*) TempElectrolytes(INT(AqEquilibriaList(I,1))), TempMolality, InParticle%AqChems(1)
                 END IF

                 !WRITE(*,*) TempMolality, I
                 !! Binary density is an exponential function of the two parameters input into the system
                 IF(TempMolality .LT. 100.0) THEN

                    !WRITE(*,*) "Before Exp, Rxn: ", I
                    !WRITE(*,*) "TempMolality: ", TempMolality
                    BinaryDensity = AqEquilibriaList(I,12) + EXP(-1.*AqEquilibriaList(I,13)*TempMolality)*		&
                         (1.*grams/cm/cm/cm-AqEquilibriaList(I,12))

                    !This prevents the binary density of NH4OH from going below 0.9 g/cm3
                    IF(AqEquilibriaList(I,13) .LT. 0. .AND. BinaryDensity .LT. 0.9) BinaryDensity = 0.9

                 ELSE IF(TempMolality .GE. 100.0 .AND. AqEquilibriaList(I,13) .GT. 0.) THEN
                    !If the molality gets very large, and exponent goes to zero
                    !(This keeps EXP from giving an overflow error)
                    BinaryDensity = AqEquilibriaList(I,12)
                 ELSE !For NH4OH at high molality, keep it above 0.9 g/cm3
                    BinaryDensity = 0.9	
                 END IF

                 !! Employ Mixing Rule of Patwardhan and Kumar (AIChE J. 39, 711-714) to Mix these ions appropriately
                 PSIij =  1000. * PatKumWeights(I,1) + AqMolecularMass(INT(AqEquilibriaList(I,1))) * TempMolality
                 !write(*,*) AqPhaseChemicalNames(INT(AqEquilibriaList(I,1)))
                 !write(*,*) 'PatKumWeights(I,1) = ', PatKumWeights(I,1)
                 !write(*,*) 'AqMolecularMass(INT(AqEquilibriaList(I,1))) = ',AqMolecularMass(INT(AqEquilibriaList(I,1)))	
                 !write(*,*) 'TempMolality = ',TempMolality

                 IF(IsNaN(PatKumWeights(I,1))) THEN
                    CALL WARN("Problem in PatKumWeights! Get NaN!")
                 END IF
                 !! Composite Density is the Numerator Sum divided by the Denomenator Sum, although we will track it as Volume
                 NumerSum = NumerSum + PSIij
                 DenomSum = DenomSum + PSIij / BinaryDensity
              END IF
           END DO

           !WRITE(*,*) "After Ion loop"
           !PAUSE
           !! Now use the composite density to add volume to the system
           Volume = Volume + (ionMass + AqMolecularMass(1)*InParticle%AqChems(1)) * DenomSum / NumerSum
           Mass   = Mass   + AqMolecularMass(1)*InParticle%AqChems(1)
           !write(*,*) 'Volume = ', Volume
           !write(*,*) 'ionMass = ',ionMass	
           !write(*,*) 'AqMolecularMass(1) = ',AqMolecularMass(1)
           !write(*,*) 'InParticle%AqChems(1) = ',InParticle%AqChems(1)
           !write(*,*) 'DenomSum = ',DenomSum
           !write(*,*) 'NumerSum = ',NumerSum

           !WRITE(*,*) "Ion solution Volume: ", Volume, DenomSum, NumerSum

           !If for some reason the above doesn't work, set inorganic solution density to 1.0
           IF(ISNaN(Volume)) THEN
              CALL WARN("Setting ion solution density to 1.0 g/cm3")
              Volume = Mass 
              !Volume = 1.0 
           END IF

           !This only includes water and inorganic ions (needed for optical property calculation)
           InorgSolnDensity = Mass/Volume
           !WRITE(*,*) "Inorganic Solution Density: ", InorgSolnDensity
        ENDIF


	!Now add in contributions from Hydrophilic Organics (assumed constant density)
	!! This assumes a fixed density for these aqueous phase chemicals
	!WRITE(*,*) "Before Aq Org loop Volume = ", Volume
	DO I = 1, HowManyAqOrgChems
		!write(*,*)
		Mass   = Mass   + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
		Volume = Volume + InParticle%AqOrgChems(I)*AqOrgMolecularMass(I)/AqOrgDensity(I) 
		!write(*,*) 'AqOrgPhaseChemicalNames(I) = ',AqOrgPhaseChemicalNames(I)
		!write(*,*) 'Volume =', Volume
		!write(*,*) 'Mass = ',Mass	
		!write(*,*) 'AqOrgMolecularMass(I) = ',AqOrgMolecularMass(I)
		!write(*,*) 'InParticle%AqOrgChems(I) = ',InParticle%AqOrgChems(I) 
		!write(*,*) 'AqOrgDensity(I)  = ',AqOrgDensity(I) 
		!write(*,*) 'DenomSum = ',DenomSum
		!write(*,*) 'NumerSum = ',NumerSum

	END DO
	
	!WRITE(*,*) "After Aq Org loop", Volume
! 28 September 2016 CMB: expanding this statement to let the user know why there is a problem
! 7 October 2016 CMB: just set them to 1 again...ugh...
!	IF(ISNaN(Volume)) STOP
!	if (IsNaN(Volume)) then
!		call error("In ParticleAttributes.h; volume is NaN after Aq Org loop, stop statement called")
!	endif
	if (isNaN(volume)) then
		call warn("Setting aq org solution density to 1.0 g/cm3")
		call flush(6)
		volume = mass
	endif

	!! This is purely the solution density (water, ions, non-solid electrolytes (e.g. NH3), aq organics)
	SolutionDensity =  Mass/Volume

	!!Calculate radius of the solution embryo
	InParticle%EmbryoRadius = (0.75 * Volume / pi)**0.33333333333
	!WRITE(*,*) "Embryo Radius: ", InParticle%EmbryoRadius
	
	!! Consider the contributions of Solid salts (except ions and water)
	InsolVolume = 0.
	InsolMass = 0.
	DO I = 2, HowManyAqChems
		SkipFlag = 0
		!Skip non-solid electrolytes (NH3, CO2, etc).
		DO C = 1, HowManyAqEqReactions
			! CMB
			if (dabs(aqequilibrialist(c,1) - I) .lt. 0.000001 .and. &
					dabs(AqEquilibriaList(c,6)) .lt. 1.0e-12) then
			!IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				SkipFlag = 1
			END IF
		END DO
			
		IF(SkipFlag .EQ. 0) THEN
			InsolMass   = InsolMass   + AqMolecularMass(I)*InParticle%AqChems(I)
			InsolVolume = InsolVolume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)

			Mass   = Mass   + AqMolecularMass(I)*InParticle%AqChems(I)
			Volume = Volume + InParticle%AqChems(I)*AqMolecularMass(I)/SolidSaltDensity(I)
		END IF
	END DO
	
	!WRITE(*,*) "After Solid Salts loop", Volume
        ! CMB: Revise to add more info to exit
		! update again: set to 1 if NaN
!	IF(ISNaN(Volume)) STOP
        if (isnan(volume)) then
            call warn("volume is NaN after solid salts loop")
		call flush(6)
        	volume = mass	! cmb add
		endif

	 ! cmb add
 !778 continue
 	 ! end cmb add
	 
	!Now add in controbutions from hydrophobic organics	
	DO I = 1, HowManyOrgChems
		!write(*,*) 'In loop'
		!write(*,*) 'InsolMass = ',InsolMass
		!write(*,*) 'OrgMolecularMass(I) = ',OrgMolecularMass(I)	

		IF(ISNaN(InParticle%OrgChems(I))) THEN
			InParticle%OrgChems(I) = 0.
		END IF

		!write(*,*) 'InParticle%OrgChems(I) = ',InParticle%OrgChems(I),OrgPhaseChemicalNames(I)
		!write(*,*) 'InsolVolume = ',InsolVolume
		!write(*,*) 'OrgDensity(I) = ',OrgDensity(I)
	
		InsolMass   = InsolMass   + OrgMolecularMass(I)*InParticle%OrgChems(I)
		InsolVolume = InsolVolume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	
		!WRITE(*,*) I, InsolVolume
		Mass   = Mass   + OrgMolecularMass(I)*InParticle%OrgChems(I)

		Volume = Volume + InParticle%OrgChems(I)*OrgMolecularMass(I)/OrgDensity(I) 
	END DO
		
!	WRITE(*,*) "After Hydrophobic Org loop", Volume
        ! CMB: Revise to add more info to exit
		! update: make NaN into 1 and don't exit'
!	IF(ISNaN(Volume)) STOP
        if (isnan(volume)) then
            !call error("In ParticleAttributes.h; volume is NaN after hydrophobic org loop")
		    call warn("volume is NaN after hydrophobic org loop")
        	volume = mass
		endif
	
	!!This is the average density of the non-aqeous part (solid salts, BC, and hydrophobic organics)
	InsolDensity = InsolMass/InsolVolume
	InParticle%InsolubleRadius = (0.75 * InsolVolume / pi)**0.33333333333

	!! This is the average density of the particle
	ParticleDensity = Mass / Volume
	!IF(Mass .LT. 0.0 ) WRITE(*,*) "Mass (<0), Volume = ", Mass, Volume

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! RESET EFFECTIVE RADIUS, too. !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	InParticle%EffectiveRadius = (0.75 * Volume / pi)**0.33333333333


	InParticle%InorgSolnDensity = InorgSolnDensity 
	InParticle%SolutionDensity = SolutionDensity
	InParticle%InsolubleDensity = InsolDensity
	InParticle%ParticleDensity = ParticleDensity

	!WRITE(*,*) "End Particle Density", Volume, Mass
	


	RETURN
END FUNCTION ParticleDensity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Particle Mass !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ParticleMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,		 ONLY : AqMolecularMass, HowManyAqChems,	&
					HowManyAqAnions, HowManyAqCations,	&
					AqCationMass, AqAnionMass, &
					HowManyOrgChems, OrgMolecularMass, &
					HowManyAqOrgChems, AqOrgMolecularMass 

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I

	ParticleMass = 0. 

	!! Consider the contributions of Aqueous Chemicals (except ions)
	!! This assumes a fixed density for these aqueous phase chemicals
	DO I = 1, HowManyAqChems
		ParticleMass   = ParticleMass + AqMolecularMass(I)*InParticle%AqChems(I)
	END DO

	DO I = 1, HowManyAqCations
		ParticleMass   = ParticleMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		ParticleMass   = ParticleMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyOrgChems
		ParticleMass   = ParticleMass + OrgMolecularMass(I)*InParticle%OrgChems(I)
	END DO

	DO I = 1, HowManyAqOrgChems
		ParticleMass   = ParticleMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO

	RETURN
END FUNCTION ParticleMass


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get a vector of weightings for the Patwardhan and Kumar	   !!
!! y_ij factors corrected for partial dissociation, incomplete !!
!! reaction sets, etc.										   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Returns 1). y_ij  2). Cation Index 3). Anion Index		   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Patwardhan and Kumar (AIChE J. 39, 711-714)				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GetPatwardhandAndKumarWeightings (InParticle)

	USE Chemistry, ONLY : HowManyAqEqReactions,	&
			      HowManyAqChems,		&
			      HowManyAqCations,		&
			      HowManyAqAnions,		&
			      AqEquilibriaList,		&
			      AqCationCharge,		&
			      AqAnionCharge,            &
			      AqPhaseChemicalNames

	USE InfrastructuralCode, ONLY : WARN, IsNaN

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: GetPatwardhandAndKumarWeightings(HowManyAqEqReactions,3)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: K, CatIndex, AnIndex
	REAL*8  :: II, JJ,														  &
			   CationicStrFr(HowManyAqCations),AnionicStrFr(HowManyAqAnions), &
			   CationChgFrac(HowManyAqCations),AnionChgFrac(HowManyAqAnions)

	!! Update the ionic strength
	!WRITE(*,*) "CRL: Before CalculateIonic Strength... "
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%IonicStr
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%AqChems
	!CALL CalculateIonicStrength(InParticle)
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%IonicStr
	!WRITE(*,*) "CRL: Ionic Strength: ", InParticle%AqChems
	IF(IsNAN(InParticle%IonicStr)) CALL WARN("Ionic Strength Error in PKweightings!")

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Charge Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Calculate Cationic First:
	II = 0.
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationMolality(K, InParticle) * AqCationCharge(K)
		II               = II + CationChgFrac(K)

	END DO
	DO K = 1, HowManyAqCations
		CationChgFrac(K) = CationChgFrac(K) /II
		IF(IsNAN(CationChgFrac(K))) CALL WARN("CationChgFrac error in PKweightings!")
		IF(IsNAN(CationChgFrac(K))) WRITE(*,*)  CationMolality(K, InParticle),AqCationCharge(K),II
	END DO

	!! And Anionic Next:
	II = 0.
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionMolality(K, InParticle) * AqAnionCharge(K)
		II              = II + AnionChgFrac(K)
	END DO
	DO K = 1, HowManyAqAnions
		AnionChgFrac(K) = AnionChgFrac(K) / II
		IF(IsNAN(AnionChgFrac(K))) CALL WARN("AnionChgFrac error in PKweightings!")
	END DO

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Make Ionic Strength Fraction Vectors !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (Inparticle%IonicStr .EQ. 0.) THEN
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = 0.
		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = 0.
		END DO
	ELSE
		DO K = 1, HowManyAqCations
			CationicStrFr(K) = CationMolality(K, InParticle) * AqCationCharge(K) * AqCationCharge(K) / 2. / InParticle%IonicStr
			IF(IsNAN(CationicStrFr(K))) CALL WARN("CationicStrFr error in PKweightings!")

		END DO
		DO K = 1, HowManyAqAnions
			AnionicStrFr(K)  = AnionMolality(K, InParticle)  * AqAnionCharge(K) * AqAnionCharge(K) /2. / InParticle%IonicStr
			IF(IsNAN(AnionicStrFr(K))) CALL WARN("AnionicStrFr error in PKweightings!")
		END DO
	END IF

	
	!! This will be the sum of all the y_ij Patwardhan and Kumar mixing rule terms
	JJ = 0.

	!! Loop over all reactions
	DO K = 1, HowManyAqEqReactions
		

		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! For partial dissociation, consider the reactions as if the partially !!
		!! dissociated form exchanges directly with the fully dissociated form, !!
		!! which is not how it would have been written.							!!
		IF(AqEquilibriaList(K,1) .LE. HowManyAqChems) THEN ! electrolyte uncharged
			CatIndex = AqEquilibriaList(K,2)
			AnIndex  = AqEquilibriaList(K,3)
		ELSE											   ! electrolyte charged

			!! Is the partially dissociated "electrolyte" a cation or an anion?
			IF (AqEquilibriaList(K,1) .LE. HowManyAqChems+HowManyAqCations) THEN  ! It's a Cation
				CatIndex = AqEquilibriaList(K,1)-HowManyAqChems
				AnIndex  = AqEquilibriaList(K,3)
			ELSE																  ! It's an Anion
				CatIndex = AqEquilibriaList(K,2)
				AnIndex  = AqEquilibriaList(K,1)-HowManyAqChems-HowManyAqCations
			END IF
		END IF
			
		!! Store the indices used to calculate this	
		GetPatwardhandAndKumarWeightings(K,2) = REAL(CatIndex)
		GetPatwardhandAndKumarWeightings(K,3) = REAL(AnIndex)

		
		!Exclude Levitocite reaction from this calculation
		! cmb
		if (dabs(AqEquilibriaList(K, 22)) .gt. 1e-40) then
		!IF(AqEquilibriaList(K,22) .NE. 0) THEN
			IF (AqPhaseChemicalNames(INT(AqEquilibriaList(K,1))) .EQ. "(NH4)3H(SO4)2") THEN
				GetPatwardhandAndKumarWeightings(K,1) = 0.
			END IF
		
		!Normal Calculation
		ELSE
			!! Calculate the Appropriate Weightings for the Patwardhan and Kumar (AIChE J. 39, 711-714)
			!! Mixed Activity Coefficient:
			!!
			!!   ln(Aw_mix) = Sum[(i,j), (Isf[i]*Cf[j]+Isf[j]*Cf[i])*ln(Aw[ij])]
			GetPatwardhandAndKumarWeightings(K,1) = AnionicStrFr(AnIndex)*CationChgFrac(CatIndex) + &
												CationicStrFr(CatIndex)*AnionChgFrac(AnIndex)
		
		END IF
		
		JJ = JJ + GetPatwardhandAndKumarWeightings(K,1)
	END DO

	!! Normalize the y_ij estimates to equal one
	DO K = 1, HowManyAqEqReactions
		GetPatwardhandAndKumarWeightings(K,1) = GetPatwardhandAndKumarWeightings(K,1) / JJ 
	END DO

	RETURN
	!write(*,*) 'GetPatwardhandAndKumarWeightings =',GetPatwardhandAndKumarWeightings 
END FUNCTION GetPatwardhandAndKumarWeightings


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Estimate Dissociated Electrolyte Concentrations !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! For some routines, it is necessary to find the hypothetical
!! electrolyte concentrations that a set of dissociated ions 
!! would make up.  This function provides a vector of such electrolyte
!! abundances, in the form of a vector matching the AqChems vector
!! in the particle structure.
!!
!! The algorithm is taken from Jacobson's Fundamentals of Atmospheric
!! Modeling, 1999 page 494
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION HypotheticalElectrolyteConcentrations (InParticle)

	USE Chemistry, ONLY : HowManyAqChems, HowManyAqCations, &
                              HowManyAqAnions, HowManyAqEqReactions, &
                              AqEquilibriaList, & 
                              AqPhaseChemicalNames, AqAnionNames, &
                              AqCationNames

	USE ModelParameters, ONLY : AqThermoNumbAllowedIterationsForHypLoop, &
                                    moles,WaterEquilibriumIndex

	! CMB: add warn
	USE InfrastructuralCode, ONLY : ERROR, REAL2STR, TRANSCRIPT, warn

	IMPLICIT NONE

	!! External Variables
	REAL*8 :: HypotheticalElectrolyteConcentrations(HowManyAqChems)
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	REAL*8  :: TempChem(HowManyAqChems+HowManyAqCations+HowManyAqAnions), DeltaMoles, Sum, ContentSum, ExcludeRxnInit
	INTEGER :: I, J, K, C, ExcludeRxn 

	!! counters and settings
	K = 0
	ExcludeRxn = 0

	!! First set all electrolytes to zero, unless they are non-solid
	!write(*,*) '(HypotheticalElectrolyteConcentrations) Particle ID: ', inparticle%particleid
10	ContentSum=0
	DO I = 1, HowManyAqChems
		TempChem(I) = 0.
		!write(*,*) 'inparticle%aqchems(', i, ') = ', inparticle%aqchems(i)
		if (inparticle%aqchems(i) .lt. 0) then
			call warn("InParticle%AqChems("//trim(real2str(real(i)))//") < 0 ("//&
					trim(real2str(inparticle%aqchems(i)))//"in HypoElectConc; ID="//trim(real2str(&
					real(inparticle%particleid))))
		endif
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			if ((abs(AqEquilibriaList(C, 1) - I) .lt. 0.000001) .and.(&
					abs(AqEquilibriaList(c, 6)) .lt. 0.000001)) then
			!IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				!if (inparticle%aqchems(i) .lt. 0) then
				!	call warn("InParticle%AqChems("//trim(real2str(real(i)))//") < 0 ("//&
				!			trim(real2str(inparticle%aqchems(i)))//"in HypoElectConc")
				!endif
				!write(*,*) 'assigning to temp chem for C index using value = ', c, i, inparticle%aqchems(i)
				TempChem(I) = InParticle%AqChems(I)
			END IF
		END DO			
		 
		ContentSum = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 1 ContentSum: ", ContentSum
	!call flush(6)

	!! Then take the Ion Concentrations as is.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		TempChem(I) = InParticle%AqChems(I)
		ContentSum  = ContentSum + InParticle%AqChems(I)
	END DO

	!WRITE(*,*) "Check 2 ContentSum: ", ContentSum

	!WRITE(*,*) "Check 3"

	!! Loop over the reactions and adjust the equations
	DO J = 1, 2
	DO I = 1, HowManyAqEqReactions
		  !WRITE(*,*) I, ExcludeRxn, AqEquilibriaList(I,1), HowManyAqChems
		IF (I .EQ. ExcludeRxn) THEN
			CYCLE

		!! Some reactions may have ions serving as electrolytes
		!! We must identify these, and then push them the other direction
		ELSE IF (AqEquilibriaList(I,1) .GT. HowManyAqChems) THEN

			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))
			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,4)
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3)), TempChem(AqEquilibriaList(I,1))*AqEquilibriaList(I,5)
		
			TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) =&
                            TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) + &
                            TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,4)
			TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) = &
			    TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) + &
			    TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,5)
			TempChem (INT(AqEquilibriaList(I,1))) = TempChem(INT(AqEquilibriaList(I,1))) - &
                            TempChem(INT(AqEquilibriaList(I,1)))*AqEquilibriaList(I,5)

			!WRITE(*,*) AqCationNames(AqEquilibriaList(I,2)), TempChem (HowManyAqChems+AqEquilibriaList(I,2))
			!WRITE(*,*) AqAnionNames(AqEquilibriaList(I,3)), TempChem (HowManyAqChems+HowManyAqCations + AqEquilibriaList(I,3))
			!WRITE(*,*) TempChem (AqEquilibriaList(I,1))


		!! Other Reactions are only water dissociations
		ELSE IF (AqEquilibriaList(I,1) .LE. 0) THEN

			!! The adjustment parameter is limited by the concentration of each ion divided
			!! by the appropriate stoicheometric coefficient.
			DeltaMoles  = MIN(TempChem(HowManyAqChems+INT(AqEquilibriaList(I,2)))/AqEquilibriaList(I,4),&
                                TempChem(HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3)))/&
                                    AqEquilibriaList(I,5))

			!! Now adjust all of the concentrations appropriately
			IF (DeltaMoles .GT. 0) THEN
				!! Don't adjust back towards water at all, just get rid of it
				TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) = &
					TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) - &
                                        DeltaMoles*AqEquilibriaList(I,4)
				TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) = &
					TempChem (HowManyAqChems+HowManyAqCations + &
                                                INT(AqEquilibriaList(I,3))) - DeltaMoles * AqEquilibriaList(I,5)
			END IF

		!! And still others are canonical
		ELSE
			! CMB: more floating-point equality fixes
			!if (abs(AqEquilibriaList(I,22)) .le. 1.0e-40) then
			IF (AqEquilibriaList(I,22) .EQ. 0) THEN
			
				!! The adjustment parameter is limited by the concentration of each ion divided
				!! by the appropriate stoicheometric coefficient.
				DeltaMoles  = MIN(TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) /&
                                        AqEquilibriaList(I,4),	& 
                                        TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) /&
                                        AqEquilibriaList(I,5))
                                !WRITE(*,*) "Delta Moles: ", DeltaMoles
				!! Now adjust all of the concentrations appropriately
				IF (DeltaMoles .GT. 0) THEN
					TempChem (INT(AqEquilibriaList(I,1))) = TempChem(INT(AqEquilibriaList(I,1))) +&
                                            DeltaMoles
                                        TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) = &
                                            TempChem (HowManyAqChems+INT(AqEquilibriaList(I,2))) - &
                                            DeltaMoles * AqEquilibriaList(I,4)
					TempChem (HowManyAqChems+HowManyAqCations + INT(AqEquilibriaList(I,3))) =&
                                            TempChem (HowManyAqChems+HowManyAqCations + &
                                                    INT(AqEquilibriaList(I,3))) - DeltaMoles * AqEquilibriaList(I,5)
				END IF
		

			ELSE
				!Levitocite case: Skip this
				!write(*,*) 'line 830 in particleattributes.h'
				CYCLE						
			END IF
		
		END IF
	END DO
	END DO

	!WRITE(*,*) "Check 4"

	!! Check to make sure have accounted for all of the ions
	SUM = 0.
	DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
		!IF(IsNan(TempChem(I))) THEN
		!	WRITE(*,*) "Error in HypElec!"
		!	STOP
		!END IF
		SUM = SUM + TempChem(I)
	END DO

	!! This is an error
	!! I changed zero limit to 1.0e-3 from 1.0e-12 to eliminate some problems (Matt Alvarado, 5/10/2006)
	!! Note that this does not error so long as the amount of ions left is less than 1/1000 of the total aqueous moles!
	!! For dry particles, just runs!
	IF (SUM .GT. ContentSum*1.e-3) THEN
	!! It should not be caught in an infinite loop!

		!WRITE(*,*) "Sum, ContentSum: ", SUM, ContentSum*1.0e-3

		!! Show us what was wrong if the calculation doesn't work right.  (Probably the ContentSum*1.e-12 is not right below.)
		CALL TRANSCRIPT("")
		CALL TRANSCRIPT("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		CALL TRANSCRIPT("+++ There are problems in Hypothetical Electrolyte Concentrations +++")
		CALL TRANSCRIPT("")

		DO I = 1, HowManyAqChems+HowManyAqCations+HowManyAqAnions

		IF (TempChem(I) .gt. 0) then
		IF (I .LE. HowManyAqChems) CALL TRANSCRIPT(TRIM(AqPhaseChemicalNames(I))//&
                        TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems .AND. I .LE. HowManyAqChems+HowManyAqCations) &
		    CALL TRANSCRIPT(TRIM(AqCationNames(I-HowManyAqChems))//TRIM(REAL2STR(TempChem(I)/moles)))
		IF (I .GT. HowManyAqChems+HowManyAqCations .AND. I .LE. HowManyAqChems+HowManyAqCations+HowManyAqAnions) &
		    CALL TRANSCRIPT(TRIM(AqAnionNames(I-HowManyAqChems-HowManyAqCations))//&
                            TRIM(REAL2STR(TempChem(I)/moles)))
		END IF

		END DO

		CALL DumpParticleContentsAtError(InParticle, "FailedHypotheticalElectrolyteConcentrationsParticle.txt",	&
						 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
                                                 InDoWaterActivity=.FALSE.,InDensity=.FALSE.,InDoSurfaceTension=.FALSE., &
                                                 InIonicStrength=.TRUE.)
		!Removed this so that routines warns of charge imbalance 
		!in transcript, but does not stop
		! CMB: replace with a warn
		! never mind, it gives a seg fault?
		!CALL warn ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")
		!CALL ERROR ("HypotheticalElectrolyteConcentrations() is stuck in what appears to be an infinite loop.")

	END IF

	!! Transfer the good stuff into a shorter vector for output
	DO I = 1, HowManyAqChems
		!WRITE(*,*) I, TempChem(I)
		HypotheticalElectrolyteConcentrations(I) = TempChem(I)
		!WRITE(*,*) 'HypotheticalElectrolyteConcentrations =',TempChem(I), I,AqPhaseChemicalNames
	END DO

	RETURN
END FUNCTION HypotheticalElectrolyteConcentrations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the TERMINAL FALL VELOCITY for particles !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION TerminalVelocity (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir

	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Vairables
	REAL*8 :: Dum

	TerminalVelocity = ReynoldsNumber (InParticle) * 0.5             &
					   * GetDynamicViscosityOfAir()  &
					   / GetAirDensity()		 &
					   / InParticle%EffectiveRadius

	RETURN
END FUNCTION TerminalVelocity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Recalculate the Radius of the Particle !!
SUBROUTINE RecalculateRadius (InParticle)

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, IFAQ, OrgPhaseChemicalNames
	!USE outputroutines, only : spillbeans
	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8  :: II

	! CMB (AER, Inc): Floating-point inequality check
	if (dabs(inparticle%numberofparticles) .le. 1.0e-40) then
	!IF (InParticle%NumberOfParticles .EQ. 0. ) THEN
		IF (InParticle%Edges(1) .GT. 0.) THEN
			InParticle%EmbryoRadius    = InParticle%Edges(1)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = InParticle%Edges(1)
		ELSE
			InParticle%EmbryoRadius    = 0.01*InParticle%Edges(2)
			InParticle%InsolubleRadius = 0.
			InParticle%EffectiveRadius = 0.01*InParticle%Edges(2)
		END IF

		InParticle%ParticleDensity = 1.
		InParticle%SolutionDensity = 1.
		InParticle%InsolubleDensity = 1.
		InParticle%InorgSolnDensity = 1.
		RETURN
	END IF

	!! Calling the Density Routine Resets the Radius, so just do that
	!Write(*,*) "In RecalculateRadius InParticle%AqChems(1) ", InParticle%AqChems(1)
	!Write(*,*) "In RecalculateRadius InParticle%OrgChems(14) ", InParticle%OrgChems(14),OrgPhaseChemicalNames(14)
	!Write(*,*) "Just before calling ParticleDensity spilling beans"
	!call spillbeans()
	II = ParticleDensity (InParticle)
	!Write(*,*) "End ParticleDensity() "
        !! Reset the Surface Tension as well
        II = SurfaceTension (InParticle)	
	!Write(*,*) "End SurfaceTension() "
	!write(*,*) 'Exiting RecalculateRadius()'
	RETURN
END SUBROUTINE RecalculateRadius

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature, used for water activity and in condensation. !!
!! It is also known as the "Kelvin Effect"						   !!
!!											   !!
!! The optional argument PerturbWater is a multiple by which to change the water content   !!
!! so we can get a sense of the derivative of the water residual.			   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION CurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens
	!WRITE(*,*) "Curvature Started", SurfaceTension(InParticle), ParticleDensity(InParticle)
	
	!Ignore this for bulk mode, no particles, or if ordered to.
	IF (ThermoBulkMode .OR. InParticle%NumberofParticles .LE. 0.0D0 .OR. IgnoreCorrectionstoDiffusivity) THEN
	   CurvatureCorrection = 1.
	   !write(*,*) 'since InParticle%NumberofParticles should be le zero, CurvatureCorrection = 1',CurvatureCorrection, InParticle%NumberofParticles 	
	ELSE
	   !! Do this in stages because ParticleDensity forces a radius + solution density update...
	   !WRITE(*,*) "Before Dens"
	   TotalDens  = InParticle%SolutionDensity
	   !WRITE(*,*) "After Dens" 
		
	   !WRITE(*,*) "Curvature Correction: ", InParticle%NumberofParticles
	   !write(*,*) 'InParticle%NumberofParticles should be gt zero', InParticle%NumberofParticles
	   CurvatureCorrection = 2.*AqMolecularMass(1)/(InParticle%Temperature*Rstar*InParticle%SolutionDensity)
	   !WRITE(*,*) CurvatureCorrection, InParticle%EmbryoRadius, SurfaceTension(InParticle)
		
	   CurvatureCorrection = EXP(CurvatureCorrection*InParticle%SurfaceTension/InParticle%EmbryoRadius)
	   !WRITE(*,*) "Curvature Correction: ", CurvatureCorrection
	END IF

	!WRITE(*,*) "Curvature Okay", CurvatureCorrection 
	RETURN
	
END FUNCTION CurvatureCorrection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Surface Tension depends on Temperature and Electrolytic Properties  !!
!! of the Aerosol's constituent species.			       !!
!!								       !!
!! Follow Li and Lu (2001): Chem. Eng. Sci., 56: 2879-2888, who use    !!
!! a Gibbs dividing surface, define surface excess, etc.  It is a 2-   !!
!! parameter-per-electrolyte system and the data is taken in when the  !!
!! electrolytes are defined.					       !!
!!								       !!
!! Many other systems seem to use the sig = sig_w + Beta I formulation !!
!! from Pruppacher and Klett (1997), but that is significantly less    !!
!! effective than this technique.				       !!
!!                                                                     !!
!! Units are dyn / cm                                                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION SurfaceTension (InParticle)

	USE GridPointFields, ONLY : SurfaceTensionOfWater, &
                                    GetRelativeHumidity

	USE Chemistry,       ONLY : AqEquilibriaList,      &
                                    HowManyAqEqReactions
	
        USE ModelParameters, ONLY : UseDonnanSurfaceTension, &
                                    Rstar,                   &
                                    HowManyAqChems,          &
                                    HowManyAqAnions,         &
                                    HowManyAqCations   
	
        USE InfrastructuralCode, ONLY : Error, REAL2STR
	IMPLICIT NONE

	!! External Variables
	TYPE (Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I
	REAL*8  :: II, Sum1, DenomSum, TotalIonMolality, Test

	! CMB: equality check
	if (inparticle%numberofparticles .le. 1.0e-40) return
	!IF (InParticle%NumberOfParticles .EQ. 0.0) RETURN

	!! CALL FindElectrolyteEquilibrium (InParticle)
	CALL CalculateIonicStrength(InParticle)
        CALL KusikMeissner(InParticle)

	IF(UseDonnanSurfaceTension) THEN !If true, Donnan; If False, Jacobson formulas
		!!Donnan's Parameterization
		!! 1. + SUM (K_i * Activity_i)
		DenomSum = 1.
		DO I = 1, HowManyAqEqReactions
			DenomSum = DenomSum + MeanActivity(I, InParticle) * AqEquilibriaList(I,15)
		END DO

!		IF (Denomsum .GT. 1.0e40) THEN
!			WRITE(*,*) DenomSum, InParticle%numberofParticles
!			CALL DumpParticleContentsAtError(InParticle, "BadActivitiesInSurfaceTensionParticle.txt", &
!                                                 InRelativeHumidity=.TRUE., InRadius=.TRUE., InpH=.TRUE., &
!                                                 InDoWaterActivity=.TRUE.,InDensity=.TRUE.,InDoSurfaceTension=.FALSE., &
!                                                 InIonicStrength=.TRUE.)
!			CALL ERROR ("In SurfaceTension(), the mean activity coefficients appear to have grown without bounds.  "// &
!					"This is likely because the water is leaching out of the system towards zero.  Perhaps the "// &
!					"relative humidity for the particle is too small?  I think it is "//trim(real2str(100*         &
!					GetRelativeHumidity()))//  &
!					"%.  Another possibility is that you are in the water equilibration routine and it is dumping "// &
!					"water from the particle to try to reduce the water activity, and it got too dry.")
!		END IF

		!! Loop over each reaction and consider the appropriate equation
		Sum1 = 0.
		DO I = 1, HowManyAqEqReactions

			!! Ignore dissociation of 
			!IF ((AqEquilibriaList(I,1) .GT. HowManyAqChems) .OR. (AqEquilibriaList(I,1) .LT. 0)) CYCLE

			!! AqIonGrid(I,14) is Surface Excess and AqIonGrid(I,15) is Equilibrium Adsorption Coefficient
			!WRITE(*,*) "Mean Act.: ", MeanActivity(I, InParticle)
			!WRITE(*,*) "Ka: ", AqEquilibriaList(I,15)
			!WRITE(*,*) "DenomSum: ", DenomSum
			Sum1 = Sum1 + AqEquilibriaList(I,14) * LOG (1. - AqEquilibriaList(I,15)	&  
				* MeanActivity(I, InParticle) / DenomSum) 

		END DO

		!! Surface Tesion of the Solution is the Surface Tension of Water + R * T * Sum1
		SurfaceTension = SurfaceTensionOfWater() &
	                 + Rstar*InParticle%Temperature*Sum1
	
		!Now, add in effect of aqueous organics using formula in
		! Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534
                !Based on Fachini et al., Nature, 401, 257-259, 1999.
 		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))
	
	ELSE !Use formulas from Jacobson, "Fund. of Atm. Modeling", 2nd ed., p. 534		
		SurfaceTension = SurfaceTensionOfWater() 
		!Ions
		TotalIonMolality = 0.
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			TotalIonMolality = TotalIonMolality + Molality(I, InParticle)
		END DO
                !Based on Pruppacher and Klett, 1997.
		SurfaceTension = SurfaceTension + 1.7*TotalIonMolality
		!Organics
                !Based on Fachini et al., Nature, 401, 257-259, 1999.
		SurfaceTension = SurfaceTension - 0.0187*InParticle%Temperature*LOG(1+628.14*MolalityAqCarbon(InParticle))

	END IF

	!!Force a minimum Surface tension
	!!Use the surface tension of Benzene at 298 K 
	!!(Also used in OrgCurvatureCorrection as the surface tension)		
	!!See Seinfeld and Pandis (1998), Table 10.2, p. 558
	IF(SurfaceTension .LT. 28.21) SurfaceTension = 28.21

	!! Units are dyn / cm
	InParticle%SurfaceTension = SurfaceTension

	RETURN
END FUNCTION SurfaceTension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the PARTICLE REYNOLDS NUMBER !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Uses a three-regime (slip, continuum, and transition) !!
!! from Beard's parameterization. (see Jacobson, 1997)   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Tested against Fig. 16.3 in Jacobson, 1999 and it     !!
!! matches well.					 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION ReynoldsNumber (InParticle)

	USE GridPointFields, ONLY : GetAirDensity,GetDynamicViscosityOfAir
	USE ModelParameters, ONLY : g, grams, cm
	IMPLICIT NONE

	!! Type input variables
	TYPE (Particle),POINTER :: InParticle

	!! Type local variables
	REAL*8 :: Dv, Cs, AirDensity, X, Y

	Dv         = GetDynamicViscosityOfAir()
	Cs         = CunninghamSlipCorrectionFactor(InParticle)
	AirDensity = GetAirDensity()

	!! Estimate the Reynolds Number and then Calculate the Real Value, 
	!! this is the estimated terminal velocity
	ReynoldsNumber = 2*InParticle%EffectiveRadius*InParticle%EffectiveRadius*(InParticle%ParticleDensity-AirDensity)*g*Cs/(9*DV)

	!! The Estimated Reynolds Number
	ReynoldsNumber = 2.*InParticle%EffectiveRadius*AirDensity*ReynoldsNumber/DV 

	!! The temporary is the Reynolds Number unless Re > 0.01
	IF (ReynoldsNumber > 0.01 .and. ReynoldsNumber < 300) THEN

		X = dLog((InParticle%ParticleDensity-AirDensity)*AirDensity*InParticle%EffectiveRadius**3.*32.*g/(3.*DV*DV))

		!! The Values here taken from Prupp and Klett, p.417 and originate from Beard's paper (1976)
		Y = -.318657d1+(0.992696d0)*X-(0.153193d-2)*X**2.-(0.987059d-3)*X**3.	&
			-(0.578878d-3)*X**4.+(0.855176d-4)*X**5.-(0.327815d-5)*X**6.

		ReynoldsNumber = Cs*dExp(Y)

	ENDIF

	IF (ReynoldsNumber .GE. 300.) THEN

		!! Use X to signify Np**1/6
		X = (InParticle%SurfaceTension**3.*AirDensity*AirDensity/(DV**4. * (InParticle%ParticleDensity-AirDensity) * g))**(1./6.)

		Y = dLog(16./3.*X*InParticle%EffectiveRadius*InParticle%EffectiveRadius*	&
				 (InParticle%ParticleDensity-AirDensity)*g/InParticle%SurfaceTension)

		ReynoldsNumber = Cs*X*dExp(-5.00015+5.23778*Y-2.04914*Y**2.+0.475294*Y**3.-0.0542819*Y**4.+0.00238449*Y**5.)

	END IF

	RETURN
END FUNCTION ReynoldsNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The correction to pressure from curvature of the organic phase.			   !!
!! It is also known as the "Kelvin Effect"						   !!
!! This uses approximate values for molecular mass and surface tension                     !!
!! of the organic phase. Since the curvature correction is mainly for                      !!
!! water above RH = 98%, this should be of minor importance.				   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION OrgCurvatureCorrection (InParticle) 

	USE Chemistry, ONLY : AqMolecularMass, HowManyAqChems, HowManyAqAnions, HowManyAqCations
	USE ModelParameters, ONLY : Rstar, micron, ThermoBulkMode, IgnoreCorrectionstoDiffusivity
	USE GridPointFields, ONLY : SurfaceTensionofWater

	IMPLICIT NONE 

	!! Type input variables
	TYPE (Particle),POINTER  :: InParticle
	REAL*8 :: TotalMoles
	INTEGER :: I

	!!Type Internal variables
	REAL*8 :: TotalDens, SurfTens, Radius
	
	!Ignore this for bulk mode or no insoluble species
	IF (ThermoBulkMode .OR. IgnoreCorrectionstoDiffusivity .OR. InParticle%NumberofParticles .EQ. 0.0) THEN
		OrgCurvatureCorrection = 1.

	ELSE
		TotalDens = InParticle%InsolubleDensity	
                !write(*,*) "In OrgCurvatureCorrection ..."
                !write(*,*) "TotalDens = ", TotalDens
		!!Assume a molecular mass of 200 g/mol
		OrgCurvatureCorrection = 2.*200.0/(InParticle%Temperature*Rstar*InParticle%InsolubleDensity)

                !write(*,*) "OrgCurvatureCorrection =",OrgCurvatureCorrection
                !write(*,*) "InParticle%Temperature = ",InParticle%Temperature
		!write(*,*) "Rstar = ",Rstar
                !write(*,*) "InParticle%InsolubleDensity = ",InParticle%InsolubleDensity 

		!Keep radius above 1 nm for calculation
		IF(InParticle%InsolubleRadius .LE. 0.001*micron) THEN
			Radius = 0.001*micron
		ELSE
			Radius = InParticle%InsolubleRadius
		END IF
		
		!!Assume surface Tension of Benzene, Seinfeld and Pandis, p. 522
		SurfTens = 28.21
		OrgCurvatureCorrection = EXP(OrgCurvatureCorrection*SurfTens/Radius)
                !write(*,*) ' SurfTens = ',SurfTens
                !write(*,*) 'Radius = ',Radius

	END IF

	!WRITE(*,*) "OrgCurvature Okay"
	RETURN
	
END FUNCTION OrgCurvatureCorrection

!Molality of carbon (mol/kg H2O) in aqueous phase for surface tension calculation (Matt Alvarado, 6/2006)
REAL*8 FUNCTION MolalityAqCarbon(InParticle)

	USE Chemistry,				ONLY :	AqOrgNumCarbon
	USE ModelParameters
	
	IMPLICIT NONE

	REAL*8  :: Total
	INTEGER :: I
	TYPE(PARTICLE),POINTER :: InParticle

	
	Total = 0.0
	
	DO I = 1, HowManyAqOrgChems
		Total = Total + InParticle%AqOrgChems(I)*AqOrgNumCarbon(I)
	END DO
	
	MolalityAqCarbon = Total / InParticle%AqChems(1) / moles * grams / WaterMolecMass * 1000.

	RETURN

END FUNCTION MolalityAqCarbon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate Mass of aqueous phase (g/particle) !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 FUNCTION AqueousSolutionMass (InParticle)

	USE ModelParameters,     ONLY : grams, cm, aqchemscale, moles, pi, avogadro
	USE Chemistry,			 ONLY : AqMolecularMass, HowManyAqChems,	&
									HowManyAqAnions, HowManyAqCations,	&
									AqCationMass, AqAnionMass, &
									HowManyAqOrgChems, AqOrgMolecularMass, &
									HowManyAqEqReactions, AqEquilibriaList

	IMPLICIT NONE

	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!! Internal Variables
	INTEGER :: I, C

	AqueousSolutionMass = 0. 

	!Water
	AqueousSolutionMass   = AqueousSolutionMass + AqMolecularMass(1)*InParticle%AqChems(1)

	!Non-solid Electrolytes
	DO I = 2, HowManyAqChems
		!Check for non-solid electrolyte (NH3, CO2, etc), and add in.
		DO C = 1, HowManyAqEqReactions
			IF(AqEquilibriaList(C,1) .EQ. I .AND. AqEquilibriaList(C,6) .EQ. 0.0) THEN
				AqueousSolutionMass = AqueousSolutionMass + AqMolecularMass(I)*InParticle%AqChems(I)
			END IF
		END DO					 
	END DO

	DO I = 1, HowManyAqCations
		AqueousSolutionMass   = AqueousSolutionMass + AqCationMass(I)*InParticle%AqChems(HowManyAqChems+I)
	END DO

	DO I = 1, HowManyAqAnions
		AqueousSolutionMass   = AqueousSolutionMass + AqAnionMass(I)*InParticle%AqChems(HowManyAqChems+HowManyAqCations+I)
	END DO

	DO I = 1, HowManyAqOrgChems
		AqueousSolutionMass   = AqueousSolutionMass + AqOrgMolecularMass(I)*InParticle%AqOrgChems(I)
	END DO
	!write(*,*) 'AqueousSolutionMass   = ',AqueousSolutionMass 
	RETURN
END FUNCTION AqueousSolutionMass

!! This subroutine calculates the volume-average refractive indices for
!! the reflective shell and then calls DMILAY for the efficienies and 
!! assymetry parameter
SUBROUTINE ShellRefIndAndRad(InParticle)
! 08/30/2010 Changed for the Photolyis Wavebands MJA
! 02/27/2012 Changed to wavelength-dependent refractive indices
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added Maxwell Garnett mixing rule
	USE ModelParameters
	USE Chemistry, ONLY : AqMolecularMass, &
			        SolidSaltDensity, AqPhaseChemicalNames, &
				OrgMolecularMass, OrgDensity, &
				AqOrgMolecularMass, AqOrgDensity, &
                                CationIonicRefraction, &
                                AnionIonicRefraction,  &
                                RefIndFlag

	USE InfrastructuralCode, ONLY : INT2STR, REAL2STR, ERROR
	
	IMPLICIT NONE
	
	!! External Variables
	TYPE(Particle),POINTER :: InParticle

	!!External Subroutines 
        EXTERNAL  DMiLay

	!! Internal Variables
	INTEGER, PARAMETER  :: MaxAng = 100
 	INTEGER	:: I, JJ, NumAng, ib, w_ind_low, sizebin
	REAL*8	:: SingleParticleVolume, VolumeConc, Factor,xxx
	REAL*8	:: Er, Ei, SumEr, SumEi, Wavelength, Sum
	REAL*8	:: SolnMolWeight, MolalRefraction, MolVol, & 
                   TotalVolConc, Rshell, ShellVolumeConc
        REAL*8  :: ErBC, EiBC, ErVA, EiVA, bcvolrat, ErMG, EiMG

	REAL*8, ALLOCATABLE	:: AqVolConc(:), &
				OrgVolConc(:), AqOrgVolConc(:), MolFrac(:), &
                                wvln_array(:)
	INTEGER, ALLOCATABLE	:: wvln_rounddown(:)
        INTEGER :: numbin
        LOGICAL :: FASTTUV

        !These are refractive indices given as data statements
        !!Wavelength-dependent refractive indices from version 3.1a
        !!of Optical Properties of Aerosols and Clouds (OPAC) 
        !! M. Hess, P. Koepke, and I. Schult (1998): 
        !!Optical Properties of Aerosols and clouds: 
        !!The software package OPAC, Bull. Am. Met. Soc., 79, 831-844.
        REAL*8 :: WAVE(61), SOOT_REAL(61), SOOT_IMAG(61)
        REAL*8 :: WATER_REAL(61), WATER_IMAG(61)
        REAL*8 :: WASO_REAL(61), WASO_IMAG(61)
        REAL*8 :: INSO_REAL(61), INSO_IMAG(61)
        REAL*8 :: SEASALT_REAL(61), SEASALT_IMAG(61)
        REAL*8 :: SULF_REAL(61), SULF_IMAG(61)
        REAL*8 :: DUST_REAL(61), DUST_IMAG(61)
        REAL*8 :: a_w, b_w, n_real, n_imag, V_molal, R_molal

	COMPLEX :: SolarShellInd, SolarCoreInd, IRShellInd, IRCoreInd, VAInd
        COMPLEX :: Eshell, EBC, Factor2, EMG, MGInd
	REAL*8	:: Angles(MaxAng), CosAng(MaxAng), M1( MaxAng, 2), &
			M2( MaxAng, 2 ), S21( MaxAng, 2 ), &
		D21( MaxAng, 2 ), SolutionRefracIndex(2), SolutionVolConc

!First specify the wavelength-dependent refractive indices of
!soot, water, organics, sulfate,
!and sea salt.
!From OPAC v3.1a, downloaded from http://www.lrz.de/~uh234an/www/radaer/opac-des.html
!Reference: M. Hess, P. Koepke, and I. Schult (1998): 
!Optical Properties of Aerosols and clouds: 
!The software package OPAC, Bull. Am. Met. Soc., 79, 831-844.
!All values at 0% RH
!Following GEOS-Chem Aerosol Optics memo of Collette Heald, (Jan. 30, 2010)
!Model BC using OPAC SOOT, organics using WASO, sulfate and nitrate salts
!as sulfate, chloride salts as sea salt.

        !Wavelength in nm
        WAVE = 1.0e3*(/ &
        2.500e-01,  3.000e-01,  3.500e-01,  4.000e-01,  4.500e-01, & 
        5.000e-01,  5.500e-01,  6.000e-01,  6.500e-01,  7.000e-01, &
        7.500e-01,  8.000e-01,  9.000e-01,  1.000e+00,  1.250e+00, &
        1.500e+00,  1.750e+00,  2.000e+00,  2.500e+00,  3.000e+00, &
        3.200e+00,  3.390e+00,  3.500e+00,  3.750e+00,  4.000e+00, &
        4.500e+00,  5.000e+00,  5.500e+00,  6.000e+00,  6.200e+00, &
        6.500e+00,  7.200e+00,  7.900e+00,  8.200e+00,  8.500e+00, &
        8.700e+00,  9.000e+00,  9.200e+00,  9.500e+00,  9.800e+00, &
        1.000e+01,  1.060e+01,  1.100e+01,  1.150e+01,  1.250e+01, &
        1.300e+01,  1.400e+01,  1.480e+01,  1.500e+01,  1.640e+01, &
        1.720e+01,  1.800e+01,  1.850e+01,  2.000e+01,  2.130e+01, &
        2.250e+01,  2.500e+01,  2.790e+01,  3.000e+01,  3.500e+01, &
        4.000e+01/)

        SOOT_REAL=(/ & !From file soot00 for black carbon
        1.620e+00,  1.740e+00,  1.750e+00,  1.750e+00,  1.750e+00, & 
        1.750e+00,  1.750e+00,  1.750e+00,  1.750e+00,  1.750e+00, &
        1.750e+00,  1.750e+00,  1.750e+00,  1.760e+00,  1.760e+00, &
        1.770e+00,  1.790e+00,  1.800e+00,  1.820e+00,  1.840e+00, &
        1.860e+00,  1.870e+00,  1.880e+00,  1.900e+00,  1.920e+00, &
        1.940e+00,  1.970e+00,  1.990e+00,  2.020e+00,  2.030e+00, &
        2.040e+00,  2.060e+00,  2.120e+00,  2.130e+00,  2.150e+00, &
        2.160e+00,  2.170e+00,  2.180e+00,  2.190e+00,  2.200e+00, &
        2.210e+00,  2.220e+00,  2.230e+00,  2.240e+00,  2.270e+00, &
        2.280e+00,  2.310e+00,  2.330e+00,  2.330e+00,  2.360e+00, &
        2.380e+00,  2.400e+00,  2.410e+00,  2.450e+00,  2.460e+00, &
        2.480e+00,  2.510e+00,  2.540e+00,  2.570e+00,  2.630e+00, &
        2.690e+00/) 
        SOOT_IMAG=(/ &
       -4.500e-01, -4.700e-01, -4.650e-01, -4.600e-01, -4.550e-01, &
       -4.500e-01, -4.400e-01, -4.350e-01, -4.350e-01, -4.300e-01, &
       -4.300e-01, -4.300e-01, -4.350e-01, -4.400e-01, -4.500e-01, &
       -4.600e-01, -4.800e-01, -4.900e-01, -5.100e-01, -5.400e-01, &
       -5.400e-01, -5.495e-01, -5.600e-01, -5.700e-01, -5.800e-01, &
       -5.900e-01, -6.000e-01, -6.100e-01, -6.200e-01, -6.250e-01, &
       -6.300e-01, -6.500e-01, -6.700e-01, -6.800e-01, -6.900e-01, &
       -6.900e-01, -7.000e-01, -7.000e-01, -7.100e-01, -7.150e-01, &
       -7.200e-01, -7.300e-01, -7.300e-01, -7.400e-01, -7.500e-01, &
       -7.600e-01, -7.750e-01, -7.900e-01, -7.900e-01, -8.100e-01, &
       -8.200e-01, -8.250e-01, -8.300e-01, -8.500e-01, -8.600e-01, &
       -8.700e-01, -8.900e-01, -9.100e-01, -9.300e-01, -9.700e-01, &
       -1.000e+00/) 

        WATER_REAL=(/ & !Refractive index from cuma00 file
        1.362e+00,  1.349e+00,  1.343e+00,  1.339e+00,  1.337e+00, &
        1.335e+00,  1.333e+00,  1.332e+00,  1.331e+00,  1.331e+00, &
        1.330e+00,  1.329e+00,  1.328e+00,  1.327e+00,  1.323e+00, &
        1.321e+00,  1.313e+00,  1.306e+00,  1.261e+00,  1.371e+00, &
        1.478e+00,  1.423e+00,  1.400e+00,  1.369e+00,  1.351e+00, &
        1.332e+00,  1.325e+00,  1.298e+00,  1.265e+00,  1.363e+00, &
        1.339e+00,  1.312e+00,  1.294e+00,  1.286e+00,  1.278e+00, &
        1.272e+00,  1.262e+00,  1.255e+00,  1.243e+00,  1.229e+00, &
        1.218e+00,  1.179e+00,  1.153e+00,  1.126e+00,  1.123e+00, &
        1.146e+00,  1.210e+00,  1.258e+00,  1.270e+00,  1.346e+00, &
        1.386e+00,  1.423e+00,  1.443e+00,  1.480e+00,  1.491e+00, &
        1.506e+00,  1.531e+00,  1.549e+00,  1.551e+00,  1.532e+00, &
        1.519e+00/) 
        WATER_IMAG=(/ &
       -3.350e-08, -1.600e-08, -6.500e-09, -1.860e-09, -1.020e-09, & 
       -1.000e-09, -1.960e-09, -1.090e-08, -1.640e-08, -3.350e-08, &
       -1.560e-07, -1.250e-07, -4.860e-07, -2.890e-06, -8.700e-06, &
       -2.000e-04, -1.000e-04, -1.100e-03, -1.740e-03, -2.720e-01, &
       -9.240e-02, -2.314e-02, -9.400e-03, -3.500e-03, -4.600e-03, &
       -1.340e-02, -1.240e-02, -1.160e-02, -1.070e-01, -8.800e-02, &
       -3.920e-02, -3.210e-02, -3.390e-02, -3.510e-02, -3.670e-02, &
       -3.790e-02, -3.990e-02, -4.150e-02, -4.440e-02, -4.790e-02, &
       -5.080e-02, -6.740e-02, -9.680e-02, -1.420e-01, -2.590e-01, &
       -3.050e-01, -3.700e-01, -3.960e-01, -4.020e-01, -4.270e-01, &
       -4.290e-01, -4.260e-01, -4.210e-01, -3.930e-01, -3.790e-01, &
       -3.700e-01, -3.560e-01, -3.390e-01, -3.280e-01, -3.360e-01, &
       -3.850e-01/) 

        !From file waso00, water-soluble sulfates, nitrates, and organic
        WASO_REAL=(/ & 
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.520e+00,  1.520e+00,  1.520e+00,  1.510e+00, &
        1.510e+00,  1.470e+00,  1.420e+00,  1.420e+00,  1.420e+00, &
        1.430e+00,  1.430e+00,  1.450e+00,  1.452e+00,  1.455e+00, &
        1.460e+00,  1.450e+00,  1.440e+00,  1.410e+00,  1.430e+00, &
        1.460e+00,  1.400e+00,  1.200e+00,  1.010e+00,  1.300e+00, &
        2.400e+00,  2.560e+00,  2.200e+00,  1.950e+00,  1.870e+00, &
        1.820e+00,  1.760e+00,  1.720e+00,  1.670e+00,  1.620e+00, &
        1.620e+00,  1.560e+00,  1.440e+00,  1.420e+00,  1.750e+00, &
        2.080e+00,  1.980e+00,  1.850e+00,  2.120e+00,  2.060e+00, &
        2.000e+00,  1.880e+00,  1.840e+00,  1.820e+00,  1.920e+00, &
        1.860e+00/) 
        WASO_IMAG=(/ &
       -3.000e-02, -8.000e-03, -5.000e-03, -5.000e-03, -5.000e-03, &
       -5.000e-03, -6.000e-03, -6.000e-03, -7.000e-03, -7.000e-03, &
       -8.500e-03, -1.000e-02, -1.300e-02, -1.550e-02, -1.900e-02, &
       -2.250e-02, -1.750e-02, -8.000e-03, -1.200e-02, -2.200e-02, &
       -8.000e-03, -7.050e-03, -5.000e-03, -4.000e-03, -5.000e-03, &
       -1.300e-02, -1.200e-02, -1.800e-02, -2.300e-02, -2.700e-02, &
       -3.300e-02, -7.000e-02, -6.500e-02, -1.000e-01, -2.150e-01, &
       -2.900e-01, -3.700e-01, -4.200e-01, -1.600e-01, -9.500e-02, &
       -9.000e-02, -7.000e-02, -5.000e-02, -4.700e-02, -5.300e-02, &
       -5.500e-02, -7.300e-02, -1.000e-01, -2.000e-01, -1.600e-01, &
       -2.420e-01, -1.800e-01, -1.700e-01, -2.200e-01, -2.300e-01, &
       -2.400e-01, -2.800e-01, -2.900e-01, -3.000e-01, -4.000e-01, &
       -5.000e-01/) 

       !From file ssam00, for sea salt
       SEASALT_REAL = (/ &
       1.510e+00,  1.510e+00,  1.510e+00,  1.500e+00,  1.500e+00, &
       1.500e+00,  1.500e+00,  1.490e+00,  1.490e+00,  1.490e+00, &
       1.490e+00,  1.480e+00,  1.480e+00,  1.470e+00,  1.470e+00, &
       1.460e+00,  1.450e+00,  1.450e+00,  1.430e+00,  1.610e+00, &
       1.490e+00,  1.480e+00,  1.480e+00,  1.470e+00,  1.480e+00, &
       1.490e+00,  1.470e+00,  1.420e+00,  1.410e+00,  1.600e+00, &
       1.460e+00,  1.420e+00,  1.400e+00,  1.420e+00,  1.480e+00, &
       1.600e+00,  1.650e+00,  1.610e+00,  1.580e+00,  1.560e+00, &
       1.540e+00,  1.500e+00,  1.480e+00,  1.480e+00,  1.420e+00, &
       1.410e+00,  1.410e+00,  1.430e+00,  1.450e+00,  1.560e+00, &
       1.740e+00,  1.780e+00,  1.770e+00,  1.760e+00,  1.760e+00, &
       1.760e+00,  1.760e+00,  1.770e+00,  1.770e+00,  1.760e+00, &
       1.740e+00/) 
       SEASALT_IMAG = (/ &
      -5.000e-06, -2.000e-06, -3.240e-07, -3.000e-08, -2.430e-08, &
      -1.550e-08, -1.000e-08, -1.600e-08, -4.240e-08, -2.000e-07, &
      -1.080e-06, -1.950e-06, -4.240e-05, -1.410e-04, -3.580e-04, &
      -5.700e-04, -7.620e-04, -1.000e-03, -4.000e-03, -1.000e-02, &
      -3.000e-03, -2.050e-03, -1.600e-03, -1.400e-03, -1.400e-03, &
      -1.400e-03, -2.500e-03, -3.600e-03, -1.100e-02, -2.200e-02, &
      -5.000e-03, -7.000e-03, -1.300e-02, -2.000e-02, -2.600e-02, &
      -3.000e-02, -2.800e-02, -2.620e-02, -1.800e-02, -1.600e-02, &
      -1.500e-02, -1.400e-02, -1.400e-02, -1.400e-02, -1.600e-02, &
      -1.800e-02, -2.300e-02, -3.000e-02, -3.500e-02, -9.000e-02, &
      -1.200e-01, -1.300e-01, -1.350e-01, -1.520e-01, -1.650e-01, &
      -1.800e-01, -2.050e-01, -2.750e-01, -3.000e-01, -5.000e-01, &
      -1.000e+00/) 

       !From file suso00, for acidic (stratospheric) H2SO4 droplets
       SULF_REAL = (/ &
        1.484e+00,  1.469e+00,  1.452e+00,  1.440e+00,  1.432e+00, &
        1.431e+00,  1.430e+00,  1.429e+00,  1.429e+00,  1.428e+00, &
        1.427e+00,  1.426e+00,  1.425e+00,  1.422e+00,  1.413e+00, &
        1.403e+00,  1.394e+00,  1.384e+00,  1.344e+00,  1.293e+00, &
        1.311e+00,  1.350e+00,  1.376e+00,  1.396e+00,  1.385e+00, &
        1.385e+00,  1.360e+00,  1.337e+00,  1.425e+00,  1.424e+00, &
        1.370e+00,  1.210e+00,  1.140e+00,  1.200e+00,  1.370e+00, &
        1.530e+00,  1.650e+00,  1.600e+00,  1.670e+00,  1.910e+00, &
        1.890e+00,  1.720e+00,  1.670e+00,  1.890e+00,  1.740e+00, &
        1.690e+00,  1.640e+00,  1.610e+00,  1.590e+00,  1.520e+00, &
        1.724e+00,  1.950e+00,  1.927e+00,  1.823e+00,  1.780e+00, &
        1.870e+00,  1.930e+00,  1.920e+00,  1.920e+00,  1.900e+00, &
        1.890e+00/) 
       SULF_IMAG = (/ &
       -1.000e-08, -1.000e-08, -1.000e-08, -1.000e-08, -1.000e-08, &
       -1.000e-08, -1.000e-08, -1.470e-08, -1.670e-08, -2.050e-08, &
       -7.170e-08, -8.630e-08, -2.550e-07, -1.530e-06, -6.940e-06, &
       -1.200e-04, -4.160e-04, -1.260e-03, -3.760e-03, -9.550e-02, &
       -1.350e-01, -1.578e-01, -1.580e-01, -1.310e-01, -1.260e-01, &
       -1.200e-01, -1.210e-01, -1.830e-01, -1.950e-01, -1.650e-01, &
       -1.280e-01, -1.760e-01, -4.880e-01, -6.450e-01, -7.550e-01, &
       -7.720e-01, -6.330e-01, -5.860e-01, -7.500e-01, -6.800e-01, &
       -4.550e-01, -3.400e-01, -4.850e-01, -3.740e-01, -1.980e-01, &
       -1.950e-01, -1.950e-01, -2.050e-01, -2.110e-01, -4.140e-01, &
       -5.900e-01, -4.100e-02, -3.025e-02, -2.352e-02, -2.925e-01, &
       -3.150e-02, -2.000e-01, -1.800e-01, -1.800e-01, -1.900e-01, &
       -2.200e-01/) 

       !From file miam00, for mineral dust (Quartz and clay particles)
       DUST_REAL = (/ &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00,  1.530e+00, &
        1.530e+00,  1.530e+00,  1.530e+00,  1.520e+00,  1.520e+00, &
        1.510e+00,  1.510e+00,  1.510e+00,  1.500e+00,  1.500e+00, &
        1.500e+00,  1.480e+00,  1.460e+00,  1.440e+00,  1.430e+00, &
        1.420e+00,  1.460e+00,  1.220e+00,  1.120e+00,  1.060e+00, &
        1.190e+00,  1.850e+00,  2.220e+00,  2.940e+00,  2.910e+00, &
        2.570e+00,  1.910e+00,  1.830e+00,  1.810e+00,  1.740e+00, &
        2.000e+00,  1.630e+00,  1.540e+00,  1.510e+00,  1.470e+00, &
        1.490e+00,  1.770e+00,  2.050e+00,  2.200e+00,  2.390e+00, &
        2.690e+00,  2.990e+00,  2.570e+00,  2.420e+00,  2.420e+00, &
        2.340e+00/) 
       DUST_IMAG = (/ &
       -3.000e-02, -2.500e-02, -1.700e-02, -1.300e-02, -8.500e-03, &
       -7.800e-03, -5.500e-03, -4.500e-03, -4.500e-03, -4.000e-03, &
       -4.000e-03, -4.000e-03, -4.000e-03, -4.000e-03, -5.000e-03, &
       -5.700e-03, -6.400e-03, -7.600e-03, -1.400e-02, -3.900e-02, &
       -2.400e-02, -1.925e-02, -1.800e-02, -1.200e-02, -6.700e-03, &
       -8.700e-03, -1.800e-02, -3.600e-02, -5.500e-02, -6.300e-02, &
       -5.200e-02, -1.300e-01, -8.900e-02, -1.200e-01, -2.100e-01, &
       -2.900e-01, -4.400e-01, -5.400e-01, -6.500e-01, -6.500e-01, &
       -5.000e-01, -2.500e-01, -2.000e-01, -3.500e-01, -5.000e-01, &
       -3.500e-01, -2.200e-01, -2.400e-01, -2.600e-01, -3.200e-01, &
       -3.700e-01, -4.700e-01, -5.700e-01, -8.200e-01, -9.400e-01, &
       -1.000e+00, -8.000e-01, -7.800e-01, -6.700e-01, -6.200e-01, &
       -7.000e-01/) 


	!! BUILD N ARRAY for AqSpecies Volume Concentrations
	ALLOCATE (AqVolConc(HowManyAqChems+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "AqVolConc at size "//	&
	     TRIM(INT2STR(HowManyAqChems+HowManyAqCations+HowManyAqAnions)))

	!! BUILD N ARRAY for Org Species Volume Concentrations
	ALLOCATE (OrgVolConc(HowManyOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate OrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyOrgChems)))

	!! BUILD N ARRAY for Aq Org Species Volume Concentrations
	ALLOCATE (AqOrgVolConc(HowManyAqOrgChems), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate AqOrgVolConc at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))

	!! BUILD N ARRAY for H2O and Ion Mol Fractions
	ALLOCATE (MolFrac(1+HowManyAqCations+HowManyAqAnions), STAT = I)
	IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate MolFrac at size "	&
						  //TRIM(INT2STR(HowManyAqOrgChems)))
	
	IF (InParticle%numberofparticles .GT. 1.0e-6) THEN
		SingleParticleVolume = (4.0*Pi/3.0)*(InParticle%effectiveradius)**3 !cm3
		
		VolumeConc = SingleParticleVolume*InParticle%numberofparticles !cm3/cm3
		
		!Convert aqueous concentrations to cm3 solution/cm3 air
		
		!Water
		AqVolConc(1)  = InParticle%AqChems(1) * &
                                InParticle%NumberOfParticles*AqMolecularMass(1)/ &
                                InParticle%SolutionDensity
		SolutionVolConc = AqVolConc(1)
		TotalVolConc = AqVolConc(1)
				
		!Solid Salts
		DO I = 2, HowManyAqChems
			AqVolConc(I)  = InParticle%AqChems(I) * &
                                        InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                        SolidSaltDensity(I)
			TotalVolConc = TotalVolConc + AqVolConc(I)
		END DO

		!Inorganic Ions
		DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
			AqVolConc(I)  = InParticle%AqChems(I) * &
                                        InParticle%NumberOfParticles*AqMolecularMass(I)/ &
                                        InParticle%SolutionDensity
			TotalVolConc = TotalVolConc + AqVolConc(I)
			SolutionVolConc = SolutionVolConc + AqVolConc(I)	
		END DO

		!Convert org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyOrgChems
			OrgVolConc(I)  = InParticle%OrgChems(I) * &
                                         InParticle%NumberOfParticles*OrgMolecularMass(I)/ &
                                         OrgDensity(I)
			TotalVolConc = TotalVolConc + OrgVolConc(I)	
		END DO
	
		!Convert hydrophilic org concentrations to cm3 solution/cm3 air
		DO I = 1, HowManyAqOrgChems
			AqOrgVolConc(I)  = InParticle%AqOrgChems(I) * &
                                           InParticle%NumberOfParticles*AqOrgMolecularMass(I) &
                                           /AqOrgDensity(I)
			TotalVolConc = TotalVolConc + AqOrgVolConc(I)	
		END DO

		!Calculate shell volume conc (Total - BC)
		ShellVolumeConc = VolumeConc - OrgVolConc(1)
			
	ELSE    !No particles, use dummy numbers
		InParticle%ExtEff = 1.0
		InParticle%ScaEff = 1.0
		InParticle%BackScaEff = 1.0
		InParticle%AssymParam = 1.0
                RETURN          
	END IF

        
        !Calculate shell and core refractive index for each wavelength bin
        FASTTUV = .TRUE.
        IF(FASTTUV) THEN !Use center wavelengths for 7 tropospheric fastTUV bins
           numbin = 7
	   ALLOCATE (wvln_array(numbin), STAT = I)
	   IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_array at size "//TRIM(INT2STR(numbin)))
	   ALLOCATE (wvln_rounddown(numbin), STAT = I)
	   IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_rounddown at size "//TRIM(INT2STR(numbin)))
           wvln_array(1) = 297.7  !in nm 
           wvln_rounddown(1) = 1  
           wvln_array(2) = 309.5       
           wvln_rounddown(2) = 2 
           wvln_array(3) = 325.5       
           wvln_rounddown(3) = 2 
           wvln_array(4) = 378.8      
           wvln_rounddown(4) = 3 
           wvln_array(5) = 447.5    
           wvln_rounddown(5) = 4 
           wvln_array(6) = 602.0
           wvln_rounddown(6) = 8 
           wvln_array(7) = 736.3 
           wvln_rounddown(7) = 10 
        ELSE !For photolysis between 250 and 700 nm, every 1 nm
           numbin = 451
	   ALLOCATE (wvln_array(numbin), STAT = I)
	   IF (I > 0) CALL ERROR("In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_array at size "//TRIM(INT2STR(numbin)))
           DO ib = 1, numbin 
              wvln_array(ib) = (250+(ib-1)*1) !in nm
           ENDDO
        ENDIF
 
        DO ib = 1, numbin

           Wavelength = wvln_array(ib) !in nm

           !Index of point on array WAVE closest to, but less than,
           !the wavelength for interpolation
           IF(FASTTUV) THEN
              w_ind_low = wvln_rounddown(ib)
           ELSE
              w_ind_low = INT(ib/50)+1
           ENDIF

           !Calculate weightings of wavelength-dependent refractive index points
           a_w = (WAVE(w_ind_low+1)-Wavelength) &
                  /(WAVE(w_ind_low+1)-WAVE(w_ind_low))
           b_w = 1-a_w

           !Calculate Core Refractive Index
           n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
           n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
	   SolarCoreInd = CMPLX(n_real, -1.0*n_imag)
           !IF (ib .eq. 419 .or. ib .eq. 420) WRITE(*,*) Wavelength, SolarShellInd, SolarCoreInd

           !We assume the core is BC, the shell everything else
           !BUT we need to handle the cases where the particle is all BC (no shell)
           ! or has no BC (all core). All core handled here.
           IF (ShellVolumeConc .LE. 0 .OR. InParticle%Dry) THEN
                !Shell and core refractive indices equal.
                SolarShellInd = SolarCoreInd
                !Calculate Core radius
	        InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
	        InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
                !Set Shell radius to above, and core a little less than shell
                InParticle%EffectiveRadius = InParticle%AbsCoreRad
                InParticle%AbsCoreRad = 0.99*InParticle%AbsCoreRad
           ELSE !Some or all shell, all shell handled at bottom
	     !Calculate shell refractive index
	     !Use the volume average dielectric constant mixing rule
	     !See Jacobson, "Fund. of Atm. Modeling", 2nd ed, p. 310
 	     SumEr = 0.
	     SumEi = 0.
 
             !Water refractive index and molal Refraction
             n_real = a_w*WATER_REAL(w_ind_low)+b_w*WATER_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*WATER_IMAG(w_ind_low+1))
             V_molal = 18.0674 !Tang, 1997: density = 0.9971 g/cm3, MW = 18.015 g/mol, V = MW/dens
	     R_molal = V_molal*(n_real**2 - 1.0)/(n_real**2 + 2) 

             !Calculate average refractive index for ionic solution
	     !Following Tang et. al., JGR V.102, D19, pp. 23269-23275, 1997. 
             !and Tang, JGR V.102, D2, pp. 1883-1893, 1997.)
	     !Step 1: Calculate mole fractions in solution
	     Sum = InParticle%AqChems(1)
	     DO I = HowManyAqChems+1, HowManyAqChems+HowManyAqCations+HowManyAqAnions
	   	Sum = Sum + InParticle%AqChems(I)
	     END DO
	     MolFrac(1) = InParticle%AqChems(1)/Sum
	     DO I = 2, 1+HowManyAqCations+HowManyAqAnions
	   	MolFrac(I) = InParticle%AqChems(HowManyAqChems+I-1)/Sum
	     END DO
		
	     !Step 2: Calculate average molecular weight	
	     SolnMolWeight = MolFrac(1)*AqMolecularMass(1)
	     DO I = 2, 1+HowManyAqCations+HowManyAqAnions
	   	SolnMolWeight = SolnMolWeight + MolFrac(I)*AqMolecularMass(HowManyAqChems+I-1)
	     END DO
		
	     !Step 3: Calculate average molal refraction
	     MolalRefraction = MolFrac(1)*R_molal !See Tang et. al, 1997
	     DO I = 2, 1+HowManyAqCations
                MolalRefraction = MolalRefraction + MolFrac(I)*CationIonicRefraction(I-1)
                !WRITE(*,*) MolFrac(I), CationIonicRefraction(I-1)
	     END DO
	     DO I = 2+HowManyAqCations, 1+HowManyAqCations+HowManyAqAnions
	   	MolalRefraction = MolalRefraction + MolFrac(I)*AnionIonicRefraction(I-HowManyAqCations-1)
                !WRITE(*,*) MolFrac(I), AnionIonicRefraction(I-HowManyAqCations-1)
	     END DO

	     !Step 4: Calculate solution molal volume
	     MolVol = SolnMolWeight/InParticle%InorgSolnDensity

	     !Step 5: Calculate Refractive Index
	     xxx = MolalRefraction/MolVol
	     SolutionRefracIndex(1) = SQRT((2.*xxx+1.)/(1.-xxx))
	     SolutionRefracIndex(2) = n_imag !Imaginary index for water
	  
	     !Step 6: Calculate Er and Ei
	     Er = SolutionRefracIndex(1)**2 - SolutionRefracIndex(2)**2
	     Ei = 2*SolutionRefracIndex(1)*SolutionRefracIndex(2)
             SumEr = SumEr + (SolutionVolConc/ShellVolumeConc)*Er
	     SumEi = SumEi + (SolutionVolConc/ShellVolumeConc)*Ei

             !Solid Salts
             DO I = 2, HowManyAqChems
                IF(RefIndFlag(I) .EQ. 0) THEN
                   n_real = a_w*WATER_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*WATER_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 1) THEN
                   n_real = a_w*SULF_REAL(w_ind_low)+b_w*SULF_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*SULF_IMAG(w_ind_low)+b_w*SULF_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 2) THEN
                   n_real = a_w*SEASALT_REAL(w_ind_low)+b_w*SEASALT_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*SEASALT_IMAG(w_ind_low)+b_w*SEASALT_IMAG(w_ind_low+1))
                ELSE IF(RefIndFlag(I) .EQ. 3) THEN
                   n_real = a_w*DUST_REAL(w_ind_low)+b_w*DUST_REAL(w_ind_low+1)
                   n_imag = -1.0*(a_w*DUST_IMAG(w_ind_low)+b_w*DUST_IMAG(w_ind_low+1))
                ENDIF
	        Er = n_real**2 - n_imag**2
	        Ei = 2*n_real*n_imag
                SumEr = SumEr + (AqVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqVolConc(I)/ShellVolumeConc)*Ei
	     END DO
         
             !Insoluble organics, skipping BC core
	     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	     Er = n_real**2 - n_imag**2
	     Ei = 2*n_real*n_imag
             DO I = 2, HowManyOrgChems
		SumEr = SumEr + (OrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (OrgVolConc(I)/ShellVolumeConc)*Ei
	     END DO

             !Hydrophilic Organics
	     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
             n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
	     Er = n_real**2 - n_imag**2
	     Ei = 2*n_real*n_imag
             DO I = 1, HowManyAqOrgChems
		SumEr = SumEr + (AqOrgVolConc(I)/ShellVolumeConc)*Er
		SumEi = SumEi + (AqOrgVolConc(I)/ShellVolumeConc)*Ei
	     END DO

            !!For core-in-shell mixing rule
	    Factor = SQRT(SumEr**2 + SumEi**2)
	    InParticle%ShellRealRefracInd(ib) = SQRT((Factor + SumEr)/2.0)
	    InParticle%ShellImagRefracInd(ib) = SQRT((Factor - SumEr)/2.0)
	    SolarShellInd = CMPLX(InParticle%ShellRealRefracInd(ib),-1.0*InParticle%ShellImagRefracInd(ib))

            !We assume the core is BC, the shell everything else
            !BUT we need to handle the cases where the particle is all BC (all core, no shell)
            ! or has no BC (all shell, no core). No core, all shell handeled here.

            IF (OrgVolConc(1) .LE. 0.0) THEN !No core
               !Set core and shell refractive indices equal
               SolarCoreInd = SolarShellInd
               !Set small core radius
               InParticle%AbsCoreRad = 0.01*InParticle%EffectiveRadius
            ELSE !Both Shell and Core
	       !Calculate Core radius
	       InParticle%AbsCoreRad = OrgVolConc(1)/InParticle%NumberofParticles
	       InParticle%AbsCoreRad = (0.75*InParticle%AbsCoreRad/Pi)**(1.0/3.0) !In cm
	       !This shouldn't happen, but be safe
               IF(InParticle%AbsCoreRad .GT. 0.9999*InParticle%EffectiveRadius) & 
                InParticle%AbsCoreRad = 0.9999*InParticle%EffectiveRadius
            ENDIF

            !!For volume average dielectric constant mixing rule
!	     n_real = a_w*SOOT_REAL(w_ind_low)+b_w*SOOT_REAL(w_ind_low+1)
!            n_imag = -1.0*(a_w*SOOT_IMAG(w_ind_low)+b_w*SOOT_IMAG(w_ind_low+1))
!            ErBC = n_real**2 - n_imag**2
!            EiBC = 2*n_real*n_imag
!            BCVolRat = (InParticle%AbsCoreRad/InParticle%EffectiveRadius)**3.0
!            ErVA = ErBC*BCVolRat+SumEr*(1.0-BCVolRat)
!            EiVA = EiBC*BCVolRat+SumEi*(1.0-BCVolRat)
!            Factor = SQRT(ErVA**2 + EiVA**2)
!	    InParticle%VARealRefracInd(ib) = SQRT((Factor + ErVA)/2.0)
!	    InParticle%VAImagRefracInd(ib) = SQRT((Factor - ErVA)/2.0)            
!            VAInd = CMPLX(InParticle%VARealRefracInd(ib),-1.0*InParticle%VAImagRefracInd(ib))

            !!For Maxwell-Garnett mixing rule
!            Eshell =  CMPLX(SumEr,-1.0*SumEi)
!            EBC =  CMPLX(ErBC,-1.0*EiBC)
!            Factor2 = (EBC-Eshell)/(EBC+2.0*Eshell)
!            EMG = Eshell*(1.0+(3.0*BCVolRat*Factor2/(1.0-BCVolRat*Factor2)))
!            ErMG = REAL(EMG)
!            EiMG = AIMAG(EMG)
!            Factor = SQRT(ErMG**2 + EiMG**2)
!	     InParticle%MGRealRefracInd(ib) = SQRT((Factor + ErMG)/2.0)
!	     InParticle%MGImagRefracInd(ib) = SQRT((Factor - ErMG)/2.0)            
!            MGInd = CMPLX(InParticle%MGRealRefracInd(ib),-1.0*InParticle%MGImagRefracInd(ib))
          ENDIF

	  NumAng = 0.

          !DMiLay needs 2*Pi divided by the wavelength
          Wavelength = 2*PI/((Wavelength)*nm)

!    DMiLay Documentation				
!    I N P U T   A R G U M E N T S

!    (Definition:  size parameter = sphere circumference / wavelength )

!      Rshell      radius of shell

!      Rcore       radius of core

!      WVNO        2*pi / wavelength

!      RindSh      COMPLEX refractive index of shell (negative
!                     imaginary part)

!      RindCo      COMPLEX refractive index of core (negative
!                     imaginary part)

!      MU          array of cosines of scattering angles (angles between
!                     directions of incident and scattered radiation).
!                     For angles between 90 and 180 degrees, use the
!                     supplement (180-angle) of the angle instead, so
!                     0.le.MU.le.1 (see comments below on M1,M2,21,D21)

!      NumAng      Number of scattering angles for which computations
!                     are required; should not exceed MaxAng
!                     (NOTE:  NumAng=0 will suppress the calculation
!                      of the scattering matrix quantitities  M1, M2,
!                      S21, D21 and save a lot of computer time)

!      MaxAng      First dimension of M1,M2,21,D21 in calling program



!    O U T P U T   A R G U M E N T S

!      (Definitions for these arguments can be found in H.C. van de
!       Hulst, Light Scattering By Small Particles, Dover Press, New
!       York, 1981 (reprint of 1957 edition); abbreviated VDH below)

!      QEXT     Efficiency factor for extinction (VDH Sec 9.32)
!               (same as corresponding quantity in MIEV)

!      Qsca     Efficiency factor for scattering (VDH Sec 9.32)
!               (same as corresponding quantity in MIEV)

!      GQS!     average(cosine theta) * Qsca (VDH Sec 9.32)
!                  (<cos theta> is usually denoted by g, hence
!                   the name of the variable)
!               (same as corresponding quantity in MIEV)

!      QBS      Backscatter cross section.
!               ( Re(SBACK)**2 + Im(SBACK)**2 ) / (Pi*XSHELL**2)
!               where the corresponding quantity from MIEV is

!               SBACK = 0.5*sum(n=1 to inf)((-1)**(n+1)(2n+1)(an-bn))

!               and an,bn are ACOE,BCOE below.

!      M1(j,k)  Element M1 of scattering matrix F' (VDH Sec 5.14);
!                  M1(j,1) refers to angle with cosine MU(j); 
!                  M1(j,2) refers to supplement of that angle.
!               (Be sure to type REAL in calling program.)

!      M2(j,k)  Element M2 of scattering matrix F' (VDH Sec 5.14);
!                  M2(j,1) refers to angle with cosine MU(j); 
!                  M2(j,2) refers to supplement of that angle.
!               (Be sure to type REAL in calling program.)

!     S21(j,k)  Element S21 of scattering matrix F' (VDH Sec 5.14);
!                  S21(j,1) refers to angle with cosine MU(j); 
!                  S21(j,2) refers to supplement of that angle.

!     D21(j,k)  Element D21 of scattering matrix F' (VDH Sec 5.14);
!                  D21(j,1) refers to angle with cosine MU(j); 
!                  D21(j,2) refers to supplement of that angle.

          !MM DEBUG - CosAng, M1, M2, S21, and D21 are all zero arrays
          !WRITE(*,*) "InParticle%AbsCoreRad: ", InParticle%AbsCoreRad
          !WRITE(*,*) "InParticle%EffectiveRadius: ", InParticle%EffectiveRadius
          !WRITE(*,*) "Wavelength: ", Wavelength
          !WRITE(*,*) "SolarShellInd: ", SolarShellInd
          !WRITE(*,*) "SolarCoreInd: ", SolarCoreInd
          !WRITE(*,*) "CosAng: ", CosAng
          !WRITE(*,*) "NumAng: ", NumAng 
          !WRITE(*,*) "InParticle%ExtEff(ib): ", InParticle%ExtEff(ib)
          !WRITE(*,*) "InParticle%ScaEff(ib): ", InParticle%ScaEff(ib)
          !WRITE(*,*) "InParticle%BackScaEff(ib): ", InParticle%BackScaEff(ib)
          !WRITE(*,*) "InParticle%AssymParam(ib): ", InParticle%AssymParam(ib)
          !WRITE(*,*) "M1: ", M1
          !WRITE(*,*) "M2: ", M2			
          !WRITE(*,*) "S21: ", S21
          !WRITE(*,*) "D21: ", D21  
          !WRITE(*,*) "MAXANG: ", MAXANG
!          WRITE(*,*) "Call 1"


!MM-DEBUG	  CALL DMiLay( InParticle%AbsCoreRad, &
!MM-DEBUG			InParticle%EffectiveRadius, &
!MM-DEBUG			Wavelength, &
!MM-DEBUG			SolarShellInd, & 
!MM-DEBUG			SolarCoreInd, &
!MM-DEBUG			CosAng, &
!MM-DEBUG			NumAng, &
!MM-DEBUG			InParticle%ExtEff(ib), &
!MM-DEBUG			InParticle%ScaEff(ib), &
!MM-DEBUG			InParticle%BackScaEff(ib), &
!MM-DEBUG			InParticle%AssymParam(ib), & 
!MM-DEBUG                        M1, M2, S21, D21, MAXANG )

          !MM-DEBUG: Replace DMiLay call with if statements to assign values
          !  to InParticle%ExtEff(ib),InParticle%ScaEff(ib),
          !  InParticle%BackScaEff(ib),InParticle%AssymParam(ib) 
          !  for each effective radius size bin (1-10) and wavelength (ib)

          !identify size bin for effective radius (cm)
          IF ((InParticle%EffectiveRadius .GT. 0.) .AND. &
              (InParticle%EffectiveRadius .LE. 2.5e-06)) THEN
            sizebin = 1
          ELSE IF ((InParticle%EffectiveRadius .GT. 2.5e-06) .AND. &
                   (InParticle%EffectiveRadius .LE. 4.0e-06)) THEN
            sizebin = 2
          ELSE IF ((InParticle%EffectiveRadius .GT. 4.0e-06) .AND. &
                   (InParticle%EffectiveRadius .LE. 6.3e-06)) THEN
            sizebin = 3
          ELSE IF ((InParticle%EffectiveRadius .GT. 6.3e-06) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.0e-05)) THEN
            sizebin = 4
          ELSE IF ((InParticle%EffectiveRadius .GT. 1.0e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.59e-05)) THEN
            sizebin = 5
          ELSE IF ((InParticle%EffectiveRadius .GT. 1.59e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 2.52e-05)) THEN
            sizebin = 6
          ELSE IF ((InParticle%EffectiveRadius .GT. 2.52e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 4.01e-05)) THEN
            sizebin = 7
          ELSE IF ((InParticle%EffectiveRadius .GT. 4.01e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 6.36e-05)) THEN
            sizebin = 8
          ELSE IF ((InParticle%EffectiveRadius .GT. 6.36e-05) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.01e-04)) THEN
            sizebin = 9
          ELSE IF ((InParticle%EffectiveRadius .GT. 1.01e-04) .AND. &
                   (InParticle%EffectiveRadius .LE. 1.0e5)) THEN
            sizebin = 10
          ELSE
            WRITE(*,*) "InParticle%EffectiveRadius = ",InParticle%EffectiveRadius
            call error("Unknown InParticle%EffectiveRadius range. STOP")
          END IF !InParticle%EffectiveRadius

	  ! DMiLay output lookup table by wavelength bin (ib) and 
	  ! effective radius size (sizebin)
          IF (ib .eq. 1) THEN                           !MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 2.121888587382035E-002
              InParticle%ScaEff(ib) = 1.307131658929511E-003
              InParticle%BackScaEff(ib) = 1.510241237540259E-004
              InParticle%AssymParam(ib) = 1.779476013214267E-005
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 0.127092654496426 
              InParticle%ScaEff(ib) = 6.046790202959535E-002
              InParticle%BackScaEff(ib) = 5.750267311267363E-003
              InParticle%AssymParam(ib) = 5.463468568745519E-003
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 0.473383790453063  
              InParticle%ScaEff(ib) = 0.333555454164078  
              InParticle%BackScaEff(ib) = 2.045947451490424E-002 
              InParticle%AssymParam(ib) = 7.862910765605213E-002
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 1.52474778252003  
              InParticle%ScaEff(ib) = 1.23767297048067 
              InParticle%BackScaEff(ib) = 8.303324626201600E-003
              InParticle%AssymParam(ib) = 0.746388408646381
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 3.16206659412624  
              InParticle%ScaEff(ib) = 2.71544421419550
              InParticle%BackScaEff(ib) = 3.200642357740086E-003
              InParticle%AssymParam(ib) = 1.86651499266298 
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 3.89127544932093  
              InParticle%ScaEff(ib) = 3.33497392058575 
              InParticle%BackScaEff(ib) = 8.155620690029193E-002
              InParticle%AssymParam(ib) = 2.36045846346754
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 2.36117476253389  
              InParticle%ScaEff(ib) = 1.75412851405968 
              InParticle%BackScaEff(ib) = 0.162787964869075
              InParticle%AssymParam(ib) = 1.01491769113054
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 2.64765819872286  
              InParticle%ScaEff(ib) = 1.92118459877779
              InParticle%BackScaEff(ib) = 0.320705458971901 
              InParticle%AssymParam(ib) = 1.52601186529154 
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 2.38669121295974 
              InParticle%ScaEff(ib) = 1.63674021138855 
              InParticle%BackScaEff(ib) = 9.046238941112521E-002
              InParticle%AssymParam(ib) = 1.30768747111401 
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.08247023748521 
              InParticle%ScaEff(ib) = 1.23183844905796
              InParticle%BackScaEff(ib) = 9.526268691288320E-004
              InParticle%AssymParam(ib) = 1.09597855140376 
            END IF !sizebin
 
          ELSE IF (ib .eq. 2) THEN                      !MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 1.940690397997060E-002
              InParticle%ScaEff(ib) = 1.119596514108643E-003
              InParticle%BackScaEff(ib) = 1.296749968046225E-004
              InParticle%AssymParam(ib) = 1.410083652156505E-005
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 0.112355100291307  
              InParticle%ScaEff(ib) = 5.192222228372283E-002
              InParticle%BackScaEff(ib) = 5.027660705562706E-003
              InParticle%AssymParam(ib) = 4.338486686936421E-003
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 0.421230415298952    
              InParticle%ScaEff(ib) = 0.294306500740802    
              InParticle%BackScaEff(ib) = 1.930914577363761E-002 
              InParticle%AssymParam(ib) = 6.350884426115723E-002
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 1.36699989379675  
              InParticle%ScaEff(ib) = 1.10586770021349 
              InParticle%BackScaEff(ib) = 8.146409454995982E-003
              InParticle%AssymParam(ib) = 0.643436122529757
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 3.02023804553219  
              InParticle%ScaEff(ib) = 2.59512177245481    
              InParticle%BackScaEff(ib) = 4.055106300621420E-003
              InParticle%AssymParam(ib) = 1.78221936537171 
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 3.89142373624037  
              InParticle%ScaEff(ib) = 3.36018385209317 
              InParticle%BackScaEff(ib) = 9.988278102159109E-002
              InParticle%AssymParam(ib) = 2.31161253263504 
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 2.67785914545257 
              InParticle%ScaEff(ib) = 2.06887595489080  
              InParticle%BackScaEff(ib) = 0.217298977443273
              InParticle%AssymParam(ib) = 1.20750055579509
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 2.65962852938415  
              InParticle%ScaEff(ib) = 1.96753057020663 
              InParticle%BackScaEff(ib) = 0.265791124410301
              InParticle%AssymParam(ib) = 1.44433188415122 
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 2.23957579055603  
              InParticle%ScaEff(ib) = 1.54243756427382 
              InParticle%BackScaEff(ib) = 8.327918065614245E-002
              InParticle%AssymParam(ib) = 1.19250234619642 
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.20230582462907 
              InParticle%ScaEff(ib) = 1.39733785469035 
              InParticle%BackScaEff(ib) = 1.908905254111742E-002
              InParticle%AssymParam(ib) = 1.23276210379463 
            END IF !sizebin

          ELSE IF (ib .eq. 3) THEN                      !MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 1.772826081778498E-002
              InParticle%ScaEff(ib) = 9.149222913046922E-004
              InParticle%BackScaEff(ib) = 1.062767186628380E-004
              InParticle%AssymParam(ib) = 1.042195272961068E-005
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 9.708241778683968E-002
              InParticle%ScaEff(ib) = 4.250553139785036E-002
              InParticle%BackScaEff(ib) = 4.203119070808736E-003
              InParticle%AssymParam(ib) = 3.212686799219858E-003
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 0.361837102343402   
              InParticle%ScaEff(ib) = 0.247862529058284  
              InParticle%BackScaEff(ib) = 1.749529294566652E-002 
              InParticle%AssymParam(ib) = 4.788040590044889E-002
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 1.18880523239633 
              InParticle%ScaEff(ib) = 0.953950114092125  
              InParticle%BackScaEff(ib) = 1.022696413910598E-002
              InParticle%AssymParam(ib) = 0.517209586039210 
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 2.74276938686630   
              InParticle%ScaEff(ib) = 2.34731786228653 
              InParticle%BackScaEff(ib) = 1.618874002603226E-002
              InParticle%AssymParam(ib) = 1.58105813178798
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 3.58937007227166 
              InParticle%ScaEff(ib) = 3.09855430248188   
              InParticle%BackScaEff(ib) = 4.939222249234157E-002
              InParticle%AssymParam(ib) = 2.08425605392468  
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 2.64615303303594  
              InParticle%ScaEff(ib) = 2.07417416079047  
              InParticle%BackScaEff(ib) = 0.213583705475907  
              InParticle%AssymParam(ib) = 1.16784221226047
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 2.31209408959056     
              InParticle%ScaEff(ib) = 1.70019209175438
              InParticle%BackScaEff(ib) = 0.193521213831159
              InParticle%AssymParam(ib) = 1.25147602437463   
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 2.12336999670260  
              InParticle%ScaEff(ib) = 1.44724291893440
              InParticle%BackScaEff(ib) = 0.143231595441244  
              InParticle%AssymParam(ib) = 1.06982952415901  
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.25793049544627 
              InParticle%ScaEff(ib) = 1.49397708874293 
              InParticle%BackScaEff(ib) = 1.395713642345309E-002
              InParticle%AssymParam(ib) = 1.31041249288240  
            END IF !sizebin

          ELSE IF (ib .eq. 4) THEN                      !MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 1.406291774772995E-002
              InParticle%ScaEff(ib) = 4.982307925897343E-004
              InParticle%BackScaEff(ib) = 5.828866856375704E-005 
              InParticle%AssymParam(ib) = 4.195622273115239E-006
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 6.505333999696981E-002
              InParticle%ScaEff(ib) = 2.317105529632476E-002
              InParticle%BackScaEff(ib) = 2.410163553939978E-003
              InParticle%AssymParam(ib) = 1.297253478060208E-003
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 0.227092530278676   
              InParticle%ScaEff(ib) = 0.142558005372046  
              InParticle%BackScaEff(ib) = 1.180038236613396E-002 
              InParticle%AssymParam(ib) = 2.001093373092828E-002
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 0.813875696970988  
              InParticle%ScaEff(ib) = 0.635700927600905 
              InParticle%BackScaEff(ib) = 2.019737939263412E-002
              InParticle%AssymParam(ib) = 0.249284427095052    
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 2.24677626027416  
              InParticle%ScaEff(ib) = 1.92008521489047  
              InParticle%BackScaEff(ib) = 3.448272124126465E-002
              InParticle%AssymParam(ib) = 1.18475642982612  
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 3.71768012114766   
              InParticle%ScaEff(ib) = 3.26697313863928  
              InParticle%BackScaEff(ib) = 4.211515392068865E-002
              InParticle%AssymParam(ib) = 2.26474184436812 
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 3.15395762049337   
              InParticle%ScaEff(ib) = 2.62845346689653  
              InParticle%BackScaEff(ib) = 9.069580730675429E-002
              InParticle%AssymParam(ib) = 1.70642861721962   
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 1.97546138502815    
              InParticle%ScaEff(ib) = 1.37829533244958 
              InParticle%BackScaEff(ib) = 0.230135277343040   
              InParticle%AssymParam(ib) = 0.732971276835079      
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 2.38073531724947  
              InParticle%ScaEff(ib) = 1.78058641865853  
              InParticle%BackScaEff(ib) = 6.271555314571442E-002
              InParticle%AssymParam(ib) = 1.37073149687549  
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.30603324113733
              InParticle%ScaEff(ib) = 1.61070243113739 
              InParticle%BackScaEff(ib) = 4.376417436252337E-002
              InParticle%AssymParam(ib) = 1.38024839802793  
            END IF !sizebin

          ELSE IF (ib .eq. 5) THEN                      !MM-DEBUG
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 1.148629920832881E-002
              InParticle%ScaEff(ib) = 2.553728806682175E-004
              InParticle%BackScaEff(ib) = 3.004715076463169E-005 
              InParticle%AssymParam(ib) = 1.542918715238243E-006
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 4.498505509724061E-002
              InParticle%ScaEff(ib) = 1.185004269625059E-002
              InParticle%BackScaEff(ib) = 1.282346107090985E-003
              InParticle%AssymParam(ib) = 4.775020967319521E-004
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 0.138416059666308   
              InParticle%ScaEff(ib) = 7.474364298821812E-002  
              InParticle%BackScaEff(ib) = 6.921461998471951E-003 
              InParticle%AssymParam(ib) = 7.491535814624178E-003
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 0.532084652621101  
              InParticle%ScaEff(ib) = 0.396180274607873  
              InParticle%BackScaEff(ib) = 2.198234421343684E-002
              InParticle%AssymParam(ib) = 0.104610322722200 
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 1.76004348814124 
              InParticle%ScaEff(ib) = 1.48017082844582   
              InParticle%BackScaEff(ib) = 1.297693602275804E-002
              InParticle%AssymParam(ib) = 0.910226004729049 
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 3.23927449303276  
              InParticle%ScaEff(ib) = 2.83434022402345  
              InParticle%BackScaEff(ib) = 1.655436664797596E-002
              InParticle%AssymParam(ib) = 1.90715576950648 
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 3.55850399641860    
              InParticle%ScaEff(ib) = 3.05895262076792  
              InParticle%BackScaEff(ib) = 5.291343608165601E-002
              InParticle%AssymParam(ib) = 2.09188528043502  
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 2.37339306979134      
              InParticle%ScaEff(ib) = 1.79231039113015  
              InParticle%BackScaEff(ib) = 0.546820892159810   
              InParticle%AssymParam(ib) = 0.901456983562644  
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 2.71943897417957   
              InParticle%ScaEff(ib) = 2.11433190822174    
              InParticle%BackScaEff(ib) = 0.153876131014613  
              InParticle%AssymParam(ib) = 1.68974923582332   
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.12879044375487 
              InParticle%ScaEff(ib) = 1.45357322190029 
              InParticle%BackScaEff(ib) = 7.306901748121303E-002
              InParticle%AssymParam(ib) = 1.16280268344626 
            END IF !sizebin

          ELSE IF (ib .eq. 6) THEN                      !MM-DEBUG    
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 8.268478254246958E-003
              InParticle%ScaEff(ib) = 7.779873607112582E-005
              InParticle%BackScaEff(ib) = 9.212939136235002E-006
              InParticle%AssymParam(ib) = 2.602070967102499E-007
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 2.645960926807461E-002
              InParticle%ScaEff(ib) = 3.591220299439265E-003
              InParticle%BackScaEff(ib) = 4.061503267999171E-004
              InParticle%AssymParam(ib) = 8.051685759276727E-005
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 6.361524696790417E-002  
              InParticle%ScaEff(ib) = 2.287929401020232E-002  
              InParticle%BackScaEff(ib) = 2.381474586401371E-003 
              InParticle%AssymParam(ib) = 1.274677295265301E-003
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 0.221955542119892   
              InParticle%ScaEff(ib) = 0.140205520571554     
              InParticle%BackScaEff(ib) = 1.164064737916311E-002
              InParticle%AssymParam(ib) = 1.953285310870732E-002
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 0.810591645799964   
              InParticle%ScaEff(ib) = 0.637047368758878    
              InParticle%BackScaEff(ib) = 2.014982220258068E-002
              InParticle%AssymParam(ib) = 0.250065646556210 
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 2.25124420745427    
              InParticle%ScaEff(ib) = 1.93013077879849    
              InParticle%BackScaEff(ib) = 3.369593126279867E-002
              InParticle%AssymParam(ib) = 1.19314667047833  
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 3.72531673802397  
              InParticle%ScaEff(ib) = 3.27184340778492   
              InParticle%BackScaEff(ib) = 4.158975971367909E-002
              InParticle%AssymParam(ib) = 2.27687227232532 
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 3.14071770978196  
              InParticle%ScaEff(ib) = 2.60574296588879 
              InParticle%BackScaEff(ib) = 9.050427764341874E-002
              InParticle%AssymParam(ib) = 1.69552287739499     
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 1.98219804722697  
              InParticle%ScaEff(ib) = 1.36568272205137
              InParticle%BackScaEff(ib) = 0.205675787979628   
              InParticle%AssymParam(ib) = 0.751780576328922   
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.19009104902301 
              InParticle%ScaEff(ib) = 1.51067313462609 
              InParticle%BackScaEff(ib) = 0.184001824287019   
              InParticle%AssymParam(ib) = 1.20377046586723  
            END IF !sizebin

          ELSE IF (ib .eq. 7) THEN                      !MM-DEBUG   
            IF (sizebin .EQ. 1) THEN
              InParticle%ExtEff(ib) = 7.036217015339360E-003
              InParticle%ScaEff(ib) = 3.474724335355355E-005
              InParticle%BackScaEff(ib) = 4.125645003152539E-006
              InParticle%AssymParam(ib) = 7.773456560967481E-008
            ELSE IF (sizebin .EQ. 2) THEN
              InParticle%ExtEff(ib) = 2.068954971605580E-002
              InParticle%ScaEff(ib) = 1.599265224338477E-003
              InParticle%BackScaEff(ib) = 1.841438916112675E-004
              InParticle%AssymParam(ib) = 2.404430602625647E-005
            ELSE IF (sizebin .EQ. 3) THEN
              InParticle%ExtEff(ib) = 4.295098702587454E-002   
              InParticle%ScaEff(ib) = 1.017877671981796E-002  
              InParticle%BackScaEff(ib) = 1.109366840617094E-003 
              InParticle%AssymParam(ib) = 3.810626939130117E-004
            ELSE IF (sizebin .EQ. 4) THEN
              InParticle%ExtEff(ib) = 0.125887270461626   
              InParticle%ScaEff(ib) = 6.403621621119997E-002 
              InParticle%BackScaEff(ib) = 6.049201522663256E-003
              InParticle%AssymParam(ib) = 5.945281814103842E-003 
            ELSE IF (sizebin .EQ. 5) THEN
              InParticle%ExtEff(ib) = 0.485885337588970  
              InParticle%ScaEff(ib) = 0.354584881948764    
              InParticle%BackScaEff(ib) = 2.110126901290764E-002
              InParticle%AssymParam(ib) = 8.665704394697525E-002
            ELSE IF (sizebin .EQ. 6) THEN
              InParticle%ExtEff(ib) = 1.60945192916300 
              InParticle%ScaEff(ib) = 1.33482386690912  
              InParticle%BackScaEff(ib) = 8.979960976123162E-003
              InParticle%AssymParam(ib) = 0.813786479371409  
            ELSE IF (sizebin .EQ. 7) THEN
              InParticle%ExtEff(ib) = 3.22854655841619   
              InParticle%ScaEff(ib) = 2.80386216469626 
              InParticle%BackScaEff(ib) = 6.900609225362380E-003
              InParticle%AssymParam(ib) = 1.92163610766603   
            ELSE IF (sizebin .EQ. 8) THEN
              InParticle%ExtEff(ib) = 3.77483958845704    
              InParticle%ScaEff(ib) = 3.24553098959117  
              InParticle%BackScaEff(ib) = 6.718209735006457E-002
              InParticle%AssymParam(ib) = 2.28403972204950      
            ELSE IF (sizebin .EQ. 9) THEN
              InParticle%ExtEff(ib) = 2.24270928344539    
              InParticle%ScaEff(ib) = 1.64470035277021 
              InParticle%BackScaEff(ib) = 0.211919973372722 
              InParticle%AssymParam(ib) = 0.877378236096786  
            ELSE IF (sizebin .EQ. 10) THEN
              InParticle%ExtEff(ib) = 2.27724359465123   
              InParticle%ScaEff(ib) = 1.56504723032691  
              InParticle%BackScaEff(ib) = 0.132656988973855   
              InParticle%AssymParam(ib) = 1.22490783048364  
            END IF !sizebin

          ELSE                                          !MM-DEBUG
            WRITE(*,*) "ib = ",ib
            call error("Unknown value for Wavelength ib. STOP")
          ENDIF !ib

	   !IF (ib .eq. 419 .or. ib .eq. 420) THEN
                !WRITE(*,*) "ExtEff: ", InParticle%ExtEff(ib)
                !WRITE(*,*) "ScaEff: ", InParticle%ScaEff(ib)
                !WRITE(*,*) "SSA: ", InParticle%ScaEff(ib)/InParticle%ExtEff(ib)
                !WRITE(*,*) "AssymParam: ", InParticle%AssymParam(ib)
           !ENDIF

	   !DMiLay gives the assymetry parameter times the scattering efficiency
	   !So divide by scattering efficiency here (MJA, 10/15/07)
	   InParticle%AssymParam(ib) = InParticle%AssymParam(ib)/InParticle%ScaEff(ib)

          !Now call for a hypothetical particle of
          !same effective size but all BC
!          WRITE(*,*) "Call 2"
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			SolarCoreInd, & 
!			SolarCoreInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffBC(ib), &
!			InParticle%ScaEffBC(ib), &
!			InParticle%BackScaEffBC(ib), &
!			InParticle%AssymParamBC(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamBC(ib) = &
!                InParticle%AssymParamBC(ib)/InParticle%ScaEffBC(ib)

          !Now call for a hypothetical particle of
          !same effective size but all shell
!          WRITE(*,*) "Call 3"
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			SolarShellInd, & 
!			SolarShellInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffShell(ib), &
!			InParticle%ScaEffShell(ib), &
!			InParticle%BackScaEffShell(ib), &
!			InParticle%AssymParamShell(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamShell(ib) = &
!                InParticle%AssymParamShell(ib)/InParticle%ScaEffShell(ib)

!          WRITE(*,*) "Call 4"
          !Now call for volume average dielectric constant mixing rule
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			VAInd, & 
!			VAInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffVA(ib), &
!			InParticle%ScaEffVA(ib), &
!			InParticle%BackScaEffVA(ib), &
!			InParticle%AssymParamVA(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamVA(ib) = &
!                InParticle%AssymParamVA(ib)/InParticle%ScaEffVA(ib)

          !Now call for Maxwell-Garnett mixing rule
!	  CALL DMiLay( InParticle%AbsCoreRad, &
!			InParticle%EffectiveRadius, &
!			Wavelength, &
!			MGInd, & 
!			MGInd, &
!			CosAng, &
!			NumAng, &
!			InParticle%ExtEffMG(ib), &
!			InParticle%ScaEffMG(ib), &
!			InParticle%BackScaEffMG(ib), &
!			InParticle%AssymParamMG(ib), & 
!                        M1, M2, S21, D21, MAXANG )

	   !DMiLay gives the assymetry parameter times the scattering eff.
	   !So divide by scattering efficiency here (MJA, 10/15/07)
!	   InParticle%AssymParamMG(ib) = &
!                InParticle%AssymParamMG(ib)/InParticle%ScaEffMG(ib)
          
	END DO

        DEALLOCATE(AqVolConc,OrgVolConc,AqOrgVolConc, MolFrac, wvln_array, STAT = I)

	RETURN
END SUBROUTINE ShellRefIndAndRad

!Calculates aerosol optical depth for a grid point at 550 nm.
!Can only be used in inversion mode!

! Modified by CMB (AER): Make Gfortran cooperate with type safety
!                        as described above                         
REAL*8 FUNCTION AerosolOpticalDepth ()
! 08/30/2010 Changed for photolysis radiation bins MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE
	
	Type(Particle),POINTER :: p_Particle

	AerosolOpticalDepth = 0.
	
	p_Particle => Particles%First

	DO WHILE (ASSOCIATED(p_Particle))

           AerosolOpticalDepth = AerosolOpticalDepth + &
                p_Particle%NumberofParticles*&
                    (PI*p_Particle%EffectiveRadius**2)*p_Particle%ExtEff(301)

           p_Particle => p_Particle%Next

	END DO
	
	AerosolOpticalDepth = AerosolOpticalDepth*InversionHeight
	
	RETURN
END FUNCTION AerosolOpticalDepth

!Calculates total aerosol extinction coefficient (m**-1) 
!average single scattering albedo, and 
!average assymetry factor 
!for a grid point at 18 radiative bands.
SUBROUTINE AerosolOptProp (ExtCoeff, SingScat, Assym, BackScatCoeff, &
                          SubExtCoeff, SubSingScat, Flag)!&
                          !ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                          !SubExtCoeffSep, SubSingScatSep)
! 08/30/2010 Changed for the photolysis wavebands MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
! 05/03/2012 Changed to calculate backscattering coefficient, units
!             for extcoeff and backscatcoeff now Mm^-1
! 08/16/2012 Fixed backscattering coefficient, added submicron extinction & SSA
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added flag to decide between mixing rules
!            Flag = 0, core-in-shell (default)
!            Flag = 1, external mixture of BC and rest
!            Flag = 2, volume average dielectric constant mixing rule
!            Flag = 3, Maxwell-Garnett mixing rule
!
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

        !Input Variables
        INTEGER :: Flag
	!Output Variables
	REAL*8, DIMENSION(451) :: ExtCoeff, SingScat, Assym, BackScatCoeff
	REAL*8, DIMENSION(451) :: SubExtCoeff, SubSingScat
        !Below are for external mixture of BC and "shell"
!	REAL*8, DIMENSION(451) :: ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep
!	REAL*8, DIMENSION(451) :: SubExtCoeffSep, SubSingScatSep

	!Internal variables
	INTEGER :: ib, numbin
        LOGICAL :: FASTTUV
	REAL*8, DIMENSION(451) :: ScatCoeff, SubScatCoeff
!	REAL*8, DIMENSION(451) :: ScatCoeffSep, SubScatCoeffSep
	REAL*8 :: AssymSum,xxx, hbs_frac, assym_sq, bc_vol_frac

        ! CMB: Gfortran doesn't like the additional definitions of the 
        !      public variables, nor does it like "Particle" being used
        !      both as a type and a variable name
	Type(Particle),POINTER :: p_Particle


        DO ib = 1, 451
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		BackScatCoeff(ib) = 0.0
		SubExtCoeff(ib) = 0.0
		SubScatCoeff(ib) = 0.0
                AssymSum = 0.0
!		ExtCoeffSep(ib) = 0.0
!		ScatCoeffSep(ib) = 0.0
!		BackScatCoeffSep(ib) = 0.0
!		SubExtCoeffSep(ib) = 0.0
!		SubScatCoeffSep(ib) = 0.0
!                AssymSumSep = 0.0
		p_Particle => Particles%First

		DO WHILE (ASSOCIATED(p_Particle))
                     IF (p_Particle%EffectiveRadius .gt. 0) THEN 
                       bc_vol_frac = (p_Particle%AbsCoreRad/&
                               p_Particle%EffectiveRadius)**3
                     ELSE
                       bc_vol_frac = 0.0
                     ENDIF
		     xxx = p_Particle%NumberofParticles*&
                        (PI*p_Particle%EffectiveRadius**2)

                     IF (FLAG .EQ. 0) THEN !Core in Shell  
			ExtCoeff(ib) = ExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
			AssymSum = AssymSum + &
                            xxx*p_Particle%ScaEff(ib)*p_Particle%AssymParam(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
                        ENDIF

                     ELSE IF (Flag .EQ. 1) THEN !External mixture
                        ExtCoeff(ib) = ExtCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ExtEffBC(ib)&
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
			AssymSum = AssymSum + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib)*&
                                p_Particle%AssymParamBC(ib) &
                                       +(1.0-bc_vol_frac)*xxx* &
                                       p_Particle%ScaEffShell(ib)*&
                                       p_Particle%AssymParamShell(ib) 
                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN	
                          SubExtCoeff(ib) = SubExtCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ExtEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
                        ENDIF
	              ELSE IF (FLAG .EQ. 2) THEN !Volume-averaged internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffVA(ib)*p_Particle%AssymParamVA(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
                        ENDIF
                      ELSE IF (FLAG .EQ. 3) THEN !Maxwell-Garnett internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffMG(ib)*&
                            p_Particle%AssymParamMG(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
                        ENDIF
                      ENDIF

                      p_Particle => p_Particle%Next
		END DO

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		SubSingScat(ib) = SubScatCoeff(ib)/SubExtCoeff(ib)                
!		SingScatSep(ib) = ScatCoeffSep(ib)/ExtCoeffSep(ib)
!		SubSingScatSep(ib) = SubScatCoeffSep(ib)/SubExtCoeffSep(ib)   
		!IF(ib .EQ. 419 .OR. ib .EQ. 420) WRITE(*,*) "SingScat: ", SingScat(ib), ScatCoeff(ib), ExtCoeff(ib) 
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
!		AssymSep(ib) = AssymSumSep/ScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff and BackScatCoeff from 1/cm to 1/Mm
		ExtCoeff(ib) = 1.0e8*ExtCoeff(ib)
		SubExtCoeff(ib) = 1.0e8*SubExtCoeff(ib)
!		ExtCoeffSep(ib) = 1.0e8*ExtCoeffSep(ib)
!		SubExtCoeffSep(ib) = 1.0e8*SubExtCoeffSep(ib)


		!hbs_frac = hemispherical backscattering fraction
                !See Wiscombe and Grams (1976) J. Atmos. Sci. 33(12), pp. 2440-2451.
                !Integrate Eq. 23 using Henyey-Greenstein phase function from Eq. 16
                assym_sq = assym(ib)**2
                hbs_frac = 0.5*((assym_sq-1)/assym(ib))*(1.0/SQRT(assym_sq+2*assym(ib)+1)-1.0/SQRT(assym_sq+1))
		BackScatCoeff(ib) = hbs_frac*ScatCoeff(ib)
                BackScatCoeff(ib) = 1.0e8*BackScatCoeff(ib)

!                assym_sq = assymSep(ib)**2
!                hbs_frac = 0.5*((assym_sq-1)/assymSep(ib))*(1.0/SQRT(assym_sq+2*assymSep(ib)+1)-1.0/SQRT(assym_sq+1))
!		BackScatCoeffSep(ib) = hbs_frac*ScatCoeffSep(ib)
!                BackScatCoeffSep(ib) = 1.0e8*BackScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
        
	RETURN
END SUBROUTINE AerosolOptProp

!Calculates total aerosol extinction coefficient (m**-1) 
!average single scattering albedo, and 
!average assymetry factor 
!for a grid point at 18 radiative bands.
SUBROUTINE AerosolOptPropFASTTUV (ExtCoeff, SingScat, Assym, BackScatCoeff, &
                          SubExtCoeff, SubSingScat, Flag)!&
                          !ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep, &
                          !SubExtCoeffSep, SubSingScatSep)
! 08/30/2010 Changed for the photolysis wavebands MJA
! 02/27/2012 Changed for photolysis radiation bins MJA	
! 05/03/2012 Changed to calculate backscattering coefficient, units
!             for extcoeff and backscatcoeff now Mm^-1
! 08/16/2012 Fixed backscattering coefficient, added submicron extinction & SSA
! 08/17/2012 Added capability for external mixtures
! 11/08/2012 Added flag to decide between mixing rules
!            Flag = 0, core-in-shell (default)
!            Flag = 1, external mixture of BC and rest
!            Flag = 2, volume average dielectric constant mixing rule
!            Flag = 3, Maxwell-Garnett mixing rule
!
	USE ModelParameters, ONLY : PI, InversionHeight

	IMPLICIT NONE

        !Input Variables
        INTEGER :: Flag
	!Output Variables
	REAL*8, DIMENSION(7) :: ExtCoeff, SingScat, Assym, BackScatCoeff
	REAL*8, DIMENSION(7) :: SubExtCoeff, SubSingScat
        !Below are for external mixture of BC and "shell"
!	REAL*8, DIMENSION(7) :: ExtCoeffSep, SingScatSep, AssymSep, BackScatCoeffSep
!	REAL*8, DIMENSION(7) :: SubExtCoeffSep, SubSingScatSep

	!Internal variables
	INTEGER :: ib, numbin
        LOGICAL :: FASTTUV
	REAL*8, DIMENSION(7) :: ScatCoeff, SubScatCoeff
!	REAL*8, DIMENSION(7) :: ScatCoeffSep, SubScatCoeffSep
	REAL*8 :: AssymSum,xxx, hbs_frac, assym_sq, bc_vol_frac

        ! CMB: Gfortran doesn't like the additional definitions of the 
        !      public variables, nor does it like "Particle" being used
        !      both as a type and a variable name
	Type(Particle),POINTER :: p_Particle


        DO ib = 1, 7
		ExtCoeff(ib) = 0.0
		ScatCoeff(ib) = 0.0
		BackScatCoeff(ib) = 0.0
		SubExtCoeff(ib) = 0.0
		SubScatCoeff(ib) = 0.0
                AssymSum = 0.0
!		ExtCoeffSep(ib) = 0.0
!		ScatCoeffSep(ib) = 0.0
!		BackScatCoeffSep(ib) = 0.0
!		SubExtCoeffSep(ib) = 0.0
!		SubScatCoeffSep(ib) = 0.0
!                AssymSumSep = 0.0
		p_Particle => Particles%First

		DO WHILE (ASSOCIATED(p_Particle))

!                     p_Particle%Temperature = TEMP
                     ! WRITE(*,*)    'cur%Temperature = ',TEMP   
                     !Recalculate optical parameters
                     !WRITE(*,*) "Call Optical, Particle #", I
                     CALL ShellRefIndAndRad(p_Particle)

                     IF (p_Particle%EffectiveRadius .gt. 0) THEN 
                       bc_vol_frac = (p_Particle%AbsCoreRad/&
                               p_Particle%EffectiveRadius)**3
                     ELSE
                       bc_vol_frac = 0.0
                     ENDIF
		     xxx = p_Particle%NumberofParticles*&
                        (PI*p_Particle%EffectiveRadius**2)

                     IF (FLAG .EQ. 0) THEN !Core in Shell  
			ExtCoeff(ib) = ExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
			AssymSum = AssymSum + &
                            xxx*p_Particle%ScaEff(ib)*p_Particle%AssymParam(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + &
                            xxx*p_Particle%ExtEff(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                            xxx*p_Particle%ScaEff(ib)
                        ENDIF

                     ELSE IF (Flag .EQ. 1) THEN !External mixture
                        ExtCoeff(ib) = ExtCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ExtEffBC(ib)&
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                            +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
			AssymSum = AssymSum + &
                            bc_vol_frac*xxx*p_Particle%ScaEffBC(ib)*&
                                p_Particle%AssymParamBC(ib) &
                                       +(1.0-bc_vol_frac)*xxx* &
                                       p_Particle%ScaEffShell(ib)*&
                                       p_Particle%AssymParamShell(ib) 
                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN	
                          SubExtCoeff(ib) = SubExtCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ExtEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ExtEffShell(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + &
                                        bc_vol_frac*xxx*p_Particle%ScaEffBC(ib) &
                                        +(1.0-bc_vol_frac)*xxx*p_Particle%ScaEffShell(ib)
                        ENDIF
	              ELSE IF (FLAG .EQ. 2) THEN !Volume-averaged internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffVA(ib)*p_Particle%AssymParamVA(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffVA(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffVA(ib)
                        ENDIF
                      ELSE IF (FLAG .EQ. 3) THEN !Maxwell-Garnett internal mixture
			ExtCoeff(ib) = ExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			ScatCoeff(ib) = ScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
			AssymSum = AssymSum + xxx*p_Particle%ScaEffMG(ib)*&
                            p_Particle%AssymParamMG(ib)

                        !Particles with diameters less than 1 micron
                        IF (2*p_Particle%EffectiveRadius .LT. 1.0e-4) THEN
 			  SubExtCoeff(ib) = SubExtCoeff(ib) + xxx*p_Particle%ExtEffMG(ib)
			  SubScatCoeff(ib) = SubScatCoeff(ib) + xxx*p_Particle%ScaEffMG(ib)
                        ENDIF
                      ENDIF

                      p_Particle => p_Particle%Next
		END DO

		!Formula from p.327, eqn 9.103 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		SingScat(ib) = ScatCoeff(ib)/ExtCoeff(ib)
		SubSingScat(ib) = SubScatCoeff(ib)/SubExtCoeff(ib)                
!		SingScatSep(ib) = ScatCoeffSep(ib)/ExtCoeffSep(ib)
!		SubSingScatSep(ib) = SubScatCoeffSep(ib)/SubExtCoeffSep(ib)   
		!IF(ib .EQ. 419 .OR. ib .EQ. 420) WRITE(*,*) "SingScat: ", SingScat(ib), ScatCoeff(ib), ExtCoeff(ib) 
		
		!Formula from p.329, eqn 9.117 of Jacobson, 2nd ed. "Fundamentals of Atm. Modeling"
		Assym(ib) = AssymSum/ScatCoeff(ib)
!		AssymSep(ib) = AssymSumSep/ScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "Assym: ", Assym(1)

		!Convert ExtCoeff and BackScatCoeff from 1/cm to 1/Mm
		ExtCoeff(ib) = 1.0e8*ExtCoeff(ib)
		SubExtCoeff(ib) = 1.0e8*SubExtCoeff(ib)
!		ExtCoeffSep(ib) = 1.0e8*ExtCoeffSep(ib)
!		SubExtCoeffSep(ib) = 1.0e8*SubExtCoeffSep(ib)


		!hbs_frac = hemispherical backscattering fraction
                !See Wiscombe and Grams (1976) J. Atmos. Sci. 33(12), pp. 2440-2451.
                !Integrate Eq. 23 using Henyey-Greenstein phase function from Eq. 16
                assym_sq = assym(ib)**2
                hbs_frac = 0.5*((assym_sq-1)/assym(ib))*(1.0/SQRT(assym_sq+2*assym(ib)+1)-1.0/SQRT(assym_sq+1))
		BackScatCoeff(ib) = hbs_frac*ScatCoeff(ib)
                BackScatCoeff(ib) = 1.0e8*BackScatCoeff(ib)

!                assym_sq = assymSep(ib)**2
!                hbs_frac = 0.5*((assym_sq-1)/assymSep(ib))*(1.0/SQRT(assym_sq+2*assymSep(ib)+1)-1.0/SQRT(assym_sq+1))
!		BackScatCoeffSep(ib) = hbs_frac*ScatCoeffSep(ib)
!                BackScatCoeffSep(ib) = 1.0e8*BackScatCoeffSep(ib)
		!IF(ib .EQ. 61) WRITE(*,*) "ExtCoeff: ", ExtCoeff(1)

	END DO
	
        
	RETURN
END SUBROUTINE AerosolOptPropFASTTUV
