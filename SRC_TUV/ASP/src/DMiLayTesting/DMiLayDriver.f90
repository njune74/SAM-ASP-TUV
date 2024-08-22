PROGRAM DMiLayDriver

  ! compiling:
  ! ifort -DIFORT -fp-stack-check -fp-model source -c -fpp -FI -O -r8 -prec_div -g -traceback -I/nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//ModuleFiles/ -module /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//ModuleFiles/ DMiLay.f
  !
  ! ifort -r8 -prec_div DMiLayDriver.f90 DMiLay.o -o DMiLayDriver
  !
  ! DMiLay.f compiled in Makfile as:
  ! ifort -DIFORT -fp-stack-check -fp-model source -c -fpp -FI -O -r8 -prec_div -g -traceback -I/nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//ModuleFiles/ -module /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//ModuleFiles/ /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//src/DMiLay.f
  !
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
  
  IMPLICIT NONE
  
  !USE DMiLay
	
  !! External Variables read in by ShellRefIndAndRad
  !TYPE(Particle),POINTER :: InParticle
  REAL*8	:: InParticleAbsCoreRad
  REAL*8	:: InParticleEffectiveRadius
  REAL*8	:: InParticleExtEffib   
  REAL*8	:: InParticleScaEffib
  REAL*8	:: InParticleBackScaEffib
  REAL*8	:: InParticleAssymParamib

  INTEGER	:: NumAng
  REAL*8	:: Wavelength
  INTEGER, PARAMETER  :: MaxAng = 100
  COMPLEX :: SolarShellInd, SolarCoreInd
  REAL*8	:: CosAng(MaxAng), M1( MaxAng, 2), M2( MaxAng, 2 ), S21( MaxAng, 2 ),D21( MaxAng, 2 )
  
  WRITE(*,*) "Starting DMiLayDriver"

  ! zero values= NumAng, CosAng, M1, M2, S21, D21
  InParticleAbsCoreRad= 1.466775439899720E-003
  InParticleEffectiveRadius=  3.602062048069621E-003 
  Wavelength= 211057.618648962 
  SolarShellInd=  CMPLX(1.14136576121514,-2.097741127316012E-003)
  SolarCoreInd=  CMPLX(1.73448000000000,-0.469080000000000)
  NumAng = 0
  InParticleExtEffib=    1.73344982058878     
  InParticleScaEffib=     1.46631605444476     
  InParticleBackScaEffib=   1.387942070974863E-003
  InParticleAssymParamib=   0.929301976265163  
 
  WRITE(*,*)
  WRITE(*,*) "InParticleAbsCoreRad: ", InParticleAbsCoreRad
  WRITE(*,*) "InParticleEffectiveRadius: ", InParticleEffectiveRadius
  WRITE(*,*) "Wavelength: ", Wavelength
  WRITE(*,*) "SolarShellInd: ", SolarShellInd
  WRITE(*,*) "SolarCoreInd: ", SolarCoreInd
  WRITE(*,*) "CosAng: ", CosAng
  WRITE(*,*) "NumAng: ", NumAng 
  WRITE(*,*) "InParticleExtEffib: ", InParticleExtEffib
  WRITE(*,*) "InParticleScaEffib: ", InParticleScaEffib
  WRITE(*,*) "InParticleBackScaEffib: ", InParticleBackScaEffib
  WRITE(*,*) "InParticleAssymParamib: ", InParticleAssymParamib
  WRITE(*,*) "M1: ", M1
  WRITE(*,*) "M2: ", M2			
  WRITE(*,*) "S21: ", S21
  WRITE(*,*) "D21: ", D21  
  WRITE(*,*) "MAXANG: ", MAXANG
  WRITE(*,*)

  !WRITE(*,*) "REAL(SolarShellInd): ", REAL(SolarShellInd)
  !WRITE(*,*) "AIMAG(SolarShellInd): ", AIMAG(SolarShellInd)
  !IF( REAL(SolarShellInd).LE.0.0 .OR. AIMAG(SolarShellInd).GT.0.0 ) THEN
  !   WRITE(*,*) "ERROR IN SolarShellInd:"
  !ELSE
  !   WRITE(*,*) "No error with SolarShellInd"
  !END IF

  !WRITE(*,*)
  !WRITE(*,*) "REAL(SolarCoreInd): ", REAL(SolarCoreInd)
  !WRITE(*,*) "AIMAG(SolarCoreInd): ", AIMAG(SolarCoreInd) 
  !IF( REAL(SolarCoreInd).LE.0.0 .OR. AIMAG(SolarCoreInd).GT.0.0 ) THEN
  !   WRITE(*,*) "ERROR IN SolarCoreInd:"
  !ELSE
  !   WRITE(*,*) "No error with SolarCoreInd"
  !END IF

  WRITE(*,*)
  WRITE(*,*) "InParticleAbsCoreRad: ", InParticleAbsCoreRad
  WRITE(*,*) "InParticleEffectiveRadius: ", InParticleEffectiveRadius
  WRITE(*,*) "Call DMiLay"

  CALL DMiLay( InParticleAbsCoreRad, &
       InParticleEffectiveRadius, &
       Wavelength, &
       SolarShellInd, & 
       SolarCoreInd, &
       CosAng, &
       NumAng, &
       InParticleExtEffib, &
       InParticleScaEffib, &
       InParticleBackScaEffib, &
       InParticleAssymParamib, & 
       M1, M2, S21, D21, MAXANG )

  WRITE(*,*) "Ending DMiLayDriver"

END PROGRAM DMiLayDriver
