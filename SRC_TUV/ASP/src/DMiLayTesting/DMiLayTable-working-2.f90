PROGRAM DMiLayTable

  ! compiling:
  ! ifort -DIFORT -fp-stack-check -fp-model source -c -fpp -FI -O -r8 -prec_div -g -traceback -I/nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//ModuleFiles/ -module /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/SRC_TUV/ASP//ModuleFiles/ DMiLay.f
  !
  ! ifort -r8 -prec_div DMiLayTable.f90 DMiLay.o -o DMiLayTable
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
  REAL*8	:: EffRadBinLim(11)
  REAL*8	:: InParticleExtEffib   
  REAL*8	:: InParticleScaEffib
  REAL*8	:: InParticleBackScaEffib
  REAL*8	:: InParticleAssymParamib

  INTEGER	:: NumAng
  REAL*8	:: Wavelength
  INTEGER, PARAMETER  :: MaxAng = 100
  COMPLEX :: SolarShellInd, SolarCoreInd
  REAL*8	:: CosAng(MaxAng), M1( MaxAng, 2), M2( MaxAng, 2 ), S21( MaxAng, 2 ),D21( MaxAng, 2 )
  
  REAL*8 :: a_w, b_w, n_real, n_imag
  REAL*8, ALLOCATABLE	:: wvln_array(:)
  INTEGER, ALLOCATABLE	:: wvln_rounddown(:)
  INTEGER :: I, ib, ief, numbin, w_ind_low
  LOGICAL :: FASTTUV
  REAL*8 :: WAVE(61), SOOT_REAL(61), SOOT_IMAG(61), WASO_REAL(61), WASO_IMAG(61)
  real*8, parameter :: LengthScale = 1.
  real*8, parameter :: nm = LengthScale * 1.d-7
  real*8, parameter  :: pi=3.1415926535898

  WRITE(*,*) "Starting DMiLayTable"

  ! Effective radius limits (cm) for bins 1 to 10
  ! NOTE: bin 10 actually has infinite upper bound, use 3.0E-4 to
  !      get ~2micrometer average
  EffRadBinLim = (/0., 2.5e-06, 4.0e-06, 6.3e-06, 1.0e-05, 1.59e-05, &
                   2.52e-05, 4.01e-05, 6.36e-05, 1.01e-04, 3.0e-04/) 


  NumAng = 0

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ParticleAttributes.h L1852, ShellRefIndAndRad subroutine           
  !Calculate shell and core refractive index for each wavelength bin
  FASTTUV = .TRUE.
  IF(FASTTUV) THEN !Use center wavelengths for 7 tropospheric fastTUV bins
     numbin = 7
     ALLOCATE (wvln_array(numbin), STAT = I)
     IF (I > 0) THEN
        WRITE(*,*) "ERROR"
        WRITE(*,*) "numbin = ", numbin
        WRITE(*,*) "In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_array at size numbin"
        STOP
     END IF

     ALLOCATE (wvln_rounddown(numbin), STAT = I)
     IF (I > 0) THEN
        WRITE(*,*) "ERROR"
        WRITE(*,*) "numbin = ", numbin
        WRITE(*,*) "In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_rounddown at size numbin"
        STOP
     END IF
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
     IF (I > 0) THEN
        WRITE(*,*) "ERROR"
        WRITE(*,*) "numbin = ", numbin
        WRITE(*,*) "In ShellRefIndAndRad(), couldn't allocate "// &
                               "wvln_array at size numbin"
        STOP
     END IF
     DO ib = 1, numbin 
        wvln_array(ib) = (250+(ib-1)*1) !in nm
     ENDDO
  ENDIF

  DO ib = 1, numbin
     WRITE(*,*) 
     WRITE(*,*) "********************** ib = ",ib," ***********************"
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
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! ParticleAttributes.h L2000, ShellRefIndAndRad subroutine
     n_real = a_w*WASO_REAL(w_ind_low)+b_w*WASO_REAL(w_ind_low+1)
     n_imag = -1.0*(a_w*WASO_IMAG(w_ind_low)+b_w*WASO_IMAG(w_ind_low+1))
     SolarShellInd = CMPLX(n_real, -1.0*n_imag)

     ! ParticleAttributes.h L2071, ShellRefIndAndRad subroutine
     !DMiLay needs 2*Pi divided by the wavelength
     Wavelength = 2*PI/((Wavelength)*nm)  

     ! loop over 10 bins assuming average effective radius for each
     ! each bin
     DO ief = 1, 10
        WRITE(*,*) 
        WRITE(*,*) "******************** ief = ",ief," *********************"
        ! set radius values (cm)
        InParticleEffectiveRadius =(EffRadBinLim(ief) + EffRadBinLim(ief+1))/2
        InParticleAbsCoreRad=0.4*InParticleEffectiveRadius

        WRITE(*,*)
        WRITE(*,*) "InParticleAbsCoreRad: R", InParticleAbsCoreRad
        WRITE(*,*) "InParticleEffectiveRadius: ", InParticleEffectiveRadius
        WRITE(*,*) "Wavelength: ", Wavelength
        WRITE(*,*) "SolarShellInd: ", SolarShellInd
        WRITE(*,*) "SolarCoreInd: ", SolarCoreInd
        !WRITE(*,*) "CosAng: ", CosAng
        WRITE(*,*) "NumAng: ", NumAng 
        !WRITE(*,*) "InParticleExtEffib: ", InParticleExtEffib
        !WRITE(*,*) "InParticleScaEffib: ", InParticleScaEffib
        !WRITE(*,*) "InParticleBackScaEffib: ", InParticleBackScaEffib
        !WRITE(*,*) "InParticleAssymParamib: ", InParticleAssymParamib
        !WRITE(*,*) "M1: ", M1
        !WRITE(*,*) "M2: ", M2			
        !WRITE(*,*) "S21: ", S21
        !WRITE(*,*) "D21: ", D21  
        WRITE(*,*) "MAXANG: ", MAXANG
        WRITE(*,*)
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

        WRITE(*,*) "Ending DMiLay"
        WRITE(*,*)
        WRITE(*,*) "InParticleExtEffib: ", InParticleExtEffib
        WRITE(*,*) "InParticleScaEffib: ", InParticleScaEffib
        WRITE(*,*) "InParticleBackScaEffib: ", InParticleBackScaEffib
        WRITE(*,*) "InParticleAssymParamib: ", InParticleAssymParamib

     ENDDO !ief      
  END DO !ib
  DEALLOCATE(wvln_array, STAT = I)
END PROGRAM DMiLayTable
