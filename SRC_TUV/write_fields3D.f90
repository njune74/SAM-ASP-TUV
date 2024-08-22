subroutine write_fields3D

use grid, only: dx, dy, dz, pres, zi	
use vars
use sgs
use rad, only: qrad
use params
use microphysics, only: nmicro_fields, micro_field, flag_number, &
     flag_micro3Dout, mkname, mklongname, mkunits, mkoutputscale, &
     index_water_vapor, GET_reffc, Get_reffi
use tracers
use GridPointFields
use ModelParameters
implicit none
 character *120 filename
 character *80 long_name
 character *8 name
 character *10 timechar
 character *4 rankchar
 character *5 sepchar
 character *6 filetype
 character *10 units
 character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,l,k,n,nfields,nfields1,binnum
integer nbins,ntrtrack

real tmp(nx,ny,nzm)
!real TotalOA_01(nx,ny,nzm)
!real TotalOA_02(nx,ny,nzm)
!real TotalOA_03(nx,ny,nzm)
!real TotalOA_04(nx,ny,nzm)
!real TotalOA_05(nx,ny,nzm)
!real TotalOA_06(nx,ny,nzm)
!real TotalOA_07(nx,ny,nzm)
!real TotalOA_08(nx,ny,nzm)
!real TotalOA_09(nx,ny,nzm)
!real TotalOA_10(nx,ny,nzm)

real :: TotOA(nx,ny,nzm)
!real TotH2O(nx,ny,nzm)
real :: TotK(nx,ny,nzm)
real :: TotNa(nx,ny,nzm)
real :: TotNH4(nx,ny,nzm)
real :: TotCl(nx,ny,nzm)
real :: TotNO3(nx,ny,nzm)
real :: TotSO4(nx,ny,nzm)
real :: TotMg(nx,ny,nzm)
real :: TotCa(nx,ny,nzm)

real*8 :: std(nx,ny,nzm) ! -CRL
real*8 :: airmols ! molar density of air (mols/m3)
real*8 :: boxmols ! mols of air in gridbox 
real*8 :: box_vol ! volume of gridbox (cm3)
real*8 :: airdens
real*8 :: airdenstm
real*8,parameter :: MwK = 0.039        ! molecular weight [kg/mol]
real*8,parameter :: MwKCl = 0.075      ! molecular weight [kg/mol]
real*8,parameter :: MwKNO3 = 0.101     ! molecular weight [kg/mol]
real*8,parameter :: MwK2SO4 = 0.174    ! molecular weight [kg/mol]
real*8,parameter :: MwKHSO4 = 0.136    ! molecular weight [kg/mol]
real*8,parameter :: MwKOH = 0.056      ! molecular weight [kg/mol]

real*8,parameter :: MwNa = 0.023       ! molecular weight [kg/mol]
real*8,parameter :: MwNaCl = 0.058     ! molecular weight [kg/mol]
real*8,parameter :: MwNaNO3 = 0.085    ! molecular weight [kg/mol]
real*8,parameter :: MwNa2SO4 = 0.142   ! molecular weight [kg/mol]
real*8,parameter :: MwNaHSO4 = 0.120   ! molecular weight [kg/mol]
real*8,parameter :: MwNaOH = 0.040     ! molecular weight [kg/mol]

real*8,parameter :: MwNH4 =0.018       ! molecular weight [kg/mol]
real*8,parameter :: MwNH3 = 0.017      ! molecular weight [kg/mol]
real*8,parameter :: MwNH4NO3 = 0.080   ! molecular weight [kg/mol]
real*8,parameter :: MwNH4Cl = 0.054    ! molecular weight [kg/mol]
real*8,parameter :: MwNH4HSO4 = 0.115  ! molecular weight [kg/mol]
real*8,parameter :: MwNH42SO4 = 0.132  ! molecular weight [kg/mol]
real*8,parameter :: MwLEV = 0.247      ! molecular weight [kg/mol]

real*8,parameter :: MwCl = 0.036       ! molecular weight [kg/mol]
real*8,parameter :: MwMgCl2 = 0.095    ! molecular weight [kg/mol]
real*8,parameter :: MwCaCl2 = 0.111    ! molecular weight [kg/mol]

real*8,parameter :: MwNO3 =0.062       ! molecular weight [kg/mol]
real*8,parameter :: MwMgNO32 =0.148    ! molecular weight [kg/mol]
real*8,parameter :: MwCaNO32 =0.164    ! molecular weight [kg/mol]

real*8,parameter :: MwSO4 =0.096       ! molecular weight [kg/mol]
real*8,parameter :: MwHSO4 =0.097      ! molecular weight [kg/mol]
!real*8,parameter :: MwCaSO4 =0.136     ! molecular weight [kg/mol]
real*8,parameter :: MwCaSO4 =0.172     ! molecular weight [kg/mol]
real*8,parameter :: MwMgSO4 =0.120     ! molecular weight [kg/mol]
real*8,parameter :: MwMgHSO42 =0.218   ! molecular weight [kg/mol]
real*8,parameter :: MwCaHSO42 =0.233   ! molecular weight [kg/mol]

real*8,parameter :: MwMg =0.024        ! molecular weight [kg/mol]
real*8,parameter :: MwMgOH2 =0.058     ! molecular weight [kg/mol]

real*8,parameter :: MwCa =0.040        ! molecular weight [kg/mol]
real*8,parameter :: MwCaOH2 =0.074     ! molecular weight [kg/mol]

Integer :: 	Kind,KNO3ind, KClind, K2SO4ind, KHSO4ind, KOHind, Naind , NaNO3ind, NaClind, &
                NaHSO4ind, Na2SO4ind, NaOHind, Mgind, MgNO32ind, MgCl2ind, MgHSO42ind, &
                MgSO4ind, MgOH2ind, Caind, CaNO32ind, CaCl2ind, CaHSO42ind, CaSO4ind, &
                CaOH2ind, NH3ind, NH4ind, NH4NO3ind, NH4Clind, NH4HSO4ind, NH42SO4ind, LEVind, &
                Clind, NO3ind, SO4ind, HSO4ind,BCind,cur_index,H2Oind,HNO3ind, IVOC1ind,IVOC2ind,&
                IVOC3ind,IVOC4ind,IVOC5ind,IVOC6ind,IVOC7ind,IVOC8ind,IVOC9ind,POA1ind,&
                POA2ind,POA3ind,POA4ind,POA5ind,POA6ind,POA7ind,POA8ind,LEVOind,&
                CBIOind,PD1ind,CPD1ind, CPD2ind,CPD3ind,HClind,H2SO4ind,AqIVOC1ind,AqIVOC2ind,&
                AqIVOC3ind,AqIVOC4ind,AqIVOC5ind,AqIVOC6ind,AqIVOC7ind,AqIVOC8ind,AqIVOC9ind,&
                AqPOA1ind,AqPOA2ind,AqPOA3ind,AqPOA4ind,AqPOA5ind,AqPOA6ind,AqPOA7ind,AqPOA8ind,&
                AqLEVOind,AqCBIOind,AqCPD1ind,AqCPD2ind,AqCPD3ind


BCind = Find_Chem("BC01")
Kind = Find_Chem("K+01")
KNO3ind = Find_Chem("KNO301")
KClind = Find_Chem("KCl01")
K2SO4ind = Find_Chem("K2SO401")
KHSO4ind = Find_Chem("KHSO401")
KOHind = Find_Chem("KOH01")
Naind = Find_Chem("Na+01")
NaNO3ind = Find_Chem("NaNO301")
NaClind = Find_Chem("NaCl01")
NaHSO4ind = Find_Chem("NaHSO401")
Na2SO4ind = Find_Chem("Na2SO401")
NaOHind = Find_Chem("NaOH01")
Mgind = Find_Chem("Mg++01")
MgNO32ind = Find_Chem("Mg(NO3)201")
MgCl2ind = Find_Chem("MgCl201")
MgHSO42ind = Find_Chem("Mg(HSO4)201")
MgSO4ind = Find_Chem("MgSO401")
MgOH2ind = Find_Chem("Mg(OH)201")
Caind = Find_Chem("Ca++01")
CaNO32ind = Find_Chem("Ca(NO3)201")
CaCl2ind = Find_Chem("CaCl201")
CaHSO42ind = Find_Chem("Ca(HSO4)201")
CaSO4ind = Find_Chem("CaSO4*2H2O01")
CaOH2ind = Find_Chem("Ca(OH)201")	
NH3ind = Find_Chem("NH301")
NH4ind = Find_Chem("NH4+01")
NH4NO3ind = Find_Chem("NH4NO301")
NH4Clind = Find_Chem("NH4Cl01")
NH4HSO4ind = Find_Chem("NH4HSO401")
NH42SO4ind = Find_Chem("(NH4)2SO401")
LEVind = Find_Chem("(NH4)3H(SO4)201")
Clind = Find_Chem("Cl-01")
NO3ind = Find_Chem("NO3-01")
SO4ind = Find_Chem("SO4--01")
HSO4ind = Find_Chem("HSO4-01")

H2Oind = Find_Chem("H2O01")
HNO3ind = Find_Chem("HNO301")

IVOC1ind = Find_Chem("IVOC101")
IVOC2ind = Find_Chem("IVOC201")
IVOC3ind = Find_Chem("IVOC301")
IVOC4ind = Find_Chem("IVOC401")
IVOC5ind = Find_Chem("IVOC501")
IVOC6ind = Find_Chem("IVOC601")
IVOC7ind = Find_Chem("IVOC701")
IVOC8ind = Find_Chem("IVOC801")
IVOC9ind = Find_Chem("IVOC901")

POA1ind = Find_Chem("POA101")
POA2ind = Find_Chem("POA201")
POA3ind = Find_Chem("POA301")
POA4ind = Find_Chem("POA401")
POA5ind = Find_Chem("POA501")
POA6ind = Find_Chem("POA601")
POA7ind = Find_Chem("POA701")
POA8ind = Find_Chem("POA801")

LEVOind = Find_Chem("LEVO01")

CBIOind = Find_Chem("CBIO01")
CPD1ind = Find_Chem("CPD101")
CPD2ind = Find_Chem("CPD201")
CPD3ind = Find_Chem("CPD301")
!AqOrg
AqIVOC1ind = Find_Chem("AqIVOC101")
AqIVOC2ind = Find_Chem("AqIVOC201")
AqIVOC3ind = Find_Chem("AqIVOC301")
AqIVOC4ind = Find_Chem("AqIVOC401")
AqIVOC5ind = Find_Chem("AqIVOC501")
AqIVOC6ind = Find_Chem("AqIVOC601")
AqIVOC7ind = Find_Chem("AqIVOC701")
AqIVOC8ind = Find_Chem("AqIVOC801")
AqIVOC9ind = Find_Chem("AqIVOC901")

AqPOA1ind = Find_Chem("AqPOA101")
AqPOA2ind = Find_Chem("AqPOA201")
AqPOA3ind = Find_Chem("AqPOA301")
AqPOA4ind = Find_Chem("AqPOA401")
AqPOA5ind = Find_Chem("AqPOA501")
AqPOA6ind = Find_Chem("AqPOA601")
AqPOA7ind = Find_Chem("AqPOA701")
AqPOA8ind = Find_Chem("AqPOA801")

AqLEVOind = Find_Chem("AqLEVO01")

AqCBIOind = Find_Chem("AqCBIO01")
AqCPD1ind = Find_Chem("AqCPD101")
AqCPD2ind = Find_Chem("AqCPD201")
AqCPD3ind = Find_Chem("AqCPD301")


!nfields=ntracers+ngas+3!! number of 3D fields to save 
!nfields=41+10+19+60!! number of 3D fields to save
nfields=130
if(.not.docloud) nfields=nfields-1
if(.not.doprecip) nfields=nfields-1
!bloss: add 3D outputs for microphysical fields specified by flag_micro3Dout
!       except for water vapor (already output as a SAM default).
if(docloud) nfields=nfields+SUM(flag_micro3Dout)-flag_micro3Dout(index_water_vapor)

if((dolongwave.or.doshortwave).and..not.doradhomo) nfields=nfields+1
if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) nfields=nfields+1

nfields1=0

if(masterproc.or.output_sep) then

  if(output_sep) then
     write(rankchar,'(i4)') rank
     sepchar="_"//rankchar(5-lenstr(rankchar):4)
  else
     sepchar=""
  end if
  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
    if(save3Dbin) then
      filetype = '.bin3D'
    else
      filetype = '.com3D'
    end if
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    open(46,file=filename,status='unknown',form='unformatted')

  else
    if(save3Dbin) then
     if(save3Dsep) then
       filetype = '.bin3D'
     else
       filetype = '.bin2D'
     end if
    else
     if(save3Dsep) then
       filetype = '.com3D'
     else
       filetype = '.com2D'
     end if
    end if
    if(save3Dsep) then
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
      open(46,file=filename,status='unknown',form='unformatted')	
    else
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//filetype//sepchar
      if(nrestart.eq.0.and.notopened3D) then
         open(46,file=filename,status='unknown',form='unformatted')	
      else
         open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
      notopened3D=.false.
    end if  

  end if

  if(masterproc) then

    if(save3Dbin) then

      write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
        write(46) real(z(k),4)
      end do
      do k=1,nzm
        write(46) real(pres(k),4)
      end do
      write(46) real(dx,4)
      write(46) real(dy,4)
      write(46) real(float(nstep)*dt/(3600.*24.)+day0,4)

    else

      write(long_name,'(8i4)') nx,ny,nzm,nsubdomains, &
                                   nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
         write(c_z(k),'(f12.3)') z(k)
      end do
      do k=1,nzm
         write(c_p(k),'(f12.3)') pres(k)
      end do
      write(c_dx,'(f12.5)') dx
      write(c_dy,'(f12.5)') dy
      write(c_time,'(f12.5)') nstep*dt/(3600.*24.)+day0
	
      write(46) long_name(1:32)
      write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

    end if ! save3Dbin

  end if ! masterproc
 
end if ! masterproc.or.output_sep

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=u(i,j,k) + ug
    end do
   end do
  end do
  name='U'
  long_name='X Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=v(i,j,k) + vg
    end do
   end do
  end do
  name='V'
  long_name='Y Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=w(i,j,k)
    end do
   end do
  end do
  name='W'
  long_name='Z Wind Component'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=p(i,j,k)
    end do
   end do
  end do
  name='PP'
  long_name='Pressure Perturbation'
  units='Pa'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


if((dolongwave.or.doshortwave).and..not.doradhomo) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qrad(i,j,k)*86400.
    end do
   end do
  end do
  name='QRAD'
  long_name='Radiative heating rate'
  units='K/day'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if
if(compute_reffc.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffc()
  name='REL'
  long_name='Effective Radius for Cloud Liquid Water'
  units='mkm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if
if(compute_reffi.and.(dolongwave.or.doshortwave).and.rad3Dout) then
  nfields1=nfields1+1
  tmp(1:nx,1:ny,1:nzm)=Get_reffi()
  name='REI'
  long_name='Effective Radius for Cloud Ice'
  units='mkm'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if


  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tabs(i,j,k)
    end do
   end do
  end do
  name='TABS'
  long_name='Absolute Temperature'
  units='K'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=qv(i,j,k)*1.e3
    end do
   end do
  end do
  name='QV'
  long_name='Water Vapor'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

if(docloud) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(qcl(i,j,k)+qci(i,j,k))*1.e3
    end do
   end do
  end do
  name='QN'
  long_name='Non-precipitating Condensate (Water+Ice)'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if


if(doprecip) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=(qpl(i,j,k)+qpi(i,j,k))*1.e3
    end do
   end do
  end do
  name='QP'
  long_name='Precipitating Water (Rain+Snow)'
  units='g/kg'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if


write(*,*) 'In write_fields3D '
do n=1,100
  nfields1=nfields1+1

   do i=1,nx
      do j=1,ny
         do l=1,nzm

            tmp(i,j,l)=tracer(i,j,l,n) ! use for values averaged of nsteps3D (in kg)

            tmp(i,j,l) = tmp(i,j,l)*airmw/GasMolecularMass(n)*1e9
            std(i,j,l) = std(i,j,l)*airmw/GasMolecularMass(n)*1e9
		
    end do
   end do
  end do
  name=tracername(n)
  long_name=tracername(n)
  units = '[ppb]'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
enddo

nfields1=nfields1+1
!write(*,*) 'In write_fields ',GasPhaseChemicalNames(167)! PANindex = 167
do i=1,nx
   do j=1,ny
      do l=1,nzm
            tmp(i,j,l)=tracer(i,j,l,167) ! use for values averaged of nsteps3D (in kg)
	
            tmp(i,j,l) = tmp(i,j,l)*airmw/GasMolecularMass(167)*1e9
            std(i,j,l) = std(i,j,l)*airmw/GasMolecularMass(167)*1e9		
      end do
   end do
end do
name=tracername(167)
long_name=tracername(167)
units = '[ppb]'
call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                               save3Dbin,dompi,rank,nsubdomains)

if(dotracers) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tkh(i,j,k)
    end do
   end do
  end do
  name='TKH'
  long_name='Eddy diffusivity'
  units='m2/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if

if(dotracers) then
  nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=tke(i,j,k)
    end do
   end do
  end do
  name='TKE'
  long_name='Turbulent kinetic energy (resolved)'
  units='m2/s2'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end if


!TotOA
nfields1=nfields1+1

    ! write(*,*)'TotOA = ', tracername(n),n
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotOA(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            do binnum=1,10
               !tmp(i,j,l)= tracer(i,j,l,(Kind+cur_index))*(airdenstm/Av*airmw)*1e12 ! kg/kg_air to ug/m3)
               !std(i,j,l)= tracer(i,j,l,(Kind+cur_index))*(airdenstm/Av*airmw)*1e12 ! kg/kg_air to ug/m3)

               cur_index = 86*(binnum-1)
                ! Org       
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA1ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               !write(*,*)'POA1 =',tracername(POA1ind+cur_index)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA2ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)

               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA3ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA4ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA5ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA6ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA7ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(POA8ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(LEVOind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(CBIOind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(CPD1ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(CPD2ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(CPD3ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC1ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC2ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC3ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC4ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC5ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC6ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC7ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC8ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               !write(*,*)'IVOC8 =',tracername(IVOC8ind+cur_index)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(IVOC9ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)

               ! AqOrg

               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA1ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA2ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA3ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               !write(*,*)'AqPOA3 =',tracername(AqPOA3ind+cur_index)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA4ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA5ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA6ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA7ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqPOA8ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqLEVOind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqCBIOind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqCPD1ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqCPD2ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqCPD3ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC1ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC2ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC3ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC4ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC5ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               !write(*,*)'AqIVOC5 =',tracername(AqIVOC5ind+cur_index)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC6ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC7ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC8ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)
               TotOA(i,j,l) = TotOA(i,j,l) + (tracer(i,j,l,(AqIVOC9ind+cur_index))*(airdenstm/Av*airmw)*1e12) ! kg/kg_air to ug/m3)

            end do ! end n

            !if (TotOA(i,j,l).gt.2) then
               !write(*,*)nstep, 'TotOA(',i,j,l,') = ', TotOA(i,j,l)
               !do binnum=1,10

               !cur_index = 86*(binnum-1)

               !write(*,*) 'Org:'        

               !write(*,*) tracername(POA1ind+cur_index),(tracer(i,j,l,(POA1ind+cur_index))*(airdenstm/Av*airmw)*1e12) 

               !write(*,*) tracername(POA2ind+cur_index), (tracer(i,j,l,(POA2ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(POA3ind+cur_index),(tracer(i,j,l,(POA3ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(POA4ind+cur_index),(tracer(i,j,l,(POA4ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(POA5ind+cur_index),(tracer(i,j,l,(POA5ind+cur_index))*(airdenstm/Av*airmw)*1e12) 

               !write(*,*) tracername(POA6ind+cur_index),(tracer(i,j,l,(POA6ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(POA7ind+cur_index),(tracer(i,j,l,(POA7ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(POA8ind+cur_index),(tracer(i,j,l,(POA8ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(LEVOind+cur_index),(tracer(i,j,l,(LEVOind+cur_index))*(airdenstm/Av*airmw)*1e12) 

               !write(*,*) tracername(CBIOind+cur_index),(tracer(i,j,l,(CBIOind+cur_index))*(airdenstm/Av*airmw)*1e12) 

               !write(*,*) tracername(CPD1ind+cur_index), (tracer(i,j,l,(CPD1ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(CPD2ind+cur_index),(tracer(i,j,l,(CPD2ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(CPD3ind+cur_index), (tracer(i,j,l,(CPD3ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC1ind+cur_index), (tracer(i,j,l,(IVOC1ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC2ind+cur_index), (tracer(i,j,l,(IVOC2ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC3ind+cur_index),(tracer(i,j,l,(IVOC3ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC4ind+cur_index),(tracer(i,j,l,(IVOC4ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC5ind+cur_index),(tracer(i,j,l,(IVOC5ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC6ind+cur_index),(tracer(i,j,l,(IVOC6ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC7ind+cur_index), (tracer(i,j,l,(IVOC7ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC8ind+cur_index),(tracer(i,j,l,(IVOC8ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(IVOC9ind+cur_index),(tracer(i,j,l,(IVOC9ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) 'AqOrg:'

               !write(*,*) tracername(AqPOA1ind+cur_index),(tracer(i,j,l,(AqPOA1ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA2ind+cur_index),(tracer(i,j,l,(AqPOA2ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA3ind+cur_index),(tracer(i,j,l,(AqPOA3ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA4ind+cur_index),(tracer(i,j,l,(AqPOA4ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA5ind+cur_index),(tracer(i,j,l,(AqPOA5ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA6ind+cur_index),(tracer(i,j,l,(AqPOA6ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA7ind+cur_index),(tracer(i,j,l,(AqPOA7ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqPOA8ind+cur_index),(tracer(i,j,l,(AqPOA8ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqLEVOind+cur_index),(tracer(i,j,l,(AqLEVOind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqCBIOind+cur_index),(tracer(i,j,l,(AqCBIOind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqCPD1ind+cur_index),(tracer(i,j,l,(AqCPD1ind+cur_index))*(airdenstm/Av*airmw)*1e12) 

               !write(*,*) tracername(AqCPD2ind+cur_index),(tracer(i,j,l, (AqCPD2ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqCPD3ind+cur_index),(tracer(i,j,l, (AqCPD3ind+cur_index))*(airdenstm/Av*airmw)*1e12) 

               !write(*,*) tracername(AqIVOC1ind+cur_index),(tracer(i,j,l,(AqIVOC1ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC2ind+cur_index),(tracer(i,j,l,(AqIVOC2ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC3ind+cur_index),(tracer(i,j,l,(AqIVOC3ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC4ind+cur_index),(tracer(i,j,l,(AqIVOC4ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC5ind+cur_index),(tracer(i,j,l,(AqIVOC5ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC6ind+cur_index),(tracer(i,j,l,(AqIVOC6ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC7ind+cur_index),(tracer(i,j,l,(AqIVOC7ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC8ind+cur_index),(tracer(i,j,l,(AqIVOC8ind+cur_index))*(airdenstm/Av*airmw)*1e12)

               !write(*,*) tracername(AqIVOC9ind+cur_index),(tracer(i,j,l,(AqIVOC9ind+cur_index))*(airdenstm/Av*airmw)*1e12)
               !end do
            !endif

         end do	! end nz	
       end do! end ny
     end do! end nx!

name='TotOA'
long_name='TotalOA'
units = 'ug/m3'
call compress3D(TotOA,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


!TotK
 nfields1=nfields1+1

     !write(*,*)'TotK = '
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotK(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) +p(i,j,l)! in Pa  
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotK(i,j,l) = TotK(i,j,l) + ((tracer(i,j,l,(Kind+cur_index))/MwK)*airmw*1e-3)! kg/kg_air to mol/mol)
               TotK(i,j,l) = TotK(i,j,l) + ((tracer(i,j,l,(KNO3ind+cur_index))/MwKNO3)*airmw*1e-3)! kg/kg_air to mol/mol)
               TotK(i,j,l) = TotK(i,j,l) + ((tracer(i,j,l,(KClind+cur_index))/MwKCl)*airmw*1e-3)! kg/kg_air to mol/mol)
               TotK(i,j,l) = TotK(i,j,l) + ((tracer(i,j,l,(KHSO4ind+cur_index))/MwKHSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               TotK(i,j,l) = TotK(i,j,l) + ((tracer(i,j,l,(KOHind+cur_index))/MwKOH)*airmw*1e-3)! kg/kg_air to mol/mol)
               TotK(i,j,l) = TotK(i,j,l) + (2*(tracer(i,j,l,(K2SO4ind+cur_index))/MwK2SO4)*airmw*1e-3)! kg/kg_air to mol/mol)
            end do ! end binnum
            
            TotK(i,j,l) = TotK(i,j,l)*MwK*(airdenstm/Av)*1e15 ! mol/mol to ug/m3
          end do
        end do
      end do

name='TotK'
long_name='TotalK'
units = 'ug/m3'
call compress3D(TotK,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


! TotNA
nfields1=nfields1+1
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotNa(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotNa(i,j,l) = TotNa(i,j,l) + ((tracer(i,j,l,(Naind+cur_index))/MwNa)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNa(i,j,l) = TotNa(i,j,l) + ((tracer(i,j,l,(NaNO3ind+cur_index))/MwNaNO3)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNa(i,j,l) = TotNa(i,j,l) + ((tracer(i,j,l,(NaClind+cur_index))/MwNaCl)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNa(i,j,l) = TotNa(i,j,l) + ((tracer(i,j,l,(NaHSO4ind+cur_index))/MwNaHSO4)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNa(i,j,l) = TotNa(i,j,l) + ((tracer(i,j,l,(NaOHind+cur_index))/MwNaOH)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNa(i,j,l) = TotNa(i,j,l) + (2*(tracer(i,j,l,(Na2SO4ind+cur_index))/MwNa2SO4)*airmw*1e-3)! kg/kg_air to mol/mol)

            end do ! end binnum
            
            TotNa(i,j,l) = TotNa(i,j,l)*MwNa*(airdenstm/Av)*1e15
          end do
        end do
      end do

name='TotNa'
long_name='TotNa'
units = 'ug/m3'
call compress3D(TotNa,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


!TotNH4
nfields1=nfields1+1
         !write(*,*) 'TotNH4 = '
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotNH4(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               !TotNH4(i,j,l) = TotNH4(i,j,l) + ((tracer(i,j,l,(NH3ind+cur_index))/MwNH3)*airmw*1e-3)! kg/kg_air to mol/mol)
               TotNH4(i,j,l) = TotNH4(i,j,l) + ((tracer(i,j,l,(NH3ind+cur_index))/MwNH3)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(NH3ind+cur_index)
               TotNH4(i,j,l) = TotNH4(i,j,l) + ((tracer(i,j,l,(NH4ind+cur_index))/MwNH4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(NH4ind+cur_index)
               TotNH4(i,j,l) = TotNH4(i,j,l) + ((tracer(i,j,l,(NH4NO3ind+cur_index))/MwNH4NO3)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(NH4NO3ind+cur_index)
               TotNH4(i,j,l) = TotNH4(i,j,l) + ((tracer(i,j,l,(NH4Clind+cur_index))/MwNH4Cl)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(NH4Clind+cur_index)
               TotNH4(i,j,l) = TotNH4(i,j,l) + ((tracer(i,j,l,(NH4HSO4ind+cur_index))/MwNH4HSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(NH4HSO4ind+cur_index)
               TotNH4(i,j,l) = TotNH4(i,j,l) + (2*(tracer(i,j,l,(NH42SO4ind+cur_index))/MwNH42SO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(NH42SO4ind+cur_index)
               TotNH4(i,j,l) = TotNH4(i,j,l) + (3*(tracer(i,j,l,(LEVind+cur_index))/MwLEV)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*)tracername(LEVind+cur_index)
            end do ! end binnum
            
            TotNH4(i,j,l) = TotNH4(i,j,l)*MwNH4*(airdenstm/Av)*1e15
          end do
        end do
      end do
name='TotNH4'
long_name='TotNH4'
units = 'ug/m3'
call compress3D(TotNH4,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

!TotCl-
nfields1=nfields1+1
     do i=1,nx
       do j=1,ny
         do l=1,nzm
         !write(*,*)'TotCl  ='
            TotCl(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotCl(i,j,l) = TotCl(i,j,l) + ((tracer(i,j,l,(Clind+cur_index))/MwCl)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCl(i,j,l) = TotCl(i,j,l) + ((tracer(i,j,l,(NH4Clind+cur_index))/MwNH4Cl)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCl(i,j,l) = TotCl(i,j,l) + ((tracer(i,j,l,(KClind+cur_index))/MwKCl)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCl(i,j,l) = TotCl(i,j,l) + ((tracer(i,j,l,(NaClind+cur_index))/MwNaCl)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCl(i,j,l) = TotCl(i,j,l) + (2*(tracer(i,j,l,(MgCl2ind+cur_index))/MwMgCl2)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCl(i,j,l) = TotCl(i,j,l) + (2*(tracer(i,j,l,(CaCl2ind+cur_index))/MwCaCl2)*airmw*1e-3)! kg/kg_air to mol/mol)

            end do ! end binnum
            
            TotCl(i,j,l) = TotCl(i,j,l)*MwCl*(airdenstm/Av)*1e15
          end do
        end do
      end do
name='TotCl'
long_name='TotCl'
units = 'ug/m3'
call compress3D(TotCl,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

!TotNO3-    
nfields1=nfields1+1
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotNO3(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotNO3(i,j,l) = TotNO3(i,j,l) + ((tracer(i,j,l,(NO3ind+cur_index))/MwNO3)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNO3(i,j,l) = TotNO3(i,j,l) + ((tracer(i,j,l,(NH4NO3ind+cur_index))/MwNH4NO3)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNO3(i,j,l) = TotNO3(i,j,l) + ((tracer(i,j,l,(KNO3ind+cur_index))/MwKNO3)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNO3(i,j,l) = TotNO3(i,j,l) + ((tracer(i,j,l,(NaNO3ind+cur_index))/MwNaNO3)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNO3(i,j,l) = TotNO3(i,j,l) + (2*(tracer(i,j,l,(MgNO32ind+cur_index))/MwMgNO32)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotNO3(i,j,l) = TotNO3(i,j,l) + (2*(tracer(i,j,l,(CaNO32ind+cur_index))/MwCaNO32)*airmw*1e-3)! kg/kg_air to mol/mol)

            end do ! end binnum
            
            TotNO3(i,j,l) = TotNO3(i,j,l)*MwNO3*(airdenstm/Av)*1e15
          end do
        end do
      end do
name='TotNO3'
long_name='TotNO3'
units = 'ug/m3'
call compress3D(TotNO3,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


!TotSO4-
nfields1=nfields1+1
     !write(*,*) 'TotSO4 = '
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotSO4(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0) + p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(SO4ind+cur_index))/MwSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(SO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(NH42SO4ind+cur_index))/MwNH42SO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(NH42SO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(K2SO4ind+cur_index))/MwK2SO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(K2SO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(Na2SO4ind+cur_index))/MwNa2SO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(Na2SO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(HSO4ind+cur_index))/MwHSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(HSO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(NH4HSO4ind+cur_index))/MwNH4HSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(NH4HSO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(KHSO4ind+cur_index))/MwKHSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(KHSO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(NaHSO4ind+cur_index))/MwNaHSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(NaHSO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(CaSO4ind+cur_index))/MwCaSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(CaSO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + ((tracer(i,j,l,(MgSO4ind+cur_index))/MwMgSO4)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(MgSO4ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + (2*(tracer(i,j,l,(LEVind+cur_index))/MwLEV)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(LEVind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + (2*(tracer(i,j,l,(MgHSO42ind+cur_index))/MwMgHSO42)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(MgHSO42ind+cur_index)
               TotSO4(i,j,l) = TotSO4(i,j,l) + (2*(tracer(i,j,l,(CaHSO42ind+cur_index))/MwCaHSO42)*airmw*1e-3)! kg/kg_air to mol/mol)
               !write(*,*) tracername(CaHSO42ind+cur_index)
            end do ! end binnum
            
            TotSO4(i,j,l) = TotSO4(i,j,l)*MwSO4*(airdenstm/Av)*1e15
          end do
        end do
      end do
name='TotSO4'
long_name='TotSO4'
units = 'ug/m3'
call compress3D(TotSO4,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

!TotMg
nfields1=nfields1+1
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotMg(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0)+p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotMg(i,j,l) = TotMg(i,j,l) + ((tracer(i,j,l,(Mgind+cur_index))/MwMg)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotMg(i,j,l) = TotMg(i,j,l) + ((tracer(i,j,l,(MgNO32ind+cur_index))/MwMgNO32)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotMg(i,j,l) = TotMg(i,j,l) + ((tracer(i,j,l,(MgCl2ind+cur_index))/MwMgCl2)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotMg(i,j,l) = TotMg(i,j,l) + ((tracer(i,j,l,(MgHSO42ind+cur_index))/MwMgHSO42)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotMg(i,j,l) = TotMg(i,j,l) + ((tracer(i,j,l,(MgOH2ind+cur_index))/MwMgOH2)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotMg(i,j,l) = TotMg(i,j,l) + ((tracer(i,j,l,(MgSO4ind+cur_index))/MwMgSO4)*airmw*1e-3)! kg/kg_air to mol/mol)

            end do ! end binnum
            
            TotMg(i,j,l) = TotMg(i,j,l)*MwMg*(airdenstm/Av)*1e15
          end do
        end do
      end do
name='TotMg'
long_name='TotMg'
units = 'ug/m3'
call compress3D(TotMg,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

!TotCa
nfields1=nfields1+1
     do i=1,nx
       do j=1,ny
         do l=1,nzm
            TotCa(i,j,l) = 0.0
            prestm = dble(pres(l)*100.d0)+p(i,j,l) ! in Pa       
            temptm = dble(tabs(i,j,l)) ! in K
	    airdenstm = prestm*Av/(gasc*temptm)/1.0d6! in molec/cm3
            !airmw g/mol
            do binnum=1,10 ! starting aerosol species after number of gas species 
               cur_index = 86*(binnum-1)
               TotCa(i,j,l) = TotCa(i,j,l) + ((tracer(i,j,l,(Caind+cur_index))/MwCa)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCa(i,j,l) = TotCa(i,j,l) + ((tracer(i,j,l,(CaNO32ind+cur_index))/MwCaNO32)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCa(i,j,l) = TotCa(i,j,l) + ((tracer(i,j,l,(CaCl2ind+cur_index))/MwCaCl2)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCa(i,j,l) = TotCa(i,j,l) + ((tracer(i,j,l,(CaHSO42ind+cur_index))/MwCaHSO42)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCa(i,j,l) = TotCa(i,j,l) + ((tracer(i,j,l,(CaOH2ind+cur_index))/MwCaOH2)*airmw*1e-3)! kg/kg_air to mol/mol)

               TotCa(i,j,l) = TotCa(i,j,l) + ((tracer(i,j,l,(CaSO4ind+cur_index))/MwCaSO4)*airmw*1e-3)! kg/kg_air to mol/mol)

            end do ! end binnum
            
            TotCa(i,j,l) = TotCa(i,j,l)*MwCa*(airdenstm/Av)*1e15
          end do
        end do
      end do
name='TotCa'
long_name='TotCa'
units = 'ug/m3'
call compress3D(TotCa,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


! Number of Particles (No...)

do n=1477,1486
  nfields1=nfields1+1
   do i=1,nx
      do j=1,ny
         do l=1,nzm!


            prestm = dble(pres(l)*100.d0)+p(i,j,l) ! in Pa   
            temptm = dble(tabs(i,j,l)) ! in K
            airdenstm = prestm*Av/(gasc*temptm)/1e6 ! in molec/cm3

            tmp(i,j,l) = tracer(i,j,l,n)*(airdenstm/Av*airmw*1e-3) ! use for values averaged of nsteps3D (in kg)
            std(i,j,l) = tracer(i,j,l,n)*(airdenstm/Av*airmw*1e-3)
            !write(*,*) 'In Write3D Aerosol', tracername(n),tmp(i,j,l)

    end do
   end do
  end do
  name=tracername(n)
  long_name=tracername(n)
  units = '#/cm3'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
end do


do n = 1,nmicro_fields
   if(docloud.AND.flag_micro3Dout(n).gt.0.AND.n.ne.index_water_vapor) then
      nfields1=nfields1+1
      do k=1,nzm
         do j=1,ny
            do i=1,nx
               tmp(i,j,k)=micro_field(i,j,k,n)*mkoutputscale(n)
            end do
         end do
         ! remove factor of rho from number, if this field is a number concentration
         if(flag_number(n).gt.0) tmp(:,:,k) = tmp(:,:,k)*rho(k)
      end do
      name=TRIM(mkname(n))
      long_name=TRIM(mklongname(n))
      units=TRIM(mkunits(n))
      call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
           save3Dbin,dompi,rank,nsubdomains)
   end if
end do

  call task_barrier()

  if(nfields.ne.nfields1) then
    if(masterproc) print*,'write_fields3D error: nfields=',nfields,'  nfields1=',nfields1
    call task_abort()
  end if
  if(masterproc) then
    close (46)
    if(RUN3D.or.save3Dsep) then
       if(dogzip3D) call systemf('gzip -f '//filename)
       print*, 'Writting 3D data. file:'//filename
    else
       print*, 'Appending 3D data. file:'//filename
    end if
  endif

 
end
