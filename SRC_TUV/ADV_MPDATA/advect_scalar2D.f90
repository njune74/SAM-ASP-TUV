
subroutine advect_scalar2D (f, u, w, rho, rhow, flux, adtype, sbin)
 	
!     positively definite monotonic advection with non-oscillatory option

use grid
use params, only: dowallx
use tracers, only: tbulk, tnum, tpmass, tdmass, fluxsave, fnumsave
implicit none


real f(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real u(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
real w(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )
real rho(nzm)
real rhow(nz)
real flux(nz)
integer adtype, sbin ! added by jrp
                     ! adtype determines if tracer is aerosol number, mass or
                     ! other
                     ! sbin is which aerosol bin (use a value of 1 if not
                     ! aerosol)
	
real mx (0:nxp1,1,nzm)
real mn (0:nxp1,1,nzm)
real uuu(-1:nxp3,1,nzm)
real www(-1:nxp2,1,nz)

real eps, dd
integer i,j,k,ic,ib,kc,kb,h
logical nonos
real iadz(nzm),irho(nzm),irhow(nzm)
real testarray1(7) ! added by jrp for TOMAS aerosol advection
real testarray2(8) ! added by jrp for TOMAS aerosol advection
integer tmpv ! temporary variable
real meanmass

real x1, x2, a, b, a1, a2, y
real andiff,across,pp,pn
andiff(x1,x2,a,b)=(abs(a)-a*a*b)*0.5*(x2-x1)
across(x1,a1,a2)=0.03125*a1*a2*x1
pp(y)= max(0.,y)
pn(y)=-min(0.,y)
	
nonos = .true.
eps = 1.e-10

j=1

www(:,:,nz)=0.
!write(*,*)'in advect_scalar2D'
! check for NaNs before advecting
!write(*,*), "f:" ,f!AD DEBUG
do k=1,nzm
  do i=1,nx
     if(isnan(f(i,j,k)))then
      print*,'NaN before advection in advect_scalar2D'
      if     (adtype.eq.tbulk)  then
         print*,'bulk tracer'
      elseif (adtype.eq.tnum)   then
         print*,'aerosol number tracer'
      elseif (adtype.eq.tpmass) then
         print*,'aerosol prognostic mass tracer'
      elseif (adtype.eq.tdmass) then
         print*,'aerosol diagnostic mass tracer'
      else
         print*,'unknown tracer'
      endif
      print*,'i,j,k,sizebin',i,j,k,sbin
      stop
     endif
  enddo
enddo

if(dowallx) then

  if(mod(rank,nsubdomains_x).eq.0) then
    do k=1,nzm
       do i=dimx1_u,1
         u(i,j,k) = 0.
       end do
    end do
  end if
  if(mod(rank,nsubdomains_x).eq.nsubdomains_x-1) then
    do k=1,nzm
       do i=nx+1,dimx2_u
         u(i,j,k) = 0.
       end do
    end do
  end if

end if

!-----------------------------------------
if(adtype.lt.tpmass)then ! do treatment for bulk and aerosol number
!WRITE(*,*) ' adtype.lt.tpmass', adtype,tpmass
 
 if(nonos) then

 do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  do i=0,nxp1
    ib=i-1
    ic=i+1
    mx(i,j,k)=max(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k))
    mn(i,j,k)=min(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k))
   end do
 end do

end if  ! nonos

  if(adtype.eq.tnum)then ! aerosol number tracer
  !WRITE(*,*) ' adtype.eq.tnum', adtype,tnum
   do k=1,nzm
     do i=-2,nxp3
      fnumsave(i,j,k,1,sbin)=f(i,j,k)
     enddo
   enddo
  endif

do k=1,nzm
  kb=max(1,k-1)
  do i=-1,nxp3
    uuu(i,j,k)=max(0.,u(i,j,k))*f(i-1,j,k)+min(0.,u(i,j,k))*f(i,j,k)
    if(adtype.eq.tnum)then
      fluxsave(i,j,k,1,1,sbin)=uuu(i,j,k)
    endif
  end do
  do i=-1,nxp2
    www(i,j,k)=max(0.,w(i,j,k))*f(i,j,kb)+min(0.,w(i,j,k))*f(i,j,k)
    if(adtype.eq.tnum)then
      fluxsave(i,j,k,1,3,sbin)=www(i,j,k)
    endif
  end do
  flux(k) = 0.
  do i=1,nx
    flux(k) = flux(k) + www(i,j,k)	
  end do
end do

do k=1,nzm
  irho(k) = 1./rho(k)
  iadz(k) = 1./adz(k)
   do i=-1,nxp2
      f(i,j,k) = f(i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) & 
                        + (www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)  
      if(adtype.eq.tnum)then
        fnumsave(i,j,k,2,sbin)=f(i,j,k)
      endif          
   end do
end do 


do k=1,nzm
  kc=min(nzm,k+1)
  kb=max(1,k-1)
  dd=2./(kc-kb)/adz(k)
  irhow(k)=1./(rhow(k)*adz(k))
  do i=0,nxp2
   ib=i-1
   uuu(i,j,k)=andiff(f(ib,j,k),f(i,j,k),u(i,j,k),irho(k)) &
      - across(dd*(f(ib,j,kc)+f(i,j,kc)-f(ib,j,kb)-f(i,j,kb)), &
              u(i,j,k), w(ib,j,k)+w(ib,j,kc)+w(i,j,k)+w(i,j,kc)) *irho(k)
   if(adtype.eq.tnum)then
      fluxsave(i,j,k,2,1,sbin)=uuu(i,j,k)
   endif

  end do
          

  do i=0,nxp1
   ib=i-1
   ic=i+1
   www(i,j,k)=andiff(f(i,j,kb),f(i,j,k),w(i,j,k),irhow(k)) &
      -across(f(ic,j,kb)+f(ic,j,k)-f(ib,j,kb)-f(ib,j,k), &
        w(i,j,k), u(i,j,kb)+u(i,j,k)+u(ic,j,k)+u(ic,j,kb)) *irho(k)
   if(adtype.eq.tnum)then
      fluxsave(i,j,k,2,3,sbin)=www(i,j,k)
   endif
  end do
end do
www(:,:,1) = 0.
!---------- non-osscilatory option ---------------

if(nonos) then

 do k=1,nzm
   kc=min(nzm,k+1)
   kb=max(1,k-1)
   do i=0,nxp1
    ib=i-1
    ic=i+1
    mx(i,j,k)=max(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mx(i,j,k))
    mn(i,j,k)=min(f(ib,j,k),f(ic,j,k),f(i,j,kb),f(i,j,kc),f(i,j,k),mn(i,j,k))
   end do
 end do

 do k=1,nzm
   kc=min(nzm,k+1)
   do i=0,nxp1
    ic=i+1
     mx(i,j,k)=rho(k)*(mx(i,j,k)-f(i,j,k))/(pn(uuu(ic,j,k)) + pp(uuu(i,j,k))+&
               iadz(k)*(pn(www(i,j,kc)) + pp(www(i,j,k)))+eps)	
     mn(i,j,k)=rho(k)*(f(i,j,k)-mn(i,j,k))/(pp(uuu(ic,j,k)) + pn(uuu(i,j,k))+&
               iadz(k)*(pp(www(i,j,kc)) + pn(www(i,j,k)))+eps)	
   end do
 end do

 do k=1,nzm
  kb=max(1,k-1)
   do i=1,nxp1
    ib=i-1
    uuu(i,j,k)= pp(uuu(i,j,k))*min(1.,mx(i,j,k), mn(ib,j,k)) &
              - pn(uuu(i,j,k))*min(1.,mx(ib,j,k),mn(i,j,k))
    if(adtype.eq.tnum)then
      fluxsave(i,j,k,2,1,sbin)=uuu(i,j,k)
    endif
   end do
   do i=1,nx
    www(i,j,k)= pp(www(i,j,k))*min(1.,mx(i,j,k), mn(i,j,kb)) &
              - pn(www(i,j,k))*min(1.,mx(i,j,kb),mn(i,j,k))
    if(adtype.eq.tnum)then
      fluxsave(i,j,k,2,3,sbin)=www(i,j,k)
    endif
    flux(k) = flux(k) + www(i,j,k)	
   end do
 end do


endif ! nonos


 do k=1,nzm
  kc=k+1
   do i=1,nx
 ! MK: added fix for very small negative values (relative to positive values) 
 !     especially  when such large numbers as
 !     hydrometeor concentrations are advected. The reason for negative values is
 !     most likely truncation error.
    f(i,j,k)= max(0., f(i,j,k) - (uuu(i+1,j,k)-uuu(i,j,k) &
                     +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k))
   end do
 end do 

else ! is an aerosol mass tracer

 ! do first pass at change in concentration
 do k=1,nzm
   do i=-1,nxp3
    if(fluxsave(i,j,k,1,1,sbin).gt.0.)then
     meanmass=f(i-1,j,k)/fnumsave(i-1,j,k,1,sbin)
    elseif(fluxsave(i,j,k,1,1,sbin).eq.0.)then
     meanmass=0.
    else
     meanmass=f(i,j,k)/fnumsave(i,j,k,1,sbin)
    endif
!    if((adtype.eq.tpmass).and.(j.ge.1).and.(j.le.ny) &
!       .and.(i.ge.1).and.(i.le.nx))then
!    if(meanmass.eq.0.)then
!    elseif(meanmass.lt.xk(sbin))then
!     print*,'uuu 1 meanmass low in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin)
!     print*,'f1',f(i-1,j,k),'fnumsave1',fnumsave(i-1,j,k,1,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,1,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,1,1,sbin)
!    elseif(meanmass.gt.xk(sbin+1))then
!     print*,'uuu 1 meanmass high in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin+1)
!     print*,'f1',f(i-1,j,k),'fnumsave1',fnumsave(i-1,j,k,1,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,1,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,1,1,sbin)
!    endif 
!    endif 
    uuu(i,j,k)=fluxsave(i,j,k,1,1,sbin)*meanmass
   enddo
 enddo

 do k=1,nzm
  kb=max(1,k-1)
   do i=-1,nxp2
    if(fluxsave(i,j,k,1,3,sbin).gt.0.)then
     meanmass=f(i,j,kb)/fnumsave(i,j,kb,1,sbin)
    elseif(fluxsave(i,j,k,1,3,sbin).eq.0.)then
     meanmass=0.
    else
     meanmass=f(i,j,k)/fnumsave(i,j,k,1,sbin)
    endif
!    if((adtype.eq.tpmass).and.(j.ge.1).and.(j.le.ny) &
!       .and.(i.ge.1).and.(i.le.nx))then
!    if(meanmass.eq.0.)then
!    elseif(meanmass.lt.xk(sbin))then
!     print*,'www 1 meanmass low in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin)
!     print*,'f1',f(i,j,kb),'fnumsave1',fnumsave(i,j,kb,1,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,1,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,1,3,sbin)
!    elseif(meanmass.gt.xk(sbin+1))then
!     print*,'www 1 meanmass high in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin+1)
!     print*,'f1',f(i,j,kb),'fnumsave1',fnumsave(i,j,kb,1,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,1,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,1,3,sbin)
!    endif 
!    endif 
    www(i,j,k)=fluxsave(i,j,k,1,3,sbin)*meanmass
   enddo
 enddo

 do k=1,nzm
  irho(k) = 1./rho(k)
  iadz(k) = 1./adz(k)
   do i=-1,nxp2
    f(i,j,k)=f(i,j,k) -(uuu(i+1,j,k)-uuu(i,j,k) &
             +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)
    if(isnan(f(i,j,k)))then
      print*,'NaN in advect_scalar3D'
      print*,'1st f mod'
      if(isnan(uuu(i+1,j,k))) print*,'uuu(i+1,j,k)'
      if(isnan(uuu(i,j,k))) print*,'uuu(i,j,k)'
      if(isnan(www(i,j,k+1))) print*,'www(i,j,k+1)'
      if(isnan(www(i,j,k))) print*,'www(i,j,k)'
      if(isnan(iadz(k))) print*,'iadz(k)'
      if(isnan(irho(k))) print*,'irho(k)'
      if     (adtype.eq.tbulk)  then
         print*,'bulk tracer'
      elseif (adtype.eq.tnum)   then
         print*,'aerosol number tracer'
      elseif (adtype.eq.tpmass) then
         print*,'aerosol prognostic mass tracer'
      elseif (adtype.eq.tdmass) then
         print*,'aerosol diagnostic mass tracer'
      else
         print*,'unknown tracer'
      endif
      print*,'i,j,k,sizebin',i,j,k,sbin
      print*,'fluxsave',fluxsave(i,j,k,1,:,sbin)
      print*,'f,fnumsave'
      print*,'in cell',f(i,j,k),fnumsave(i,j,k,1,sbin)
      print*,'i-1',f(i-1,j,k),fnumsave(i-1,j,k,1,sbin)
      if (k.gt.1) print*,'k-1',f(i,j,k-1),fnumsave(i,j,kb,1,sbin)
      stop
    endif
   end do
 end do

 ! do second pass at changing concentrations

 do k=1,nzm
   do i=1,nxp1
   ib=i-1
    if(fluxsave(i,j,k,2,1,sbin).gt.0.)then
     meanmass=f(ib,j,k)/fnumsave(ib,j,k,2,sbin)
    elseif(fluxsave(i,j,k,2,1,sbin).eq.0.)then
     meanmass=0.
    else
     meanmass=f(i,j,k)/fnumsave(i,j,k,2,sbin)
    endif
!    if(adtype.eq.tpmass)then
!    if(meanmass.eq.0.)then
!    elseif(meanmass.lt.xk(sbin))then
!     print*,'uuu 2 meanmass low in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin)
!     print*,'f1',f(i-1,j,k),'fnumsave1',fnumsave(i-1,j,k,2,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,2,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,2,1,sbin)
!    elseif(meanmass.gt.xk(sbin+1))then
!     print*,'uuu 2 meanmass high in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin+1)
!     print*,'f1',f(i-1,j,k),'fnumsave1',fnumsave(i-1,j,k,2,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,2,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,2,1,sbin)
!    endif 
!    endif 
    uuu(i,j,k)=fluxsave(i,j,k,2,1,sbin)*meanmass
   enddo
 enddo
 
do k=1,nzm
  kb=max(1,k-1)
   do i=1,nx
    if(fluxsave(i,j,k,2,3,sbin).gt.0.)then
     meanmass=f(i,j,kb)/fnumsave(i,j,kb,2,sbin)
    elseif(fluxsave(i,j,k,2,3,sbin).eq.0.)then
     meanmass=0.
    else
     meanmass=f(i,j,k)/fnumsave(i,j,k,2,sbin)
    endif
!    if(adtype.eq.tpmass)then
!    if(meanmass.eq.0.)then
!    elseif(meanmass.lt.xk(sbin))then
!     print*,'www 2 meanmass low in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin)
!     print*,'f1',f(i,j,kb),'fnumsave1',fnumsave(i,j,kb,2,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,2,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,2,3,sbin)
!    elseif(meanmass.gt.xk(sbin+1))then
!     print*,'www 2 meanmass high in bin',sbin
!     print*,'i',i,'j',j,'k',k
!     print*,'meanmass',meanmass,'xk(sbin)',xk(sbin+1)
!     print*,'f1',f(i,j,kb),'fnumsave1',fnumsave(i,j,kb,2,sbin)
!     print*,'f2',f(i,j,k),'fnumsave2',fnumsave(i,j,k,2,sbin)
!     print*,'fluxsave',fluxsave(i,j,k,2,3,sbin)
!    endif 
!    endif 
    www(i,j,k)=fluxsave(i,j,k,2,3,sbin)*meanmass
   enddo
 enddo

 do k=1,nzm
  kc=k+1
   do i=1,nx
    f(i,j,k)=f(i,j,k) -(uuu(i+1,j,k)-uuu(i,j,k) &
                      +(www(i,j,k+1)-www(i,j,k))*iadz(k))*irho(k)

    if(isnan(f(i,j,k)))then
      print*,'NaN in advect_scalar3D'
      print*,'2nd f mod'
      if(isnan(uuu(i+1,j,k))) print*,'uuu(i+1,j,k)'
      if(isnan(uuu(i,j,k))) print*,'uuu(i,j,k)'
      if(isnan(www(i,j,k+1))) print*,'www(i,j,k+1)'
      if(isnan(www(i,j,k))) print*,'www(i,j,k)'
      if(isnan(iadz(k))) print*,'iadz(k)'
      if(isnan(irho(k))) print*,'irho(k)'
      if     (adtype.eq.tbulk)  then
         print*,'bulk tracer'
      elseif (adtype.eq.tnum)   then
         print*,'aerosol number tracer'
      elseif (adtype.eq.tpmass) then
         print*,'aerosol prognostic mass tracer'
      elseif (adtype.eq.tdmass) then
         print*,'aerosol diagnostic mass tracer'
      else
         print*,'unknown tracer'
      endif
      print*,'i,j,k,sizebin',i,j,k,sbin
      print*,'fluxsave',fluxsave(i,j,k,2,:,sbin)
      print*,'f,fnumsave'
      print*,'in cell',f(i,j,k),fnumsave(i,j,k,2,sbin)
      print*,'i-1',f(i-1,j,k),fnumsave(i-1,j,k,2,sbin)
      kb=max(1,k-1)
      print*,'k-1',f(i,j,kb),fnumsave(i,j,kb,2,sbin)
      stop
    endif
   end do
 end do

endif

end subroutine advect_scalar2D


