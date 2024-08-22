subroutine tracer_bound_noperiodic(f,dimx1,dimx2,dimy1,dimy2,dimz,i_1, i_2, j_1, j_2, id)
	
! if tracerperiodic is set to false, the tracers are not wrapped using periodic
! boundary conditions. Instead the tracers at the boundary are set equal
! to the tracers just before the boundary.

use grid
use tracers
implicit none
	
integer dimx1, dimx2, dimy1, dimy2, dimz
integer i_1, i_2, j_1, j_2
real f(dimx1:dimx2, dimy1:dimy2, dimz)
integer id   ! tracer number - no longer used
	
real buffer((nx+ny)*3*nz)	! buffer for sending data
	
integer i, j, k, n
integer i1, i2, j1, j2
	
i1 = i_1 - 1
i2 = i_2 - 1
j1 = j_1 - 1
j2 = j_2 - 1

!----------------------------------------------------------------------
!  Send buffers to neighbors
!----------------------------------------------------------------------


	if(RUN3D) then

! "South":	

	     do k=1,dimz
	       do j=-j1,0
	         do i=1,nx
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(i,1,k)
	         end do
	       end do
	     end do

! "South-West":	
	
	     do k=1,dimz
	       do j=-j1,0
	         do i=-i1,0
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(1,1,k)
	         end do
	       end do
	     end do

! "North-West":

	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=-i1,0
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(1,ny,k)
	         end do
	       end do
	     end do

! "North":

	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=1,nx
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(i,ny,k)
	         end do
	       end do
	     end do

! "North-East":
	  
	     do k=1,dimz
	       do j=nyp1,nyp1+j2
	         do i=nxp1,nxp1+i2
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(nx,ny,k)
	         end do
	       end do
	     end do


! "South-East":
	  	  
	     do k=1,dimz
	       do j=-j1,0
	         do i=nxp1,nxp1+i2
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(nx,1,k)
	         end do
	       end do
	     end do
	     

	endif

! "West":
 	     !write(*,*) 'In tracer_bound_noperiodic.f90'
	     do k=1,dimz
	       do j=1,ny
	         do i=-i1,0
	           !f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(1,j,k)
	         end do
	       end do
	     end do

! "East":

	     do k=1,dimz
	       do j=1,ny
	         do i=nxp1,nxp1+i2
!	           f(i,j,k) = tracerbound(k,id)
	           f(i,j,k) = f(nx,j,k)
	         end do
	       end do
	     end do


end subroutine tracer_bound_noperiodic
	     
	     
