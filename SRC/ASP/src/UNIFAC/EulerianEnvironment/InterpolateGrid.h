!! ASP (c), 2004-2012, Matt Alvarado (mjalvara@mit.edu)
!! Based on MELAM of H.D.Steele (c) 2000-2004
!!
!! File Description:
!! InterpolateGrid.h
!! This file interpolates values between grid points,
!! finds the physical location closest to a grid point, and
!! finds the grid point closest to a physical location.
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! UPDATE HISTORY							     !!
!!									     !!
!! Month  Year   Name              Description				     !!
!! 07     2006   Matt Alvarado     Began Update History			     !!
!! 02/15  2012   Matt Alvarado     Removing most of these grids, making ASP  !!
!!                                 a one-box model or subroutine. This       !!
!!                                 header file no longer needed              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This file contains the following functions and subroutines:
!! 1. FUNCTION InterpolateGrid (XP, YP, ZP, Grid, Index) RESULT (OutValue)
!! 2. FUNCTION LocateGridPoint(X, Y, Z, Index) RESULT (Location)
!! 3. SUBROUTINE ClosestGridPoint(XX, YY, ZZ, Index, X, Y, Z) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This subroutine interpolates a data value from one of the Eulerian Grids.	!!
!! In the first draft (July 2001), the interpolation is bilinear as per Press	!!
!! et al's Numerical Recipes (1992).  Which weights each of the 8 nearest grid	!!
!! points according to the volume not bounded (opposing corner) by each and the!!
!! free point.  See Press et al. Section 3.6 for further discussion.	!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL*8 FUNCTION InterpolateGrid (XP, YP, ZP, Grid, Index) RESULT (OutValue)
		
		use ModelParameters
		use InfrastructuralCode
		implicit none

		!! In-Out Variables
		real*8  :: XP, YP, ZP, Grid(:,:,:)
		integer :: Index					!! Index is tells which grid spacing is used

		!! Internal Variables
		real*8  :: widthX, widthY, widthZ	!! the distance between grid points in each direction
		integer :: ix, iy, iz				!! the closest grid point in each direction (lower)
		real*8  :: dx, dy, dz				!! the cartesian distance between each lower gridpoint and the free point
		integer :: dim						!! number of dimensions of interpolation
		logical :: doXYZ(3)					!! which of the dimensions to do

		!! SCSCSCSC
		logical :: Scaffolding = .FALSE. ! .TRUE.

		dim		 = 0
		doXYZ(1) = .FALSE.
		doXYZ(2) = .FALSE.
		doXYZ(3) = .FALSE.

		!! See if the Value is out of Range and stop the 
		!! integration if it is.
		if (XP > DomainX .or. XP < 0. .or.  &
			YP > DomainY .or. YP < 0. .or.  & 
			ZP > DomainZ .or. ZP < 0.) then

			CALL ERROR("Location (x,y,z)=("//TRIM(REAL2STR(XP/lengthscale))//","//TRIM(REAL2STR(YP/lengthscale))//","//   &
				   TRIM(REAL2STR(ZP/lengthscale))//") is out of range.  (Subroutine InterpolateGrid()).")
end if
		!! If this dimension is to be interpolated over, 
		!! then there is a bit of pre-processing
		IF (XGridPoints(Index) > 1.) THEN

			!! Determine the cell grid spacing
			widthX = DomainX / (XGridPoints(Index) - 1.) 

			!! Locate the free point within the grid
			!! The +1 enters because the first gridpoint on the border of the domain
			ix = FLOOR (XP / widthX) + 1.

			!! If the point is exactly on a gridpoint, then no need to interpolate there
			IF (MODULO(XP,widthX) .ne. 0) THEN
				dim      = dim + 1
				doXYZ(1) = .TRUE.

				!! Determine the distance to each of the closest gridpoints
				dx = XP - (ix - 1.) * widthX
			END IF
		ELSE
			ix = 1
		END IF

		IF (YGridPoints(Index) > 1.) THEN
			widthY = DomainY / (YGridPoints(Index) - 1.)
			iy = FLOOR (YP / widthY) + 1.
			IF (MODULO(YP,widthY) .ne. 0) THEN
				dim      = dim + 1
				doXYZ(2) = .TRUE.
				dy       = YP - (iy - 1.) * widthY
			END IF
		ELSE 
			iy = 1
		END IF

		IF (ZGridPoints(Index) > 1.) THEN
			widthZ = DomainZ / (ZGridPoints(Index) - 1.)
			iz = FLOOR (ZP / widthZ) + 1.
			IF (MODULO(ZP,widthZ) .ne. 0) THEN
				dim = dim + 1
				doXYZ(3) = .TRUE.
				dz = ZP - (iz - 1.) * widthZ
			END IF
		ELSE
			iz = 1
		END IF

		!! Now use these values to interpolate from the nearest gridpoints
		!! Each of the values is weighted by the volume defined (at the corners) 
		!! by the free point and the opposing gridpoint
		
		!!!!!!!!!!!!!!!!!!!!!!!
		!! 3-D INTERPOLATION !!
		IF (dim == 3) THEN
			OutValue =														&
			(Grid(  ix,  iy,  iz)*(widthX-dx)*(widthY-dy)*(widthZ-dz)+		&
			 Grid(  ix,  iy,iz+1)*(widthX-dx)*(widthY-dy)*dz+				&
			 Grid(  ix,iy+1,  iz)*(widthX-dx)*dy*(widthZ-dz)+				&
			 Grid(  ix,iy+1,iz+1)*(widthX-dx)*dy*dz+						&
			 Grid(ix+1,  iy,  iz)*dx*(widthY-dy)*(widthZ-dz)+				&
			 Grid(ix+1,  iy,iz+1)*dx*(widthY-dy)*dz+						&
			 Grid(ix+1,iy+1,  iz)*dx*dy*(widthZ-dz)+						&
			 Grid(ix+1,iy+1,iz+1)*dx*dy*dz) / widthX / widthY / widthZ
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!
		!! 2-D INTERPOLATION !!
		IF (dim == 2) THEN
		IF (.NOT.doXYZ(1) .AND. doXYZ(2) .AND. doXYZ(3)) THEN
			OutValue =										&
			(Grid(ix,  iy,  iz)*(widthY-dy)*(widthZ-dz)+	&
			 Grid(ix,  iy,iz+1)*(widthY-dy)*dz+				&
			 Grid(ix,iy+1,  iz)*dy*(widthZ-dz)+				&
			 Grid(ix,iy+1,iz+1)*dy*dz) / widthY / widthZ
		END IF
		IF (doXYZ(1) .AND. .NOT.doXYZ(2) .AND. doXYZ(3)) THEN
			OutValue =										&
			(Grid(  ix,iy,  iz)*(widthX-dx)*(widthZ-dz)+	&
			 Grid(  ix,iy,iz+1)*(widthX-dx)*dz+				&
			 Grid(ix+1,iy,  iz)*dx*(widthZ-dz)+				&
			 Grid(ix+1,iy,iz+1)*dx*dz) / widthX / widthZ
		END IF
		IF (doXYZ(1) .AND. doXYZ(2) .AND. .NOT.doXYZ(3)) THEN
			OutValue =										&
			(Grid(  ix,  iy,iz)*(widthX-dx)*(widthY-dy)+	&
			 Grid(  ix,iy+1,iz)*(widthX-dx)*dy+				&
			 Grid(ix+1,  iy,iz)*dx*(widthY-dy)+				&
			 Grid(ix+1,iy+1,iz)*dx*dy) / widthX / widthY
		END IF
		END IF

		!!!!!!!!!!!!!!!!!!!!!!!
		!! 1-D INTERPOLATION !!
		IF (dim == 1) THEN
		IF (doXYZ(1) .AND. .NOT.doXYZ(2) .AND. .NOT.doXYZ(3)) &
			OutValue = (Grid(ix,iy,iz)*(widthX-dx)+Grid(ix+1,iy,iz)*dx) / widthX	

		IF (.NOT.doXYZ(1) .AND. doXYZ(2) .AND. .NOT.doXYZ(3)) &
			OutValue = (Grid(ix,iy,iz)*(widthY-dy)+Grid(ix,iy+1,iz)*dy) / widthY

		IF (.NOT.doXYZ(1) .AND. .NOT.doXYZ(2) .AND. doXYZ(3)) &
			OutValue = (Grid(ix,iy,iz)*(widthZ-dz)+Grid(ix,iy,iz+1)*dz) / widthZ

		END IF

		!!!!!!!!!!!!!!!!!!!!!!!
		!! 0-D INTERPOLATION !!
		IF (dim == 0) THEN
			OutValue = Grid(ix,iy,iz)
		END IF

	END FUNCTION InterpolateGrid

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Reverse engineer the physical location of a given !!
	!! gridpoint.  Take indices and a grid number and	 !!
	!! return the 3-space coordinates.					 !!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	FUNCTION LocateGridPoint(X, Y, Z, Index) RESULT (Location)

		USE ModelParameters
		IMPLICIT NONE

		!! In-Out Variables
		INTEGER :: X, Y, Z, Index	!! Index is tells which grid spacing is used
		REAL*8  :: Location(3)
	
		!! DON'T SWITCH THIS WITHOUT CHANGING IT IN CLOSESTGRIDPOINT() AS WELL
		!! (but don't set it true or it fucks up the volume calculation)
		LOGICAL, PARAMETER :: PointsOnBoundaries = .FALSE.

		!! Determine the cell grid spacing
		IF (XGridPoints(Index) > 1) THEN 
			IF (PointsOnBoundaries) THEN
				Location(1) = (X - 1)   * DomainX / (XGridPoints(Index) - 1.)
			ELSE
				Location(1) = (X - 0.5) * DomainX / XGridPoints(Index)
			END IF
		ELSE
			Location(1) = DomainX / 2.
		END IF

		IF (YGridPoints(Index) > 1) THEN 
			IF (PointsOnBoundaries) THEN
				Location(2) = (Y - 1)   * DomainY / (YGridPoints(Index) - 1.)
			ELSE
				Location(2) = (Y - 0.5) * DomainY / YGridPoints(Index)
			END IF
		ELSE
			Location(2) = DomainY / 2.
		END IF

		IF (ZGridPoints(Index) > 1) THEN
			IF (PointsOnBoundaries) THEN
				Location(3) = (Z - 1)   * DomainZ / (ZGridPoints(Index) - 1.)
			ELSE
				Location(3) = (Z - 0.5) * DomainZ / ZGridPoints(Index)
			END IF
		ELSE
			Location(3) = DomainZ / 2.
		END IF

	END FUNCTION LocateGridPoint

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!! Identifies the closest gridpoint in a given spacing	!!
	!! to a given physical location.						!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ClosestGridPoint(XX, YY, ZZ, Index, X, Y, Z) 

		USE ModelParameters
		IMPLICIT NONE

		!! In-Out Variables
		INTEGER :: X, Y, Z, Index	!! Index is tells which grid spacing is used
		REAL*8  :: XX, YY, ZZ

		!! DON'T SWITCH THIS WITHOUT CHANGING IT IN CLOSESTGRIDPOINT() AS WELL
		!! (but don't set it true or it fucks up the volume calculation)
		LOGICAL, PARAMETER :: PointsOnBoundaries = .FALSE.

		!! Determine the cell grid spacing
		IF (XGridPoints(Index) > 1) THEN 
			IF (PointsOnBoundaries) THEN
				X = ANINT(((XGridPoints(Index) - 1.) * XX) / DomainX + 1)
				X = MAX(1,X)
				X = MIN(X,XGridPoints(Index))
			ELSE 
				X = CEILING(XGridPoints(Index) * XX / DomainX)
			END IF
		ELSE
			X = 1.
		END IF

		IF (YGridPoints(Index) > 1) THEN 
			IF (PointsOnBoundaries) THEN
				Y = ANINT(((YGridPoints(Index) - 1.) * YY) / DomainY + 1)
				Y = MAX(1,Y)
				Y = MIN(Y,YGridPoints(Index))
			ELSE 
				Y = CEILING(YGridPoints(Index) * YY / DomainY)
			END IF
		ELSE
			Y = 1.
		END IF

		IF (ZGridPoints(Index) > 1) THEN 
			IF (PointsOnBoundaries) THEN
				Z = ANINT(((ZGridPoints(Index) - 1.) * ZZ) / DomainZ + 1)
				Z = MAX(1,Z)
				Z = MIN(Z,ZGridPoints(Index))
			ELSE
				Z = CEILING(ZGridPoints(Index) * ZZ / DomainZ)
			END IF
		ELSE
			Z = 1.
		END IF

	END SUBROUTINE ClosestGridPoint
