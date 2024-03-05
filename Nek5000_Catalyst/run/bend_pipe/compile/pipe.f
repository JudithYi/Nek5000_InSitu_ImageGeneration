! pipe dimensions: radius and length
#define RAD 0.5
#define ZLENPIPE 15.0
!======================================================================
      subroutine uservp(ix,iy,iz,ieg)
      return
      end subroutine
!======================================================================
      subroutine userf(ix,iy,iz,ieg)
      return
      end subroutine
!======================================================================
      subroutine userq(ix,iy,iz,ieg)
      return
      end subroutine
!======================================================================
      subroutine userchk()
      implicit none

      include 'SIZE'
      include 'GEOM'            ! xm1, ym1, zm1
      include 'SOLN'            ! T
      include 'MASS'            ! BM1
      include 'TSTEP'           ! ISTEP, PI
      include 'INPUT'           ! PARAM(12) (DT)

!-----------------------------------------------------------------------
      if (ISTEP.eq.0) then
         call frame_start       ! start framework
         call outpost(vx,vy,vz,pr,t,'ini')
      endif

      call frame_monitor        ! monitor simulation

      call chkpt_main           ! save/load files for full-restart
      call stat_avg
      
      call gsem_main            ! synthetic eddy method

      call lambda2(t(1,1,1,1,1))

      call trunc_main()
  
      if (ISTEP.eq.NSTEPS.or.LASTEP.eq.1) then
         call frame_end         ! finalise framework
      endif
      
      return
      end subroutine
!======================================================================
      subroutine userbc (ix,iy,iz,iside,ieg)
      implicit none

      include 'SIZE'
      include 'NEKUSE'
      include 'PARALLEL'

      integer ix,iy,iz,iside,ieg,e
      real delta, U0
      real dongOutflow ! functions declaration

      e = gllel(ieg)
      pa = 0.0
      if (cbu.eq.'o  ') then       !Dong outflow BC
         U0 = 1.0                  ! characteristic velocity
         delta = 0.1               ! small positive constant
         pa = dongOutflow(ix,iy,iz,e,iside,U0,delta)         
      endif

      return
      end subroutine
!======================================================================
      subroutine useric (ix,iy,iz,ieg)
      implicit none

      include 'SIZE'  
      include 'NEKUSE'

      integer :: ix,iy,iz,ieg

      Uy = 0.
      Ux = 0.
      Uz = 0.
                    
      return
      end subroutine
!======================================================================
      subroutine usrdat
      implicit none

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

      call setbc(1,1,'v  ') ! inlet
      call setbc(2,1,'o  ') ! outlet
      call setbc(3,1,'W  ') ! wall

      return
      end subroutine
!======================================================================
 
      subroutine usrdat2
      implicit none
      include 'SIZE' ! nelv
      include 'GEOM' 
      include 'NEKUSE' ! XM1
      include 'INPUT' 
      include 'PARALLEL' 
      include 'SOLN' 

      integer i, iel
      real angle, circumf, bent_radius, bent_phi, z_offset, pi

      integer ii, jj, kk
      ! use some reals to force ratios
      real nL, nelF, nelUp, nelBent, nelDown1, nelDown2
      real lPipe, lUp, lBent, lDown1, lDown2

      pi=4.0*atan(1.d0)

c		Rescale incoming pipe radius to be r0
c
      call rescale_x(xm1,-RAD,RAD)
      call rescale_x(ym1,-RAD,RAD)

c		Rescale incoming pipe length to [0,1]
c
      call rescale_x(zm1,0.,1.)

c		Mesh dilatation
c
      nL          = 510. ! total # elements streamwise direction
      lPipe       = 1.
      z_offset    = -10.
      bent_radius = 1.5
      bent_phi    = 0.5*pi

      lUp    = abs(z_offset)
      lBent  = bent_radius*bent_phi
      lDown1 = 4.
      lDown2 = 16.

      nelF     = 900.
      nelUp    = 150. ! These values are derived from MATLAB script mesh_resolution
      nelBent  = 53.  ! imposing deltaz+_max = 9 in the outer bend
      nelDown1 = 67.
      nelDown2 = 240.

      if (nid.eq.0) then
         write(*,*) 'nL = ', nL
         write(*,*) 'L  = ', lPipe
         write(*,*) ''
         write(*,*) 'lUp   = ', lUp
         write(*,*) 'lBent = ', lBent
         write(*,*) 'lDown1 = ', lDown1
         write(*,*) 'lDown2 = ', lDown2
         write(*,*) ''
         write(*,*) 'nelF     = ', nelF
         write(*,*) 'nelUp    = ', nelUp
         write(*,*) 'nelBent  = ', nelBent
         write(*,*) 'nelDown1 = ', nelDown1
         write(*,*) 'nelDown2 = ', nelDown2
      endif

      do i=1,nx1*ny1*nz1*nelv
        if (zm1(i,1,1,1).gt.0.and.zm1(i,1,1,1).le.lPipe/nL*nelUp) then
           zm1(i,1,1,1) = zm1(i,1,1,1)*lUp/(nelUp/nL*lPipe)
          
        elseif (zm1(i,1,1,1).gt.lPipe/nL*nelUp .and.
     $          zm1(i,1,1,1).le.lPipe/nL*(nelUp+nelBent)) then

                zm1(i,1,1,1) = (zm1(i,1,1,1) - nelUp/nL*lPipe)*
     $                          lBent/(nelBent/nL*lPipe) + lUp

        elseif (zm1(i,1,1,1).gt.lPipe/nL*(nelUp+nelBent) .and. 
     $          zm1(i,1,1,1).le.lPipe/nL*(nelUp+nelBent+nelDown1)) then
                zm1(i,1,1,1) = (zm1(i,1,1,1) 
     $                       - (nelUp+nelBent)/nL*lPipe)*lDown1/
     $                         (nelDown1/nL*lPipe) + (lUp+lBent)

        elseif (zm1(i,1,1,1).gt.lPipe/nL*(nelUp+nelBent+nelDown1)) then
                zm1(i,1,1,1) = (zm1(i,1,1,1) - (nelUp+nelBent+nelDown1)/
     $                         nL*lPipe)*lDown2/(nelDown2/nL*lPipe) + 
     $                         (lUp+lBent+lDown1)
        endif
      enddo
      param(59) = 1. ! 1 = deformed mesh

      if (abs(bent_phi).gt.1e-10) then

c     Bend and translate the straight pipe
c     The sweep of the bent arc is bent_phi
c     The pipe inlet is moved at (0,r1,0) so that the center of the bend
c     is in (0,0,0)

      circumf = bent_radius*bent_phi
c
      do i=1,nx1*ny1*nz1*nelv
        zm1(i,1,1,1) = zm1(i,1,1,1) + z_offset
        xm1(i,1,1,1) = xm1(i,1,1,1) + bent_radius

        if (zm1(i,1,1,1).gt.0.and.zm1(i,1,1,1).le.circumf) then
          angle=zm1(i,1,1,1)/bent_radius
          zm1(i,1,1,1) = xm1(i,1,1,1)*sin(angle)
          xm1(i,1,1,1) = xm1(i,1,1,1)*cos(angle)
        elseif (zm1(i,1,1,1).gt.circumf) then

          angle = zm1(i,1,1,1) - circumf

          zm1(i,1,1,1) = xm1(i,1,1,1)*sin(bent_phi)
     $    + angle*cos(bent_phi)
          xm1(i,1,1,1) = xm1(i,1,1,1)*cos(bent_phi) 
     $    - angle*sin(bent_phi)
        endif
      enddo
      param(59) = 1. ! 1 = deformed mesh

      endif
      call outpost(vx,vy,vz,pr,t,'ms0')
      
      return
      end subroutine

!======================================================================
      subroutine usrdat3
      implicit none

      return
      end subroutine
!======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     register modules
      call io_register
      call chkpt_register
      call stat_register
      call gsem_register
      call trunc_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     initialise modules
      call chkpt_init
      call stat_init
      call gsem_init
      call trunc_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------

      
      return
      end subroutine
!======================================================================
!> @brief Generate single eddy for a given family
!! @param[out] epos       eddy position
!! @param[out] eps        eddy orientation
!! @param[in]  nfam       family number
!! @param[in]  ifinit     intial distribution
      subroutine usr_gen_eddy(epos,eps,nfam,ifinit)
      implicit none

      include 'SIZE'
      include 'TSTEP'           ! pi
      include 'GSYEMD'

      ! argument list
      real epos(ldim)
      integer eps(ldim)
      integer nfam
      logical ifinit

      ! local variables
      real rho, theta, vrtmp(ldim)
      integer il

      real yp_cut
      parameter (yp_cut=0.4862)
      
      ! functions
      real mth_ran_rng
!-----------------------------------------------------------------------
      ! get random position with respect to the begining of coordinate system
      ! this must be adapted to considered amily inflow shape (in current example circle)
      rho = yp_cut*sqrt(mth_ran_rng(0.0,1.0))  
      theta = mth_ran_rng(0.,2.0*pi)

      vrtmp(1) = rho*cos(theta)
      vrtmp(2) = rho*sin(theta)
      if (ifinit) then
         vrtmp(ldim) = mth_ran_rng(-gsem_bext(nfam),gsem_bext(nfam))
      else
         vrtmp(ldim) = -gsem_bext(nfam)
      endif

      ! rotate coordinates with respect to family normal
      call mth_rot3Da(epos,vrtmp,
     $        gsem_raxs(1,nfam),gsem_rang(nfam))

      ! shift vertex position with respect to family centre
      ! notice; rotation must be done first
      do il=1,ldim
        epos(il) = epos(il) + gsem_bcrd(il,nfam)
      enddo

      ! get random orientation
      do il=1,ldim
         rho = mth_ran_rng(0.0,1.0)
         if (rho.gt.0.5) then
            eps(il) = 1
         else
            eps(il) = -1
         endif
      enddo

      return
      end subroutine
!=======================================================================
      function dongOutflow(ix,iy,iz,iel,iside,u0,delta)
      include 'SIZE'
      include 'SOLN'
      include 'GEOM'
      real sn(3)
      ux = vx(ix,iy,iz,iel)
      uy = vy(ix,iy,iz,iel)
      uz = vz(ix,iy,iz,iel)
      call getSnormal(sn,ix,iy,iz,iside,iel)
      vn = ux*sn(1) + uy*sn(2) + uz*sn(3) 
      S0 = 0.5*(1.0 - tanh(vn/u0/delta))
      dongOutflow = -0.5*(ux*ux+uy*uy+uz*uz)*S0
      return
      end
C-----------------------------------------------------------------------
!> @brief Provide element coordinates and local numbers (user interface)
!! @param[out]  idir              mapping (uniform) direction
!! @param[out]  ctrs              2D element centres
!! @param[out]  cell              local element numberring
!! @param[in]   lctrs1,lctrs2     array sizes
!! @param[out]  nelsort           number of local 3D elements to sort
!! @param[out]  map_xm1, map_ym1  2D coordinates of mapped elements
!! @param[out]  ierr              error flag
      subroutine user_map2d_get(idir,ctrs,cell,lctrs1,lctrs2,nelsort,
     $     map_xm1,map_ym1,ierr)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! [XYZ]C
      include 'GEOM'            ! [XYZ]M1

!     argument list
      integer idir
      integer lctrs1,lctrs2
      real ctrs(lctrs1,lctrs2)  ! 2D element centres  and diagonals 
      integer cell(lctrs2)      ! local element numberring
      integer nelsort           ! number of local 3D elements to sort
      real map_xm1(lx1,lz1,lelt), map_ym1(lx1,lz1,lelt)
      integer ierr              ! error flag

!     local variables
      integer ntot              ! tmp array size for copying
      integer el ,il ,jl        ! loop indexes
      integer nvert             ! vertex number
      real rnvert               ! 1/nvert
      real xmid,ymid            ! 2D element centre
      real xmin,xmax,ymin,ymax  ! to get approximate element diagonal
      integer ifc               ! face number

!     dummy arrays
      real xcoord(8,LELT), ycoord(8,LELT) ! tmp vertex coordinates

#ifdef DEBUG
!     for testing
      character*3 str1, str2
      integer iunit, ierrl
      ! call number
      integer icalldl
      save icalldl
      data icalldl /0/
#endif

!-----------------------------------------------------------------------
!     initial error flag
      ierr = 0
!     set important parameters
!     uniform direction; should be taken as input parameter
!     x-> 1, y-> 2, z-> 3
      idir = 3
      
!     get element midpoints
!     vertex number
      nvert = 2**NDIM
      rnvert= 1.0/real(nvert)

!     eliminate uniform direction
      ntot = 8*NELV
      if (idir.EQ.1) then  ! uniform X
         call copy(xcoord,YC,ntot) ! copy y
         call copy(ycoord,ZC,ntot) ! copy z
      elseif (idir.EQ.2) then  ! uniform Y
         call copy(xcoord,XC,ntot) ! copy x
         call copy(ycoord,ZC,ntot) ! copy z
      elseif (idir.EQ.3) then  ! uniform Z
         call copy(xcoord,XC,ntot) ! copy x
         call copy(ycoord,YC,ntot) ! copy y
      endif

!     set initial number of elements to sort
      nelsort = 0
      call izero(cell,NELT)

!     for every element
      do el=1,NELV
!     element centre
         xmid = xcoord(1,el)
         ymid = ycoord(1,el)
!     element diagonal
         xmin = xmid
         xmax = xmid
         ymin = ymid
         ymax = ymid
         do il=2,nvert
            xmid=xmid+xcoord(il,el)
            ymid=ymid+ycoord(il,el)
            xmin = min(xmin,xcoord(il,el))
            xmax = max(xmax,xcoord(il,el))
            ymin = min(ymin,ycoord(il,el))
            ymax = max(ymax,ycoord(il,el))
         enddo
         xmid = xmid*rnvert
         ymid = ymid*rnvert

!     count elements to sort
            nelsort = nelsort + 1
!     2D position
!     in general this coud involve some curvilinear transform
            ctrs(1,nelsort)=xmid
            ctrs(2,nelsort)=ymid
!     reference distance
            ctrs(3,nelsort)=sqrt((xmax-xmin)**2 + (ymax-ymin)**2)
            if (ctrs(3,nelsort).eq.0.0) then
               ierr = 1
               return
            endif
!     element index
            cell(nelsort) = el
      enddo

!     provide 2D mesh
!     in general this coud involve some curvilinear transform
      if (idir.EQ.1) then  ! uniform X
         ifc = 4
         do el=1,NELV
            call ftovec(map_xm1(1,1,el),ym1,el,ifc,nx1,ny1,nz1)
            call ftovec(map_ym1(1,1,el),zm1,el,ifc,nx1,ny1,nz1)
         enddo
      elseif (idir.eq.2) then  ! uniform y
         ifc = 1
         do el=1,nelv
            call ftovec(map_xm1(1,1,el),xm1,el,ifc,nx1,ny1,nz1)
            call ftovec(map_ym1(1,1,el),zm1,el,ifc,nx1,ny1,nz1)
         enddo
      elseif (idir.eq.3) then  ! uniform z
         ifc = 5
         do el=1,nelv
            call ftovec(map_xm1(1,1,el),xm1,el,ifc,nx1,ny1,nz1)
            call ftovec(map_ym1(1,1,el),ym1,el,ifc,nx1,ny1,nz1)
         enddo
      endif

#ifdef DEBUG
!     testing
      ! to output refinement
      icalldl = icalldl+1
      call io_file_freeid(iunit, ierrl)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalldl
      open(unit=iunit,file='map2d_usr.txt'//str1//'i'//str2)
      
      write(iunit,*) idir, NELV, nelsort
      write(iunit,*) 'Centre coordinates and cells'
      do el=1,nelsort
         write(iunit,*) el, ctrs(:,el), cell(el)
      enddo
      write(iunit,*) 'GLL coordinates'
      do el=1,nelsort
         write(iunit,*) 'Element ', el
         write(iunit,*) 'XM1'
         do il=1,nz1
            write(iunit,*) (map_xm1(jl,il,el),jl=1,nx1)
         enddo
         write(iunit,*) 'YM1'
         do il=1,nz1
            write(iunit,*) (map_ym1(jl,il,el),jl=1,nx1)
         enddo
      enddo
      close(iunit)
#endif

      return
      end subroutine
!=======================================================================
!> @brief Provide velocity, deriv. and vort. in required coordinates and normalise pressure
!! @param[out]   lvel             velocity
!! @param[out]   dudx,dvdx,dwdx   velocity derivatives
!! @param[out]   vort             vorticity
!! @param[inout] pres             pressure
      subroutine user_stat_trnsv(lvel,dudx,dvdx,dwdx,vort,pres)
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'               ! if3d
      include 'GEOM'

      ! argument list
      real lvel(LX1,LY1,LZ1,LELT,3) ! velocity array
      real dudx(LX1,LY1,LZ1,LELT,3) ! velocity derivatives; U
      real dvdx(LX1,LY1,LZ1,LELT,3) ! V
      real dwdx(LX1,LY1,LZ1,LELT,3) ! W
      real vort(LX1,LY1,LZ1,LELT,3) ! vorticity
      real pres(LX1,LY1,LZ1,LELT)   ! pressure

      ! local variables
      integer itmp              ! dummy variable
      integer il, jl            ! loop index
      integer ifll              ! field number for object definition
      real vrtmp(lx1*lz1)       ! work array for face
      real vrtmp2(2)            ! work array
      
      ! functions
      real vlsum
!-----------------------------------------------------------------------
      ! Velocity transformation; simple copy
      itmp = NX1*NY1*NZ1*NELV
      call copy(lvel(1,1,1,1,1),VX,itmp)
      call copy(lvel(1,1,1,1,2),VY,itmp)
      call copy(lvel(1,1,1,1,3),VZ,itmp)

      ! Derivative transformation
      ! No transformation
      call gradm1(dudx(1,1,1,1,1),dudx(1,1,1,1,2),dudx(1,1,1,1,3),
     $      lvel(1,1,1,1,1))
      call gradm1(dvdx(1,1,1,1,1),dvdx(1,1,1,1,2),dvdx(1,1,1,1,3),
     $      lvel(1,1,1,1,2))
      call gradm1(dwdx(1,1,1,1,1),dwdx(1,1,1,1,2),dwdx(1,1,1,1,3),
     $      lvel(1,1,1,1,3))

      ! get vorticity
      if (IF3D) then
         ! curlx
         call sub3(vort(1,1,1,1,1),dwdx(1,1,1,1,2),
     $        dvdx(1,1,1,1,3),itmp)
         ! curly
         call sub3(vort(1,1,1,1,2),dudx(1,1,1,1,3),
     $        dwdx(1,1,1,1,1),itmp)
      endif
      ! curlz
      call sub3(vort(1,1,1,1,3),dvdx(1,1,1,1,1),dudx(1,1,1,1,2),itmp)
      
      ! normalise pressure
      ! in this example I integrate pressure over top faces marked "W"
      ifll = 1     ! I'm interested in velocity bc
      ! relying on mesh structure given by genbox set face number
      jl = 3
      call rzero(vrtmp2,2)  ! zero work array
      itmp = LX1*LZ1
      do il=1,nelv   ! element loop
         if (cbc(jl,il,ifll).eq.'W  ') then
            vrtmp2(1) = vrtmp2(1) + vlsum(area(1,1,jl,il),itmp)
            call ftovec(vrtmp,pres,il,jl,lx1,ly1,lz1)
            call col2(vrtmp,area(1,1,jl,il),itmp)
            vrtmp2(2) = vrtmp2(2) + vlsum(vrtmp,itmp)
         endif
      enddo
      ! global communication
      call gop(vrtmp2,vrtmp,'+  ',2)
      ! missing error check vrtmp2(1) == 0
      vrtmp2(2) = -vrtmp2(2)/vrtmp2(1)
      ! remove mean pressure
      itmp = LX1*LY1*LZ1*NELV
      call cadd(pres,vrtmp2(2),itmp)

      return
      end subroutine
!======================================================================

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
