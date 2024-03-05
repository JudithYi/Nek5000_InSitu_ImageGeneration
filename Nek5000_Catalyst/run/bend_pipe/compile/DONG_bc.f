c-----------------------------------------------------------------------
c USER-DEFINED FUNCTIONS 
c    Flow over an ellipse
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      return
      end
C-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      ffx = 0.0
      ffy = 0.0
      ffz = 0.0   ! you can fix the flow rate through param(54), param(55)
      return
      end
      subroutine userq  (ix,iy,iz,ieg)
      return
      end
!--------------------------------------------------
      subroutine userchk
      !implicit none   commented by saleh
      include 'SIZE'
      include 'GEOM'                    ! xm1, ym1, zm1
      include 'SOLN'                    ! T
      include 'MASS'                    !BM1 for lambda2
      include 'TSTEP'                   ! ISTEP
      include 'INPUT'           ! PARAM(12) (DT)
      include 'CHKPOINT'
!srr!      include 'USERPAR'
      !include 'TOTAL'
      real ffx_new,ffy_new,ffz_new  !saleh
      common /cforce/ ffx_new,ffy_new,ffz_new
c   for torque calculations
      real x0(3)
      save x0
      common /ctorq/ dragx(0:maxobj),dragpx(0:maxobj),dragvx(0:maxobj)
     $             , dragy(0:maxobj),dragpy(0:maxobj),dragvy(0:maxobj)
     $             , dragz(0:maxobj),dragpz(0:maxobj),dragvz(0:maxobj)
c
     $             , torqx(0:maxobj),torqpx(0:maxobj),torqvx(0:maxobj)
     $             , torqy(0:maxobj),torqpy(0:maxobj),torqvy(0:maxobj)
     $             , torqz(0:maxobj),torqpz(0:maxobj),torqvz(0:maxobj)
c
     $             , dpdx_mean,dpdy_mean,dpdz_mean
     $             , dgtq(3,4)
      real e2
      integer n
      real ubar
!-------------------------------------------------- 
c      ! Restart code
      if (ISTEP.eq.0) then
         CHKPTSTEP=uparam(4)
         if (uparam(3).eq.1) then
            IFCHKPTRST=.TRUE.
         else
            IFCHKPTRST=.FALSE.
         endif
      endif
      call checkpoint           ! Restart check
!--------------------------------------------------
      ! Stats code
      !!call stat_avg
      n=nx1*ny1*nz1*nelv
      !computing tauW
      if (istep.eq.0) then
         call set_obj                   ! objects for surface integrals
         call rzero(x0,3)               ! define x0=0, note: torque w.r.t. x0
      endif
      call torque_calc(1.0,x0,.false.,.false.) ! wall shear
      wall_area= ZLENPIPE*RAD
      wall_area=2.*PI*wall_area
      !write(*,*) 'wall_area=',wall_area,ZLENPIPE
!      tauw=0.5*(dragz(0)+dragz(1))/wall_area   !wall shear stress 
      !Note: we have 2 objects starting from 0, see /ctorq/. So we do averaging. To see that, make the 3rd arg in torque_calc .true. 
      rho=1.0
!      if (nid.eq.0) write(*,*) 'time, uTau: ',time,sqrt(tauw/rho)
!      write(*,*) "dragOnEllipse: ",dragz(0),dragz(1)
      return
      end
!--------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      !implicit none
      integer ix,iy,iz,iside,eg,e
      include 'SIZE' 
      include 'PARALLEL'  ! GLLEL
      include 'NEKUSE'
      e = gllel(ieg)
!      pa = 0.0
      if (cbu.eq.'o  ') then       !Dong outflow BC
!         U0 = 1.0                  ! characteristic velocity
!         delta = 0.1               ! small positive constant
!         pa = dongOutflow(ix,iy,iz,e,iside,U0,delta)         
      else   !'v  ' & 'on '
         ux= 1.0 
         uy= 0.0 
         uz= 0.0 
      endif
      temp = 0.0
      return
      end
!--------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      implicit none
      integer ix,iy,iz,ieg
      integer meanProfType   !saleh
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      ux=1.0
      uy=0.0
      uz=0.0
      temp=0
      return
      end
!--------------------------------------------------
      subroutine usrdat
      return
      end
!--------------------------------------------------
      subroutine usrdat2
      !implicit none
      include 'SIZE'
      include 'TOTAL'
      !set all non-periodic BCs here. This is required due to generating mesh by gmsh and converting it by gmsh2nek
      !Here are the IDs according to ellipse.msh:
      ! 1 "inlet"
      ! 2 "outlet"
      ! 3 "wall"
      ! 4 "freestreamUp"
      ! 5 "freestreamLo"
      do iel=1,nelv
         do ifc=1,2*ndim
            id_face = bc(5,ifc,iel,1)
            if (id_face.eq.1) then        ! inlet 
               cbc(ifc,iel,1) = 'v  '
            elseif (id_face.eq.2) then    ! outlet
               cbc(ifc,iel,1) = 'o  '     ! use either 'O  ' with sponge or 'o  ' with Dong BC
            elseif (id_face.eq.3) then    ! wall
               cbc(ifc,iel,1) = 'W  '
            elseif (id_face.eq.4 .OR. id_face.eq.5) then    ! surface 4/5 for free-stream
               cbc(ifc,iel,1) = 'ON '
            endif
         enddo
      enddo
      return
      end
!--------------------------------------------------
      subroutine usrdat3
c      implicit none
      return
      end
c-----------------------------------------------------------------------
      subroutine set_obj  ! define objects for surface integrals
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e,f
c
c     Define new objects
c
      nobj = 1
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F' 
         hcode( 2,ii) = 'F' 
         hcode( 3,ii) = 'F' 
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj
c
      if (maxobj.lt.nobj) write(6,*) 'increase maxobj in SIZEu. rm *.o'
      if (maxobj.lt.nobj) call exitt
c
      nxyz = nx1*ny1*nz1
      nface=2*ndim
      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,1).eq.'W  ') then
c            write(*,*) 'wallElements=',e,f   !saleh, to check wall elements
            iobj=1   !wall
            if (iobj.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               ieg = lglel(e)
               object(iobj,mem,1) = ieg
               object(iobj,mem,2) = f
c              write(6,1) iobj,mem,f,ieg,e,nid,' OBJ'
    1          format(6i9,a4)
            endif
c
         endif
      enddo
      enddo
c     write(6,*) 'number',(nmember(k),k=1,4)
c
      return
      end
C-----------------------------------------------------------------------
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