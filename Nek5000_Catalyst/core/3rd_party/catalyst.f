      subroutine catalyst_init()
#ifdef CATALYST
      include "mpif.h"
      CHARACTER*4 core
      integer par_rank, err
      real timer_catalyst
      common /PARALLELS/ par_rank, core
      save /PARALLELS/
      common /timerCatalyst/ timer_catalyst
c      implicit none
      call MPI_COMM_RANK(MPI_COMM_WORLD, par_rank, err)
      write(core, 10) par_rank
 10   format (I4)
      core = adjustl(trim(core))
      open(unit=44, file='perf/'//'core_'//adjustl(trim(core))//'.csv')
#endif
      call coprocessorinitialize()
      ! Add user defined pipelines
      call catalyst_usrpipe()
      timer_catalyst=0.0
      end

      subroutine catalyst_end
      common /nekmpi/ nid_,np_,nekcomm
      common /timerCatalyst/ timer_catalyst
      
#ifdef CATALYST
      close(unit=44)
#endif
      call coprocessorfinalize()
      write(*,*) nid_, " : ",timer_catalyst
      end

      subroutine catalyst_process()
      include "mpif.h"
      include 'SIZE'
      include 'TSTEP'
      common /timerCatalyst/ timer_catalyst
      real start_catalyst
      start_catalyst=mpi_wtime()
      if(ifoutfld) call catalyst_update()
      timer_catalyst=timer_catalyst+(mpi_wtime()-start_catalyst)
      end subroutine

      subroutine catalyst_update()
      include "mpif.h"
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
#ifdef CATALYST
      double precision before, after, lafter, cat, sim
      common /PARALLELS/ before, after, cat , sim
      save /PARALLELS/
      common /nekmpi/ nid_,np_,nekcomm
      integer flag, dim
      real t_dum (lx1,ly1,lz1,nelt)
      integer i, j, k, l
      do l = 1, nelt 
         do k = 1, lz1 
            do j = 1, ly1 
               do i = 1, lx1 
                  t_dum(i,j,k,l)=0.01*nid_
               enddo
            enddo
        enddo
      enddo

      before = MPI_Wtime()
#endif
c      lafter=0.0
      call requestdatadescription(istep, time, flag)
      if (flag .ne. 0) then
         call needtocreategrid(flag)
         dim = 2
         if (IF3D) dim = 3
         call creategrid(xm1, ym1, zm1, lx1, ly1, lz1, nelt, dim)
         call add_scalar_field(pr, "pressure"//char(0))
         call add_vector_field(vx, vy, vz, dim, "velocity"//char(0))
         call add_scalar_field(t_dum, "temperature"//char(0))
         call coprocess()
      end if
#ifdef CATALYST
      lafter = after
      after = MPI_WTIME()

      cat = after-before
      sim = before-lafter

c      lafter = after

      WRITE (44, *)  before, ',', after, ',', cat, ',', sim
#endif

      end

      subroutine catalyst_update2(
     &vx_2,vy_2,vz_2,pr_2)
      include "mpif.h"
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
#ifdef CATALYST
      real vx_2 (lx1,ly1,lz1,nelt)
      real vy_2 (lx1,ly1,lz1,nelt)
      real vz_2 (lx1,ly1,lz1,nelt)
      real pr_2 (lx1,ly1,lz1,nelt)
      double precision before, after, lafter, cat, sim
      common /PARALLELS/ before, after, cat , sim
      save /PARALLELS/
      common /nekmpi/ nid_,np_,nekcomm
      integer flag, dim
      real t_dum (lx1,ly1,lz1,nelt)
      integer i, j, k, l
      do l = 1, nelt 
         do k = 1, lz1 
            do j = 1, ly1 
               do i = 1, lx1 
                  t_dum(i,j,k,l)=0.01*nid_
               enddo
            enddo
        enddo
      enddo

      before = MPI_Wtime()
#endif
c      lafter=0.0
      call requestdatadescription(istep+10000, time+100, flag)
      if (flag .ne. 0) then
         call needtocreategrid(flag)
         dim = 2
         if (IF3D) dim = 3
         call creategrid(xm1, ym1, zm1, lx1, ly1, lz1, nelt, dim)
         call add_scalar_field(pr_2, "pressure"//char(0))
         call add_vector_field(vx_2,vy_2,vz_2,dim,"velocity"//char(0))
         call add_scalar_field(t_dum, "temperature"//char(0))
         call coprocess()
      end if
#ifdef CATALYST
      lafter = after
      after = MPI_WTIME()

      cat = after-before
      sim = before-lafter

c      lafter = after

      WRITE (44, *)  before, ',', after, ',', cat, ',', sim
#endif

      end
