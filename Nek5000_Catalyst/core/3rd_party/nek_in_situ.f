c-----------------------------------------------------------------------
      subroutine in_situ_init(type_index)
      integer type_index
#ifdef VISIT
      call visit_init()
#elif CATALYST
      call catalyst_init()
#elif ADIOS2
      call adios2_init()
#elif MALLEABLE
      call adios2_init_malleable(type_index)
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_check()
#ifdef VISIT
      call visit_check()
#elif CATALYST
      call catalyst_process()
#elif ADIOS2
      call adios2_write()
#elif MALLEABLE
      call adios2_write()
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_end()
#ifdef VISIT
      call visit_end()
#elif CATALYST
      call catalyst_end()
#elif ADIOS2
      call adios2_end()
#elif MALLEABLE
      call adios2_end()
#endif
      end
c-----------------------------------------------------------------------

