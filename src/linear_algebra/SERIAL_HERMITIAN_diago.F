!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
! headers
!
#include<y_memory.h>
!
subroutine SERIAL_HERMITIAN_diago(n,M,E)
 !
 use pars,           ONLY: SP
 use drivers,        ONLY: l_nl_optics,l_real_time
 use linear_algebra, ONLY: magma_init_done,magma_setup,la_xheev
#ifdef _MAGMA
 use magma,          ONLY: magmaf_cheevd_m, magmaf_zheevd_m
#endif
 use timing_m,       ONLY: timing
 use com,            ONLY: msg
#ifdef _OPENMP
 use omp_lib
#endif
 use y_memory_alloc
 !
 implicit none
 !
 integer     :: n
 complex(SP) :: M(n,n)
 real(SP)    :: E(n)
 !
 character(64)     :: subname="SERIAL_HERMITIAN_diago"
 integer           :: lwork, lrwork, liwork, i_fail
 integer,     allocatable :: v_int(:)
 real(SP),    allocatable :: v_real(:)
 complex(SP), allocatable :: v_cmplx(:)
 !
 integer     :: nthreads_save

#if defined _DOUBLE
#  define magmaf_xheevd_m magmaf_zheevd_m
#else
#  define magmaf_xheevd_m magmaf_cheevd_m
#endif

 if(.not.(l_nl_optics.or.l_real_time)) call timing('SER_HERM_diago',OPR='start')
 !
 ! Thread safety
#if defined _OPENMP
 nthreads_save=omp_get_max_threads()
 call omp_set_num_threads(1)
#endif

 !
 ! Workspace
 !
 lwork=-1
 allocate(v_cmplx(1))
 !
 lrwork=max(1,3*n-2)
#if defined _MAGMA
 lrwork=1
#endif
 !
 allocate(v_real(lrwork))
 !
#if defined _MAGMA
 !
 lrwork=-1
 liwork=-1
 allocate(v_int(1))
 !
 if (.not.magma_init_done) then
   call msg('sr','Initializing magma library in geev')
   call magma_setup()
 endif
 !   
 call magmaf_xheevd_m(1,'v','u', n,M,size(M,1),E,v_cmplx,lwork,v_real,&
 &                       lrwork,v_int,liwork,i_fail)
#else
 call la_xheev('V','U',n,M,size(M,1),E,v_cmplx,lwork,v_real,i_fail)
#endif
 !
 if(i_fail /=0 ) call error(' _XHEEV_ v_cmplx size query failed')
 !
 lwork=nint(real(v_cmplx(1)))
 deallocate(v_cmplx)
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   YAMBO_ALLOC(v_cmplx,(lwork))
 else
   allocate(v_cmplx(lwork))
 endif
 !
#if defined _MAGMA
 !
 lrwork=nint(v_real(1))
 deallocate(v_real)
 liwork=v_int(1)
 deallocate(v_int)
 !
 if(.not.(l_nl_optics.or.l_real_time)) then
   YAMBO_ALLOC(v_real,(lrwork))
   YAMBO_ALLOC(v_int,(liwork))
 else
   allocate(v_real(lrwork))
   allocate(v_int(liwork))
 endif
 !   
 ! use magma (if available)
 !   
 ! 
 call magmaf_xheevd_m(1,'v','u',n, M,size(M,1),E,v_cmplx,lwork,v_real,&
&                     lrwork,v_int,liwork,i_fail)
 !
#else
 !
 call la_xheev('V','U',n,M,size(M,1),E,v_cmplx,lwork,v_real,i_fail)
 !
#endif
 !
 if(i_fail /= 0) call error(' XHEEV failed')
 !
 ! Thread safety
#if defined _OPENMP
 call omp_set_num_threads(nthreads_save)
#endif
 !
 if(.not.(l_nl_optics.or.l_real_time)) call timing('SER_HERM_diago',OPR='stop')
 !
end subroutine

