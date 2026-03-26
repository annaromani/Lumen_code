!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
! headers
!
#include<y_memory.h>
!
subroutine SERIAL_diago(n,M,E,V_left,V_right)
 !
 !  The right eigenvector V_right(j) of M satisfies
 !                   M * V_right(j) = E(j) * V_right(j)
 !  where E(j) is its eigenvalue.
 !
 !  The left eigenvector V_left(j) of M satisfies
 !                V_left(j)**H * M = E(j) * V_left(j)**H
 !  where V_left(j)**H denotes the conjugate transpose of V_left(j).
 !
 use pars,            ONLY: SP,schlen
 use linear_algebra,  ONLY: LINEAR_ALGEBRA_error,LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS,&
&                           magma_init_done,magma_setup,la_xgeev
#ifdef _MAGMA
 use magma,           ONLY: magmaf_zgeev_m,magmaf_cgeev_m
#endif
 use timing_m,        ONLY: timing
#ifdef _OPENMP
 use omp_lib
#endif
 use y_memory_alloc
 !
 implicit none
 !
 character(schlen) :: subname
 integer           :: n
 complex(SP)       :: M(n,n),E(n),V_left(n,n),V_right(n,n)
 type(LALGEBRA_WS) :: WS
 !
 integer :: nthreads_save

#if defined _MAGMA
# if defined _DOUBLE
#   define magmaf_xgeev_m magmaf_zgeev_m
# else
#   define magmaf_xgeev_m magmaf_cgeev_m
# endif
#endif

#if defined _MAGMA
 subname="SER_GENERAL_magma"
#else
 subname="SER_GENERAL_lapack"
#endif

 call timing(trim(subname),OPR='start')
 !
 ! Thread safety
#if defined _OPENMP
 nthreads_save=omp_get_max_threads()
 call omp_set_num_threads(1)
#endif

 !
 ! Workspace
 !
 allocate(WS%v_cmplx(1))
 YAMBO_ALLOC(WS%v_real,(2*n))
 !
#if defined _MAGMA
 if (magma_init_done) then
   call magmaf_xgeev_m('V','V', n,M,n,E,V_left,n,V_right,n,&
&                      WS%v_cmplx,-1,WS%v_real,WS%i_fail)
 else
#endif
   call la_xgeev('V','V',n,M,n,E,V_left,n,V_right,n,WS%v_cmplx,-1,WS%v_real,WS%i_fail)
#if defined _MAGMA
 endif
#endif
 !
 if(WS%i_fail /=0 ) call LINEAR_ALGEBRA_error(' _XGEEV_ ','WS%v_cmplx size query failed')
 !
 WS%dim=int(real(WS%v_cmplx(1)))
 deallocate(WS%v_cmplx)
 YAMBO_ALLOC(WS%v_cmplx,(WS%dim))

#if defined _MAGMA
 if (magma_init_done) then
   !   
   ! use magma (if available)
   call magmaf_xgeev_m('V','V', n,M,n,E,V_left,n,V_right,n,&
   &                    WS%v_cmplx,WS%dim,WS%v_real,WS%i_fail)
   !
 else   
#endif
   !   
   ! use lapack (otherwise)
   call la_xgeev('V','V',n,M,n,E,V_left,n,V_right,n,WS%v_cmplx,WS%dim,WS%v_real,WS%i_fail)
   !
#if defined _MAGMA
 endif
#endif
 !
 if(WS%i_fail /= 0) call LINEAR_ALGEBRA_error(' _XGEEV_ ','failed')
 !
 call LINEAR_ALGEBRA_WS_reset(WS)
 !
 ! Thread safety
#if defined _OPENMP
 call omp_set_num_threads(nthreads_save)
#endif
 !   
 call timing(trim(subname),OPR='stop')
 !
end subroutine
