!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF DS
!
! headers
!
#include<dev_defs.h>
!
subroutine INVERSION_driver(M_slk)
 !
 ! * INV:   M_slk/M => M_slk/M
 !
 use pars,          ONLY:SP,cONE,cZERO
 use linear_algebra,ONLY:USE_LK,USE_SLK
 use matrix,        ONLY:PAR_matrix
 use parallel_m,    ONLY:linear_algebra_is_parallel
 use timing_m,      ONLY:timing
 !
 implicit none
 !
 type(PAR_matrix) :: M_slk
 !
 ! Work Spaces ...
 !
 complex(SP)         :: det_
 integer             :: N,lib_in_use
 !
 call timing("INVERSION",OPR="start")
 !
 N=M_slk%N
 !
 lib_in_use=USE_LK
#if defined _SCALAPACK 
 !if (M_slk%kind=="SLK") lib_in_use=USE_SLK
 if (linear_algebra_is_parallel) lib_in_use=USE_SLK
 !
 if(lib_in_use==USE_SLK) then 
   call PARALLEL_inversion( M_slk, N ) 
 endif
#endif
 !
 if (lib_in_use==USE_LK) then
   call SERIAL_inversion(N,M_slk%blc(:,:,M_slk%I),det_,.false.)
 endif
 !
 call timing("INVERSION",OPR="stop")
 !
end subroutine INVERSION_driver
!
! 
subroutine LIN_SYSTEM_driver(M_slk,B_slk)
 !
 ! * LIN_SYSTEM: M_slk + B_slk => B_slk
 !
 use pars,          ONLY:SP
 use linear_algebra,ONLY:USE_LK,USE_SLK
 use matrix,        ONLY:PAR_matrix
 use parallel_m,    ONLY:linear_algebra_is_parallel
 use timing_m,      ONLY:timing
 use gpu_m,         ONLY:have_gpu
 use devxlib,       ONLY:devxlib_mapped
 !
 implicit none
 !
 type(PAR_matrix) :: M_slk
 type(PAR_matrix) :: B_slk
 !
 ! Workspace
 !
 integer  :: N,lib_in_use
 logical  :: USE_GPU
 !
 call timing("LIN_SYSTEM",OPR="start")
 !
 N=M_slk%N
 !
 USE_GPU=.false.
 if (have_gpu.and.devxlib_mapped(DEV_VAR(M_slk%blc))) USE_GPU=.true.
 !
 lib_in_use=USE_LK
#if defined _SCALAPACK 
 !if (M_slk%kind=="SLK") lib_in_use=USE_SLK
 if (linear_algebra_is_parallel.and..not.USE_GPU) lib_in_use=USE_SLK
 !
 if(lib_in_use==USE_SLK) then
   call PARALLEL_lin_system( M_slk, B_slk, N, N)
 endif
#endif
   !
 if(lib_in_use==USE_LK) then
   !
   if (USE_GPU) then
#if defined _GPU && defined _PGI
     call SERIAL_lin_system_gpu(N,N,DEV_VAR(M_slk%blc)(:,:,M_slk%I),DEV_VAR(B_slk%blc)(:,:,B_slk%I))
#else
     call error('[LA] GPU support required but not present')
#endif
   else
     call SERIAL_lin_system(N,N,M_slk%blc(:,:,M_slk%I),B_slk%blc(:,:,B_slk%I))
   endif
   !
 endif
 !
 call timing("LIN_SYSTEM",OPR="stop")
 !
end subroutine LIN_SYSTEM_driver
!
!
subroutine MATMUL_driver(M_slk,B_slk,C_slk)
 !
 ! ####> MAT_MUL
 ! C = M * B
 !
 use pars,          ONLY:SP,cONE,cZERO
 use linear_algebra,ONLY:USE_LK,USE_SLK
 use matrix,        ONLY:PAR_matrix
 use wrapper,       ONLY:M_by_M
 use parallel_m,    ONLY:linear_algebra_is_parallel
 use timing_m,      ONLY:timing
 !
 implicit none
 !
 type(PAR_matrix) :: M_slk
 type(PAR_matrix) :: B_slk
 type(PAR_matrix) :: C_slk
 !
 ! Workspace
 !
 integer  :: N,lib_in_use
 !
 call timing("MATMUL",OPR="start")
 !
 N=M_slk%N
 !
 lib_in_use=USE_LK
#if defined _SCALAPACK 
 !if (M_slk%kind=="SLK") lib_in_use=USE_SLK
 if (linear_algebra_is_parallel) lib_in_use=USE_SLK
 !
 if(lib_in_use==USE_SLK) then
   call PARALLEL_M_by_M( M_slk, B_slk, C_slk, 'N', 'N', N, N , N )
 endif
#endif
 !
 if(lib_in_use==USE_LK) then
   call M_by_M('N','N',N,N,N,cONE,M_slk%blc(:,:,M_slk%I),N,&
   &           B_slk%blc(:,:,B_slk%I),N,cZERO,C_slk%blc(:,:,C_slk%I),N)
 endif
 !
 call timing("MATMUL",OPR="stop")
 !
end subroutine MATMUL_driver
