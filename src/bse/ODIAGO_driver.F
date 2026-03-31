!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
! headers
!
#include<y_memory.h>
!
subroutine ODIAGO_driver(N,neig,ncpu_blacs,M_slk,V_right,V_left,E_cmpl,E_real)
 !
 use pars,          ONLY:SP,cZERO
 use matrix,        ONLY:PAR_matrix,MATRIX_reset
 use linear_algebra,ONLY:USE_LK,USE_SLK
 use parallel_m,    ONLY:ncpu
 use interfaces,    ONLY:MATRIX_transfer
 use parallel_int,  ONLY:PP_redux_wait,PP_wait
 use BS_solvers,    ONLY:BSS_ldiago_solver
#if defined _SCALAPACK
 use SLK_m,         ONLY:SLK_ORTHO
#endif
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(in)           :: N
 integer, intent(in)           :: neig
 integer, intent(in)           :: ncpu_blacs
 type(PAR_matrix), intent(in)  :: M_slk
 complex(SP),      intent(out) :: V_right(N,neig)
 complex(SP),      intent(out), optional :: V_left(N,neig)
 real(SP),         intent(out), optional :: E_real(neig)
 complex(SP),      intent(out), optional :: E_cmpl(neig)
 !
 ! Workspace
 !
 integer   :: rows(2),cols(2)
 type (PAR_matrix)   :: V_loc
 !
 if (present(E_cmpl).and.ncpu_blacs>1) call error("Parallel non-Hermitian diago not available")
 !
#if defined _MAGMA
 if (ncpu_blacs==1 .and. .not.magma_init_done) then
   call msg('sr','Initializing magma library in geev')
   call magma_setup()
 endif
#endif
 !
 if (ncpu_blacs==1) then
   !
   if (present(E_cmpl)) then
     !
     !... non-Hermitian case
     call SERIAL_diago(N,M_slk%blc(:,:,1),E_cmpl,V_left,V_right)
     !
   else
     !
     !... Hermitian case
     if (neig==N) call SERIAL_HERMITIAN_diago(N,M_slk%blc(:,:,1),E_real)
     if (neig< N) call SERIAL_HERMITIAN_diagx(N,M_slk%blc(:,:,1),E_real,neig)
     V_right=M_slk%blc(:,:neig,1)
     !
   endif
   !
   !
 else
   !
   if (present(E_cmpl)) then
     !
     !... non-Hermitian case
     !call PARALLEL_PSEUDOHERM_diag?( M_slk, V_loc, N, E_real, neig )
     !
   else
     !
     call MATRIX_init("SLK",  V_loc, N, 1 )
     V_loc%kind="SLK"
#if defined _SCALAPACK
     if (SLK_ORTHO%coordinate(1)==-1) then
       rows=(/1,1/)
       cols=(/1,1/)
     else
#endif
       rows=V_loc%BLCrows
       cols=V_loc%BLCcols
#if defined _SCALAPACK
     endif
#endif
     YAMBO_ALLOC(V_loc%blc,(rows(1):rows(2),cols(1):cols(2),1))
     V_loc%blc=cZERO
     !
     !... Hermitian case, QR diago
     if(neig==N) then
       if (BSS_ldiago_solver=="o") then
         ! QR Algorithm, less memory, slower (pcheev)
         call PARALLEL_HERMITIAN_diago( M_slk, V_loc, N, E_real )
       else
         ! Divide and conquer more memory, faster, ideal for HPC (pcheevd), O(n3)
         call PARALLEL_HERMITIAN_diagd( M_slk, V_loc, N, E_real )
       endif
     endif
     if(neig< N) then
       if (BSS_ldiago_solver=="x") then
         ! Bisection and Inverse Iteration, ok for very few eigenvalues (pcheevx), O(n3)
         call PARALLEL_HERMITIAN_diagx( M_slk, V_loc, N, E_real, neig )
       else
         ! Multiple Relatively Robust Representations, more scalable, better for
         ! more eigenvalues or anyway for larger matrices (pcheevr), O(n2)
         call PARALLEL_HERMITIAN_diagr( M_slk, V_loc, N, E_real, neig )
       endif
     endif
     !
     call MATRIX_transfer( M_in = V_loc , M=V_right )
     YAMBO_FREE(V_loc%blc)
     call MATRIX_reset(V_loc)
     !
   endif
   !
 endif
 !
end subroutine ODIAGO_driver
