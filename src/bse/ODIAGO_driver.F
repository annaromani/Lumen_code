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
subroutine ODIAGO_driver(N,neig,M_slk,V_right,V_left,E_cmpl,E_real)
 !
 use pars,          ONLY:SP,cZERO
 use matrix,        ONLY:PAR_matrix,MATRIX_reset
 use linear_algebra,ONLY:USE_LK,USE_SLK
 use parallel_m,    ONLY:ncpu
 use interfaces,    ONLY:MATRIX_transfer
 use parallel_int,  ONLY:PP_redux_wait,PP_wait
#if defined _SCALAPACK
 use SLK_m,         ONLY:SLK_ORTHO
#endif
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(in)           :: N
 integer, intent(in)           :: neig
 type(PAR_matrix), intent(in)  :: M_slk
 complex(SP),      intent(out) :: V_right(N,neig)
 complex(SP),      intent(out), optional :: V_left(N,neig)
 real(SP),         intent(out), optional :: E_real(neig)
 complex(SP),      intent(out), optional :: E_cmpl(neig)
 !
 ! Workspace
 !
 integer   :: lib_in_use,rows(2),cols(2)
 type (PAR_matrix)   :: V_loc
 !
 lib_in_use=USE_LK
#if defined _SCALAPACK 
 if(ncpu>1) lib_in_use=USE_SLK
#endif
!
 if (present(E_cmpl).and.ncpu>1) then
   call warning('Parallel non-herm diag not yet implemented. Switching to LK')
   lib_in_use=USE_LK
 endif
 !
#if defined _MAGMA
 if (lib_in_use==USE_LK .and. .not.magma_init_done) then
   call msg('sr','Initializing magma library in geev')
   call magma_setup()
 endif
#endif
 !
 if (lib_in_use==USE_LK) then
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
     !
   endif
   !
   V_right=M_slk%blc(:,:neig,1)
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
     !... Hermitian case
     if(neig==N) call PARALLEL_HERMITIAN_diago( M_slk, V_loc, N, E_real )
     if(neig< N) call PARALLEL_HERMITIAN_diagx( M_slk, V_loc, N, E_real, neig )
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
