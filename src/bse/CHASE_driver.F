!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM MG DS BM NM
!
! headers
!
#include<y_memory.h>
!
subroutine CHASE_driver(N,M,M_slk,V_left,V_right,E_cmpl,E_real,neig)
 !
 use pars,          ONLY:SP
 use matrix,        ONLY:PAR_matrix,MATRIX_reset
 use linear_algebra,ONLY:USE_LK,USE_SLK
 use parallel_m,    ONLY:ncpu
 use interfaces,    ONLY:MATRIX_transfer
 use y_memory_alloc
 !
 implicit none
 !
 integer,          optional :: N
 complex(SP),      optional :: M(:,:)
 type(PAR_matrix), optional :: M_slk
 complex(SP),      optional :: V_left(:,:)
 complex(SP),      optional :: V_right(:,:)
 real(SP),         optional :: E_real(:)
 complex(SP),      optional :: E_cmpl(:)
 integer,          optional :: neig
 !
 ! Workspace
 !
 integer   :: lib_in_use
 real(SP), allocatable :: E_raux(:)
 type (PAR_matrix)   :: M_loc,V_loc
 complex(SP), allocatable :: M_tmp(:,:)
 !
 lib_in_use=USE_LK
#if defined _SCALAPACK 
 lib_in_use=USE_SLK
#endif
!
 if (present(E_cmpl).and.ncpu>1) then
   !
   call warning('Parallel non-herm diag not yet implemented. Switching to LK')
   lib_in_use=USE_LK
   !
 endif
 !
 if (lib_in_use==USE_LK) then
   !
   if (present(E_cmpl)) then
     !
     !... non-Hermitian case
     YAMBO_ALLOC( E_raux,(neig) )
     !
     if (present(M_slk)) then
       call SERIAL_PSEUDOHERM_chase(N,M_slk%blc(:,:,M_slk%I),E_raux,V_left,V_right,neig)
     else
       call SERIAL_PSEUDOHERM_chase(N,M,E_raux,V_left,V_right,neig)
     endif
     !
     E_cmpl(1:neig)=cmplx(E_raux(1:neig),kind=SP)
     YAMBO_FREE( E_raux )
     !
   else
     !
     !... Hermitian case
     if (present(M_slk)) then
       call SERIAL_HERMITIAN_chase(N,M_slk%blc(:,:,M_slk%I),E_real,V_right,neig)
     else
       call SERIAL_HERMITIAN_chase(N,M,E_real,V_right,neig)
     endif
     !
   endif
   !
 else
   !
   call MATRIX_init("SLK",  V_loc, N, 1 )
   YAMBO_ALLOC(V_loc%blc,(V_loc%BLCrows(1):V_loc%BLCrows(2),V_loc%BLCcols(1):V_loc%BLCcols(2),V_loc%Nb))
   !
   if (present(E_cmpl)) then
     !
     !... non-Hermitian case
     !call PARALLEL_PSEUDOHERM_chase( M_slk, V_loc, N, E_real, neig )
     !
   else
     !
     !... Hermitian case
     call PARALLEL_HERMITIAN_chase( M_slk, V_loc, N, E_real, neig )
     !
   endif
   !
   YAMBO_ALLOC(M_tmp,(N,N))
   call MATRIX_transfer( M_in = V_loc , M=M_tmp )
   YAMBO_FREE(V_loc%blc)
   call MATRIX_reset(V_loc)
   V_right(:,:)=M(:,:neig)
   !
 endif
 !
end subroutine CHASE_driver
