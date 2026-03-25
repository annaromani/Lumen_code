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
subroutine CHASE_driver(N,neig,M_slk,V_right,V_left,E_cmpl,E_real)
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
 integer   :: lib_in_use
 real(SP), allocatable :: E_raux(:)
 !
 lib_in_use=USE_LK
#if defined _SCALAPACK 
 if(ncpu>1) lib_in_use=USE_SLK
 !lib_in_use=USE_SLK
#endif
!
 if (present(E_cmpl).and.ncpu>1) then
   !
   call warning('Parallel non-herm chase not yet implemented. Switching to LK')
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
     call SERIAL_PSEUDOHERM_chase(N,M_slk%blc(:,:,M_slk%I),E_raux,V_left,V_right,neig)
     !
     E_cmpl(1:neig)=cmplx(E_raux(1:neig),kind=SP)
     YAMBO_FREE( E_raux )
     !
   else
     !
     !... Hermitian case
     call SERIAL_HERMITIAN_chase(N,M_slk%blc(:,:,1),E_real,V_right,neig)
     !
   endif
   !
 else
   !
   if (present(E_cmpl)) then
     !
     !... non-Hermitian case
     !call PARALLEL_PSEUDOHERM_chase( M_slk, V_loc, N, E_real, neig )
     !
   else
     !
     !... Hermitian case
     call PARALLEL_HERMITIAN_chase( M_slk, V_right, N, E_real, neig)
     !
     call PP_redux_wait(V_right)
     !
     if (ncpu>SLK_ORTHO%n_CPU) then
       call PP_redux_wait(E_real)
       E_real=E_real/real(SLK_ORTHO%n_CPU,SP)
     endif
     !
   endif
   !
 endif
 !
end subroutine CHASE_driver
