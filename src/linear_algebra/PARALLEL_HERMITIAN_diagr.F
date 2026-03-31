!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
! Copyright (C) 2004 WanT Group
!
! Authors (see AUTHORS file for details): CR AF DS
!
! headers
!
#include<y_memory.h>
!
subroutine PARALLEL_HERMITIAN_diagr(M, V, N, E, neig)
 !
 ! perform the diagonalization by using scalapack pheevr
 ! and conputes only neig eigenpairs
 !
 use pars,           ONLY:SP
 use parallel_int,   ONLY:PP_redux_wait
 use SLK_m,          ONLY:SLK_POOL,SLK_ORTHO
 use matrix,         ONLY:PAR_matrix
 use timing_m,       ONLY:timing
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(in)             :: N,neig
 type(PAR_matrix), intent(in)    :: M
 type(PAR_matrix), intent(inout) :: V
 real(SP),         intent(out)   :: E(neig)
 !
 character(64)      :: subname="PARALLEL_HERMITIAN_diagr"
 integer            :: lwork,lrwork,liwork,neigsl,neigsV,info,nprow,npcol,ii
 integer            :: ivar(1)
 real(SP)           :: rvar(1)
 complex(SP)        :: cvar(1)
 real(SP)           :: abstol,orfac
 real(SP),    allocatable :: Etemp(:)
 real(SP),    allocatable :: v_real(:)
 complex(SP), allocatable :: v_cmplx(:)
 integer,     allocatable :: iwork(:)
 !
 E=0.0_SP
 !
 nprow = SLK_ORTHO%grid(1)
 npcol = SLK_ORTHO%grid(2)
 orfac  = 0.001
 !
 if (SLK_ORTHO%coordinate(1)==-1) return
 !
 call timing('PAR_HERM_diagr',OPR='start')
 !
 allocate(Etemp(N))
 !
 lwork =-1
 lrwork=-1
 liwork=-1
 !
#if defined _DOUBLE
 call PZHEEVR&
#else
 call PCHEEVR&
#endif
 & ('V','I','U',N, M%blc(:,:,1), 1,1, M%desc, 1._SP, 1._SP, 1, neig, &
 & neigsl, neigsV, Etemp, V%blc(:,:,1), 1,1, V%desc, cvar, lwork, rvar, lrwork, &
 & ivar, liwork, info)

 !
 lwork = nint(real(cvar(1))*1.0001_SP)
 lrwork = nint(real(rvar(1))*1.0001_SP)
 liwork = ivar(1)
 !
 YAMBO_ALLOC(v_cmplx,(lwork))
 YAMBO_ALLOC(v_real,(lrwork))
 YAMBO_ALLOC(iwork,(liwork))
 !
 Etemp = 0.0_SP
 V%blc(:,:,1) = 0.0
 !
#if defined _DOUBLE
 call PZHEEVR&
#else
 call PCHEEVR&
#endif
 & ('V','I','U',N, M%blc(:,:,1), 1,1, M%desc, 1._SP, 1._SP, 1, neig, &
 & neigsl, neigsV, Etemp, V%blc(:,:,1), 1,1, V%desc, v_cmplx, lwork, v_real, lrwork, &
 & iwork, liwork, info)
 !
 do ii = 1,neigsl
   E(ii) = Etemp(ii)
 enddo
 !
 if(info/=0) call error(trim(subname)//': P(Z/C)HEEVR failed')
 !
 if (SLK_ORTHO%task/=0) E=0.0_SP
 !
 call PP_redux_wait(E,COMM=SLK_POOL%INTRA_comm)
 !
 YAMBO_FREE(v_cmplx)
 YAMBO_FREE(v_real)
 YAMBO_FREE(iwork)
 !
 call timing('PAR_HERM_diagr',OPR='stop')
 !
end subroutine PARALLEL_HERMITIAN_diagr
