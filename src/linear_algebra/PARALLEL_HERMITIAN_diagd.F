!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
! Copyright (C) 2004 WanT Group
!
! Authors (see AUTHORS file for details): AM AF DS
!
! headers
!
#include<y_memory.h>
!
subroutine PARALLEL_HERMITIAN_diagd(M, V, N, E)
 !
 ! perform the diagonalization by using scalapack
 !
 use pars,           ONLY:SP,cZERO
 use stderr,         ONLY:intc
 use parallel_int,   ONLY:PP_redux_wait
 use SLK_m,          ONLY:SLK_POOL,SLK_ORTHO
 use matrix,         ONLY:PAR_matrix
 use timing_m,       ONLY:timing
 use y_memory_alloc
 !
 implicit none
 !
 integer,          intent(in)    :: N
 type(PAR_matrix), intent(in)    :: M
 type(PAR_matrix), intent(inout) ::  V
 real(SP),         intent(out)   :: E(N)
 !
 character(64)     :: subname="PARALLEL_HERMITIAN_diagd"
 integer           :: lwork,lrwork,liwork,i_fail
 integer           :: ivar(1)
 real(SP)          :: rvar(1)
 complex(SP)       :: cvar(1)
 integer,     allocatable :: v_int(:)
 real(SP),    allocatable :: v_real(:)
 complex(SP), allocatable :: v_cmplx(:)
 !
 E=0.0_SP
 !
 if (SLK_ORTHO%coordinate(1)==-1) return
 !
 call timing('PAR_HERM_diagd',OPR='start')
 !
 lwork=-1
 lrwork=-1
 liwork=-1
 rvar=0._SP
 cvar=cZERO
 !
#if defined _DOUBLE
 call PZHEEVD&
#else
 call PCHEEVD&
#endif
 & ('V','U',N,M%blc(:,:,1),1,1,M%desc,E,V%blc(:,:,1),1,1,V%desc,&
 &          cvar,lwork,rvar,lrwork,ivar,liwork,i_fail)
 !
 if(i_fail/=0) call error(trim(subname)//': P(Z/C)HEED init failed with error '//trim(intc(i_fail)))
 !
 lwork = nint(real(cvar(1)))
 lrwork = nint(real(rvar(1)))
 liwork = ivar(1)
 YAMBO_ALLOC(v_cmplx,(lwork))
 YAMBO_ALLOC(v_real,(lrwork))
 YAMBO_ALLOC(v_int,(liwork))
 !
#if defined _DOUBLE
 call PZHEEVD&
#else
 call PCHEEVD&
#endif
 & ('V','U',N,M%blc(:,:,1),1,1,M%desc,E,V%blc(:,:,1),1,1,V%desc,&
 &          v_cmplx,lwork,v_real,lrwork,v_int,liwork,i_fail)
 !
 if(i_fail/=0) call error(trim(subname)//': P(Z/C)HEEVD diago failed with error '//trim(intc(i_fail)))
 !
 if (SLK_ORTHO%task/=0) E=0.0_SP
 !
 call PP_redux_wait(E,COMM=SLK_POOL%INTRA_comm)
 !
 call timing('PAR_HERM_diagd',OPR='stop')
 !
end subroutine PARALLEL_HERMITIAN_diagd

