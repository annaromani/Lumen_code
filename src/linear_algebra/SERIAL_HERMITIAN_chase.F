!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF DS
!
! #CHASE STAMP      
!
! headers
!
#include<y_memory.h>
!
subroutine SERIAL_HERMITIAN_chase(n,M,E,Vout,neig)
 !
 use pars,           ONLY:SP,cZERO,cONE
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,&
&                         LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use timing_m,       ONLY:timing
 use chase_diag
 use y_memory_alloc
 !
 implicit none
 !
 integer,     intent(in)  :: n,neig
 complex(SP), intent(in)  :: M(n,n)
 complex(SP), intent(out) :: Vout(n,neig)
 real(SP),    intent(out) :: E(neig)

 integer     :: nex, deg, ierr
 real(SP)    :: tol
 character   :: mode, opt, qr
 !
 character(64)     :: subname="SERIAL_HERMITIAN_chase"
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP),    allocatable :: Etemp(:)

 call timing('SERIAL_HERMITIAN_chase',OPR='start')
 !
 E=0._SP
 !
 nex = min(n-neig,neig)
 !
 deg = 20
 tol = 1e-6_SP
 mode = 'R'
 opt = 'S'
 qr = 'C'
 !
 YAMBO_ALLOC(Vtemp,(n,neig+nex))
 YAMBO_ALLOC(Etemp,(neig+nex))
 Vtemp = cZERO
 Etemp = 0._SP
 !
 call timing('ChASE_INIT',OPR='start')
 !
#if defined _DOUBLE
 call zchase_init(n, neig, nex, M, n, Vtemp, Etemp, ierr)
#else
 call cchase_init(n, neig, nex, M, n, Vtemp, Etemp, ierr)
#endif
 !
 call timing('ChASE_INIT',OPR='stop')
 !
 if (ierr/=1) call error("ChASE did not init properly")
 !
#if defined _DOUBLE
 call zchase(deg,tol,mode,opt,qr)
#else
 call cchase(deg,tol,mode,opt,qr)
#endif
 ! 
#if defined _DOUBLE
 call zchase_finalize(ierr)
#else
 call cchase_finalize(ierr)
#endif
 !
 if (ierr/=0) call error("ChASE did not finalize properly")
 !
 E(1:neig) = real(Etemp(1:neig),SP)
 Vout(:,1:neig) = Vtemp(:,1:neig)
 !
 YAMBO_FREE(Etemp)
 YAMBO_FREE(Vtemp)
 !
 call timing('SERIAL_HERMITIAN_chase',OPR='stop')
 !
end subroutine
