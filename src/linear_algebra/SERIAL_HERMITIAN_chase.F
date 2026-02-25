!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): AM AF
!
! #CHASE STAMP      
!
! headers
!
#include<y_memory.h>
!
subroutine SERIAL_HERMITIAN_chase(n,M,E,neig)
 !
 use pars,           ONLY:SP,cZERO,cONE
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,&
&                         LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use timing_m,       ONLY:timing
#if defined _DOUBLE
 use chase_diag,     ONLY:zchase_init,zchase,zchase_finalize
#else
 use chase_diag,     ONLY:cchase_init,cchase,cchase_finalize
#endif
 use y_memory_alloc
 !
 implicit none
 !
 integer     :: n,neig
 complex(SP) :: M(n,n)
 real(SP)    :: E(neig)

 integer     :: nev, nex, deg, ierr
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
 nev = min(neig,n)
 nex = min(n-nev,nev)
 !
 deg = 20
 tol = 1e-5_SP
 mode = 'R'
 opt = 'S'
 qr = 'C'
 !
 YAMBO_ALLOC(Vtemp,(n,nev+nex))
 YAMBO_ALLOC(Etemp,(nev+nex))
 Vtemp = cZERO
 Etemp = 0._SP
 !
 call timing('ChASE_INIT',OPR='start')
 !
#if defined _DOUBLE
 call zchase_init(n, nev, nex, M, n, Vtemp, Etemp, ierr)
#else
 call cchase_init(n=n, nev=nev, nex=nex, h=M, ldh=n, v=Vtemp, ritzv=Etemp, init=ierr)
#endif
 !
 call timing('ChASE_INIT',OPR='stop')
 !
 if (ierr/=0) call error("ChASE did not init properly")
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
 E(1:nev) = real(Etemp(1:nev),SP)
 M(:,1:nev) = Vtemp(:,1:nev)
 if (nev<n) M(:,nev+1:) = cZERO
 !
 YAMBO_FREE(Etemp)
 YAMBO_FREE(Vtemp)
 !
 call timing('SERIAL_HERMITIAN_chase',OPR='stop')
 !
end subroutine
