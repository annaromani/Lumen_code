!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): CR AF DS
!
! headers
!
#include<y_memory.h>
!
subroutine PARALLEL_HERMITIAN_chase(Mat, V_blc, N, E, neig)
 !
 use pars,           ONLY:SP,cZERO
 use parallel_int,   ONLY:PP_redux_wait
 use SLK_m,          ONLY:SLK_POOL,SLK_ORTHO
 use linear_algebra, ONLY:LINEAR_ALGEBRA_error,LINEAR_ALGEBRA_WS_reset,LALGEBRA_WS
 use BS_solvers,     ONLY:BSS_slepc_tol
 use parallel_m,     ONLY:PAR_COM_WORLD,MPI_COMM_WORLD
 use parallel_int,   ONLY:PP_redux_wait
 use matrix,         ONLY:PAR_matrix
 use timing_m,       ONLY:timing
 use chase_diag
 use y_memory_alloc
 !
 implicit none
 !
 integer,          intent(in)   :: N,neig
 type(PAR_matrix), intent(in)   :: Mat
 complex(SP),      intent(out)  :: V_blc(N,neig)
 real(SP),         intent(out)  :: E(neig)

 integer   :: i,j,k,nr,nb,nprow,npcol,m,myrow,mycol,nl,rows(2)
 integer   :: nex, ierr, deg
 character :: mode, opt, qr, grid_major
 !
 character(64)     :: subname="PARALLEL_HERMITIAN_chase"
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP),    allocatable :: Etemp(:)
 !
 integer :: myrank_intra, myrank_inter, intra_size, inter_size,&
&           myrank_slk_ortho_comm, slk_ortho_comm_size

 call timing('PARALLEL_HERMITIAN_chase',OPR='start')
 !
 nex = neig
 !
 deg = 20
 mode = 'R'
 opt = 'S'
 qr = 'C'
 grid_major = 'R' !----- Check grid major /!\. -- Checked in SLK_ORTHO_init
 !
 nprow = SLK_ORTHO%grid(1)
 npcol = SLK_ORTHO%grid(2)
 myrow = SLK_ORTHO%coordinate(1)
 mycol = SLK_ORTHO%coordinate(2)
 nb    = Mat%BLCnrows
 m     = Mat%nrows
 nl    = Mat%ncols
 rows  = Mat%BLCrows
 !
 allocate(Vtemp(m,neig+nex))
 allocate(Etemp(neig+nex))
 !
 Vtemp = 0.0
 Etemp = 0.0
 ! 
 call mpi_comm_rank(MPI_COMM_WORLD,myrank_slk_ortho_comm,ierr)
 call mpi_comm_size(MPI_COMM_WORLD,slk_ortho_comm_size  ,ierr)
 ! 
 call mpi_barrier(MPI_COMM_WORLD,ierr)
 !
 call timing('ChASE_INIT',OPR='start')
 !
#if defined _DOUBLE
 !
 ! call pzchase_init(N,neig,nex,m,nl,Mat%blc(:,:,Mat%I),m,Vtemp,Etemp,&
 !&                  nprow,npcol,grid_major,comm,init)
 !
 call pzchase_init_blockcyclic(N, neig, nex, nb, nb, Mat%blc(:,:,1), m, Vtemp, Etemp,&
&                              nprow, npcol, grid_major, 0, 0, MPI_COMM_WORLD, ierr)
#else
 call pcchase_init_blockcyclic(N, neig, nex, nb, nb, Mat%blc(:,:,1), m, Vtemp, Etemp,&
&                              nprow, npcol, grid_major, 0, 0, MPI_COMM_WORLD, ierr)
#endif
 call timing('ChASE_INIT',OPR='stop')
 !
 if (ierr/=1) call warning("ChASE did not init properly")

#if defined _DOUBLE
 call pzchase(deg,BSS_slepc_tol,mode,opt,qr)
#else
 call pcchase(deg,BSS_slepc_tol,mode,opt,qr)
#endif
 ! 
#if defined _DOUBLE
 call pzchase_finalize(ierr)
#else
 call pcchase_finalize(ierr)
#endif
 !
 if (ierr/=0) call warning("ChASE did not finalize properly")
 !
 call timing('ChASE_Reorganize',OPR='start')
 !
 E(1:neig) = Etemp(1:neig)
 !
 V_blc=cZERO
 if (neig > nl * mycol) then
   do i = 1,min(n,neig-n*mycol)
     V_blc(rows(1):rows(2),i) = Vtemp(:, nl * mycol + i)
   enddo
 endif
 call PP_redux_wait(V_blc)
 ! 
 call timing('ChASE_Reorganize',OPR='stop')
 !
 deallocate(Etemp)
 deallocate(Vtemp)
 call timing('PARALLEL_HERMITIAN_chase',OPR='stop')
 !
end subroutine PARALLEL_HERMITIAN_chase
