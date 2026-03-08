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
 use parallel_m,     ONLY:myid,MPI_COMM_WORLD
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

 integer   :: i,j,k,nr,nbr,nbc,nprow,npcol,myrow,mycol,ncols,nrows,rows(2),cols(2)
 integer   :: nex, ierr, deg, nmax, nelem
 integer   :: comm
 character :: mode, opt, qr, grid_major
 !
 character(64)     :: subname="PARALLEL_HERMITIAN_chase"
 complex(SP), allocatable :: Vtemp(:,:)
 real(SP),    allocatable :: Etemp(:)
 !
 integer :: myrank_intra, myrank_inter, intra_size, inter_size,&
&           myrank_slk_ortho_comm, slk_ortho_comm_size

 V_blc=cZERO
 if (SLK_ORTHO%coordinate(1)==-1) return

 call timing('PARALLEL_HERMITIAN_chase',OPR='start')
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
 nbr   = Mat%BLCnrows
 nbc   = Mat%BLCncols
 nrows = Mat%nrows
 ncols = Mat%ncols
 rows  = Mat%rows
 cols  = Mat%cols
 comm  = SLK_ORTHO%COMM
 !
 nex = neig
 !
 allocate(Vtemp(nrows,neig+nex))
 allocate(Etemp(neig+nex))
 !
 Vtemp = 0.0
 Etemp = 0.0
 ! 
 call mpi_comm_rank(comm,myrank_slk_ortho_comm,ierr)
 call mpi_comm_size(comm,slk_ortho_comm_size  ,ierr)
 !write(*,*) "mpi comm out ",myrank_slk_ortho_comm,slk_ortho_comm_size
 ! 
 call mpi_barrier(comm,ierr)
 !write(*,*) "mpi barrier overcome ",myrank_slk_ortho_comm,myid
 !write(*,*) myid,comm,MPI_COMM_WORLD,rows,cols,nprow,npcol
 !write(100+myid,*) rows,cols,comm
 !write(100+myid,*) N, neig, nex, nbr, nbc, nprow, npcol
 !do i=cols(1),cols(2)
 !  write(100+myid,*) Mat%blc(rows(1):rows(2),i,1)
 !enddo
 !
 call timing('ChASE init',OPR='start')
 !
#if defined _DOUBLE
 !
 ! call pzchase_init(N,neig,nex,m,nl,Mat%blc(:,:,Mat%I),m,Vtemp,Etemp,&
 !&                  nprow,npcol,grid_major,comm,init)
 !
 call pzchase_init_blockcyclic(N, neig, nex, nbr, nbc, Mat%blc(:,:,1), nbr, Vtemp, Etemp,&
 &                              nprow, npcol, grid_major, 0, 0, comm, ierr)
#else
 !call pcchase_init_blockcyclic(N, neig, nex, nrows, ncols, &
 !&                  Mat%blc(rows(1):rows(2),cols(1):cols(2),1), nrows, Vtemp, Etemp,&
 call pcchase_init_blockcyclic(N, neig, nex, nbr, nbc, Mat%blc(:,:,1), nbr, Vtemp, Etemp,&
 &                              nprow, npcol, grid_major, 0, 0, comm, ierr)
#endif
 call timing('ChASE init',OPR='stop')
 !
 !write(*,*) "init block_cyclic done "
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
 call timing('ChASE redux',OPR='start')
 !
 E(1:neig) = Etemp(1:neig)
 !
 if (neig >= cols(1)) then
   nmax=min(neig,cols(2))
   nelem=nmax-cols(1)+1
   V_blc(rows(1):rows(2),cols(1):nmax) = Vtemp(:nrows,cols(1):nmax)
 endif
 ! 
 call timing('ChASE redux',OPR='stop')
 !
 deallocate(Etemp)
 deallocate(Vtemp)
 call timing('PARALLEL_HERMITIAN_chase',OPR='stop')
 !
end subroutine PARALLEL_HERMITIAN_chase
