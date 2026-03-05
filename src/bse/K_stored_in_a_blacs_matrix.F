!
! License-Identifier: GPL
!
! Copyright (C) 2024 The Yambo Team
!
! Authors (see AUTHORS file for details): NM DS
!
! headers
!
#include <y_memory.h>
!
#if defined _LDIAGO
subroutine K_stored_in_a_blacs_matrix_ldiago(i_BS_mat, diago_mat, mpicxt, blacs_blk_size )
  !
  ! Here we store the BSE kernel into a blacs matrix using K_fill_block_cyclic
  !
  use pars,           ONLY:cI,cONE, SP
  use BS,             ONLY:BS_K_dim,BS_H_dim,BS_blk,n_BS_blks,BS_K_coupling,&
    &                        BS_res_ares_n_mat,l_BS_ares_from_res
  use parallel_m,     ONLY:MPI_COMM_WORLD
  use com,            ONLY:msg
  use ldiago_interface
#if defined _MPI
  USE mpi
#endif
  !
  implicit none
  !
  integer,             intent(in)   :: i_BS_mat
  type(c_ptr),         intent(out)  :: mpicxt, diago_mat
  integer(YDIAGO_INT), intent(out)  :: blacs_blk_size 
  !
  ! FIX ME : Move to a better place (modules)
  integer, parameter                 :: blacs_block_size = 64
  ! blacs block size
  ! 64 is good default but give a variable so advanced expert users can play with it
  !
  ! Local variables
  integer                            :: nProcs, my_rank, neig_tmp
  integer(YDIAGO_INT)                :: SL_H_dim, ProcX, ProcY
  integer                            :: ierr !! mpi error code 
  !
  ! Start the function
  !
  if(     BS_K_coupling) SL_H_dim=BS_H_dim
  if(.not.BS_K_coupling) SL_H_dim=BS_K_dim(i_BS_mat)
  !
  ! Allocate the 2D block cyclic matrix
  !
#if defined _MPI
  ! First create a blacs grid
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)

  ProcX = int(sqrt(real(nProcs)),kind=YDIAGO_INT)
  ProcY = ProcX
  ! Maximize the number of cpus participating
  if ( ProcX * (ProcX + 1) .le. nProcs) ProcY = ProcX + 1
  if ( ProcX * (ProcX + 2) .le. nProcs) ProcY = ProcX + 2
  ! Geev block size must be smaller than mat dim and >=8
  ! Start with the default 64
  blacs_blk_size = blacs_block_size
  ! if not try 32
  if (blacs_blk_size .gt. SL_H_dim/ProcY) blacs_blk_size = blacs_blk_size/2
  ! if not try 16
  if (blacs_blk_size .gt. SL_H_dim/ProcY) blacs_blk_size = blacs_blk_size/2
  ! Last try.
  if (blacs_blk_size .gt. SL_H_dim/ProcY) then
    blacs_blk_size = 8 !SL_H_dim
    ProcX = 1
    ProcY = 1
  endif

  call msg("s","BLACS grid",(/ProcX,ProcY/))
#else
  ProcX=1
  ProcY=1
  blacs_blk_size = 1
#endif


  mpicxt = BLACScxtInit_Fortran('R', MPI_COMM_WORLD, ProcX, ProcY)
  if (.not. c_associated(mpicxt)) then
    call error("Failed to initiate BLACS context")
  end if

  diago_mat = init_D_Matrix(SL_H_dim, SL_H_dim, blacs_blk_size, blacs_blk_size, mpicxt)
  if (.not. c_associated(diago_mat)) then
    call error("Failed to initiate block cyclic BSE matrix")
  end if
  !
  ! Fill the block cyclic matrix
  call K_fill_block_cyclic(i_BS_mat, diago_mat)
  !
end subroutine K_stored_in_a_blacs_matrix_ldiago
#endif
!
!
#if defined _SCALAPACK
subroutine K_stored_in_a_blacs_matrix_yambo(M_blacs)
 !
 use pars,           ONLY:cZERO
 use matrix,         ONLY:PAR_matrix
 use linear_algebra, ONLY:diago
 use interfaces,     ONLY:MATRIX_transfer
 use timing_m,       ONLY:timing
 use parallel_m,     ONLY:ncpu
 use BS,             ONLY:BS_blk,BS_H_dim,BS_Blocks_symmetrize_K
 use y_memory_alloc
 !
 implicit none
 !
 type(PAR_matrix) :: M_blacs
 !
 ! workspace
 !
 integer :: rows(2)
 integer :: cols(2)
 !
 call timing('K_stored_in_a_blacs_matrix',OPR='start')
 !
 call PARALLEL_assign_LA_COMMs("Response_T_space",DIAGO,ncpu)
 !
 call BS_Blocks_symmetrize_K()
 ! 
 ! SLK matrix initialization before the Kernel loading
 !
 ! WORKS always
 call MATRIX_init( "SLK", M_blacs, BS_H_dim, 1 )
 M_blacs%kind="SLK"
 ! WORKS only for ncpu>1
 !call SLK_MATRIX_init( M_blacs, BS_H_dim, 1 )
 rows=M_blacs%BLCrows
 cols=M_blacs%BLCcols
 YAMBO_ALLOC(M_blacs%blc,(rows(1):rows(2),cols(1):cols(2),1))
 M_blacs%blc=cZERO
 ! 
 ! WORKS always
 call MATRIX_transfer(M_out=M_blacs, M_bse= BS_blk )
 ! WORKS only for ncpu>1
 !call MATRIX_BSE_to_PAR(BS_blk,M_blacs)
 !
 call K_add_E_diagonal_slk(M_blacs)
 !
 call timing('K_stored_in_a_blacs_matrix',OPR='stop')
 !
end subroutine K_stored_in_a_blacs_matrix_yambo
#endif
