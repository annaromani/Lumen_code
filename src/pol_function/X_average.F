!
! License-Identifier: GPL
!
! Copyright (C) 2022 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS MN
!
! headers
!
#include<dev_defs.h>
!
subroutine X_average_init(q0,syms_average,dirs_average)
 !
 use pars,          ONLY:SP,rZERO
 use R_lattice,     ONLY:q0_def_norm
 use vec_operate,   ONLY:v_norm
 use D_lattice,     ONLY:i_space_inv
 use parser_m,      ONLY:parser
 use X_m,           ONLY:OPTICAL_do_average,n_OPTICAL_dir_to_eval,OPTICAL_averaged_dirs_string
 !
 implicit none
 !
 real(SP),intent(inout) :: q0(3)
 logical, intent(out)   :: syms_average,dirs_average
 !
 ! Work Space
 !
 integer :: ic
 !
 OPTICAL_averaged_dirs_string="none"
 !
 call parser("OptDipAverage",dirs_average)
 call parser("OptSymAverage",syms_average)
 !
 ! The following is not really needed since X_irredux does again the same operation
 ! However it is kept for back compatibility of old databases, since I've removed
 ! this operation from INIT.F
 !
 q0(:)=q0/v_norm(q0)*q0_def_norm
 !
 if (.not.dirs_average) then
   n_OPTICAL_dir_to_eval=1
   return
 endif
 !
 n_OPTICAL_dir_to_eval=0
 OPTICAL_averaged_dirs_string=""
 do ic=1,3
   if (abs(q0(ic))<=0._SP) cycle
   if (ic==1) OPTICAL_averaged_dirs_string=trim(OPTICAL_averaged_dirs_string)//'1'
   if (ic==2) OPTICAL_averaged_dirs_string=trim(OPTICAL_averaged_dirs_string)//'2'
   if (ic==3) OPTICAL_averaged_dirs_string=trim(OPTICAL_averaged_dirs_string)//'3'
   n_OPTICAL_dir_to_eval=n_OPTICAL_dir_to_eval+1
   OPTICAL_do_average(n_OPTICAL_dir_to_eval)=ic
 enddo
 !
 if (i_space_inv==0) then
   do ic=1,3
     if (abs(q0(ic))<=0._SP) cycle
     n_OPTICAL_dir_to_eval=n_OPTICAL_dir_to_eval+1
     OPTICAL_do_average(n_OPTICAL_dir_to_eval)=-ic
   enddo
 endif
 !
end subroutine X_average_init
!
!
subroutine X_average_q0_def(iq_dir,q0,X_what)
 !
 use pars,          ONLY:SP,schlen,rZERO
 use stderr,        ONLY:intc
 use vec_operate,   ONLY:v_norm,c2a
 use X_m,           ONLY:OPTICAL_do_average,n_OPTICAL_dir_to_eval
 !
 implicit none
 !
 integer,           intent(in)  :: iq_dir
 real(SP),          intent(out) :: q0(3)
 character(schlen), intent(out) :: X_what
 !
 integer :: ic,idir
 !
 X_what="("//trim(intc(iq_dir))//" of "//trim(intc(n_OPTICAL_dir_to_eval))//")"
 !
 q0=rZERO
 ic = OPTICAL_do_average(iq_dir)
 if(ic>0) idir = +1
 if(ic<0) idir = -1
 q0( abs(ic) )= idir*1._SP
 call c2a(v_in=q0,mode='ka2c')
 q0=q0/v_norm(q0)
 !
end subroutine X_average_q0_def
!
!
subroutine X_average_q0(X_par,Xen,Xk,Xw,X,Dip)
 !
 use pars,          ONLY:SP,rZERO,cZERO,schlen
 use drivers,       ONLY:l_rpa_IP 
 use matrix,        ONLY:PAR_matrix
 use electrons,     ONLY:levels
 use DIPOLES,       ONLY:DIPOLE_t
 use frequency,     ONLY:w_samp
 use R_lattice,     ONLY:bz_samp
 !use IO_m,          ONLY:io_RESPONSE
 use matrix,        ONLY:MATRIX_copy!,MATRIX_reset
 use X_m,           ONLY:OPTICAL_do_average,n_OPTICAL_dir_to_eval,X_t
 !
 implicit none
 !
 type(PAR_matrix), intent(inout) :: X_par
 type(levels)  :: Xen 
 type(bz_samp) :: Xk,q
 type(X_t)     :: X
 type(DIPOLE_t):: Dip
 type(w_samp)  :: Xw
 !
 ! Workspace
 !
 type(PAR_matrix)  :: X_par_average
 integer           :: iq_dir
 character(schlen) :: X_what
 !
 if (n_OPTICAL_dir_to_eval==1) return
 !
 call MATRIX_copy(X_par,X_par_average,.TRUE.)
 !
 X_par_average%blc=X_par%blc
 X_par%blc=cZERO
 !
 !if (io_RESPONSE) call MATRIX_reset(X_par)
 !
 do iq_dir=2,n_OPTICAL_dir_to_eval
   !
   call X_average_q0_def(iq_dir,X%q0,X_what)
   !
   call X_irredux(1,"Xo"//trim(X_what),X_par,Xen,Xk,Xw,X,Dip)
   if (.not.l_rpa_IP) call X_redux(1,"X"//trim(X_what),X_par,Xw,X)
   !
   ! AF: even with GPU support, accumulation is done on the host
   !     X_par%blc is assumed to be up-to-date
   !
   X_par_average%blc=X_par_average%blc+X_par%blc
   X_par%blc=cZERO
   !
   !if (io_RESPONSE.and.iq_dir<n_OPTICAL_dir_to_eval) call MATRIX_reset(X_par)
   !
 enddo
 !
 X_par%blc=X_par_average%blc/real(n_OPTICAL_dir_to_eval,SP)
 !
end subroutine X_average_q0
!
!
subroutine X_symmetrize_q0(X_par, NG, Xw)
 !
 use pars,          ONLY: SP, DP, cZERO
 use matrix,        ONLY: PAR_matrix
 use frequency,     ONLY: w_samp
 use D_lattice,     ONLY: nsym, i_time_rev
 use R_lattice,     ONLY: g_rot
 use parallel_int,  ONLY: PP_redux_wait
 use com,           ONLY: msg

 implicit none

 type(PAR_matrix), intent(inout) :: X_par
 integer,          intent(in)    :: NG
 type(w_samp),     intent(in)    :: Xw

 integer :: ig1, ig2, ig1_sym, ig2_sym, isym, iw
 integer :: nsym_loop, ig1_loc
 complex(SP), allocatable :: X_full(:,:)
 !
 if (nsym <= 1) return
 !
 call msg('r', 'Applying q->0 symmetries to the macroscopic dielectric matrix')
 !
 allocate(X_full(NG, NG))
 !
 do iw = 1, Xw%n_freqs
    X_full = cZERO
    if (X_par%nrows>0) X_full(X_par%rows(1):X_par%rows(2), 1:NG) = X_par%blc(1:X_par%nrows, 1:NG, iw)
    call PP_redux_wait(X_full, COMM=X_par%INTER_comm%COMM)
    if (X_par%nrows>0) X_par%blc(1:X_par%nrows, 1:NG, iw) = cZERO
    ! Note :
    ! Time reversal maps (X(q,...,E) to X(q,...,-E)
    ! so we only apply to E = 0 case
    if (abs(Xw%p(iw)) < 1.e-5_SP) then
       nsym_loop = nsym
    else
       nsym_loop = nsym / (i_time_rev + 1)
    endif
    do ig2 = 1, NG
       do ig1_loc = 1, X_par%nrows
          ig1 = ig1_loc + X_par%rows(1) - 1
          do isym = 1, nsym_loop
             ig1_sym = g_rot(ig1, isym)
             ig2_sym = g_rot(ig2, isym)
             !
             if (ig1_sym > 0 .and. ig1_sym <= NG .and. ig2_sym > 0 .and. ig2_sym <= NG) then
                if (isym > nsym / (i_time_rev + 1)) then
                   X_par%blc(ig1_loc, ig2, iw) = X_par%blc(ig1_loc, ig2, iw) + &
                        conjg(X_full(ig1_sym, ig2_sym)) / real(nsym_loop, SP)
                else
                   X_par%blc(ig1_loc, ig2, iw) = X_par%blc(ig1_loc, ig2, iw) + &
                        X_full(ig1_sym, ig2_sym) / real(nsym_loop, SP)
                endif
             endif
          enddo
       enddo
    enddo
 enddo
 !
 deallocate(X_full)
 !
end subroutine X_symmetrize_q0
