!
! License-Identifier: GPL
!
! Copyright (C) 2022 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DS
!
subroutine X_AVERAGE_setup(what,i_counter,q0)
 !
 use pars,          ONLY:SP,rZERO
 use vec_operate,   ONLY:v_norm,c2a
 use R_lattice,     ONLY:q0_def_norm
 use D_lattice,     ONLY:i_space_inv
 use parser_m,      ONLY:parser
 use X_m,           ONLY:OPTICAL_do_average,n_OPTICAL_dir_to_eval,OPTICAL_averaged_dirs_string
 !
 implicit none
 !
 character(*)           :: what
 integer                :: i_counter
 real(SP),intent(inout) :: q0(3)
 !
 ! Work Space
 !
 logical :: is_there_an_average_to_do
 integer :: ic,idir
 real(SP), SAVE :: q0_save(3)
 !
 if (what=="INIT") then
   !
   OPTICAL_averaged_dirs_string="none"
   !
   call parser("OptDipAverage",is_there_an_average_to_do)
   !
   ! The following is not really needed since X_irredux does again the same operation
   ! However it is kept for back compatibility of old databases, since I've removed
   ! this operation from INIT.F
   !
   q0(:)=q0/v_norm(q0)*q0_def_norm
   q0_save=q0
   !
   if (.not.is_there_an_average_to_do) then
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
   return
   !
 endif
 !
 if (what=="Q-LOOP") then
   n_OPTICAL_dir_to_eval=1
   if (i_counter>1) return
   n_OPTICAL_dir_to_eval=max(count(OPTICAL_do_average/=0),1)
   return
 endif
 !
 if (what=="DIR-LOOP") then
   q0=rZERO
   ic = OPTICAL_do_average(i_counter)
   if(ic>0) idir = +1
   if(ic<0) idir = -1
   q0( abs(ic) )= idir*1._SP
   call c2a(v_in=q0,mode='ka2c')
   q0(:)=q0/v_norm(q0)
   return
 endif
 !
 if (what=="CLOSE") then
   q0=q0_save
   return
 endif
 !
end subroutine X_AVERAGE_setup
