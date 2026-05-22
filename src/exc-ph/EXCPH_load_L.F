!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM FP
!
! headers
!
#include<y_memory.h>
!
! --- Theory ---
! (5) Figure out rotations/symmetry of exciton eigenvectors and EPC m.e.
! - FP
!
subroutine EXCPH_load_L(iq_bse,X,what,Ltype)
 !! **EXCPH_load_L - Load BSE L matrix from ndb.BSS_diago**
 !!
 !! Wrapper around io_BSS_diago that loads exciton eigenvectors and
 !! eigenvalues for a given q-point. Used in excitonic Raman and
 !! exciton-phonon coupling to provide the exciton states |S> needed
 !! to evaluate matrix elements <S|g_nu|S'> and dipoles <0|v|S>.
 !!
 !! Arguments:
 !!   iq_bse  - BSE q-point index
 !!   X       - response function descriptor (passed through to io_BSS_diago)
 !!   what    - controls what is loaded:
 !!               'check'       : open DB, validate dimensions against EXCPH_states/EXCPH_sum,
 !!                               set L_kind_in or L_kind_out from BSE_L_kind
 !!               'eigenvalues' : read eigenvalues BS_E (and residuals if 'residuals' in what)
 !!               'load'        : read full eigenvectors BS_mat + eigenvalues BS_E
 !!                               (and residuals BS_R_left/BS_R_right if 'residuals' in what)
 !!   Ltype   - which L matrix is being loaded:
 !!               'Lin'  : incoming photon leg; sets L_kind_in
 !!               'Lout' : outgoing photon leg; sets L_kind_out, optionally reads
 !!                        from a different job directory (Lout_path)
 !!
 !! Validation (on 'check'):
 !!   - EXCPH_states(2) must be <= BSS_n_eig  (enough eigenstates stored)
 !!   - EXCPH_sum(2)    must be <= BS_H_dim   (virtual states within Hamiltonian size)
 !!   - BSS_write_eig_2_db must be .TRUE.     (eigenvectors present in DB)
 !!
 use pars,           ONLY:lchlen
 use BS_solvers,     ONLY:io_BSS_diago,BS_mat,BSS_n_eig,BSS_free,BSS_write_eig_2_db
 use IO_m,           ONLY:NONE,OP_APP_CL,OP_RD_CL,DUMP,OP_RD,RD_CL_IF_END,OP_WR_CL,deliver_IO_error_message
 use R_lattice,      ONLY:bz_samp
 use BS,             ONLY:BS_H_dim,BSE_L_kind,BS_res_ares_n_mat,BS_K_dim
 use EXCPH,          ONLY:BS_E,BS_R_left,BS_R_right,EXCPH_states,EXCPH_sum,Lout_path,L_kind_in,L_kind_out
 use IO_int,         ONLY:io_control
 use stderr,         ONLY:intc
 use X_m,            ONLY:X_t
 use com,            ONLY:jobstr,msg
 !
 use y_memory_alloc
 !
 implicit none
 ! 
 type(X_t)     ::X
 type(bz_samp) ::k,q
 integer, intent(in) :: iq_bse     
 character(*)        :: what,Ltype
 !
 ! Work space
 !
 integer :: ID_BS,io_err
 character(lchlen) :: jobstr_save
 logical :: l_read_residuals
 !
 if(index(Ltype,'Lout')>0.and.Lout_path/='none') then
   jobstr_save=jobstr
   jobstr     =trim(Lout_path)
 endif
 !
 l_read_residuals=.FALSE.
 if(index(what,'residuals')>0) l_read_residuals=.TRUE.
 !
 if(index(what,'check')>0) then
   !
   call io_control(ACTION=OP_RD_CL,COM=NONE,MODE=DUMP,SEC=(/1/),ID=ID_BS)
   io_err=io_BSS_diago(iq_bse,1,ID_BS,X)
   if(io_err/=0) call error(" Error checking "//trim(Ltype)//" for q = "//trim(intc(iq_bse)))
   !
   !if (BS_res_ares_n_mat/=1) call error("Excitonic Gkkp: anti-resonant part not implemented")
   !
   if(.not.BSS_write_eig_2_db) call error("BSS diago DB does not contain the eigenstates")
   !
   BS_H_dim=BS_K_dim(1)
   !
   if(EXCPH_states(2)>BSS_n_eig)  &
&     call error("Number of exc. states in the excphon coupling must         <= number of BSE eigenstates")
   !
   if(EXCPH_sum(2)>BS_H_dim)      &
&     call error("Number of virtual exc. states in the excphon coupling must <= size of the exc. Hamiltonian")
   !
   if(iq_bse==1) then
     call msg('rs','Found '//trim(Ltype)//' with dimension: '//trim(intc(BSS_n_eig))//' x '//trim(intc(BS_H_dim)))
     if(index(Ltype,'Lout')>0) then
       L_kind_out=BSE_L_kind
       call msg('rs','L kind out: '//trim(L_kind_out))
     elseif(index(Ltype,'Lin')>0) then
       L_kind_in =BSE_L_kind
       call msg('rs','L kind  in: '//trim(L_kind_in))
     endif
   endif
   !
 elseif(index(what,'eigenvalues')>0) then
   !
   ! Read eigenvalues and optionally residuals
   !
   call io_control(ACTION=OP_RD_CL,COM=NONE,SEC=(/2/),ID=ID_BS)
   if(l_read_residuals) then
      io_err=io_BSS_diago(iq_bse,1,ID_BS,X,bsE=BS_E,bsRr=BS_R_right,bsRl=BS_R_left)
   else
      io_err=io_BSS_diago(iq_bse,1,ID_BS,X,bsE=BS_E)
   endif
   if(io_err/=0) call error(" Error loading "//trim(Ltype)//" for q = "//trim(intc(iq_bse)))
   !
 elseif(index(what,'load')>0) then
   !      
   ! Read BSE eigenvectors and eigenvalues and optionally residuals
   !
   call io_control(ACTION=OP_RD_CL,COM=NONE,SEC=(/2,3/),ID=ID_BS)
   if(l_read_residuals) then
      io_err=io_BSS_diago(iq_bse,1,ID_BS,X,bsE=BS_E,BS_mat=BS_mat,bsRr=BS_R_right,bsRl=BS_R_left)
   else
      io_err=io_BSS_diago(iq_bse,1,ID_BS,X,bsE=BS_E,BS_mat=BS_mat)
   endif
   if(io_err/=0) call error(" Error loading "//trim(Ltype)//" for q = "//trim(intc(iq_bse)))
   !
 else
   call error(" Wrong call to the EXCPH_load_L function! ")
 endif
 !
 if(index(Ltype,'Lout')>0.and.Lout_path/='none') then
   jobstr     =jobstr_save
 endif
 !
end subroutine EXCPH_load_L
