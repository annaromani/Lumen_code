! 
! License-Identifier: GPL
!
! Copyright (C) 2023 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
! headers
!
#include<y_memory.h>
!
subroutine excitons_symmetries(iq,Xen,Xk,qpt,BS_E,Res)
 !
 use pars,           ONLY:SP,DP,cZERO,schlen
 use units,          ONLY:HA2EV
 use parser_m,       ONLY:parser
 use YPPm,           ONLY:deg_energy,EXCITONS_n_user_states,EXCITONS_user_indexes, &
 &                        mix_thresh,symm_thresh
 use R_lattice,      ONLY:bz_samp!,coll_scattering
 use electrons,      ONLY:levels,n_sp_pol
 use BS_solvers,     ONLY:BSS_n_eig
 use BS,             ONLY:BS_bands
 use wave_func,      ONLY:io_WF_phases,WF_phases_b_map,WF_phases
 use stderr,         ONLY:intc
 use D_lattice,      ONLY:i_time_rev,nsym,nrot,dl_sop,dl_tra
 use vec_operate,    ONLY:sort,degeneration_finder
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:DUMP,REP,OP_RD_CL,OP_RD,RD_CL,RD
 use com,            ONLY:msg,of_open_close
 !
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(in) :: iq
 type(bz_samp)       :: Xk,qpt
 type(levels)        :: Xen
 complex(SP)         :: BS_E(BSS_n_eig)
 real(SP)            :: Res(BSS_n_eig)
 !
 real(SP) :: BS_E_real(EXCITONS_n_user_states),Res_sum
 logical  :: l_min_mem,l_mix_thresh_set
 integer  :: S_index(EXCITONS_n_user_states)
 integer  :: first_exc(EXCITONS_n_user_states),n_of_exc(EXCITONS_n_user_states),n_deg_grp,n_of_exc_max,iexc_ref
 integer  :: ik,i_sp_pol,iqbz,nsymq,nrotq,nsz(5),iexc,il1,i_grp,i_sub_grp,io_err,ID,n_exc_sbgrp
 !
 complex(SP), allocatable :: exc_trace(:,:,:)
 complex(SP), allocatable :: exc_reps(:,:,:,:)
 !
 integer, allocatable  :: exc_sbgr(:,:)
 integer, allocatable  :: n_sub_grps(:)
 integer, allocatable  :: ref_exc_grp(:,:)
 integer, allocatable  :: ref_exc_sbgr(:,:)
 !
 character(schlen)    :: filebse,titles(5),loc_msg,loc_fmt
 character(schlen), allocatable :: exc_character(:,:)
 !
 integer, external :: io_QINDX
 !
 ! Electrons part (pre-requisite)
 !
 call k_build_up_BZ_tables(Xk)
 call k_ibz2bz(Xk,'iku',.true.)
 call k_small_group(Xk)
 !
 call WF_phase_matrices(Xen,Xk,BS_bands,(/1,Xk%nibz/))
 !
 nsz=0
 YAMBO_ALLOC(WF_phases_b_map,(BS_bands(1):BS_bands(2),Xk%nibz,n_sp_pol,2))
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1,2/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,0,ID,nsz)
 !
 call parser("MinimizeMem",l_min_mem)
 call parser("Mix_thresh",l_mix_thresh_set)
 if (.not.l_mix_thresh_set) mix_thresh=symm_thresh 
 !
 if (.not.l_min_mem) then
   nsz(1:2)=maxval(WF_phases_b_map(:,:,:,1))
   nsz(3)  =maxval(WF_phases_b_map(:,:,:,2))
   nsz(4)  =nsym
   nsz(5)  =Xk%nbz
   YAMBO_ALLOC(WF_phases,(nsz(1),nsz(2),nsz(3),nsz(4),nsz(5)))
   !
   do i_sp_pol=1,n_sp_pol
     do ik=1,Xk%nibz
       nsz(1:2)=maxval(WF_phases_b_map(:,ik,i_sp_pol,1))
       nsz(3)  =maxval(WF_phases_b_map(:,ik,i_sp_pol,2))
       nsz(4)  =nsym
       nsz(5)  =Xk%nstar(ik)
       il1=1
       if(ik>0) il1= sum(Xk%nstar(:ik-1))+1
       call io_control(ACTION=RD,COM=REP,SEC=(/3/),MODE=DUMP,ID=ID)
       io_err=io_WF_phases( BS_bands,ik,0,i_sp_pol,ID,nsz, &
       & WF_phases(:nsz(1),:nsz(2),:nsz(3),:,il1:il1+nsz(5)-1) )
     enddo
   enddo
 endif
 !
 call build_inverse_BS_eh_table_from_BSS_table(Xk%nbz)
 !
 ! Exciton part
 !
 call k_build_up_BZ_tables(qpt)
 call k_ibz2bz(qpt,'iku',.true.)
 call k_small_group(qpt)
 !
 !coll_scattering=.true.
 !call io_control(ACTION=OP_RD_CL,COM=REP,SEC=(/6/),ID=ID)
 !io_err=io_QINDX(Xk,qpt,ID,'minus_q')
 !if (io_err/=0) then
 !  call msg('s',"Missing k/q scattering database")
 !  call error('Please (re)run the setup K_grids="C"')
 !endif
 !
 ! The q point for which I want to compute the irrep.
 !isq =q%star(iq,1)
 iqbz=qpt%k_table(iq,1)
 !iGq =qpt%g_table(iq,isq)
 !
 iexc_ref=EXCITONS_user_indexes(1)
 do iexc=1,EXCITONS_n_user_states
   il1=EXCITONS_user_indexes(iexc)
   BS_E_real(iexc)=real(BS_E(il1))
 enddo
 !
 call sort(BS_E_real,indx=S_index)
 !
 call degeneration_finder(EXCITONS_n_user_states,first_exc,n_of_exc,n_deg_grp,Er=BS_E_real,&
 &                                     deg_accuracy=deg_energy,Include_single_values=.true.)
 !
 n_of_exc_max=maxval(n_of_exc)
 !n_of_exc_max=3
 !
 call msg('s',"Scan of eneriges gave number of reps ",n_deg_grp)
 !
 nsymq=qpt%grp_nsym(iqbz)
 nrotq=nsymq
 if ( any( qpt%grp_table(:nrotq,iqbz)>nrot ) ) nrotq=nrotq/(1+i_time_rev)
 !
 !write(*,*) qpt%grp_table(:nrotq,iqbz)
 !
 ! Detect the point group
 !
 call detect_point_group(nrotq,real(dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),DP))
 !
 filebse="exc_qpt"//trim(intc(iq))//"_E+syms_sorted"
 !
 call of_open_close(filebse,'ot')
 call msg("o "//trim(filebse),"# Symmetries at qpt ",qpt%pt(iq,:),INDENT=0)
 call msg("o "//trim(filebse),"#")
 !
 call print_point_group("o "//trim(filebse),nrotq,nsymq,dl_sop(:,:,qpt%grp_table(:nsymq,iqbz)),&
 &                                                  dl_tra(:,qpt%grp_table(:nrotq,iqbz)),.true.)
 !
 ! Compute the reps
 YAMBO_ALLOC(exc_reps,(12,n_of_exc_max,n_of_exc_max,n_deg_grp))
 YAMBO_ALLOC(exc_trace,(12,n_of_exc_max,n_deg_grp))
 allocate(exc_character(n_of_exc_max,n_deg_grp))
 exc_reps=0._SP
 exc_trace=0._SP
 exc_character=""
 !
 YAMBO_ALLOC(exc_sbgr,(n_of_exc_max,n_deg_grp))
 YAMBO_ALLOC(n_sub_grps,(n_deg_grp))
 YAMBO_ALLOC(ref_exc_grp,(n_of_exc_max,n_deg_grp))
 YAMBO_ALLOC(ref_exc_sbgr,(n_of_exc_max,n_deg_grp))
 !
 call compute_exc_reps(EXCITONS_n_user_states,n_deg_grp,n_of_exc_max,&
 &          n_of_exc,first_exc,S_index,iqbz,l_min_mem,Xk,qpt,exc_reps)
 !
 call io_control(ACTION=RD_CL,COM=REP,SEC=(/0/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,0,ID,(/0,0,0,0,0/))
 !
 ! Now I refine the number of irreps using the reps.
 ! The first guess was obtained with the degeneration finder
 !
 call refine_groups(EXCITONS_n_user_states,n_deg_grp,n_of_exc_max,n_of_exc,first_exc, &
                        &   n_sub_grps,ref_exc_sbgr,exc_sbgr,ref_exc_grp,exc_reps)
 !
 ! Finally I compute the exciton character
 !
 call get_character_exc(EXCITONS_n_user_states,n_deg_grp,n_of_exc_max,n_of_exc,first_exc, &
                        &   n_sub_grps,ref_exc_sbgr,exc_sbgr,iexc_ref,&
                        &   exc_trace,exc_reps,exc_character,.true.,iq)
 !
 YAMBO_FREE(exc_reps)
 !
 titles=(/"-   E[eV]   -","-   Int.    -","-  i_lambda -","-   n_deg   -","-  exc-symm -"/)
 call msg("o "//trim(filebse),"#")
 call msg("o "//trim(filebse),'#',titles,INDENT=0,USE_TABS=.false.)
 call msg("o "//trim(filebse),"#")
 !
 loc_fmt="(F10.4,4x,E10.4,4x,i8,8x,i4,8x,a)"
 !
 do i_grp=1,n_deg_grp
   do i_sub_grp=1,n_sub_grps(i_grp)
     !
     il1=first_exc(i_grp)+ref_exc_sbgr(i_sub_grp,i_grp)-1
     n_exc_sbgrp=nint(real(exc_trace(1,i_sub_grp,i_grp),SP))
     Res_sum=sum(Res(il1+iexc_ref-1:il1+iexc_ref-1+n_exc_sbgrp-1))
     !write(100,*) BS_E_real(il1)*HA2EV,real(BS_E(il1+iexc_ref-1:il1+iexc_ref-1+n_exc_sbgrp-1))*HA2EV
     !
     write(loc_msg,loc_fmt) BS_E_real(il1)*HA2EV,Res_sum,il1+iexc_ref-1,n_exc_sbgrp,&
     &              trim(exc_character(i_sub_grp,i_grp))
     call msg("o "//trim(filebse),'',loc_msg,INDENT=-2)
     !
   enddo
 enddo
 !
 YAMBO_FREE(exc_trace)
 deallocate(exc_character)
 !
 YAMBO_FREE(exc_sbgr)
 YAMBO_FREE(n_sub_grps)
 YAMBO_FREE(ref_exc_grp)
 YAMBO_FREE(ref_exc_sbgr)
 !
 call of_open_close(filebse)
 !
end subroutine excitons_symmetries
!
!
subroutine compute_exc_reps(nstates,n_deg_grp,n_of_exc_max,n_of_exc,&
&                   first_exc,S_index,iqbz,l_min_mem,Xk,qpt,loc_reps)
 !
 use pars,           ONLY:SP
 use LIVE_t,         ONLY:live_timing
 use interfaces,     ONLY:EXC_WF_symm_qpoint
 use R_lattice,      ONLY:bz_samp
 use BS,             ONLY:BS_H_dim
 use wrapper,        ONLY:Vstar_dot_V
 use D_lattice,      ONLY:nclass,elem
 !
 implicit none
 !
 integer,       intent(in)  :: nstates,n_deg_grp
 integer,       intent(in)  :: n_of_exc_max,n_of_exc(nstates)
 integer,       intent(in)  :: first_exc(nstates)
 integer,       intent(in)  :: S_index(nstates)
 integer,       intent(in)  :: iqbz
 logical,       intent(in)  :: l_min_mem
 complex(SP),   intent(out) :: loc_reps(12,n_of_exc_max,n_of_exc_max,n_deg_grp)
 !
 type(bz_samp), intent(inout) :: qpt,Xk
 !
 ! Workspace
 !
 complex(SP)  :: EXC_WF_rot1(BS_H_dim,n_of_exc_max)
 complex(SP)  :: EXC_WF_rot2(BS_H_dim,n_of_exc_max)
 !
 integer  :: iq,iRqbz,isqc,isqy,isq_small_grp,iGqc,iGqy,nrotq,nsz(5),&
 &           iclass,ielem,iexc,il1,il2,i_grp,i_sub_grp,io_err,ID!,n_elem_used
 !
 !write(*,*) "point group detected ",nclass,n_deg_grp
 !
 ! Compute excitonic reps
 !
 iq=qpt%sstar(iqbz,1)
 !
 !call live_timing("Exciton symmetries:",sum(nelem(:))*n_deg_grp)
 call live_timing("Excitonic reps:",nclass*n_deg_grp)
 !
 do iclass=1,nclass
   ! The following loop is useful for debugging only
   ! Any element of the class should be ok 
   !n_elem_used=nelem(iclass)
   do ielem=1,1
     !
     isq_small_grp=elem(ielem,iclass)
     !
     ! The rotated q-point
     isqc=qpt%grp_table(isq_small_grp,iqbz)       ! the symmetry of the loop
     isqy=qpt%s_table(iq,isqc)                    ! the symmetry used by yambo
     !
     iRqbz=qpt%k_table(iq,isqc)                   ! the rotated q-point
     !
     iGqc=qpt%k_table(iq,isqc)                    ! the G shift
     iGqy=qpt%k_table(iq,isqy)                    ! the G shift
     !
     do i_grp=1,n_deg_grp
       !
       ! Get A(l,R*q)= Ry * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,qpt,first_exc(i_grp),n_of_exc(i_grp),&
       &       iq,isqy,iGqy,ID,S_index,EXC_WF_rot1(:,:n_of_exc(i_grp)),l_min_mem)
       !
       ! Get R * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,qpt,first_exc(i_grp),n_of_exc(i_grp),&
       &       iq,isqc,iGqc,ID,S_index,EXC_WF_rot2(:,:n_of_exc(i_grp)),l_min_mem)
       !
       ! Compute the trace of the Irrep
       !
       do il1=1,n_of_exc(i_grp)
         do il2=1,n_of_exc(i_grp)
           loc_reps(iclass,il2,il1,i_grp) = &
           & Vstar_dot_V(BS_H_dim,EXC_WF_rot1(:,il2),EXC_WF_rot2(:,il1))
         enddo
       enddo
       !
       call live_timing(steps=1)
       !
     enddo
   enddo
 enddo
 !
 call live_timing()
 !
end subroutine compute_exc_reps
!
!
subroutine refine_groups(nstates,n_deg_grp,n_of_exc_max,n_of_exc,first_exc, &
                        &   n_sub_grps,ref_exc_sbgr,exc_sbgr,ref_exc_grp,&
                        &   loc_reps   )
 !
 use pars,           ONLY:SP
 use com,            ONLY:msg
 use D_lattice,      ONLY:nclass
 use YPPm,           ONLY:mix_thresh
 use stderr,         ONLY:intc
 !
 implicit none
 !
 integer,      intent(in)  :: nstates,n_deg_grp,n_of_exc_max
 integer,      intent(in)  :: first_exc(nstates),n_of_exc(nstates)
 integer,      intent(out) :: n_sub_grps(n_deg_grp),ref_exc_sbgr(n_of_exc_max,n_deg_grp)
 integer,      intent(out) :: exc_sbgr(n_of_exc_max,n_deg_grp),ref_exc_grp(n_of_exc_max,n_deg_grp)
 complex(SP),  intent(in)  :: loc_reps(12,n_of_exc_max,n_of_exc_max,n_deg_grp)
 !
 ! Workspace
 !
 integer  :: iexc,jexc,il1,il2,i_grp
 !
 ! From each rep grep the size of the irreps
 !
 ref_exc_grp=-1
 n_sub_grps=0
 !
 do i_grp=1,n_deg_grp
   do il1=1,n_of_exc(i_grp)
     if ( ref_exc_grp(il1,i_grp) > 0 ) cycle
     ref_exc_grp(il1,i_grp)=il1
     iexc=first_exc(i_grp)+il1-1
     do il2=il1+1,n_of_exc(i_grp)
       ref_exc_grp(il2,i_grp)=il1
       jexc=first_exc(i_grp)+il2-1
       !write(*,*) i_grp,il1,il2,loc_reps(:,il2,il1,i_grp)
       call msg('s',"Rep "//trim(intc(i_grp))//" ,Exc. mix "//trim(intc(iexc))//","&
       &                  //trim(intc(jexc)),sum(abs(loc_reps(:,il2,il1,i_grp))))
       if ( sum(abs(loc_reps(:,il2,il1,i_grp)))>=nclass*mix_thresh) cycle
       ref_exc_grp(il2,i_grp)=-1
     enddo
   enddo
   !
   do il1=1,n_of_exc(i_grp)
     if(il1/=ref_exc_grp(il1,i_grp)) cycle
     n_sub_grps(i_grp)=n_sub_grps(i_grp)+1
     exc_sbgr(il1,i_grp)=n_sub_grps(i_grp)
     ref_exc_sbgr(n_sub_grps(i_grp),i_grp)=il1
   enddo
   !
   do il1=1,n_of_exc(i_grp)
     exc_sbgr(il1,i_grp)=exc_sbgr(ref_exc_grp(il1,i_grp),i_grp)
   enddo
   !
   if (n_sub_grps(i_grp)>1) then
     call msg('s',"Rep "//trim(intc(i_grp))//" split into "&
     &              //trim(intc(n_sub_grps(i_grp)))//" irreps")
   endif
   !
 enddo
 !
 call msg('s',"Scan of rep matrices gave number of irreps ",sum(n_sub_grps(:n_deg_grp)))
 !
end subroutine refine_groups
!
!
subroutine get_character_exc(nstates,n_deg_grp,n_of_exc_max,n_of_exc,first_exc, &
                        &   n_sub_grps,ref_exc_sbgr,exc_sbgr,iexc_ref,&
                        &   loc_trace,loc_reps,loc_character,l_prt_debug,iq)
 !
 use pars,           ONLY:SP,schlen,cI,cZERO
 use units,          ONLY:HA2EV
 use YPPm,           ONLY:symm_thresh
 use stderr,         ONLY:intc,real2ch
 use com,            ONLY:msg,of_open_close
 use D_lattice,      ONLY:which_irr,char_mat,nclass,& !nelem,&
 &                        name_class,name_rap
 !
 implicit none
 !
 integer,      intent(in)  :: nstates,n_deg_grp,n_of_exc_max,iexc_ref
 integer,      intent(in)  :: first_exc(nstates),n_of_exc(nstates)
 integer,      intent(in)  :: n_sub_grps(n_deg_grp),ref_exc_sbgr(n_of_exc_max,n_deg_grp)
 integer,      intent(in)  :: exc_sbgr(n_of_exc_max,n_deg_grp)
 complex(SP),  intent(out) :: loc_trace(12,n_of_exc_max,n_deg_grp)
 complex(SP),  intent(in)  :: loc_reps(12,n_of_exc_max,n_of_exc_max,n_deg_grp)
 character(*), intent(out) :: loc_character(n_of_exc_max,n_deg_grp)
 logical,      intent(in)  :: l_prt_debug
 integer,      intent(in)  :: iq
 !
 ! Workspace
 !
 logical     :: debug_file_to_open
 logical     :: l_mixed_states,l_partial_space,l_char_table_imag,l_char_not_found,case_simple,case_cmplx
 complex(SP) :: tmp_trace(12),char_tmp(12)
 !
 real(SP)             :: tmp_real
 character(9)         :: vec_ch(12)
 character(schlen)    :: filedebug,titles(3),tmp_ch,loc_msg,loc_fmt
 !
 integer  :: ich,jrap,irap,krap,iclass,ielem,iexc,il1,il2,i_grp,i_sub_grp,max_deg
 !
 integer  :: which_class(12),ncmplx
 !
 ! Compute the trace of the irreps
 !
 loc_trace=cZERO
 !
 do i_grp=1,n_deg_grp
   do i_sub_grp=1,n_sub_grps(i_grp)
     !write(*,*) i_grp," sub grp ",i_sub_grp," with ref ",ref_exc_sbgr(i_sub_grp,i_grp)
     do il1=1,n_of_exc(i_grp)
       if( exc_sbgr(il1,i_grp)/=i_sub_grp ) cycle
       !write(*,*) i_sub_grp,"exc in sbgr ",il1
       do iclass=1,nclass
         loc_trace(iclass,i_sub_grp,i_grp)=loc_trace(iclass,i_sub_grp,i_grp) + &
         &                      loc_reps(iclass,il1,il1,i_grp)
       enddo
     enddo
     !write(*,*) i_grp," sub grp ",i_sub_grp," with ref ",loc_trace(:nclass,i_sub_grp,i_grp)
   enddo
 enddo
 !
 ! DEBUG <
 !write(*,*) ""
 !do irap=1,nclass
 !  write(*,*) 0.0_SP,nelem(irap),real(char_mat(irap,which_irr(:nclass)),SP)
 !enddo
 ! DEBUG >
 !
 l_char_table_imag=maxval(abs(real(cI*char_mat,SP)))>0.01_SP
 !
 max_deg=nint(maxval(real(abs(char_mat(:,1)),SP)))
 ! DS TEST: removing this for excitons
 !if (l_char_table_imag) max_deg=max_deg*2
 !
 debug_file_to_open=.true.
 l_char_not_found=.false.
 !
 do iclass=1,nclass
   irap=which_irr(iclass)
   which_class(irap)=iclass
 enddo
 !
 do i_grp=1,n_deg_grp
   do i_sub_grp=1,n_sub_grps(i_grp)
   !
   !il=first_exc(i_grp)
   il1=first_exc(i_grp)+ref_exc_sbgr(i_sub_grp,i_grp)-1
   ! DEBUG <
   !write(*,*) BS_E_real(il)*HA2EV,i_grp,n_of_exc(i_grp),real(loc_trace(:nclass,i_grp))
   ! DEBUG >
   !
   l_mixed_states = nint(abs(loc_trace(1,i_sub_grp,i_grp))) > max_deg 
   !
   l_partial_space=.false.
   do iclass=1,nclass
     tmp_real=abs(loc_trace(iclass,i_sub_grp,i_grp))
     tmp_real=abs(tmp_real-real(nint(tmp_real),SP))
     l_partial_space = tmp_real > symm_thresh
     if(l_partial_space) exit
   enddo
   if (l_partial_space) then
     tmp_ch="partial space "//trim(real2ch(tmp_real))
     l_char_not_found=.true.
   endif
   !
   if ( .not. (l_mixed_states.or.l_partial_space) ) then
     do irap=1,nclass
       char_tmp(:nclass)=cmplx(char_mat(irap,which_irr(:nclass)),kind=SP)
       tmp_trace(:nclass)=loc_trace(:nclass,i_sub_grp,i_grp)
       ! Point groups without imaginary part
       case_simple=all(abs(     tmp_trace(:nclass)       -     char_tmp(:nclass)) <=symm_thresh)
       ! Point groups with imaginary part: 2 states with conjg character are summed together 
       case_cmplx =all(abs(real(tmp_trace(:nclass))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
       if ( .not. (case_simple .or. case_cmplx) ) cycle
       tmp_ch=trim(name_rap(irap))
       l_char_not_found=.false.
       exit
     enddo
     ! Try again with combination of 2 states
     if(irap==nclass+1) l_mixed_states=.true.
   endif
   !
   if (l_mixed_states) then
     main_loop2: do irap=1,nclass
       do jrap=1,nclass
         char_tmp(:nclass)=cmplx(char_mat(irap,which_irr(:nclass))+char_mat(jrap,which_irr(:nclass)),kind=SP)
         tmp_trace(:nclass)=loc_trace(:nclass,i_sub_grp,i_grp)
         case_simple=all(abs(     tmp_trace(:nclass)       -     char_tmp(:nclass) )<=symm_thresh)
         case_cmplx =all(abs(real(tmp_trace(:nclass))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
         if ( .not. (case_simple .or. case_cmplx) ) cycle
         tmp_ch=trim(name_rap(irap)(1:5))//" + "//trim(name_rap(jrap)(1:5))
         l_char_not_found=.false.
         exit main_loop2
       enddo
     enddo main_loop2
     ! Try again with combination of 3 states
     if (l_char_not_found) then
       main_loop3: do irap=1,nclass
         do jrap=1,nclass
         do krap=1,nclass
           char_tmp(:nclass)=cmplx(char_mat(irap,which_irr(:nclass))&
           &                      +char_mat(jrap,which_irr(:nclass))&
           &                      +char_mat(krap,which_irr(:nclass)),kind=SP)
           tmp_trace(:nclass)=loc_trace(:nclass,i_sub_grp,i_grp)
           case_simple=all(abs(     tmp_trace(:nclass)       -     char_tmp(:nclass) )<=symm_thresh)
           case_cmplx =all(abs(real(tmp_trace(:nclass))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
           if ( .not. (case_simple .or. case_cmplx) ) cycle
           tmp_ch=trim(name_rap(irap)(1:5))//" + "//trim(name_rap(jrap)(1:5))//" + "//trim(name_rap(krap)(1:5))
           l_char_not_found=.false.
           exit main_loop3
         enddo
         enddo
       enddo main_loop3
     endif
     if(irap==nclass+1) then
       tmp_ch="mixed states (not sum)"
       l_char_not_found=.true.
     endif
   endif
   !
   loc_character(i_sub_grp,i_grp)=trim(tmp_ch)
   !
   if ((l_char_not_found.or.l_mixed_states).and.l_prt_debug) then
     !
     if (debug_file_to_open) then
       filedebug="exc_qpt"//trim(intc(iq))//"_E+syms_sorted_dbg"
       call of_open_close(filedebug,'ot')
       loc_fmt="(i8,12x,F4.2,12x,a)"
       debug_file_to_open=.false.
       !
       call msg("o "//trim(filedebug),"# Symmetries debug at qpt ",iq,INDENT=0)
       call msg("o "//trim(filedebug),"#")
       !
       titles=(/"i_lambda"," n_deg  "," symm.  "/)
       call msg("o "//trim(filedebug),"#")
       call msg("o "//trim(filedebug),'#',titles,INDENT=0,USE_TABS=.true.)
       call msg("o "//trim(filedebug),"#")
       !
     endif
     !
     write(loc_msg,loc_fmt) il1+iexc_ref-1,real(loc_trace(1,i_sub_grp,i_grp)),trim(tmp_ch)
     call msg("o "//trim(filedebug),'',loc_msg,INDENT=-2)
     !
     ncmplx=1
     if ( any(aimag(loc_trace(:,i_sub_grp,i_grp))>1.E-4_SP) ) ncmplx=2
     do ielem=1,ncmplx
       !
       if (ielem==1) loc_msg="Real part analysis "
       if (ielem==2) loc_msg="Img. part analysis"
       call msg("o "//trim(filedebug),'',loc_msg,INDENT=-2)
       vec_ch=""
       do iclass=1,nclass
         vec_ch(iclass)=trim(name_class(iclass))
         do ich=len_trim(vec_ch(iclass)),5
           vec_ch(iclass)=" "//vec_ch(iclass)
         enddo
       enddo
       call msg("o "//trim(filedebug),"Class",vec_ch(:nclass),INDENT=0)
       !
       loc_msg=" Sum"
       do irap=1,nclass
         iclass=which_class(irap)
         if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(loc_trace(iclass,i_sub_grp,i_grp))
         if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(loc_trace(iclass,i_sub_grp,i_grp))
         loc_msg=trim(loc_msg)//trim(tmp_ch)
       enddo
       !
       call msg("o "//trim(filedebug),'',loc_msg,INDENT=-2)
       !
       do il2=1,n_of_exc(i_grp)
         if( exc_sbgr(il2,i_grp)/=i_sub_grp ) cycle
         iexc=first_exc(i_grp)+il2-1
         loc_msg="  "//trim(intc(iexc))
         do irap=1,nclass
           iclass=which_class(irap)
           if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(loc_reps(iclass,il2,il2,i_grp))
           if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(loc_reps(iclass,il2,il2,i_grp))
           loc_msg=trim(loc_msg)//trim(tmp_ch)
         enddo
         call msg("o "//trim(filedebug),'',loc_msg,INDENT=-2)
       enddo
       !
     enddo
     call msg("o "//trim(filedebug),'')
   endif
   !
   enddo
 enddo
 !
 if (.not.debug_file_to_open) call of_open_close(filedebug)
 !
end subroutine get_character_exc
