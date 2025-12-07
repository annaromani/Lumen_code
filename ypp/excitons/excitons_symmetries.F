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
subroutine excitons_symmetries(iq,Xen,Xk,qpt,BS_E)
 !
 use pars,           ONLY:SP,DP,cZERO,cI,schlen
 use units,          ONLY:HA2EV
 use parser_m,       ONLY:parser
 use YPPm,           ONLY:deg_energy,EXCITONS_n_user_states,&
 &                        EXCITONS_user_indexes,symm_thresh
 use R_lattice,      ONLY:bz_samp,coll_scattering
 use electrons,      ONLY:levels,n_sp_pol
 use BS_solvers,     ONLY:BS_mat,BSS_n_eig
 use BS,             ONLY:BS_bands,BS_H_dim,BS_K_dim
 use wave_func,      ONLY:io_WF_phases,WF_phases_b_map,WF_phases
 use stderr,         ONLY:intc,real2ch
 use D_lattice,      ONLY:which_irr,char_mat,nclass,elem,nelem,name_class,&
 &                        name_rap,i_time_rev,nsym,nrot,dl_sop,dl_tra
 use vec_operate,    ONLY:sort,degeneration_finder
 use wrapper,        ONLY:Vstar_dot_V
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:DUMP,REP,OP_RD_CL,OP_RD,RD_CL,RD
 use com,            ONLY:msg,of_open_close
 use LIVE_t,         ONLY:live_timing
 !
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(in) :: iq
 type(bz_samp)       :: Xk,qpt
 type(levels)        :: Xen
 complex(SP)         :: BS_E(BSS_n_eig)
 !
 real(SP) :: BS_E_real(EXCITONS_n_user_states),tmp_real,ref_val
 logical  :: l_min_mem,l_mixed_states,l_partial_space,l_char_table_imag,l_char_not_found,case_simple,case_cmplx
 integer  :: S_index(EXCITONS_n_user_states)
 integer  :: first_exc(EXCITONS_n_user_states),n_of_exc(EXCITONS_n_user_states),n_deg_grp,n_of_exc_max
 integer  :: ik,i_sp_pol,iqbz,iRqbz,isqc,isqy,isq_small_grp,iGqc,iGqy,nrotq,ncmplx,nsz(5),&
 &           ich,jrap,irap,iclass,ielem,iexc,jexc,il1,il2,i_grp,i_sub_grp,io_err,ID,n_elem_used,max_deg
 !
 complex(SP), allocatable :: EXC_WF_rot1(:,:),EXC_WF_rot2(:,:)
 complex(SP) :: trace(12,EXCITONS_n_user_states,EXCITONS_n_user_states),trace_loc(12),char_tmp(12)
 complex(SP) :: exc_reps(12,EXCITONS_n_user_states,EXCITONS_n_user_states,EXCITONS_n_user_states)
 !
 integer  :: exc_sbgr(EXCITONS_n_user_states,EXCITONS_n_user_states),n_sub_grps(EXCITONS_n_user_states),&
 &           ref_exc_grp(EXCITONS_n_user_states,EXCITONS_n_user_states),&
 &           ref_exc_sbgr(EXCITONS_n_user_states,EXCITONS_n_user_states),which_class(12)
 !
 character(9)         :: vec_ch(12)
 character(schlen)    :: filename,titles(4),tmp_ch,loc_msg,loc_fmt
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
 YAMBO_ALLOC(WF_phases_b_map,(BS_bands(1):BS_bands(2),Xk%nibz,n_sp_pol,2))
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1,2/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,0,ID,(/0,0,0,0,0/))
 !
 call parser("MinimizeMem",l_min_mem)
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
 YAMBO_ALLOC(EXC_WF_rot1,(BS_H_dim,n_of_exc_max))
 YAMBO_ALLOC(EXC_WF_rot2,(BS_H_dim,n_of_exc_max))
 !
 filename="exc_qpt"//trim(intc(iq))//"_E+syms_sorted"
 !
 call of_open_close(filename,'ot')
 call msg("o "//trim(filename),"# Symmetries at qpt ",qpt%pt(iq,:),INDENT=0)
 call msg("o "//trim(filename),"#")
 !
 nrotq=qpt%grp_nsym(iqbz)
 if ( any( qpt%grp_table(:nrotq,iqbz)>nrot ) ) nrotq=nrotq/(1+i_time_rev)
 !
 !write(*,*) qpt%grp_table(:nrotq,iqbz)
 !
 call detect_point_group(nrotq,real(dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),DP))
 !
 call print_point_group("o "//trim(filename),nrotq,nsym,dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),&
 &                                                  dl_tra(:,qpt%grp_table(:nrotq,iqbz)),.true.)
 !
 !write(*,*) "point group detected ",nclass,n_deg_grp
 !
 titles=(/" E[eV]  "," n_deg  "," symm.  ","i_lambda"/)
 call msg("o "//trim(filename),"#")
 call msg("o "//trim(filename),'#',titles,INDENT=0,USE_TABS=.true.)
 call msg("o "//trim(filename),"#")
 !
 ! Compute excitonic reps
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
       & iq,isqy,iGqy,ID,S_index,EXC_WF_rot1(:,:n_of_exc(i_grp)),l_min_mem)
       !
       ! Get R * A(l,q)
       !
       call EXC_WF_symm_qpoint(Xk,qpt,first_exc(i_grp),n_of_exc(i_grp),&
       & iq,isqc,iGqc,ID,S_index,EXC_WF_rot2(:,:n_of_exc(i_grp)),l_min_mem)
       !
       ! Compute the trace of the Irrep
       !
       do il1=1,n_of_exc(i_grp)
         do il2=1,n_of_exc(i_grp)
           exc_reps(iclass,il2,il1,i_grp) = &
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
 call io_control(ACTION=RD_CL,COM=REP,SEC=(/0/),MODE=DUMP,ID=ID)
 io_err=io_WF_phases(BS_bands,0,0,0,ID,(/0,0,0,0,0/))
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
       call msg('s',"Rep "//trim(intc(i_grp))//" ,Exc. mix "//trim(intc(iexc))//","&
       &                  //trim(intc(jexc)),sum(abs(exc_reps(:,il2,il1,i_grp))))
       if ( sum(abs(exc_reps(:,il2,il1,i_grp)))>=nclass*symm_thresh) cycle
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
 ! Compute the trace of the irreps
 !
 trace=cZERO
 !
 do i_grp=1,n_deg_grp
   do i_sub_grp=1,n_sub_grps(i_grp)
     !write(*,*) i_grp," sub grp ",i_sub_grp," with ref ",ref_exc_sbgr(i_sub_grp,i_grp)
     do il1=1,n_of_exc(i_grp)
       if( exc_sbgr(il1,i_grp)/=i_sub_grp ) cycle
       !write(*,*) i_sub_grp,"exc in sbgrp ",il1
       do iclass=1,nclass
         trace(iclass,i_sub_grp,i_grp)=trace(iclass,i_sub_grp,i_grp) + &
         &                      exc_reps(iclass,il1,il1,i_grp)
       enddo
     enddo
   enddo
 enddo
 !
 ! Finally I compute the exciton character
 !
 loc_fmt="(F10.4,12x,F4.2,12x,a,12x,i4)"
 l_char_table_imag=maxval(abs(real(cI*char_mat)))>0.01_SP
 !
 max_deg=nint(maxval(real(abs(char_mat(:,1)))))
 if (l_char_table_imag) max_deg=max_deg*2
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
   !write(*,*) BS_E_real(il)*HA2EV,i_grp,n_of_exc(i_grp),real(trace(:nclass,i_grp))
   ! DEBUG >
   !
   l_mixed_states = nint(abs(trace(1,i_sub_grp,i_grp))) > max_deg 
   !
   l_partial_space=.false.
   do iclass=1,nclass
     tmp_real=abs(trace(iclass,i_sub_grp,i_grp))
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
       char_tmp(:nclass)=char_mat(irap,which_irr(:nclass))
       trace_loc(:nclass)=trace(:nclass,i_sub_grp,i_grp)
       ! Point groups without imaginary part
       case_simple=all(abs(     trace_loc(:nclass)       -     char_tmp(:nclass)) <=symm_thresh)
       ! Point groups with imaginary part: 2 states with conjg character are summed together 
       case_cmplx =all(abs(real(trace_loc(:nclass))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
       if ( .not. (case_simple .or. case_cmplx) ) cycle
       tmp_ch=trim(name_rap(irap))
       l_char_not_found=.false.
       exit
     enddo
     if(irap==nclass+1) then
       tmp_ch="char not found"
       l_char_not_found=.true.
     endif
   else if (l_mixed_states) then
     main_loop: do irap=1,nclass
       do jrap=1,nclass
         char_tmp(:nclass)=char_mat(irap,which_irr(:nclass))+char_mat(jrap,which_irr(:nclass))
         trace_loc(:nclass)=trace(:nclass,i_sub_grp,i_grp)
         case_simple=all(abs(     trace_loc(:nclass)       -     char_tmp(:nclass) )<=symm_thresh)
         case_cmplx =all(abs(real(trace_loc(:nclass))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
         if ( .not. (case_simple .or. case_cmplx) ) cycle
         tmp_ch=trim(name_rap(irap))//"+"//trim(name_rap(jrap))
         l_char_not_found=.false.
         exit main_loop
       enddo
     enddo main_loop
     if(irap==nclass+1) then
       tmp_ch="mixed states (not sum)"
       l_char_not_found=.true.
     endif
   endif
   !
   write(loc_msg,loc_fmt) BS_E_real(il1)*HA2EV,real(trace(1,i_sub_grp,i_grp)),trim(tmp_ch),il1
   call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
   if (l_char_not_found) then
     ncmplx=1
     if ( any(aimag(trace(:,i_sub_grp,i_grp))>1.E-4_SP) ) ncmplx=2
     do ielem=1,ncmplx
       !
       if (ielem==1) loc_msg="Real part analysis "
       if (ielem==2) loc_msg="Img. part analysis"
       call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
       vec_ch=""
       do iclass=1,nclass
         vec_ch(iclass)=trim(name_class(iclass))
         do ich=len_trim(vec_ch(iclass)),5
           vec_ch(iclass)=" "//vec_ch(iclass)
         enddo
       enddo
       call msg("o "//trim(filename),"Class",vec_ch(:nclass),INDENT=0)
       !
       loc_msg=" Sum"
       do irap=1,nclass
         iclass=which_class(irap)
         if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(trace(iclass,i_sub_grp,i_grp))
         if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(trace(iclass,i_sub_grp,i_grp))
         loc_msg=trim(loc_msg)//trim(tmp_ch)
       enddo
       !
       call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
       !
       do il2=1,n_of_exc(i_grp)
         if( exc_sbgr(il2,i_grp)/=i_sub_grp ) cycle
         iexc=first_exc(i_grp)+il2-1
         loc_msg="  "//trim(intc(iexc))
         do irap=1,nclass
           iclass=which_class(irap)
           if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(exc_reps(iclass,il2,il2,i_grp))
           if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(exc_reps(iclass,il2,il2,i_grp))
           loc_msg=trim(loc_msg)//trim(tmp_ch)
         enddo
         call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
       enddo
       !
     enddo
     call msg("o "//trim(filename),'')
   endif
   !
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
 call of_open_close(filename)
 !
end subroutine excitons_symmetries
