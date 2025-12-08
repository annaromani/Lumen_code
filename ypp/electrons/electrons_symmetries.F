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
subroutine electrons_symmetries(Xk,Xen,blim,klim)
 !
 use pars,           ONLY:DP,SP,cZERO,schlen,cI
 use units,          ONLY:HA2EV
 use electrons,      ONLY:levels,n_sp_pol,n_max_deg
 use R_lattice,      ONLY:bz_samp
 use stderr,         ONLY:intc,real2ch
 use D_lattice,      ONLY:elem,which_irr,char_mat,nclass,elem,nelem,&
&                         i_time_rev,nsym,nrot,dl_sop,dl_tra
 use vec_operate,    ONLY:degeneration_finder
 use wave_func,      ONLY:io_WF_phases,WF_phases_b_map
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:VERIFY,DUMP,REP,OP_RD,RD,RD_CL
 use com,            ONLY:msg,of_open_close
 use LIVE_t,         ONLY:live_timing
 use y_memory_alloc
 !
 implicit none
 !
 type(levels) :: Xen
 type(bz_samp):: Xk
 integer, intent(in) :: blim(2),klim(2)
 !
 complex(SP)       :: trace(12,blim(2)-blim(1)+1)
 character(schlen) :: character(blim(2)-blim(1)+1)
 !
 integer  :: iclass,irap,jrap,i_sp_pol,ikbz,ipbz,ik,isk,ibk,istark,ib1,iel,blim_n(2),ich,       &
 &           ielem,isp,k_count,s_count,n_deg_grp,i_grp,nb,nrotk,max_deg,nsz(5),ncmplx,list(48)
 !
 integer  :: io_err,ID_WF_phases
 !
 integer, allocatable :: first_el(:),n_of_el(:)
 !
 character(schlen)    :: filename,titles(4),tmp_ch,loc_msg,ref_msg,loc_fmt
 !
 nb=blim(2)-blim(1)+1
 !
 do ielem=1,48
   list(ielem)=ielem
 enddo
 !
 allocate(first_el(nb))
 allocate( n_of_el(nb))
 !
 call k_build_up_BZ_tables(Xk)
 call k_ibz2bz(Xk,'iku',.true.)
 call k_small_group(Xk)
 !
 call WF_phase_matrices(Xen,Xk,blim,klim)
 !
 nsz=0
 YAMBO_ALLOC(WF_phases_b_map,(blim(1):blim(2),klim(1):klim(2),n_sp_pol,2))
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1,2/),MODE=DUMP,ID=ID_WF_phases)
 io_err=io_WF_phases(blim,0,0,0,ID_WF_phases,nsz)
 !
 call live_timing("Electron symmetries:",n_sp_pol*(klim(2)-klim(1)+1))
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=klim(1),klim(2)
     !
     call degeneration_finder(nb,first_el,n_of_el,n_deg_grp,Er=Xen%E(blim(1):blim(2),ik,i_sp_pol),&
     &                                      deg_accuracy=0.0001/HA2EV,Include_single_values=.true.)
     !
     ! This loop can be used for debug
     do istark=1,1 !Xk%nstar(ik)
       !
       !if ( .not.PAR_IND_WF_k%element_1D(ik) ) cycle
       !
       isk =Xk%star(ik,istark)
       ikbz=Xk%k_table(ik,isk)
       !
       nrotk=Xk%grp_nsym(ikbz)
       if ( any( Xk%grp_table(:nrotk,ikbz)>nrot ) ) nrotk=nrotk/(1+i_time_rev)
       !
       filename="elec_k"//trim(intc(ik))//"_sp"//trim(intc(i_sp_pol))//"_syms"
       !
       ref_msg=""
       loc_msg=""
       do isp=1,nrotk
         write(tmp_ch,"(2x,i2)") isp
         ref_msg=trim(ref_msg)//trim(tmp_ch)
         write(tmp_ch,"(2x,i2)") Xk%grp_table(isp,ikbz)
         loc_msg=trim(loc_msg)//trim(tmp_ch)
       enddo
       !
       call of_open_close(filename,'ot')
       !call msg("o "//trim(filename),"# Symmetries at k ",Xk%ptbz(ikbz,:),INDENT=0)
       !call msg("o "//trim(filename),"#"," Symmetries at ikbz "//trim(intc(ikbz)),INDENT=0)
       call msg("o "//trim(filename),"# Symmetries at k ",Xk%pt(ik,:),INDENT=0)
       call msg("o "//trim(filename),"#")
       call msg("o "//trim(filename),"# Syms smallk grp ",ref_msg,INDENT=0)
       call msg("o "//trim(filename),"# Syms gamma  grp ",loc_msg,INDENT=0)
       call msg("o "//trim(filename),"#")
       !
       call detect_point_group(nrotk,real(dl_sop(:,:,Xk%grp_table(:nrotk,ikbz)),DP))
       !
       call print_point_group("o "//trim(filename),nrotk,nsym,dl_sop(:,:,Xk%grp_table(:nrotk,ikbz)),&
       &                                                  dl_tra(:,Xk%grp_table(:nrotk,ikbz)),.true.)
       !
       call get_character_el(nb,n_deg_grp,n_of_el,first_el, &
       &           blim,ID_WF_phases,Xk,                    &
       &           trace,character,.true.,ik,i_sp_pol)
       !
       loc_fmt="(F10.4,12x,i4,16x,F4.2,16x,a)"
       !
       titles=(/" E[eV]  ","   ib   "," n_deg  "," symm.  "/)
       call msg("o "//trim(filename),"#")
       call msg("o "//trim(filename),'#',titles,INDENT=0,USE_TABS=.true.)
       call msg("o "//trim(filename),"#")
       !
       do i_grp=1,n_deg_grp
         !
         ib1=first_el(i_grp)
         !
         write(loc_msg,loc_fmt) Xen%E(first_el(i_grp),ik,i_sp_pol)*HA2EV,ib1,&
         &                   real(trace(1,i_grp)),trim(character(i_grp))
         call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
         !
       enddo
       !
       call of_open_close(filename)
       !
     enddo
     !
     call live_timing(steps=1)
     !
   enddo 
 enddo
 !
 call live_timing()
 !
 nsz=0
 call io_control(ACTION=RD_CL,COM=REP,SEC=(/0/),MODE=DUMP,ID=ID_WF_phases)
 io_err=io_WF_phases(blim,0,0,0,ID_WF_phases,nsz)
 !
end subroutine electrons_symmetries
!
!
subroutine get_character_el(nb,n_deg_grp,n_of_el,first_el, &
                        &   blim,ID_WF_phases,Xk,          &
                        &   loc_trace,loc_character,l_prt_debug,ik,i_sp_pol)
 !
 use pars,           ONLY:SP,cI,cZERO,schlen
 use D_lattice,      ONLY:elem,which_irr,char_mat,nclass,elem,nelem,&
&                         nsym,sop_tab,which_class,name_rap,name_class
 use wave_func,      ONLY:io_WF_phases
 use com,            ONLY:msg,of_open_close
 use YPPm,           ONLY:symm_thresh
 use R_lattice,      ONLY:bz_samp
 use stderr,         ONLY:intc,real2ch
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:DUMP,REP,OP_RD,RD,RD_CL
 use y_memory_alloc
 !
 implicit none
 !
 integer,       intent(in)    :: nb,blim(2)
 integer,       intent(inout) :: ID_WF_phases
 type(bz_samp), intent(in)    :: Xk
 integer,       intent(in)    :: n_deg_grp
 integer,       intent(in)    :: first_el(nb),n_of_el(nb)
 complex(SP),   intent(out)   :: loc_trace(12,nb)
 !complex(SP),   intent(in)    :: loc_reps(12,nb,nb) ! WF_phases
 character(*),  intent(out)   :: loc_character(nb)
 logical,       intent(in)    :: l_prt_debug
 integer,       intent(in)    :: ik,i_sp_pol
 !
 ! Workspace
 !
 character(9)      :: vec_ch(12)
 complex(SP)       :: char_tmp(12)
 real(SP)          :: tmp_real
 !
 character(schlen) :: filedebug,titles(3),tmp_ch,loc_msg,ref_msg,loc_fmt
 !
 complex(SP), allocatable :: WF_phases(:,:,:,:,:)
 !
 logical  :: l_mixed_states,l_partial_space,l_char_table_imag,l_char_not_found,&
 &           case_simple,case_cmplx,debug_file_to_open
 integer  :: iclass,irap,jrap,ikbz,ipbz,isk,ibk,istark,ib1,iel,blim_n(2),ich,       &
 &           ielem,isp,k_count,s_count,i_grp,max_deg,nsz(5),ncmplx,list(48)
 !
 integer  :: io_err
 !
 istark=1
 isk =Xk%star(ik,istark)
 ikbz=Xk%k_table(ik,isk)
 !
 l_char_table_imag=maxval(abs(real(cI*char_mat)))>0.01_SP
 !
 max_deg=nint(maxval(real(abs(char_mat(:,1)))))
 if (l_char_table_imag) max_deg=max_deg*2
 !
 nsz=(/maxval(n_of_el),maxval(n_of_el),n_deg_grp,nsym,Xk%nstar(ik)/)
 YAMBO_ALLOC(WF_phases,(nsz(1),nsz(2),nsz(3),nsz(4),1))
 !
 call io_control(ACTION=RD,COM=REP,SEC=(/3/),MODE=DUMP,ID=ID_WF_phases)
 io_err=io_WF_phases(blim,ik,istark,i_sp_pol,ID_WF_phases,nsz,WF_phases)
 !
 do i_grp=1,n_deg_grp
   !
   ib1=first_el(i_grp)
   !
   loc_trace(:,i_grp)=cZERO
   do iclass=1,nclass
     ! For debugging, an averages over the elements can be performed
     do ielem=1,1 !nelem(iclass)
       isp=Xk%grp_table(elem(ielem,iclass),ikbz)
       !ipbz=Xk%k_table(ik,sop_tab(isp,isk))
       do ibk=1,n_of_el(i_grp)
         loc_trace(iclass,i_grp)=loc_trace(iclass,i_grp)+WF_phases(ibk,ibk,i_grp,isp,1)
       enddo
     enddo
     !loc_trace(iclass,i_grp)=loc_trace(iclass,i_grp)/real(nelem(iclass),SP)
   enddo
   !
   l_mixed_states = nint(abs(loc_trace(1,i_grp))) > max_deg 
   !
   l_partial_space=.false.
   do iclass=1,nclass
     tmp_real=abs(loc_trace(iclass,i_grp))
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
     ! Detection
     do irap=1,nclass
       char_tmp(:nclass)=char_mat(irap,which_irr(:nclass))
       ! Point groups without imaginary part
       case_simple=all(abs(     loc_trace(:nclass,i_grp)       -     char_tmp(:nclass)) <=symm_thresh)
       ! Point groups with imaginary part: 2 states with conjg character are summed together 
       case_cmplx =all(abs(real(loc_trace(:nclass,i_grp))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
       if ( .not. (case_simple .or. case_cmplx) ) cycle
       tmp_ch=trim(name_rap(irap))
       l_char_not_found=.false.
       exit
     enddo
     if(irap==nclass+1) l_char_not_found=.true.
   endif
   if (l_mixed_states.or.l_char_not_found) then
     main_loop: do irap=1,nclass
       do jrap=1,nclass
         char_tmp(:nclass)=char_mat(irap,which_irr(:nclass))+char_mat(jrap,which_irr(:nclass))
         case_simple=all(abs(     loc_trace(:nclass,i_grp)       -     char_tmp(:nclass) )<=symm_thresh)
         case_cmplx =all(abs(real(loc_trace(:nclass,i_grp))/2._SP-real(char_tmp(:nclass)))<=symm_thresh)
         if ( .not. (case_simple .or. case_cmplx) ) cycle
         tmp_ch=trim(name_rap(irap))//"+"//trim(name_rap(jrap))
         l_char_not_found=.false.
         exit main_loop
       enddo
     enddo main_loop
     if(irap==nclass+1) then
       if(l_mixed_states)   tmp_ch="mixed states (not sum)"
       if(l_char_not_found) tmp_ch="char not found (not sum)"
       l_char_not_found=.true.
     endif
   endif
   !
   loc_character(i_grp)=trim(tmp_ch)
   !
   if (l_char_not_found.and.l_prt_debug) then
     !
     if (debug_file_to_open) then
       filedebug="elec_k"//trim(intc(ik))//"_sp"//trim(intc(i_sp_pol))//"_syms_dbg"
       call of_open_close(filedebug,'ot')
       loc_fmt="(i4,12x,F4.2,12x,a)"
       debug_file_to_open=.false.
     endif
     !
     call msg("o "//trim(filedebug),"# Symmetries debug at kpt ",ik,INDENT=0)
     call msg("o "//trim(filedebug),"#")
     !
     titles=(/"   ib   "," n_deg  "," symm.  "/)
     call msg("o "//trim(filedebug),"#")
     call msg("o "//trim(filedebug),'#',titles,INDENT=0,USE_TABS=.true.)
     call msg("o "//trim(filedebug),"#")
     !
     ncmplx=1
     if ( any(aimag(loc_trace(:,i_grp))>1.E-4_SP) ) ncmplx=2
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
         if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(loc_trace(iclass,i_grp))
         if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(loc_trace(iclass,i_grp))
         loc_msg=trim(loc_msg)//trim(tmp_ch)
       enddo
       !
       call msg("o "//trim(filedebug),'',loc_msg,INDENT=-2)
       !
       do ibk=1,n_of_el(i_grp)
         iel=first_el(i_grp)+ibk-1
         loc_msg="  "//trim(intc(iel))
         do irap=1,nclass
           iclass=which_class(irap)
           isp=Xk%grp_table(elem(1,iclass),ikbz)
           if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(WF_phases(ibk,ibk,i_grp,isp,1))
           if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(WF_phases(ibk,ibk,i_grp,isp,1))
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
 !
 if (.not.debug_file_to_open) call of_open_close(filedebug)
 !
 YAMBO_FREE(WF_phases)
 !
end subroutine get_character_el
