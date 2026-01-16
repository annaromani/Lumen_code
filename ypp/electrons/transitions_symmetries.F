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
subroutine transitions_symmetries(iq,Xk,Xen,qpt,blim,klim)
 !
 use pars,           ONLY:DP,SP,cZERO,schlen,cI
 use units,          ONLY:HA2EV
 use electrons,      ONLY:levels,n_sp_pol,n_max_deg
 use R_lattice,      ONLY:bz_samp,qindx_X
 use stderr,         ONLY:intc,real2ch
 use D_lattice,      ONLY:elem,which_irr,char_mat,nclass,elem,nelem,&
&                         i_time_rev,nsym,nrot,dl_sop,dl_tra
 use vec_operate,    ONLY:sort,degeneration_finder
 use wave_func,      ONLY:io_WF_phases,WF_phases_b_map
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:VERIFY,DUMP,REP,OP_RD,RD,RD_CL
 use com,            ONLY:msg,of_open_close
 use LIVE_t,         ONLY:live_timing
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(in) :: iq
 type(levels) :: Xen
 type(bz_samp):: Xk,qpt
 integer, intent(in) :: blim(2),klim(2)
 !
 complex(SP), allocatable      :: trace(:,:)
 character(len=:),allocatable :: chartable(:)
 !
 integer  :: iclass,irap,jrap,i_sp_pol,ikbz,ik,isk,ibk,istark,ib1,iel,blim_n(2),ich,i_tr1,i_trk,&
 &           ielem,isp,k_count,s_count,n_deg_grp,i_grp,nb,nrotk,max_deg,nsz(5),ncmplx,list(48), &
 &           iv,ic,ikmq,ikmqbz,ikmqg0,istarkmq,iqbz,n_trans_k,nsymq,nrotq
 !
 integer  :: io_err,ID_WF_phases
 !
 integer,  allocatable :: first_tr(:),n_of_tr(:),table_tr(:,:),S_index(:)
 real(SP), allocatable :: E_trans(:)
 !
 character(schlen)    :: filename,titles(4),tmp_ch,loc_msg,ref_msg,loc_fmt
 !
 nb=blim(2)-blim(1)+1
 !
 do ielem=1,48
   list(ielem)=ielem
 enddo
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
 ! Transitions part
 !
 call k_build_up_BZ_tables(qpt)
 call k_ibz2bz(qpt,'iku',.true.)
 call k_small_group(qpt)
 !
 ! The q point for which I want to compute the irrep.
 !isq =qpr%star(iq,1)
 iqbz=qpt%k_table(iq,1)
 !iGq =qpt%g_table(iq,isq)
 !
 nsymq=qpt%grp_nsym(iqbz)
 nrotq=nsymq
 if ( any( qpt%grp_table(:nrotq,iqbz)>nrot ) ) nrotq=nrotq/(1+i_time_rev)
 !
 ! Detect the point group
 !
 call detect_point_group(nrotq,real(dl_sop(:,:,qpt%grp_table(:nrotq,iqbz)),DP))
 !
 call live_timing("Transitions symmetries:",n_sp_pol*(klim(2)-klim(1)+1))
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ik=klim(1),klim(2)
     !
     do istark=1,1 !Xk%nstar(ik)
       !
       isk =Xk%star(ik,istark)
       ikbz=Xk%k_table(ik,isk)
       !
       ! (R)^-1 k_BZ - q_IBZ = p_BZ - q_IBZ
       ikmqbz=qindx_X(iq,ikbz,1)
       ikmqg0=qindx_X(iq,ikbz,2)
       !
       ikmq    =Xk%sstar(ikmqbz,1)
       istarkmq=Xk%sstar(ikmqbz,3)
       !
       n_trans_k=0
       do iv=blim(1),Xen%nbm(i_sp_pol)
         do ic=Xen%nbf(i_sp_pol)+1,blim(2)
           n_trans_k=n_trans_k+1
         enddo
       enddo
       !
       allocate(first_tr(n_trans_k))
       allocate( n_of_tr(n_trans_k))
       allocate(table_tr(n_trans_k,2))
       allocate( E_trans(n_trans_k))
       allocate( S_index(n_trans_k))
       !
       n_trans_k=0
       do iv=blim(1),Xen%nbm(i_sp_pol)
         do ic=Xen%nbf(i_sp_pol)+1,blim(2)
           n_trans_k=n_trans_k+1
           table_tr(n_trans_k,:)=(/iv,ic/)
           E_trans(n_trans_k)=Xen%E(ic,ik,i_sp_pol)-Xen%E(iv,ikmq,i_sp_pol)
         enddo
       enddo
       !
       call sort(E_trans,indx=S_index)
       !
       call degeneration_finder(n_trans_k,first_tr,n_of_tr,n_deg_grp,Er=E_trans,&
       &                   deg_accuracy=0.0001/HA2EV,Include_single_values=.true.)
       !
       filename="trans_k"//trim(intc(ik))//"_sp"//trim(intc(i_sp_pol))//"_syms"
       !
       call of_open_close(filename,'ot')
       call msg("o "//trim(filename),"# Symmetries at k ",Xk%pt(ik,:),INDENT=0)
       call msg("o "//trim(filename),"#")
       !
       !call print_point_group("o "//trim(filename),nrotq,nsymq,dl_sop(:,:,qpt%grp_table(:nsymq,iqbz)),&
       !&                                                  dl_tra(:,qpt%grp_table(:nrotq,iqbz)),.true.)
       !
       allocate(trace(12,n_deg_grp))
       allocate(character(len=schlen) :: chartable(n_deg_grp))

       !
       call get_character_trans(n_trans_k,n_deg_grp,n_of_tr,first_tr, &
       &           blim,ID_WF_phases,Xk,qpt,table_tr,S_index,         &
       &           trace,chartable,.true.,ik,ikmq,iq,i_sp_pol)
       !
       loc_fmt="(F10.4,12x,a,12x,F4.2,12x,a)"
       !
       titles=(/" E[eV]  ","iv -> ic"," n_deg  "," symm.  "/)
       call msg("o "//trim(filename),"#")
       call msg("o "//trim(filename),'#',titles,INDENT=0,USE_TABS=.true.)
       call msg("o "//trim(filename),"#")
       !
       do i_grp=1,n_deg_grp
         !
         i_tr1=first_tr(i_grp)
	 i_trk=S_index(i_tr1)
         !
         iv=table_tr(i_trk,1)
         ic=table_tr(i_trk,2)
	 tmp_ch=trim(intc(iv))//" --> "//trim(intc(ic))
         write(loc_msg,loc_fmt) E_trans(i_tr1)*HA2EV,trim(tmp_ch),&
         &                   real(trace(1,i_grp)),trim(chartable(i_grp))
         call msg("o "//trim(filename),'',loc_msg,INDENT=-2)
         !
       enddo
       !
       call of_open_close(filename)
       !
       deallocate(trace)
       deallocate(chartable)
       !
     enddo
     !
     deallocate(first_tr)
     deallocate( n_of_tr)
     deallocate(table_tr)
     deallocate( E_trans)
     deallocate( S_index)
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
 YAMBO_FREE(WF_phases_b_map)
 !
end subroutine transitions_symmetries
!
!
subroutine get_character_trans(n_trans_k,n_deg_grp,n_of_tr,first_tr,  &
                        &   blim,ID_WF_phases,Xk,qpt,table_tr,S_index,&
                        &   loc_trace,loc_character,l_prt_debug,ik,ikmq,iq,i_sp_pol)
 !
 use pars,           ONLY:SP,cI,cZERO,schlen
 use D_lattice,      ONLY:elem,which_irr,char_mat,nclass,elem,nelem,&
&                         nsym,sop_tab,which_class,name_rap,name_class
 use wave_func,      ONLY:io_WF_phases,WF_phases_b_map
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
 integer,       intent(in)    :: n_trans_k,blim(2)
 integer,       intent(inout) :: ID_WF_phases,S_index(n_trans_k)
 type(bz_samp), intent(in)    :: Xk,qpt
 integer,       intent(in)    :: n_deg_grp,table_tr(n_trans_k,2)
 integer,       intent(in)    :: first_tr(n_trans_k),n_of_tr(n_trans_k)
 complex(SP),   intent(out)   :: loc_trace(12,n_deg_grp)
 !complex(SP),   intent(in)    :: loc_reps(12,nb,nb) ! WF_phases
 character(*),  intent(out)   :: loc_character(n_deg_grp)
 logical,       intent(in)    :: l_prt_debug
 integer,       intent(in)    :: ik,ikmq,iq,i_sp_pol
 !
 ! Workspace
 !
 character(9)      :: vec_ch(12)
 complex(SP)       :: char_tmp(12),phase_vv,phase_cc
 real(SP)          :: tmp_real
 !
 character(schlen) :: filedebug,titles(3),tmp_ch,loc_msg,ref_msg,loc_fmt
 !
 complex(SP), allocatable :: WF_phases_vv(:,:,:,:,:)
 complex(SP), allocatable :: WF_phases_cc(:,:,:,:,:)
 !
 logical  :: l_mixed_states,l_partial_space,l_char_table_imag,l_char_not_found,&
 &           case_simple,case_cmplx,debug_file_to_open
 integer  :: iclass,irap,jrap,ikbz,isk,ibk,istark,ib1,iel,blim_n(2),ich,       &
 &           ielem,isp,k_count,s_count,i_grp,max_deg,nsz(5),ncmplx,list(48)
 !
 integer  :: i_trk,i_trk_all,i_tr1,iv,ic,ic_bnd,ic_grp,iv_bnd,iv_grp,isq,iqbz
 !
 integer  :: io_err
 !
 istark=1
 isk =Xk%star(ik,istark)
 ikbz=Xk%k_table(ik,isk)
 !
 debug_file_to_open=.true.
 !
 iqbz=qpt%k_table(iq,1)
 !
 l_char_table_imag=maxval(abs(real(cI*char_mat)))>0.01_SP
 !
 max_deg=nint(maxval(real(abs(char_mat(:,1)))))
 if (l_char_table_imag) max_deg=max_deg*2
 !
 nsz(1:2)=maxval(WF_phases_b_map(:,ik,i_sp_pol,1))
 nsz(3)  =maxval(WF_phases_b_map(:,ik,i_sp_pol,2))
 nsz(4)  =nsym
 nsz(5)  =Xk%nstar(ik)
 YAMBO_ALLOC(WF_phases_cc,(nsz(1),nsz(2),nsz(3),nsz(4),nsz(5)))
 !
 call io_control(ACTION=RD,COM=REP,SEC=(/3/),MODE=DUMP,ID=ID_WF_phases)
 io_err=io_WF_phases(blim,ik,0,i_sp_pol,ID_WF_phases,nsz,WF_phases_cc)
 !
 !
 nsz(1:2)=maxval(WF_phases_b_map(:,ikmq,i_sp_pol,1))
 nsz(3)  =maxval(WF_phases_b_map(:,ikmq,i_sp_pol,2))
 nsz(5)  =Xk%nstar(ikmq)
 YAMBO_ALLOC(WF_phases_vv,(nsz(1),nsz(2),nsz(3),nsz(4),nsz(5)))
 !
 call io_control(ACTION=RD,COM=REP,SEC=(/3/),MODE=DUMP,ID=ID_WF_phases)
 io_err=io_WF_phases(blim,ikmq,0,i_sp_pol,ID_WF_phases,nsz,WF_phases_vv)
 !
 do i_grp=1,n_deg_grp
   !
   i_tr1=first_tr(i_grp)
   !
   loc_trace(:,i_grp)=cZERO
   do iclass=1,nclass
     ! For debugging, an averages over the elements can be performed
     do ielem=1,1 !nelem(iclass)
      isq=qpt%grp_table(elem(ielem,iclass),iqbz)
      ! DS: to devide what to do here ...
      !do isq=1,qpt%nstar(iq)
       !ipbz=Xk%k_table(ik,sop_tab(isp,isk))
       do i_trk=1,n_of_tr(i_grp)
         iv=table_tr(S_index(i_trk+i_tr1-1),1)
         ic=table_tr(S_index(i_trk+i_tr1-1),2)
         ic_bnd=WF_phases_b_map(ic,ik,i_sp_pol,1)
         ic_grp=WF_phases_b_map(ic,ik,i_sp_pol,2)
         iv_bnd=WF_phases_b_map(iv,ikmq,i_sp_pol,1)
         iv_grp=WF_phases_b_map(iv,ikmq,i_sp_pol,2)
         phase_vv=conjg(WF_phases_vv(iv_bnd,iv_bnd,iv_grp,isq,1))
         phase_cc=      WF_phases_cc(ic_bnd,ic_bnd,ic_grp,isq,1)
         loc_trace(iclass,i_grp)=loc_trace(iclass,i_grp)+phase_vv*phase_cc
       enddo
      !enddo
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
       !
       filedebug="trans_k"//trim(intc(ik))//"_sp"//trim(intc(i_sp_pol))//"_syms_dbg"
       call of_open_close(filedebug,'ot')
       loc_fmt="(i4,12x,F4.2,12x,a)"
       debug_file_to_open=.false.
       call msg("o "//trim(filedebug),"# Symmetries debug at kpt ",ik,INDENT=0)
       call msg("o "//trim(filedebug),"#")
       !
       titles=(/"  itr   "," n_deg  "," symm.  "/)
       call msg("o "//trim(filedebug),"#")
       call msg("o "//trim(filedebug),'#',titles,INDENT=0,USE_TABS=.true.)
       call msg("o "//trim(filedebug),"#")
       !
     endif
     !
     write(loc_msg,loc_fmt) i_tr1,real(loc_trace(1,i_grp)),trim(tmp_ch)
     call msg("o "//trim(filedebug),'',loc_msg,INDENT=-2)
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
       do i_trk=1,n_of_tr(i_grp)
         i_trk_all=first_tr(i_grp)+i_trk-1
         loc_msg="  "//trim(intc(iel))
         do irap=1,nclass
           iclass=which_class(irap)
           isp=Xk%grp_table(elem(1,iclass),ikbz)
           ! DS debug to be fixed
           !if(ielem==1) write(tmp_ch,"(2x,F5.2,2x)")  real(WF_phases(ibk,ibk,i_grp,isp,1))
           !if(ielem==2) write(tmp_ch,"(2x,F5.2,2x)") aimag(WF_phases(ibk,ibk,i_grp,isp,1))
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
 YAMBO_FREE(WF_phases_vv)
 YAMBO_FREE(WF_phases_cc)
 !
end subroutine get_character_trans
