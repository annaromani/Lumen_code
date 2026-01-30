!
! License-Identifier: GPL
!
! Copyright (C) 2022 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
! headers
!
#include<y_memory.h>
!
subroutine RT_transient_ARPES(Xen,Xk,Xq)
 !
 use parser_m,        ONLY:parser
 use stderr,          ONLY:intc
 use vec_operate,     ONLY:sort
 use units,           ONLY:HA2EV
 use pars,            ONLY:SP,cONE,cI,cZERO,schlen
 use R_lattice,       ONLY:bz_samp,bz_samp_duplicate,bz_samp_duplicate_nosymm,bz_samp_reset,qindx_X
 use D_lattice,       ONLY:lattice,nsym
 use BS,              ONLY:BS_K_dim,BS_bands
 use BS_solvers,      ONLY:BSS_eh_table,BS_mat,BSS_n_eig
 use YPPm,            ONLY:BS_E,DIPs,exc_kind,EXCITONS_user_indexes   !,EXCITONS_n_user_states
 use YPPm,            ONLY:BANDS_path,BANDS_steps,CIRCUIT_k,BSiq  !,USER_k
 use YPP_interfaces,  ONLY:k_special
 use YPP_real_time,   ONLY:RT_conf,GreenF_T_and_W,TRabs_Exc
 use parallel_m,      ONLY:PAR_IND_DIPk_ibz,PAR_IND_VAL_BANDS_DIP,PAR_IND_CON_BANDS_DIP,myid
 use parallel_int,    ONLY:PP_redux_wait,PARALLEL_global_indexes
 use electrons,       ONLY:levels,spin,n_sp_pol
 use com,             ONLY:msg,of_open_close
 use functions,       ONLY:bose_f,boltzman_f !,Lorentzian_func
 use interpolate,     ONLY:INTERPOLATE_use_symm,INTERP_obj, &
 &                         INTERPOLATION_driver_seed,INTERPOLATION_driver_do
 use LIVE_t,          ONLY:live_timing
 !
 use y_memory_alloc
 !
 implicit none
 !
 type(levels),  intent(in)    :: Xen
 type(bz_samp), intent(inout) :: Xk,Xq
 !
 type(bz_samp)         :: k_nosymm,CIRCUIT_kmq
 logical               :: l_full_expansion,l_single_exciton,l_exciton_distribution,warn
 integer               :: iq,iq_bz,i_l,i_l_ref,i_c,i_v,i_cvk,i_kibz,i_kbz,i_kmq_ibz,i_kmq_bz,i_karpes,i_ks,i_sp_pol,i_w,i_l_v, &
 &                        io_err,ID_ARPES,ID_ARPES_int,ID_bands,NV,NK,q_range(2),l_range(2),                               &
 &                        EXC_table(1000,2),N_exc,j_exc,EXC_indx(1000)
 real(SP)              :: dw,Ev,Z,min_exc_E,N_lq,EXC_N_lq(1000),EXC_E_lq(1000)
 complex(SP)           :: F_tmp(RT_conf%W_n_steps),PES(RT_conf%W_n_steps)
 real(SP),allocatable  :: ARPES_int(:,:,:)
 integer, allocatable  :: BS_E_degs(:)
 real(SP),allocatable  :: Res(:)
 !
 character(schlen)     :: dist_head(4),headings(3),file_name,title
 real(SP)              :: dist_data(4),odata(3),DIP_dir(3,3)
 !
 complex(SP),external      ::Lorentzian_FT
 !
 warn=.false.
 exc_kind="exciton"
 !
 call excitons_read(Xk,Xk,Xen,Xen,1,"check",warn,io_err)
 if (io_err==0) then
   call excitons_read(Xk,Xk,Xen,Xen,1,"eigenvalues",warn,io_err)
   call excitons_read(Xk,Xk,Xen,Xen,1,"eigenvectors",warn,io_err)
 else
   call error("BSE information not found or not correct")
 endif
 !
 !
 ! Start with bands interpolation
 !
 call section('=',"Loading special Points for the "//trim(lattice)//" lattice")
 !============================================================================= 
 !
 call k_special()
 !
 call k_circuit_driver( BANDS_path, BANDS_steps )
 !
 ! Now move to excitonic ARPES
 !==============================
 !
 call parser("FullExpansion", l_full_expansion)
 call parser("EXCTemp", l_exciton_distribution)
 l_single_exciton=.not.l_exciton_distribution
 !
 if (.not. (l_single_exciton .or. l_exciton_distribution ) ) call error(" Select either Exc states or Exc temp in input")
 !
 ! Par Setup
 !============
 DIPs%ib = BS_bands
 call PARALLEL_global_indexes(Xen,Xk,Xq," ",RESET=.TRUE.)
 call PARALLEL_global_indexes(Xen,Xk,Xq,"DIPOLES",Dip=DIPs)
 !
 if (l_single_exciton) then
   ! ... I/O
   call DIPOLE_IO(Xk,Xen,DIPs,'read',io_err,'DIP')
   if (io_err/=0) return
 endif
 !
 ! Frequency setup
 !============
 !call RT_1D_Fourier_setup(GreenF_T_and_W,RT_conf,1)
 !AMBO_FREE(GreenF_T_and_W%F_w)
 !AMBO_FREE(GreenF_T_and_W%F_t)
 !
 YAMBO_ALLOC(GreenF_T_and_W%W,(RT_conf%W_n_steps))
 GreenF_T_and_W%W=cZERO
 !
 dw=(RT_conf%X_W_range(2) - RT_conf%X_W_range(1))/(RT_conf%W_n_steps-1)
 do i_w = 1, RT_conf%W_n_steps
   GreenF_T_and_W%W(i_w)=cmplx((i_w-1)*dw+RT_conf%X_W_range(1),RT_conf%damp_factor)
 enddo
 !
 if (l_single_exciton) then
   !
   ! Set the exciton
   !=================
   allocate(BS_E_degs(BSS_n_eig))
   allocate(Res(BSS_n_eig))
   call excitons_get_user_states(BS_E_degs,Res,BSS_n_eig,.true.)
   deallocate(BS_E_degs)
   deallocate(Res)
   !
   call RT_TRabs_rotate_exciton("DIP_iR",DIP_dir,cONE,Xen,Xk)
   i_l_ref=EXCITONS_user_indexes(1)
   deallocate(EXCITONS_user_indexes)
   !
   q_range = (/BSiq,BSiq/)
   l_range= (/i_l_ref,i_l_ref/)
   !
 else if (l_exciton_distribution) then
   !
   min_exc_E=minval(real(BS_E)) ! I suppose the BSE eigenvalues are real
   call excitons_read(Xk,Xk,Xen,Xen,1,"clean",warn,io_err)
   !
   do iq=2,Xq%nibz
     call excitons_read(Xk,Xk,Xen,Xen,iq,"check",warn,io_err)
     if (io_err/=0) then
       call excitons_read(Xk,Xk,Xen,Xen,iq,"clean",warn,io_err)
       cycle
     endif
     call excitons_read(Xk,Xk,Xen,Xen,iq,"eigenvalues",warn,io_err)   
     min_exc_E=min(min_exc_E,minval(real(BS_E))) ! I suppose the BSE eigenvalues are real
     call excitons_read(Xk,Xk,Xen,Xen,iq,"clean",warn,io_err)
   enddo
   !
   if(TRabs_Exc(1)<=0) TRabs_Exc(1)=1
   if(TRabs_Exc(2)==0) TRabs_Exc(2)=BSS_n_eig
   if(TRabs_Exc(2)> BSS_n_eig) then
     TRabs_Exc(2)=BSS_n_eig
     call warning(" Number of excitons in the sum reduced to"//trim(intc(BSS_n_eig)))
   endif
   !
   q_range = (/1,Xq%nibz/)
   l_range= (/TRabs_Exc(1),TRabs_Exc(2)/)
   !
 endif
 !
 call section('=','ARPES Interpolation ')
 !============================================================================= 
 !
 NV=Xen%nbf(1)-BS_bands(1)+1
 NK=Xk%nibz
 if (l_full_expansion) NK=Xk%nbz
 !
 ! Last dimension is for number of excitons
 ! Since I use it inside the loop over excitons, it is always 1
 YAMBO_ALLOC(ARPES_int,(BS_bands(1):Xen%nbf(1),NK,1))
 !
 ID_ARPES=3
 YAMBO_ALLOC(INTERP_obj(ID_ARPES)%INTERP,(RT_conf%W_n_steps,CIRCUIT_K%nibz,2))
 INTERP_obj(ID_ARPES)%INTERP=0._SP
 ! 
 PES = cZERO
 Z   = 0._SP
 N_exc = 0
 !
 if(l_exciton_distribution) call live_timing("ARPES Interpolation",(q_range(2)-q_range(1)+1)*(l_range(2)-l_range(1)+1))
 !
 call INTERPOLATION_driver_seed("E_K_PATH",ID_bands,BS_bands,Xk,SEED_E=Xen,NEAREST=.false.,QUIET=l_exciton_distribution)
 !
 !do iq_bz=q_range(1),q_range(2)
 do iq=q_range(1),q_range(2)
   !
   ! Bands interpolation
   !
   call bz_samp_reset(CIRCUIT_kmq)
   call bz_samp_duplicate(CIRCUIT_k,CIRCUIT_kmq)
   !
   do i_kibz=1,CIRCUIT_K%nibz
     CIRCUIT_kmq%pt(i_kibz,:)=CIRCUIT_k%pt(i_kibz,:)-Xq%pt(iq,:)
   enddo
   !
   call INTERPOLATION_driver_do("E_K_PATH",ID_bands,CIRCUIT_kmq)
   !
   if (iq==q_range(1)) then
     call plot_interpolated_values &
     &   (Xk,Xen,0._SP,CIRCUIT_k,(/ID_bands/),BS_bands,1,(/n_sp_pol/),"interpolated","bands",0)
   endif
   !
   ! ARPES interpolation
   !
   if (l_exciton_distribution) then
     call excitons_read(Xk,Xk,Xen,Xen,iq,"check",warn,io_err)
     if (io_err/=0) then
       call excitons_read(Xk,Xk,Xen,Xen,iq,"clean",warn,io_err)
       call live_timing(steps=(l_range(2)-l_range(1)+1))
       cycle
     endif
     call excitons_read(Xk,Xk,Xen,Xen,iq,"eigenvalues",warn,io_err)   
     call excitons_read(Xk,Xk,Xen,Xen,iq,"eigenvectors",warn,io_err)
   endif
   !
   do i_l=l_range(1),l_range(2)
     !
     ARPES_int = cZERO
     !
     if (l_single_exciton)       N_lq=1._SP
     if (l_exciton_distribution) N_lq=boltzman_f(real(BS_E(i_l))-min_exc_E)
     !
     if (N_lq<1.E-3) then
       if (l_exciton_distribution) call live_timing(steps=1)
       cycle
     endif
     !
     N_exc=N_exc+1
     Z=Z+N_lq
     !
     if (l_exciton_distribution) then
       j_exc=N_exc
       if(N_exc>1000) j_exc=minloc(EXC_N_lq,1)
       !
       EXC_table(j_exc,:)=(/i_l,iq/)
       EXC_N_lq(j_exc)   =N_lq
       EXC_E_lq(j_exc)   =BS_E(i_l)
     endif
     !
     !if (l_single_exciton)       iq=iq_bz
     !if (l_exciton_distribution) iq=                  ! ********* To Fix ************
     !
     do i_cvk = 1,BS_K_dim(1)
       !
       i_kbz    = BSS_eh_table(i_cvk,1)
       i_v      = BSS_eh_table(i_cvk,2)
       i_c      = BSS_eh_table(i_cvk,3)
       i_sp_pol = spin(BSS_eh_table(i_cvk,:))
       !
       i_kibz   = Xk%sstar(i_kbz,1)
       i_ks     = Xk%sstar(i_kbz,2)
       !
       i_kmq_bz = qindx_X(iq,i_kbz,1)
       i_kmq_ibz= Xk%sstar(i_kmq_bz,1) 
       !
       if ( (.not.l_full_expansion) .and. i_ks>1 ) cycle
       if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
       if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
       if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
       !
       if(     l_full_expansion) i_karpes=i_kbz
       if(.not.l_full_expansion) i_karpes=i_kibz
       !
       ! Quantity to be interpolated.
       ! Associate to every pole eps_(v,k-q)+E_l(q) an intensity
       ! Here we are inside the loop over excitons, so the last index is always one
       ! Instead I sum over the conduction band index
       !
       ARPES_int(i_v,i_karpes,1)=ARPES_int(i_v,i_karpes,1)+abs(BS_mat(i_cvk,i_l))**2
       !
       ! Constract the PES integrating over k
       !
       do i_w=1,RT_conf%W_n_steps
         F_tmp(i_w)=Lorentzian_FT(GreenF_T_and_W%W(i_w),(Xen%E(i_v,i_kmq_ibz,i_sp_pol)+BS_E(i_l)),"r")
       enddo
       !
       F_tmp=F_tmp*abs(BS_mat(i_cvk,i_l))**2*N_lq
       !
       PES(:)=PES(:)+F_tmp(:)
       !
     enddo
     !
     call PP_redux_wait(ARPES_int)
     !
     if (l_full_expansion) then
       !
       INTERPOLATE_use_symm=.false.
       call k_ibz2bz(Xk,'i',.true.)
       call bz_samp_duplicate_nosymm(Xk,k_nosymm)
       !
       call INTERPOLATION_driver_seed("ARPES",ID_ARPES_int,(/BS_bands(1),Xen%nbf(1)/),k_nosymm,&
       &                                             SEED_R3=ARPES_int,QUIET=l_exciton_distribution)
       !
     else
       !
       call INTERPOLATION_driver_seed("ARPES",ID_ARPES_int,(/BS_bands(1),Xen%nbf(1)/),Xk,SEED_R3=ARPES_int,&
       &                                                                           QUIET=l_exciton_distribution)
       !
     endif
     !
     call INTERPOLATION_driver_do("ARPES",ID_ARPES_int,CIRCUIT_k)
     !
     i_karpes=0
     !
     do i_kibz=1,CIRCUIT_K%nibz
       !
       if (all(abs(CIRCUIT_K%pt(i_kibz,:))<1.E-5_SP)) i_karpes=i_kibz
       !
       do i_v=BS_bands(1),Xen%nbf(1)
         !
         Ev=INTERP_obj(ID_bands)%INTERP(i_v,i_kibz,1)
         !
         do i_w=1,RT_conf%W_n_steps
           F_tmp(i_w)=Lorentzian_FT(GreenF_T_and_W%W(i_w),(Ev+BS_E(i_l)),"r")
         enddo
         !
	 ! Here I have the k-interpolated pole, times the k-interpolated intensity,
         ! times the exciton occupation. This is at fixed exciton and q-index
         F_tmp=F_tmp*INTERP_obj(ID_ARPES_int)%INTERP(i_v,i_kibz,1)*N_lq
         !
	 ! F_tmp is frequency dependent. The valence and exciton index are not
	 ! needed anymore. I can accumulate
         INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,1)=INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,1)+ real(F_tmp,SP)
         INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,2)=INTERP_obj(ID_ARPES)%INTERP(:,i_kibz,2)-aimag(F_tmp)
         !
       enddo
     enddo !i_k circuit
     !
     call INTERPOLATION_driver_end(ID_ARPES_int)
     !
     if (l_exciton_distribution) call live_timing(steps=1)
     !
   enddo ! i_lambda
   !
   call excitons_read(Xk,Xk,Xen,Xen,iq,"clean",warn,io_err)
   !
   YAMBO_FREE(INTERP_obj(ID_bands)%INTERP)
   !
 enddo  ! iq
 !
 call INTERPOLATION_driver_end(ID_bands)
 !
 if (l_exciton_distribution) call live_timing()
 !
 YAMBO_FREE(ARPES_int)
 !
 call PP_redux_wait(PES)
 !
 PES=PES/Z
 INTERP_obj(ID_ARPES)%INTERP=INTERP_obj(ID_ARPES)%INTERP/Z
 !
 if (l_exciton_distribution) then
   !
   ! Print the excitonic weights
   !
   call sort(arrin=EXC_N_lq(:N_exc),indx=EXC_indx(:N_exc))
   !
   file_name="Exc_distribution"
   call of_open_close(trim(file_name),'ot')
   dist_head=(/"E_lq","N_lq","i_l "," iq "/)
   call msg('o '//trim(file_name),"#",dist_head,INDENT=0,USE_TABS=.true.)
   call msg('o '//trim(file_name),"#")
   !
   do j_exc=N_exc,1,-1
     dist_data=(/EXC_E_lq(EXC_indx(j_exc))*HA2EV,EXC_N_lq(j_exc)/Z,&
     & real(EXC_table(EXC_indx(j_exc),1),SP),real(EXC_table(EXC_indx(j_exc),2),SP)/)
     call msg('o '//trim(file_name),'',dist_data,INDENT=-2,USE_TABS=.true.)
   enddo
   !
   call of_open_close(trim(file_name))
   !
 endif
 !
 ! Print the integrated PES
 !
 file_name="G_w_integrated"
 call of_open_close(trim(file_name),'ot')
 headings=(/" E[eV]  ","Im[G(w)]","Re[G(w)]"/)
 call msg('o '//trim(file_name),"#",headings,INDENT=0,USE_TABS=.true.)
 !
 do i_w=1,RT_conf%W_n_steps
   odata=(/real(GreenF_T_and_W%W(i_w),SP)*HA2EV,real(cI*PES(i_w),SP),real(PES(i_w),SP)/)
   call msg('o '//trim(file_name),'',odata,INDENT=-2,USE_TABS=.true.)
 enddo
 !
 call of_open_close(trim(file_name))
 !
 ! DEBUG <
 !!
 !! Print G(w,k) for specific k  (only gamma for now)
 !!
 !file_name="G_k1_w"
 !call of_open_close(trim(file_name),'ot')
 !headings=(/" E[eV]  ","Im[G(w)]","Re[G(w)]"/)
 !call msg('o '//trim(file_name),"#",headings,INDENT=0,USE_TABS=.true.)
 !call msg('o '//trim(file_name),"#")
 !!
 !do i_w=1,RT_conf%W_n_steps
 !  odata=(/real(GreenF_T_and_W%W(i_w),SP)*HA2EV,      &
 !  &      INTERP_obj(ID_ARPES)%INTERP(i_w,i_karpes,2), &
 !  &      INTERP_obj(ID_ARPES)%INTERP(i_w,i_karpes,1) /)
 !  call msg('o '//trim(file_name),'',odata,INDENT=-2,USE_TABS=.true.)
 !enddo
 !!
 !call of_open_close(trim(file_name))
 !!
 ! DEBUG >
 !
 call section('=','Plotting')
 !============================================================================= 
 !
 if(CIRCUIT_k%nbz>0) call plot_interpolated_values(Xk,Xen,0._SP,CIRCUIT_k,&
 &                                                  (/ID_ARPES/),(/1,RT_conf%W_n_steps/),1,(/2/),"interpolated","G_k_w",0)
 !if(   USER_k%nbz>0) call plot_interpolated_values(Xk,Xen,0._SP,USER_k,   &
 !&                                                  (/ID_ARPES/),(/1,NBW/),1,(/2/),"built_in",    "G_k_w",0)
 !
 !
 call INTERPOLATION_driver_end(0)
 !
end subroutine RT_transient_ARPES
