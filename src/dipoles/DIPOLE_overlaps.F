!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA
!
! headers
!
#include<dev_defs.h>
#include<y_memory.h>
!      
!> @brief Calculate overlaps between wave-function <\psi_k | \psi_{k+q} > for q = 1,2
!
!! @param[it]       Xk                  k-points sampling for the xhi
!! @param[it]       Dip                 Dip dipoles range
!! @param[out]      DIP_S               Overlap matrices
!
!
subroutine DIPOLE_overlaps(Xk,Dip)
 !    
 use pars,              ONLY:SP,cZERO,pi
 use LIVE_t,            ONLY:live_timing
 use R_lattice,         ONLY:k_map,bz_samp,WF_shifts,G_m_G,minus_G
 use wave_func,         ONLY:WF,WAVEs,WF_load_mode
 use vec_operate,       ONLY:k_periodic_idx
 use DIPOLES,           ONLY:DIPOLE_t,DIP_S
 use ALLOC,             ONLY:DIPOLE_ALLOC_elemental
 use electrons,         ONLY:n_spinor,n_sp_pol
 use wrapper,           ONLY:Vstar_dot_V_gpu
 use wave_func,         ONLY:wf_ng_1st_BZ,wf_ng_overlaps,wf_ng,WF_copy
 use parallel_m,        ONLY:PAR_IND_DIPk_bz,PAR_IND_DIPk_bz_ID,PAR_COM_DIPk_ibz_A2A, &
&                            PAR_IND_CON_BANDS_OVLP,PAR_IND_VAL_BANDS_OVLP,           &
&                            PAR_IND_OVLPk_ibz,PAR_COM_DIPk_ibz_INDEX
 use parallel_int,      ONLY:PP_redux_wait,PARALLEL_WF_index,PARALLEL_WF_distribute
 use interfaces,        ONLY:WF_symm_kpoint_g2,eval_G_minus_G,WF_load,WF_free
 use gpu_m,             ONLY:have_gpu
 use timing_m,          ONLY:timing
 use y_memory_alloc
 !
 implicit none
 !
 type(bz_samp), intent(in)    :: Xk 
 type(DIPOLE_t),intent(in)    :: Dip
 !
 ! Work Space
 !
 integer           :: ikbz,ikbzp,is,isp,ib,ibp,i_sp_pol,ik,ikp,istep,max_step
 integer           :: iG0,ig,igp,ig_shift,ng_tmp,wf_ng_save,wf_ng_1st_BZ_save
 integer           :: id,idp,idx_kp(3),idx_k(3),shift(3),g0_idx(3,2),i_wf,i_wfp
 real(SP)          :: g0_length(3)
 !
 !WF_skg1
 !type (WAVEs) :: WF_ik,WF_ikp
 !WF_skg2
 complex(SP), allocatable DEV_ATTR :: WF_ik(:,:,:)
 complex(SP), allocatable DEV_ATTR :: WF_ikp(:,:,:)
 !
 logical           :: USE_shifted_wf
 !
 call timing("DIPOLE_overlaps",OPR="start")
 !
 if(.not.allocated(k_map%k_map_dir)) call k_map_nearest_by_dir(Xk,k_map)
 if(.not.allocated(Xk%k_table))      call k_build_up_BZ_tables(Xk)
 !
 DIP_S=cZERO
 !
 call k_find_smallest_g0(g0_idx,g0_length)
 iG0=maxval(g0_idx)
 !
 ! Generate G_m_G of the needed size
 !
 YAMBO_FREE(G_m_G)
 ng_tmp=eval_G_minus_G(iG0,iG0,opr_in="+")
 !
 ! New iG0 for composed shifts
 iG0=maxval(G_m_G)
 !
 YAMBO_FREE(G_m_G)
 ng_tmp=eval_G_minus_G(wf_ng_1st_BZ,iG0)
 !
 wf_ng_save       =wf_ng
 wf_ng_1st_BZ_save=wf_ng_1st_BZ
 !
 if (maxval(G_m_G(:,1)) < wf_ng_1st_BZ) then
   !
   call warning(" Forced lowering on cutoff on WFs")
   wf_ng_1st_BZ  =maxval(G_m_G(:,1))
   wf_ng         =maxval(G_m_G(:,1))
   !
 endif
 !
 wf_ng_overlaps=maxval(G_m_G)
 !
 call PARALLEL_WF_distribute(K_index=PAR_IND_OVLPk_ibz,CLEAN_UP=.TRUE.)
 !
 call PARALLEL_WF_index( )
 !
 if (trim(WF_load_mode)=="all") then
   !WF_skg1
   !call error("not implemented with loading of all WFs")
   !call WF_load(WF,0,1,(/1,Dip%ib(2)/),(/1,Xk%nibz/),&
   !&            wf_ng_in=wf_ng_overlaps,space='G',title='-Overlaps')
   !WF_skg2
   call WF_load(WF,0,1,(/1,Dip%ib(2)/),(/1,Xk%nibz/),&
   &            space='G',title='-Overlaps')
 endif
 !
 YAMBO_ALLOC_GPU(WF_ikp,(wf_ng_overlaps,n_spinor,Dip%ib(2)))
 YAMBO_ALLOC_GPU(WF_ik, (wf_ng_overlaps,n_spinor,Dip%ib(2)))
 !
!$OMP WORKSHARE
 DIP_S=cZERO
!$OMP END WORKSHARE
 !
 call live_timing("Overlaps",PAR_IND_DIPk_bz%n_of_elements(PAR_IND_DIPk_bz_ID+1)*n_sp_pol)
 !
 do i_sp_pol=1,n_sp_pol
   !
   do ikbz=1,Xk%nbz
     !
     if(.not.PAR_IND_DIPk_bz%element_1D(ikbz)) cycle
     !
     idx_k=k_map%k_map_inv(ikbz,:)
     !
     ik = Xk%sstar(ikbz,1)
     is = Xk%sstar(ikbz,2)
     !
     if (trim(WF_load_mode)=="on-the-fly") then
       !WF_skg1
       !call WF_load(WF,0,1,(/1,Dip%ib(2)/),(/ik,ik/),(/i_sp_pol,i_sp_pol/),&
       !&            wf_ng_in=wf_ng_overlaps,space='G',title='-Oscill/G space/Overlaps',quiet=.true.)
       !call WF_copy(WF,WF_ik)
       !WF_skg2
       call WF_load(WF,0,1,(/1,Dip%ib(2)/),(/ik,ik/),(/i_sp_pol,i_sp_pol/),&
       &            space='G',title='-Oscill/G space/Overlaps',quiet=.true.)
     endif
     !
     shift=WF_shifts(ikbz,:)
     ig_shift=1
     do id=1,3
       if(shift(id)==0) cycle
       ig=k_map%g0_idx(id,shift(id))
       igp=minus_G(ig_shift)
       ig_shift=G_m_G(ig,igp)
     enddo 
     !
     !WF_skg1
     !call WF_symm_kpoint_g((/1,Dip%ib(2)/),ik,i_sp_pol,is,ig_shift,WF_ik)
     !WF_skg2
     call WF_symm_kpoint_g2((/1,Dip%ib(2)/),ik,i_sp_pol,is,ig_shift,WF,WF_ik)
     !
     if (trim(WF_load_mode)=="on-the-fly") call WF_free(WF)
     !
     do id=1,3  ! Loop on the three crystal directions in the BZ
       !
       if(k_map%max_kdir(id)==1) cycle ! Non-periodic direction, skip calculation
       !
       max_step=1 ! Only the first neighbor
       if(k_map%max_kdir(id)>=6) max_step=2  ! First and second neighbors
       !
       do istep=1,max_step
         !
         idx_kp=idx_k
         idx_kp(id)=idx_kp(id)+istep   ! Neighbor on the right
         !
         USE_shifted_wf=.false.
         if(idx_kp(id)>k_map%max_kdir(id)) USE_shifted_wf=.true. ! Check if the neighbor is out of the BZ
         !
         idx_kp=k_periodic_idx(idx_kp,k_map)  ! Bring it back by periodicity
         ikbzp=k_map%k_map_dir(idx_kp(1),idx_kp(2),idx_kp(3))  ! Find the corresponding index in the BZ
         !
         ! Find the index in the IBZ and the symmetry
         ikp = Xk%sstar(ikbzp,1)
         isp = Xk%sstar(ikbzp,2)
         !
         if (trim(WF_load_mode)=="on-the-fly") then
           call WF_load(WF,0,1,(/1,Dip%ib(2)/),(/ikp,ikp/),(/i_sp_pol,i_sp_pol/),&
           &              space='G',title='-Oscill/G space/Overlaps',quiet=.true.)
           !WF_skg1
           !call WF_copy(WF,WF_ikp)
         !WF_skg1
         !else
         !  ! to be fixed
         !  !AMBO_ALLOC(WF_ikp%c,(wf_ng_overlaps,n_spinor,WF%N))
         !  WF_ikp%c=cZERO
         endif
         !
         ! Shift the wave-function by a G-vector if the neighbor is out of the BZ (USE_shifed_wf=.true.)
         shift=WF_shifts(ikbzp,:)
         if(USE_shifted_wf) then
           !
           ! In this case a points k+q is out of the BZ
           ! now I have different possibilities:
           ! 1) the original point ikbz that corresponds to (k+q) was
           !    between 0-1, therefore I have to shift it to generate
           !    the k+q => case WF_shifts==0
           ! 2) the orignal points ikbz that corresponds to (k+q) was 
           !    shifted in the opposite direction the two
           !    shifts cancel out and I get shift==0 ==> WF_shifts==-1
           ! 3) the orignal points ikbz that corresponds to (k+q) was
           !    already shifted in the right direction I keep this shift
           !    case WF_shifts==1    
           ! 
           ! see also comments in src/bz_ops/k_map_nearest_by_dir.F 
           ! 
           ! WF_shifts==1 --> kpt(id) < 0
           if(WF_shifts(ikbzp,id)==1) shift(id)=1
           ! WF_shifts==0 --> 0 < kpt(id) < 1
           if(WF_shifts(ikbzp,id)==0) shift(id)=1
           ! WF_shifts==-1 --> kpt(id) > 1
           if(WF_shifts(ikbzp,id)==-1) shift(id)=0
           !
         endif
         !
         ig_shift=1
         do idp=1,3
           if(shift(idp)==0) cycle
           ig=k_map%g0_idx(idp,shift(idp))
           igp=minus_G(ig_shift)
           ig_shift=G_m_G(ig,igp)
         enddo 
	 !
	 !WF_skg1
	 !call WF_symm_kpoint_g((/1,Dip%ib(2)/),ikp,i_sp_pol,isp,ig_shift,WF_ikp)
         !WF_skg2
	 call WF_symm_kpoint_g2((/1,Dip%ib(2)/),ikp,i_sp_pol,isp,ig_shift,WF,WF_ikp)
         !
         if (trim(WF_load_mode)=="on-the-fly") call WF_free(WF)
         !
         do ibp=1,Dip%ib(2)
           if(.not.PAR_IND_CON_BANDS_OVLP%element_1D(ibp)) cycle
	   !WF_skg1
	   !i_wfp=WF_ik%index(ibp,ikp,i_sp_pol)
           do ib=1,Dip%ib(2)
             if(.not.PAR_IND_VAL_BANDS_OVLP%element_1D(ib)) cycle
	     !WF_skg1
	     !i_wf=WF_ik%index(ib,ik,i_sp_pol)
             !DIP_S(ib,ibp,id+(istep-1)*3,ikbz,i_sp_pol)= &
             !& Vstar_dot_V_gpu(wf_ng_overlaps*n_spinor,WF_ik%c(:,:,i_wf),WF_ikp%c(:,:,i_wfp))
	     !WF_skg2
             DIP_S(ib,ibp,id+(istep-1)*3,ikbz,i_sp_pol)= &
             & Vstar_dot_V_gpu(wf_ng_overlaps*n_spinor,WF_ik(:,:,ib),WF_ikp(:,:,ibp))
           enddo ! ibp
         enddo ! ib
         !
	 !WF_skg1
         !call WF_free(WF_ikp)
         !
       enddo ! istep
     enddo ! idir
     !
     call PP_redux_wait(DIP_S(:,:,:,ikbz,i_sp_pol),COMM=PAR_COM_DIPk_ibz_A2A%COMM)
     !
     call live_timing(steps=1)
     !
   enddo ! ikbz
   !
   !WF_skg1
   !call WF_free(WF_ik)
   !
   ! I need to distributed overlaps among all processors because
   ! the construction of covariant dipoles require the first and second 
   ! neighboar of each k-point
   !
   call PP_redux_wait(DIP_S(:,:,:,:,i_sp_pol),COMM=PAR_COM_DIPk_ibz_INDEX%COMM)
   !
 enddo ! i_sp_pol
 !
 call live_timing()
 !
 wf_ng       =wf_ng_save
 wf_ng_1st_BZ=wf_ng_1st_BZ_save
 !
 !WF_skg2
 if (have_gpu) then
   YAMBO_FREE_GPU(WF_ik)
   YAMBO_FREE_GPU(WF_ikp)
 endif
 !
 YAMBO_FREE(WF_ik)
 YAMBO_FREE(WF_ikp)
 !
 if (trim(WF_load_mode)=="all") call WF_free(WF)
 !
 call timing("DIPOLE_overlaps",OPR="stop")
 !
end subroutine Dipole_Overlaps
