!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA DS AF
!
! headers
!
#include<dev_defs.h>                                                    
!
!
subroutine WF_symm_kpoint_g_tmp(b_to_load,ik,i_sp_pol,isymm,ig0,WF,WF_k_out)
 !
 ! Subroutine needed by dipole overlaps
 !
 use pars,          ONLY:SP,cZERO
 use electrons,     ONLY:n_spinor
 use wave_func,     ONLY:WAVEs,wf_nc_k,DEV_VAR(wf_igk)
 use wave_func,     ONLY:wf_ng,wf_ng_overlaps
 use vec_operate,   ONLY:c2a
 use R_lattice,     ONLY:DEV_VAR(g_rot),DEV_VAR(g_phs),DEV_VAR(G_m_G),&
 &                       k_pt,rl_sop,g_vec
 use D_lattice,     ONLY:idt_index,DEV_VAR(spin_sop),sop_inv,dl_tra, &
 &                       dl_sop,S_contains_TR,nsym
#if defined _GPU
 use devxlib,       ONLY:devxlib_memcpy_d2d,devxlib_conjg_d
#endif
 !
 implicit none
 !
 integer,       intent(in) :: b_to_load(2)
 integer,       intent(in) :: ik,i_sp_pol,isymm,ig0
 type (WAVEs), target, intent(inout) :: WF
 complex(SP) DEV_ATTR, target, intent(out) :: WF_k_out(wf_ng_overlaps,n_spinor,b_to_load(1):b_to_load(2))
 !
 ! Work space
 !
 complex(SP)              DEV_ATTR :: WF_k_tmp(wf_ng_overlaps,n_spinor)
 complex(SP), pointer     DEV_ATTR :: WF_p(:,:)
 complex(SP), pointer     DEV_ATTR :: WF_p_out(:,:)
 complex(SP)              DEV_ATTR :: ssop_loc(n_spinor,n_spinor)
 integer                  DEV_ATTR :: wf_tab_loc(wf_ng,3)
 !
 real(SP)    :: krot(3),G_cc(3)
 complex(SP) :: gphase,kphase
 !
 integer :: ib,i_wf,isymm_m1
 integer :: ic,ic_max,ig,igp,igr,igs,i_spinor
 !
 ic_max=wf_nc_k(ik)
 do ic=1,wf_nc_k(ik)
   ig=DEV_VAR(wf_igk)(ic,ik)
   if(ig<=wf_ng) cycle
   ic_max=ic-1
   exit
 enddo
 !
 if(isymm==idt_index.and.ig0==1) then
   WF_k_out=cZERO
   do ib=b_to_load(1),b_to_load(2)
     i_wf=WF%index(ib,ik,i_sp_pol)
     WF_p => DEV_VAR(WF%c)(:,:,i_wf)
     WF_p_out => WF_k_out(:,:,ib)
     !
     if (n_spinor==1) then
       !
       !DEV_ACC_DEBUG data present(WF_p_out,WF_p)
       !DEV_ACC parallel loop private(ig)
       !DEV_CUF kernel do(1) <<<*,*>>>
       !DEV_OMPGPU target map(present,alloc:WF_p_out,WF_p)
       !DEV_OMPGPU teams loop
       do ic=1,ic_max
         ig=DEV_VAR(wf_igk)(ic,ik)
         WF_p_out(ig,1)=WF_p(ig,1)
       enddo
       !DEV_OMPGPU end target
       !DEV_ACC_DEBUG end data
       !
     else
       !
       !DEV_ACC_DEBUG data present(WF_p_out,WF_p)
       !DEV_ACC parallel loop private(ig)
       !DEV_CUF kernel do(1) <<<*,*>>>
       !DEV_OMPGPU target map(present,alloc:WF_p_out,WF_p)
       !DEV_OMPGPU teams loop
       do ic=1,ic_max
         ig=DEV_VAR(wf_igk)(ic,ik)
         WF_p_out(ig,1)=WF_p(ig,1)
         WF_p_out(ig,2)=WF_p(ig,2)
       enddo
       !DEV_OMPGPU end target
       !DEV_ACC_DEBUG end data
       !
     endif
     !
     nullify(WF_p)
     nullify(WF_p_out)
   enddo
   return
 endif
 !
 isymm_m1=sop_inv(isymm)
 !  
 !krot=k_pt(ik,:)
 !call c2a(v_in=krot,mode='ki2c')
 !krot=matmul(dl_sop(:,:,isymm),krot)
 krot=matmul(rl_sop(:,:,isymm),k_pt(ik,:))
 call c2a(v_in=krot,mode='ki2c')
 kphase=cmplx(cos(dot_product(dl_tra(:,isymm),krot)),sin(dot_product(dl_tra(:,isymm),krot)))
 !
 G_cc=g_vec(ig0,:)
 call c2a(v_in=G_cc,mode='ki2c')
 gphase=cmplx(cos(dot_product(dl_tra(:,isymm),G_cc)),sin(dot_product(dl_tra(:,isymm),G_cc)))
 !
 kphase=kphase*gphase
 !
 !if(ik==13) then
 !  call c2a(v_in=krot,mode='kc2i')
 !  call c2a(v_in=G_cc,mode='kc2i')
 !  write(*,*) isymm,isymm_m1,ig0,kphase,gphase,dot_product(dl_tra(:,isymm),G_cc),G_cc,krot
 !endif
 !
 if (n_spinor==2) ssop_loc=DEV_VAR(spin_sop)(:,:,isymm)
 !
 do ic=1,ic_max
   ig=DEV_VAR(wf_igk)(ic,ik)
   igr=DEV_VAR(g_rot)(ig,isymm)
   wf_tab_loc(ic,1)=ig
   wf_tab_loc(ic,2)=igr
   wf_tab_loc(ic,3)=DEV_VAR(G_m_G)(igr,ig0)
 enddo
 !
 do ib=b_to_load(1),b_to_load(2)
   !
   WF_k_tmp=cZERO
   !
   i_wf=WF%index(ib,ik,i_sp_pol)
   WF_p => DEV_VAR(WF%c)(:,:,i_wf)
   WF_p_out => WF_k_out(:,:,ib)
   !
   if (n_spinor==1) then
     !
     !DEV_ACC_DEBUG data present(WF_k_tmp,g_rot,WF_p)
     !DEV_ACC parallel loop private(ig,igr,igs)
     !DEV_CUF kernel do(1) <<<*,*>>>
     !DEV_OMPGPU target map(present,alloc:WF_k_tmp,g_rot,WF_p)
     !DEV_OMPGPU teams loop
     do ic=1,ic_max
       ig =wf_tab_loc(ic,1)
       igr=wf_tab_loc(ic,2)
       igs=wf_tab_loc(ic,3)
       gphase=cmplx(cos(DEV_VAR(g_phs)(igr,isymm)),sin(DEV_VAR(g_phs)(igr,isymm)))
       WF_k_tmp(igs,1)=WF_p(ig,1)*conjg(gphase)*conjg(kphase)
     enddo
     !DEV_OMPGPU end target
     !DEV_ACC_DEBUG end data
     !
   else
     !
     !DEV_ACC_DEBUG data present(WF_k_tmp,g_rot,WF_p)
     !DEV_ACC parallel loop private(ig,igs)
     !DEV_CUF kernel do(1) <<<*,*>>>
     !DEV_OMPGPU target map(present,alloc:WF_k_tmp,g_rot,WF_p)
     !DEV_OMPGPU teams loop
     do ic=1,ic_max
       ig =wf_tab_loc(ic,1)
       igr=wf_tab_loc(ic,2)
       igs=wf_tab_loc(ic,3)
       gphase=cmplx(cos(DEV_VAR(g_phs)(igr,isymm)),sin(DEV_VAR(g_phs)(igr,isymm)))
       WF_k_tmp(igs,1)=ssop_loc(1,1)*WF_p(ig,1)+ssop_loc(1,2)*WF_p(ig,2)
       WF_k_tmp(igs,2)=ssop_loc(2,1)*WF_p(ig,1)+ssop_loc(2,2)*WF_p(ig,2)
       WF_k_tmp(igs,:)=WF_k_tmp(igs,:)*conjg(gphase)*conjg(kphase)
     enddo
     !DEV_OMPGPU end target
     !DEV_ACC_DEBUG end data
   endif
   !
   ! take care of time reversal if needed
#if defined _GPU
   if(S_contains_TR(isymm)) call devxlib_conjg_d(WF_k_tmp)
   call devxlib_memcpy_d2d(WF_p_out,WF_k_tmp)
#else
   if(S_contains_TR(isymm)) WF_k_tmp=conjg(WF_k_tmp)
   WF_p_out=WF_k_tmp
#endif
   !
   nullify(WF_p)
   nullify(WF_p_out)
   !
 enddo
 !
end subroutine WF_symm_kpoint_g_tmp
