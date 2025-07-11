!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): MG CA DS AF
!
#include<dev_defs.h>                                                    
!
subroutine WF_symm_kpoint_g(b_to_load,ik,i_sp_pol,isymm,ig0,WF,ssop,isref,igref)
 !
 use pars,          ONLY:SP,cZERO
 use electrons,     ONLY:n_spinor
 use wave_func,     ONLY:WAVEs,wf_nc_k,DEV_VAR(wf_igk)
 use wave_func,     ONLY:wf_ng,wf_ng_1st_BZ
 use vec_operate,   ONLY:c2a
 use R_lattice,     ONLY:DEV_VAR(g_rot),DEV_VAR(g_phs),DEV_VAR(G_m_G),&
 &                       k_pt,rl_sop,g_vec
 use D_lattice,     ONLY:idt_index,DEV_VAR(spin_sop),sop_inv,dl_tra, &
 &                       dl_sop,i_time_rev,nsym
 !
 implicit none
 !
 integer,       intent(in) :: b_to_load(2)
 integer,       intent(in) :: ik,i_sp_pol,isymm,ig0
 type (WAVEs), target, intent(inout) :: WF
 complex(SP) DEV_ATTR, optional, intent(in) :: ssop(2,2)
 integer,              optional, intent(in) :: isref,igref
 !
 ! Work space
 !
 complex(SP)              DEV_ATTR :: WF_k_tmp(wf_ng_1st_BZ,n_spinor)
 complex(SP), pointer     DEV_ATTR :: WF_p(:,:)
 complex(SP)              DEV_ATTR :: ssop_loc(2,2)
 integer                  DEV_ATTR :: wf_tab_loc(wf_ng,2)
 !
 real(SP)    :: krot(3),G_cc(3)
 complex(SP) :: gphase,kphase
 !
 integer :: ib,i_wf,isymm_m1
 integer :: ic,ig,igp,igr,igs,i_spinor
 !
#if defined _GPU
 call error("subroutine WF_symm_kpoint_g is not gpu ported")
#endif
 !
 !if(isymm==idt_index) write(*,*) k_pt(ik,:)
 if(isymm==idt_index) return
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
 if (n_spinor==2) then
   if(present(ssop)) then
     ssop_loc=ssop
   else
     ssop_loc=DEV_VAR(spin_sop)(:,:,isymm)
   endif
 endif
 !
 if (present(igref).or.present(isref)) then
   if (.not.present(isref)) call error("WF_symm_kpoint_g: igref provided, but not isref")
   if (.not.present(igref)) call error("WF_symm_kpoint_g: isref provided, but not igref")
   do ic=1,wf_nc_k(ik)
     ig=DEV_VAR(wf_igk)(ic,ik)
     igp=DEV_VAR(g_rot)(ig,isref)
     igs=DEV_VAR(G_m_G)(igp,igref)
     igr=DEV_VAR(g_rot)(igs,isymm)
     wf_tab_loc(ic,1)=igr
     wf_tab_loc(ic,2)=DEV_VAR(G_m_G)(igr,ig0)
   enddo
 else
   do ic=1,wf_nc_k(ik)
     ig=DEV_VAR(wf_igk)(ic,ik)
     igr=DEV_VAR(g_rot)(ig,isymm)
     wf_tab_loc(ic,1)=igr
     wf_tab_loc(ic,2)=DEV_VAR(G_m_G)(igr,ig0)
   enddo
 endif
 !
 do ib=b_to_load(1),b_to_load(2)
   !
   WF_k_tmp=cZERO
   !
   i_wf=WF%index(ib,ik,i_sp_pol)
   WF_p => DEV_VAR(WF%c)(:,:,i_wf)
   !
   if (n_spinor==1) then
     !
     !DEV_ACC_DEBUG data present(WF_k_tmp,g_rot,WF_p)
     !DEV_ACC parallel loop private(ig,igr,igs)
     !DEV_CUF kernel do(1) <<<*,*>>>
     !DEV_OMPGPU target map(present,alloc:WF_k_tmp,g_rot,WF_p)
     !DEV_OMPGPU teams loop
     do ic=1,wf_nc_k(ik)
       igr=wf_tab_loc(ic,1)
       igs=wf_tab_loc(ic,2)
       gphase=cmplx(cos(DEV_VAR(g_phs)(igr,isymm)),sin(DEV_VAR(g_phs)(igr,isymm)))
       WF_k_tmp(igs,1)=WF_p(ig,1)*conjg(gphase)
     enddo
     !
   else
     !
     !DEV_ACC_DEBUG data present(WF_k_tmp,g_rot,WF_p)
     !DEV_ACC parallel loop private(ig,igr,igs)
     !DEV_CUF kernel do(1) <<<*,*>>>
     !DEV_OMPGPU target map(present,alloc:WF_k_tmp,g_rot,WF_p)
     !DEV_OMPGPU teams loop
     do ic=1,wf_nc_k(ik)
       igr=wf_tab_loc(ic,1)
       igs=wf_tab_loc(ic,2)
       gphase=cmplx(cos(DEV_VAR(g_phs)(igr,isymm)),sin(DEV_VAR(g_phs)(igr,isymm)))
       WF_k_tmp(igs,1)=ssop_loc(1,1)*WF_p(ig,1)+ssop_loc(1,2)*WF_p(ig,2)
       WF_k_tmp(igs,2)=ssop_loc(2,1)*WF_p(ig,1)+ssop_loc(2,2)*WF_p(ig,2)
       WF_k_tmp(igs,:)=WF_k_tmp(igs,:)*conjg(gphase)
     enddo
   endif
   !
   ! Take care of translations
   WF_k_tmp=WF_k_tmp*conjg(kphase)
   !
   ! take care of time reversal if needed
   !if(isymm>nsym/(1+i_time_rev)) call dev_conjg(WF_k_tmp)
   if(isymm>nsym/(1+i_time_rev)) WF_k_tmp=conjg(WF_k_tmp)
   !
   !call dev_memcpy(WF_p,WF_k_tmp)
   WF_p=WF_k_tmp
   !
   nullify(WF_p)
   !
 enddo
 !
end subroutine WF_symm_kpoint_g
