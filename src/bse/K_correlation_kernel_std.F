!
! License-Identifier: GPL
!
! Copyright (C) 2018 The Yambo Team
!
! Authors (see AUTHORS file for details): AM, DS
!
! headers
!
#include<dev_defs.h>
#include<y_memory.h>
!
!************************* WARNING **************************
!* This subroutine is not correct in presence of            *
!* non symmorphyc symmetri operations                       * 
!* An alternative subroutine that applies symmetries in a   *
!* different way exist: K_correlation_kernel_dir.F          *
!* This second subroutine is automatically called in the    *
!* such case                                                *
!************************************************************
!
function K_correlation_kernel_std(i_block,i_p,i_pmq,                               &
&                                 i_k_s,i_kp_s,i_n_k,i_n_p,i_kmq_s,i_kp_mq_s,i_m_k,i_m_p,  &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_n,i_p_sp_pol_n,i_k_sp_pol_m,i_p_sp_pol_m, &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1,&
&                                 BS_n_g_W,O1,O2,O_times_W)
 !
 ! Standard implementation. This is the default algorithm.
 !
 ! Alternative implementation to K_correlation_kernel_dir
 ! Numerically more efficient,
 ! although more complex than the direct implementation
 !
 use pars,         ONLY:SP,cZERO,cONE,pi
 use BS,           ONLY:BS_blk,BS_W,BS_W_is_diagonal,WF_phase
 use D_lattice,    ONLY:nsym,S_contains_TR
 use R_lattice,    ONLY:DEV_VAR(G_m_G),DEV_VAR(g_rot),DEV_VAR(g_phs)
 use wrapper,      ONLY:V_by_V_pwise_gpu,Vstar_dot_V_gpu
 use devxlib,      ONLY:devxlib_conjg_d,devxlib_xgemv_gpu
 use y_memory_alloc
 !
 implicit none
 !
 integer,     intent(in) :: i_block,i_p,i_pmq,                                            &
&                                 i_k_s,i_kp_s,i_n_k,i_n_p,i_kmq_s,i_kp_mq_s,i_m_k,i_m_p, &
&                                 i_kmq_t,i_pmq_t,i_k_sp_pol_n,i_p_sp_pol_n,i_k_sp_pol_m,i_p_sp_pol_m, &
&                                 iq_W,iq_W_s,ig_W,i_k_s_m1,iq_W_s_mq,ig_W_mq,i_kmq_s_m1
 integer,     intent(in) :: BS_n_g_W
 complex(SP), intent(inout) DEV_ATTR  :: O1(BS_n_g_W),O2(BS_n_g_W),O_times_W(BS_n_g_W)
 !
 complex(SP)  :: K_correlation_kernel_std
 !
 integer      :: iO1,i_block_1,iO2,i_block_2,i_g1,i_g2,i_g3,i_g2p,i_g3p
 complex(SP)  :: PHASE_1,PHASE_2,gphase1,gphase2
 complex(SP), pointer DEV_ATTR :: O_c1_p(:,:), O_c2_p(:,:)
 !
 if (i_k_sp_pol_n/=i_p_sp_pol_n .or. i_k_sp_pol_m/=i_p_sp_pol_m) then
   K_correlation_kernel_std=cZERO
   return
 endif
 !
 !
 iO1       = BS_blk(i_block)%O_table(1,i_kp_s,1,1,i_n_k,i_n_p,i_k_sp_pol_n)
 i_block_1 = BS_blk(i_block)%O_table(2,i_kp_s,1,1,i_n_k,i_n_p,i_k_sp_pol_n)
 PHASE_1=WF_phase(i_p,i_kp_s,i_n_p,i_p_sp_pol_n)
 if (PHASE_1==cZERO  ) PHASE_1=1._SP
 if (PHASE_1==-99._SP) PHASE_1=1._SP
 !
 iO2       = BS_blk(i_block)%O_table(1,i_kp_mq_s,i_kmq_t,i_pmq_t,i_m_k,i_m_p,i_k_sp_pol_m)
 i_block_2 = BS_blk(i_block)%O_table(2,i_kp_mq_s,i_kmq_t,i_pmq_t,i_m_k,i_m_p,i_k_sp_pol_m)
 PHASE_2=WF_phase(i_pmq,i_kp_mq_s,i_m_p,i_p_sp_pol_m)
 if (PHASE_2==cZERO  ) PHASE_2=1._SP
 if (PHASE_2==-99._SP) PHASE_2=1._SP
 !
 !
 O_c1_p => DEV_VAR(BS_blk(i_block_1)%O_c)(:,:)
 O_c2_p => DEV_VAR(BS_blk(i_block_2)%O_c)(:,:)
 !
 !DEV_ACC_DEBUG data present(O_c1_p,O_c2_p,G_m_G,g_rot,g_phs,O1,O2)
 !DEV_ACC parallel loop private(i_g2,i_g3,i_g2p,i_g3p,gphase1,gphase2)
 !DEV_CUF kernel do(1) <<<*,*>>>
 !DEV_OMP parallel do default(shared), private(i_g1,i_g2,i_g3,i_g2p,i_g3p,gphase1,gphase2)
 do i_g1=1,BS_n_g_W
   i_g2=DEV_VAR(G_m_G)(DEV_VAR(g_rot)(i_g1,iq_W_s)   ,ig_W   )
   i_g3=DEV_VAR(G_m_G)(DEV_VAR(g_rot)(i_g1,iq_W_s_mq),ig_W_mq)
   gphase1=cmplx(cos(DEV_VAR(g_phs)(i_g2,iq_W_s)   ),sin(DEV_VAR(g_phs)(i_g2,iq_W_s)   ))
   gphase2=cmplx(cos(DEV_VAR(g_phs)(i_g3,iq_W_s_mq)),sin(DEV_VAR(g_phs)(i_g3,iq_W_s_mq)))
   i_g2p=DEV_VAR(g_rot)(i_g2,i_k_s_m1  )
   i_g3p=DEV_VAR(g_rot)(i_g3,i_kmq_s_m1)
   gphase1=cmplx(cos(DEV_VAR(g_phs)(i_g2p,i_k_s_m1)  ),sin(DEV_VAR(g_phs)(i_g2p,i_k_s_m1)  ))*gphase1
   gphase2=cmplx(cos(DEV_VAR(g_phs)(i_g3p,i_kmq_s_m1)),sin(DEV_VAR(g_phs)(i_g3p,i_kmq_s_m1)))*gphase2
   O1(i_g1)=O_c1_p(i_g2p,iO1)*PHASE_1*gphase1
   O2(i_g1)=O_c2_p(i_g3p,iO2)*PHASE_2*gphase2
 enddo
 !DEV_OMP end parallel do
 !DEV_ACC_DEBUG end data
 !
 !
 if (S_contains_TR(i_k_s)  ) call devxlib_conjg_d(O1,range1=(/1,BS_n_g_W/))
 if (S_contains_TR(i_kmq_s)) call devxlib_conjg_d(O2,range1=(/1,BS_n_g_W/))
 !
 if (     BS_W_is_diagonal) call V_by_V_pwise_gpu(BS_n_g_W,O_times_W,O1,BS_W(:,1,iq_W))
 if (.not.BS_W_is_diagonal) call devxlib_xgemv_gpu('T',BS_n_g_W,BS_n_g_W,cONE, &
&                                        BS_W(:,:,iq_W),BS_n_g_W,O1,1,cZERO,O_times_W,1)
 !
 K_correlation_kernel_std=Vstar_dot_V_gpu(BS_n_g_W,O2,O_times_W)*4._SP*pi
 !
end function K_correlation_kernel_std
