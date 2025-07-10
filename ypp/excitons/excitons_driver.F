!
! License-Identifier: GPL
!
! Copyright (C) 2015 The Yambo Team
!
! Authors (see AUTHORS file for details): AM DV DS
!
! headers
!
#include<y_memory.h>
!
subroutine excitons_driver(k,Xk,en,Xen,q)
 !
 use pars,          ONLY:SP,pi,schlen
 use parser_m,      ONLY:parser
 use R_lattice,     ONLY:bz_samp
 use stderr,        ONLY:intc
 use electrons,     ONLY:levels,spin
 use drivers,       ONLY:l_exc_rad_life
 use electrons,     ONLY:levels,spin
 use YPP_interfaces,ONLY:excitons_sort_and_report
 use YPPm,          ONLY:l_sort,l_exc_wf,l_spin,l_exc_symm,l_amplitude,EXCITONS_user_indexes,BSiq, &
&                        BS_R_right,BS_E,BS_E_SOC_corr,l_interp,l_dipoles
 use BS_solvers,    ONLY:BSS_n_eig
#if defined _YPP_ELPH
 use YPP_ELPH,      ONLY:l_eliashberg,l_gkkp,l_ph_ass_dos
#endif
#if defined _YPP_RT
 use YPP_real_time, ONLY:l_RT_abs
#endif
 use y_memory_alloc
 !
 implicit none
 !
 type(bz_samp) ::Xk,k,q
 type(levels)  ::Xen,en
 !
 ! Work Space
 !
 logical              :: l_flag
 integer ,allocatable :: BS_E_degs(:)
 !... I/0
 integer              :: iq,io_err
 !
 l_flag=l_sort .or. l_exc_wf .or.l_amplitude .or. l_spin .or. l_exc_symm .or. l_exc_rad_life .or. l_interp .or. l_dipoles
#if defined _YPP_ELPH
 l_flag=l_flag .or. l_eliashberg .or. l_gkkp .or. l_ph_ass_dos
#endif
#if defined _YPP_RT
 l_flag=l_flag .or. l_RT_abs
#endif
 !
 if(.not.l_flag) return
 !
 if(.not.l_interp) call section('*','Excitonic Properties @ Q-index #'//trim(intc(BSiq)))
 if(     l_interp) call section('*','Excitonic Properties')
 !
 iq = BSiq
 !
 call                excitons_read(k,Xk,en,Xen,BSiq,"check",.true.,io_err)
 !
 if (io_err==0) call excitons_read(k,Xk,en,Xen,BSiq,"eigenvalues",.true.,io_err)
 !
 if (io_err/=0) then
   call warning("BSE databases not found.")
#if defined _YPP_RT
   ! IP Transient absorption 
   if (l_RT_abs)   call RT_transient_absorption(Xen,Xk,q)
#endif
   return
 endif
 !
 ! Sort energies and write to output
 !===================================
 if (l_sort) call excitons_sort_and_report(iq,BS_R_right,BS_E,BS_E_SOC_corr=BS_E_SOC_corr)
 !
 call parser("States",l_flag)
 !
 if (l_flag) then
   !
   ! Sort energies to find degenerate states
   !========================================
   YAMBO_ALLOC(BS_E_degs,(BSS_n_eig))
   call excitons_find_degeneracies(BS_E,BS_E_degs)
   !
   ! Define the USER set of excitonic states
   !=========================================
   if(.not.l_sort) call excitons_get_user_states(BS_E_degs)
   !
 endif
 !
 ! Exciton dispersion interpolation
 ! =======================================
 if(l_interp) then
   call excitons_bands(k,Xk,en,Xen,q)
   return
 endif
 !
 ! Phonon assisted exciton density of states
 ! =========================================
#if defined _YPP_ELPH
 if(l_ph_ass_dos) then
   call excitons_ph_ass_dos(k,Xk,en,Xen,q)
   return
 endif
#endif
 !
 if(.not.l_sort)   call excitons_read(k,Xk,en,Xen,BSiq,"eigenvectors",.true.,io_err)
 ! 
 ! Exciton's symmetries
 !=====================
 if (l_exc_symm)  call excitons_symmetries(iq,Xen,Xk,q,BS_E)
 ! 
 ! Exciton's spin
 !=====================
 if (l_spin)       call excitons_spin(Xk,BS_R_right,BS_E,BS_E_degs)
 !
 ! Exciton's Amplitude
 !=====================
 if (l_amplitude.or.l_dipoles)  call excitons_amplitudes(Xk,Xen,q,BS_E_degs,iq)
 !
 ! Exciton's lifetimes from BSE
 !=======================================
 if (l_exc_rad_life) call excitons_lifetimes(Xen,Xk,q,BS_R_right,BS_E,BS_E_degs)
 !
 ! Exciton's WFs
 !=====================
 if (l_exc_wf)     call exciton_WFs(Xk,BS_E_degs,iq)
 !
 ! Excitonic Eliashberg function 
 !================================
#if defined _YPP_ELPH
 !
 if (l_eliashberg) call ELPH_general_gFsq(k,en,Xk,Xen,q,BS_E_degs)
 !
 if (l_gkkp)       call ELPH_excitonic_gkkp(Xk)
 !
#endif
 !
#if defined _YPP_RT
 !
 ! Transient absorption 
 !======================
 if (l_RT_abs)     call RT_transient_absorption(Xen,Xk,q)
 !
#endif
 !
 ! CLEAN
 !
 call excitons_read(k,Xk,en,Xen,BSiq,"clean",.true.,io_err)
 !
 if((.not.l_sort).and.l_flag) deallocate(EXCITONS_user_indexes)
 !
end subroutine
