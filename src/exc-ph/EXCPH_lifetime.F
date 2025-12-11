!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM FP AC
!
! headers
!
#include<y_memory.h>
!
!> @brief Calculate exciton non-radiative life-time induced by phonon-scattering
!!
subroutine EXCPH_lifetime(k,q,X)
 !
 use units,          ONLY:HA2EV,HA2CMm1
 use pars,           ONLY:SP,cZERO,rZERO,schlen,pi
 use electrons,      ONLY:levels
 use R_lattice,      ONLY:bz_samp,bz_map
 use parser_m,       ONLY:parser
 use com,            ONLY:msg,of_open_close
 use X_m,            ONLY:X_t
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:manage_action,OP_RD,REP,RD_CL_IF_END,DUMP
 use EXCPH,          ONLY:EXCPH_sum,EXCPH_states,BS_E,BS_all_E,BS_E_in,L_damping,n_exc_in, &
&                         EXC_LifeTime,min_EXC_E,EXCPH_Gkkp_sq,DB_Q_map,l_no_matrix_elements, &
&                         RESIDUALS_and_EIGENVALUES_alloc,RESIDUALS_and_EIGENVALUES_free,alphaQ, &
&                         l_EXCPH_offdiago,EXCPH_Gkkp,EXC_LineWidth,l_interp_G2 
 use ELPH,           ONLY:use_PH_DbGd,elph_branches,FAN_deltaE_treshold,Q_io_map
 use BS_solvers,     ONLY:deg_exc_thrshld
 use vec_operate,    ONLY:degeneration_finder,iku_v_norm,iku_v_2D_norm,degeneration_average
 use frequency,      ONLY:w_samp
 use stderr,         ONLY:intc
 use ELPH,           ONLY:ph_modes
 use parallel_m,     ONLY:PP_indexes,myid,PP_indexes_reset
 use parallel_int,   ONLY:PP_redux_wait,PARALLEL_index
 use LIVE_t,         ONLY:live_timing
 !
 use y_memory_alloc
 !
 implicit none
 !
 type(levels)  ::E
 type(bz_samp) ::k,q
 type(X_t)     ::X
 !
 ! Work space
 !
 integer  :: ID_exc_gkkp,ID_gkkp,io_err,IO_ACT
 integer, external :: io_EXCPH_gkkp,io_ELPH
 integer  :: iq,iq_s,iq_bz,i_alpha,j_alpha,i_beta,il,i_star,iq_db,igrp
 integer  :: elph_branches_save(2)
 integer  :: min_pos(1),first_exc_grp(n_exc_in,q%nibz),n_of_exc_grp(n_exc_in,q%nibz),n_deg_grp(q%nibz)
 logical  :: l_keep_deg,l_extend_output
 real(SP) :: ExcPhM,EXC_life_average,EXC_line_average,EXC_energy_average,Pole
 complex(SP) :: ExcPhW
 real(SP) , allocatable :: EXC_life_of_q(:,:),EXC_pole_of_q(:,:)
 type(PP_indexes)::px
 type(bz_map) :: q_map
 real(SP), allocatable :: EXCPH_Gkkp_sq_pm(:,:,:,:)
 !
 ! Initialization
 call parser('NoMatrxEl',l_no_matrix_elements)
 call parser('KeepDeg',  l_keep_deg)
 call parser('ExtendOutput', l_extend_output)
 call parser('ExcPhOffDiago', l_EXCPH_offdiago)
 !
 call section('*','EXCPH-LifeTime')
 !
 call k_build_up_BZ_tables(k)
 call k_build_up_BZ_tables(q)
 !
 call k_ibz2bz(q,"i",.false.)
 !
 ! Load excitonic energies 
 !
 call EXCPH_load_L(1,X,'check','Lout')
 !
 ! Allocations 
 !
 call RESIDUALS_and_EIGENVALUES_alloc(q)

 !
 ! Load excitonic energies 
 !
 call EXCPH_load_L(1,X,'eigenvalues','Lin')
 BS_E_in(:)=BS_E(EXCPH_states(1):EXCPH_states(2)) ! L_in can be different from L_out 
 !
 do iq=1,q%nibz
   call EXCPH_load_L(iq,X,'eigenvalues','Lout')
   BS_all_E(:,iq)=real(BS_E(:)) ! I suppose the BSE eigenvalues are real
 enddo
 !
 ! Apply 2D stretching if requested
 if(alphaQ/=rZERO) then
   call msg('rs','Excitonic bands 2D stretching :',alphaQ)
   do iq=1,q%nibz
     BS_all_E(:,iq)=BS_all_E(:,iq)+alphaQ*iku_v_2D_norm(q%pt(iq,:3))
   enddo
 endif
 !
 ! Here we find the degenerations of the first n_exc_in values of BS_all_E at all q-points
 ! and we set their energies to exactly the same value (average). This is to
 ! improve the quality of the interpolation in the case of double grid 
 ! (try to avoid unphysical splitting of degenerate states) 
 !
 do iq=1,q%nibz
   call degeneration_finder(n_exc_in,first_exc_grp(:,iq),n_of_exc_grp(:,iq),n_deg_grp(iq),&
   &                      Er=BS_all_E(EXCPH_states(1):EXCPH_states(2),iq),&
   &                 deg_accuracy=deg_exc_thrshld,Include_single_values=.true.)
   !
   call degeneration_average(n_exc_in,first_exc_grp(:,iq),n_of_exc_grp(:,iq),n_deg_grp(iq),&
   &                      Er=BS_all_E(EXCPH_states(1):EXCPH_states(2),iq))
 enddo
 !
 ! Load the Exciton-phonon elements and phonon energies
 !
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1/),MODE=DUMP,ID=ID_exc_gkkp)
 io_err=io_EXCPH_gkkp(ID_exc_gkkp)
 if(io_err/=0) call error('excph_gkkp databases not found!')
 !
 elph_branches_save = elph_branches
 call io_control(ACTION=OP_RD,COM=REP,SEC=(/1/),MODE=DUMP,ID=ID_gkkp)
 io_err=io_ELPH(ID_gkkp,'gkkp_expanded no_matrix_elements')
 elph_branches = elph_branches_save
 if(io_err/=0) call error('gkkp_expanded databases not found!')
 !
 ! Info on excitons
 ! 
 min_pos  =minloc(BS_all_E(1,:))
 min_EXC_E=minval(BS_all_E(1,:))
 !
 call msg('rs','Minimum excitonic energy : ',min_EXC_E*HA2EV,' eV, at q-point :'//intc(min_pos(1)))
 call msg('rs','FAN treshold: ',FAN_deltaE_treshold*HA2CMm1,' cm-1')
 call msg('rs','Self-energy broadening : ',L_damping*HA2EV*1000._SP,'meV')
 if (l_EXCPH_offdiago) then
   call msg('rs','Interband transitions (off-diagonal self-energy) included.')
 endif
 !
 ! Turn on double-grid if presents
 call EXCPH_double_grid_setup(q)
 !
 ! Calculate the life-times
 !
 if(l_interp_G2) then
   !
   if(.not.use_PH_DbGd) &
&     call error("G2 interpolation possible only with double-grid")
   if(l_EXCPH_offdiago) &
&     call error("G2 interpolation not compatible with EXCPH off-diagoanl")
   !
   ! Map q-points on a regular grid for linear interpolation of |G^2|
   !
   call k_map_nearest_by_dir(q,q_map)
   !
 endif
 !
 ! Allocation
 call EXCPH_Ltime_alloc()
 !
 call PP_indexes_reset(px)
 call PARALLEL_index(px,(/q%nbz/))
 !
 call live_timing('Lifetimes',px%n_of_elements(myid+1))
 !
 do iq=1,q%nibz
 do i_star=1,q%nstar(iq)
   !
   iq_s=q%star(iq,i_star)
   iq_bz = q%k_table(iq, iq_s)
   !
   if (.not.px%element_1D(iq_bz)) cycle
   !
   if(l_interp_G2) call eval_G2_grad(iq_bz)
   !
   call load_EXCPH(iq_bz)
   !
   do i_alpha=EXCPH_states(1),EXCPH_states(2)
     !
     do i_beta=EXCPH_sum(1),EXCPH_sum(2)
       !
       do il=elph_branches(1),elph_branches(2)
         !
         ! Intraband terms only (lifetimes)
         ExcPhM=EXCPH_Gkkp_sq(il,i_beta,i_alpha)
         ! Including interband terms (linewidths)
         if (l_EXCPH_offdiago) then
           !
           ExcPhW=cZERO
           do j_alpha=EXCPH_states(1),EXCPH_states(2)
             !
             ExcPhW=ExcPhW+EXCPH_Gkkp(il,i_beta,j_alpha)*conjg(EXCPH_Gkkp(il,i_beta,i_alpha))
             !
           enddo
           !
         endif
         !
         if(l_no_matrix_elements) ExcPhM=1._SP 
         !
         if(.not.use_PH_DbGd) Pole=Simple_Pole()
         if(use_PH_DbGd)      Pole=DbGrid_Pole(l_interp_G2)
         !
         if(l_interp_G2) then
           !      
           ! If we interpoalte the numberator it will automatically
           ! included in the Pole variable
           !
           EXC_LifeTime(i_alpha)=EXC_LifeTime(i_alpha)+Pole
         else
           EXC_LifeTime(i_alpha)=EXC_LifeTime(i_alpha)+ExcPhM*Pole
         endif
         if (l_EXCPH_offdiago) then
           EXC_LineWidth(i_alpha)=EXC_LineWidth(i_alpha)+REAL(ExcPhW,SP)*Pole
         endif
         if (l_extend_output) then
           EXC_life_of_q(i_alpha,iq_bz)=EXC_life_of_q(i_alpha,iq_bz)+ExcPhM*Pole
           EXC_pole_of_q(i_alpha,iq_bz)=EXC_pole_of_q(i_alpha,iq_bz)+Pole
         endif
         ! DS Debug <
         ! write(100+i_alpha,*) i_alpha,EXC_LifeTime(i_alpha)
         ! DS Debug >
         !
       enddo ! i_l
       !
     enddo ! i_beta
     !
   enddo ! i_alpha
   !
   call live_timing(steps=1)
   !
 enddo ! i_start
 enddo ! iq
 !
 call live_timing( )
 !
 call PP_redux_wait(EXC_LifeTime)
 if (l_EXCPH_offdiago) then
   call PP_redux_wait(EXC_LineWidth)
 endif
 if (l_extend_output) then
   call PP_redux_wait(EXC_life_of_q)
   call PP_redux_wait(EXC_pole_of_q)
 endif
 !
 call PP_indexes_reset(px)
 !
 ! Average lifetimes over degenerate groups (use input energies BS_E_in)
 !
 call degeneration_finder(n_exc_in,first_exc_grp(:,1),n_of_exc_grp(:,1),n_deg_grp(1),&
 &                    Er=real(BS_E_in(EXCPH_states(1):EXCPH_states(2))),&
 &                 deg_accuracy=0.001/HA2EV,Include_single_values=.true.)
 !
 
 call degeneration_average(n_exc_in,first_exc_grp(:,1),n_of_exc_grp(:,1),n_deg_grp(1),&
 &                        Er=EXC_LifeTime(EXCPH_states(1):EXCPH_states(2)))
 if (l_EXCPH_offdiago)  call degeneration_average(n_exc_in,first_exc_grp(:,1),n_of_exc_grp(:,1),n_deg_grp(1), &
 &                        Er=EXC_LineWidth(EXCPH_states(1):EXCPH_states(2)))
 !
 ! Convert to final units          
 EXC_LifeTime=EXC_LifeTime*2._SP*pi/real(q%nbz,SP)
 if (l_EXCPH_offdiago) EXC_LineWidth=EXC_LineWidth*2._SP*pi/real(q%nbz,SP)
 !
 ! Print output and write on file
 !
 call print_and_write()
 !
 if (l_extend_output) then
   call print_LifeTimes_of_q("L")
   call print_LifeTimes_of_q("P")
 endif
 !
 ! Deallocations
 call EXCPH_Ltime_free()
 !
 contains
  !
  subroutine EXCPH_Ltime_alloc()
    implicit none
    !
    if (l_extend_output) then
      allocate(EXC_life_of_q(EXCPH_states(1):EXCPH_states(2),q%nbz))
      allocate(EXC_pole_of_q(EXCPH_states(1):EXCPH_states(2),q%nbz))
      EXC_life_of_q=rZERO
      EXC_pole_of_q=rZERO
    endif
    if(l_interp_G2) then
      YAMBO_ALLOC(EXCPH_Gkkp_sq_pm,(ph_modes,EXCPH_sum(1):EXCPH_sum(2),EXCPH_states(1):EXCPH_states(2),6))
    endif
    !
    YAMBO_ALLOC(EXC_LifeTime,(EXCPH_states(1):EXCPH_states(2)))
    EXC_LifeTime=rZERO
    if (l_EXCPH_offdiago) then
      YAMBO_ALLOC(EXC_LineWidth,(EXCPH_states(1):EXCPH_states(2)))
      EXC_LineWidth=rZERO
    endif
    !  
  end subroutine EXCPH_Ltime_alloc
  !
  subroutine EXCPH_Ltime_free()
    implicit none
    !
    YAMBO_FREE(EXC_LifeTime)
    if (l_EXCPH_offdiago) then
      YAMBO_FREE(EXC_LineWidth)
    endif
    call RESIDUALS_and_EIGENVALUES_free()
    if (l_extend_output)  then
       deallocate(EXC_life_of_q)
       deallocate(EXC_pole_of_q)
    endif
    if(l_interp_G2) then
      YAMBO_FREE(EXCPH_Gkkp_sq_pm)
    endif
    ! 
  end subroutine EXCPH_Ltime_free
  !
  subroutine eval_G2_grad(iq_bz_in)
    use pars,           ONLY:rZERO
    use vec_operate,    ONLY:k_periodic_idx
    implicit none
    integer, intent(in) :: iq_bz_in
    !
    ! Working space
    !
    integer  :: id_red,iq_plus,iq_minus,idx(3),istep
    real(SP) :: delta_q
    !
    EXCPH_Gkkp_sq_pm=rZERO
    !
    do id_red=1,3  ! loop on reciprocal lattice
      !
      if(q_map%max_kdir(id_red)==1) then ! non-periodic directions
        EXCPH_Gkkp_sq_pm(:,:,:,id_red)  =rZERO
        EXCPH_Gkkp_sq_pm(:,:,:,id_red+3)=rZERO
	cycle
      endif
      !
      ! Generate q_minus=Q-\delta_q
      !
      idx=q_map%k_map_inv(iq_bz_in,:)
      idx(id_red)=idx(id_red)-istep
      idx=k_periodic_idx(idx,q_map)
      iq_minus=q_map%k_map_dir(idx(1),idx(2),idx(3))
      !
      call load_EXCPH(iq_minus)
      !
      EXCPH_Gkkp_sq_pm(:,:,:,id_red)=EXCPH_Gkkp_sq
      !
      ! Generate q_plus=Q+\delta_q
      !
      idx=q_map%k_map_inv(iq_bz_in,:)
      idx(id_red)=idx(id_red)+istep
      idx=k_periodic_idx(idx,q_map)
      iq_plus=q_map%k_map_dir(idx(1),idx(2),idx(3))
      !
      call load_EXCPH(iq_plus)
      !
      EXCPH_Gkkp_sq_pm(:,:,:,id_red+3)=EXCPH_Gkkp_sq
      !
    enddo ! loop on id_red 
    !
    return
    !  
  end subroutine eval_G2_grad
  !
  subroutine load_EXCPH(iq_bz_in)
    implicit none
    integer, intent(in) :: iq_bz_in
    !
    IO_ACT=manage_action(RD_CL_IF_END,iq_bz_in,1,q%nbz,PAR_INDEX=px)
    !
    call io_control(ACTION=IO_ACT,COM=REP,SEC=(/1+iq_bz_in/),MODE=DUMP,ID=ID_exc_gkkp)
    io_err=io_EXCPH_gkkp(ID_exc_gkkp)
    if(io_err/=0) call error('excph_gkkp databases not found!')
    !
    iq_db=Q_io_map(iq_bz_in)
    call io_control(ACTION=IO_ACT,COM=REP,SEC=(/1+iq_db/),MODE=DUMP,ID=ID_gkkp)
    io_err=io_ELPH(ID_gkkp,'gkkp_expanded no_matrix_elements')
    if(io_err/=0) call error('gkkp_expanded databases not found!')
    !
  end subroutine load_EXCPH 
  !
  function Simple_pole()
    use EXCPH,          ONLY:min_EXC_E,DB_Q_map
    use ELPH,           ONLY:FAN_deltaE_treshold,PH_freqs_sq,Q_io_map
    use functions,      ONLY:bose_f,boltzman_f,Lorentzian_func
    use units,          ONLY:HA2EV
    implicit none
    !
    real(SP) :: Simple_pole
    real(SP) :: E_alpha,E_beta
    real(SP) :: ph_E,N_bose,F_boltz
    integer  :: iq_db
    !
    Simple_pole=rZERO
    !
    E_alpha=real(BS_E_in(i_alpha),SP)
    E_beta =BS_all_E(i_beta,iq)
    !
    iq_db=Q_io_map(iq_bz)
    ph_E=sqrt(abs(PH_freqs_sq(iq_db,il)))
    !
    if(ph_E<FAN_deltaE_treshold) return
    !
    ! Remove degenerate states E(alpha,iq_in)=E(beta,iq_in) [currently iq_in=1]
    !
    if(iq_db==1.and.abs(E_alpha-E_beta)<FAN_deltaE_treshold .and. .not.l_keep_deg) return
    ! even if keeping deg states, still remove the one where i_alpha=i_beta
    if(iq_db==1.and.i_alpha==i_beta .and. l_keep_deg) return
    ! [FP] This removes resonances, I think it is not correct
    !if(abs(E_alpha-E_beta-ph_E)<FAN_deltaE_treshold.or.abs(E_alpha-E_beta+ph_E)<FAN_deltaE_treshold) return
    !
    N_bose  =bose_f(ph_E)
    !
    F_boltz=boltzman_f(E_beta-min_EXC_E)  ! Bose function? Boltzman function?
    !
    Simple_pole=(1._SP+N_bose+F_boltz)*Lorentzian_func(E_alpha-E_beta-ph_E,L_damping) &
&                    +(N_bose-F_boltz)*Lorentzian_func(E_alpha-E_beta+ph_E,L_damping)
    !
    ! DEBUG <
    ! write(200+i_alpha,*) E_alpha*HA2EV,E_beta*HA2EV,ph_E*HA2EV
    ! DEBUG >
    !
  end function Simple_Pole
  !
  function DbGrid_Pole(l_interp_EXCPH)
    use EXCPH,          ONLY:min_EXC_E_DbGd,EXC_E_DbGd
    use ELPH,           ONLY:FAN_deltaE_treshold,PH_freqs,PH_freqs_sq,Q_io_map
    use functions,      ONLY:bose_f,boltzman_f,Lorentzian_func
    use vec_operate,    ONLY:c2a
    implicit none
    !
    logical, intent(in) :: l_interp_EXCPH
    real(SP) :: DbGrid_Pole
    real(SP) :: E_alpha,E_beta
    real(SP) :: ph_E,N_bose,F_boltz
    integer  :: nq_around,iq_fg,iE_fg,id_red
    real(SP) :: weight_dbgd
    real(SP), external :: get_ph_E_DbGd
    real(SP) :: q_fg(3),q_cg(3),G2_var,delta_q(3)
    !
    DbGrid_Pole=rZERO
    !
    E_alpha=real(BS_E_in(i_alpha),SP)
    !
    call c2a(v_in=q%ptbz(iq_bz,:),v_out=q_cg,mode='ki2a') ! To be converted in crystal 
    !
    nq_around  =q%FGbz%k_range(iq_bz,2)-q%FGbz%k_range(iq_bz,1) 
    weight_dbgd=1._SP/real(nq_around+1._SP,SP)                  
    !
    do iq_fg=q%FGbz%k_range(iq_bz,1),q%FGbz%k_range(iq_bz,2)
      !
      !
      if(iq_fg==q%FGbz%k_range(iq_bz,1)) then
        !
        ! Here we get original points for coarse grids
        ! (partially counters inaccuracies in the interpolation in the exc.
        ! case)
        !
        iq_db=Q_io_map(iq_bz) ! This mapping is because sometimes the phonon q_bz differs from the exciton one in yambo
        ph_E=sqrt(abs(PH_freqs_sq(iq_db,il)))
        !
        E_beta=BS_all_E(i_beta,iq)  ! Getting this in the IBZ for the equivalent iq to iq_bz
        !
        q_fg =q%FGbz%pt(iq_fg,:)                 ! q-point for the G2 interpolation
        !
        ! If there is the double-grid I read the gamma phonons from there
        ! because they include the LO-TO splitting and the acustic sum rule
        !
        if(iq_bz==1) ph_E=PH_freqs%FG%E(il,1,1)
        !
      else
        iE_fg=q%FGbz%E_map(iq_fg)                ! Map w(q) from BZ to IBZ
        !
        call c2a(v_in=q%FGbz%pt(iq_fg,:),v_out=q_fg,mode='ki2a')      ! q-point for the G2 interpolation
        !
        ph_E =PH_freqs%FG%E(il,iE_fg,1)          ! w(q_db) phonon energy in the D-Grid
        !
        E_beta=EXC_E_DbGd(i_beta,iq_fg)          ! E(q_bz) exciton energy in the D-Grid
      endif
      !
      if(ph_E<FAN_deltaE_treshold) cycle
      !
      ! Remove degenerate states E(alpha,iq_in)=E(beta,iq_fg(iq_in)) [currently iq_in=1]
      !
      if(iq==1.and.abs(E_alpha-E_beta)<FAN_deltaE_treshold .and. .not.l_keep_deg) return
      ! even if keeping deg states, still remove the one where i_alpha=i_beta
      if(iq==1.and.i_alpha==i_beta .and. l_keep_deg) return
      ! [FP] This removes resonances, I think it is not correct
      !if(abs(E_alpha-E_beta-ph_E)<FAN_deltaE_treshold.or.abs(E_alpha-E_beta+ph_E)<FAN_deltaE_treshold) return
      !
      N_bose  =bose_f(ph_E)
      F_boltz=boltzman_f(E_beta-min_EXC_E_DbGd)  ! Bose function? Boltzman function?
      !
      if(l_interp_EXCPH) then
        !      
        ! EXCPH_gkkp(q_fg)=EXCPH_gkkp(q(cg)+\Delta_q(:)*\Grad EXCPH_gkkp(:)
	!
	delta_q=q_fg-q_cg
	delta_q=delta_q-nint(delta_q) ! Bring in the BZ between [-0.5:0.5]
        !
	G2_var=EXCPH_Gkkp_sq(il,i_beta,i_alpha)
	do id_red=1,3
          if(q_map%max_kdir(id_red)==1) cycle ! non-periodic directions
	  if(delta_q(id_red)<0._SP) then 
	    ! Left interpolation
     	    G2_var=G2_var+delta_q(id_red)*(EXCPH_Gkkp_sq_pm(il,i_beta,i_alpha,id_red) & 
&                                         -EXCPH_Gkkp_sq(il,i_beta,i_alpha))/q_map%q_step(id_red)
	  else
	    ! Right interpolation
     	    G2_var=G2_var+delta_q(id_red)*(EXCPH_Gkkp_sq_pm(il,i_beta,i_alpha,id_red+3) &
&                                         -EXCPH_Gkkp_sq(il,i_beta,i_alpha))/q_map%q_step(id_red)
	  endif
	enddo
        !
      else
        G2_var=1._SP      
      endif
      !
      DbGrid_Pole=DbGrid_Pole+G2_var*(1._SP+N_bose+F_boltz)*Lorentzian_func(E_alpha-E_beta-ph_E,L_damping) &
&                                   +(N_bose-F_boltz)*Lorentzian_func(E_alpha-E_beta+ph_E,L_damping)
      !
    enddo
    !
    DbGrid_pole=DbGrid_pole*weight_dbgd
    !
  end function DbGrid_pole
  !
  subroutine print_and_write()
    use pars,           ONLY:schlen
    use OUTPUT,         ONLY:OUTPUT_driver
    use stderr,         ONLY:real2ch,intc,log2ch
    !
    implicit none
    !
    character(schlen) :: what='QP'
    !
    call OUTPUT_driver(what,action="reset")
    call OUTPUT_driver(what,action="open")
    !
    call OUTPUT_driver(what,ORDER=(/"Eo"/))
    call OUTPUT_driver(what,ORDER=(/"Width [me"/))
    !
    call msg('rs','Exciton lifetimes:')
    !
    call msg('o '//what,'#')    
    call msg('o '//what,'# === Exciton LifeTimes === ')    
    call msg('o '//what,'#')    
    call msg('o '//what,'# Exciton Sum                  : '//trim(intc(EXCPH_sum(1)))//' - '//trim(intc(EXCPH_sum(2))))    
    call msg('o '//what,'# Exciton States               : '//trim(intc(EXCPH_states(1)))//' - '//trim(intc(EXCPH_states(2))))    
    call msg('o '//what,'# Self-Energy broadening       :'//trim(real2ch(L_damping*HA2EV*1000._SP))//' meV')    
    call msg('o '//what,'# FAN threshold                 :'//trim(real2ch(FAN_deltaE_treshold*HA2EV*1000._SP))//' meV')    
    call msg('o '//what,'# PH Double-Grid                : '//log2ch(use_PH_DbGd))    
    !
    do i_alpha=EXCPH_states(1),EXCPH_states(2)
      !
      call OUTPUT_driver(what,TITLES=(/"Band   "/),I_VALUES=(/i_alpha/))
      call OUTPUT_driver(what,TITLES=(/"E"/),R_VALUES=(/real(BS_E_in(i_alpha))/),UNIT="eV")
      call OUTPUT_driver(what,TITLES=(/"Life^-1"/),R_VALUES=(/EXC_LifeTime(i_alpha)/),UNIT="meV")
      if (l_EXCPH_offdiago) then
        call OUTPUT_driver(what,TITLES=(/"Width"/),R_VALUES=(/EXC_LineWidth(i_alpha)/),UNIT="meV")
      endif
     
      call OUTPUT_driver(what,action="write")
      !
      call msg('rs','Exciton ['//trim(intc(i_alpha))//'] Lifetime: ',EXC_LifeTime(i_alpha)*HA2EV*1000._SP,'meV')
      !
    enddo  
    !
    call OUTPUT_driver(what,action="close")
    !
  end subroutine 
  !
  subroutine print_LifeTimes_of_q(MODE)
    ! MODE="L"[ife]:
    !   - print q-map of lifetimes for each exc_in state for q-space plot
    !     We print l(q)/sum_q l(q)
    ! MODE="P"[ole]:
    !   - Print value of the Pole [F_occ*delta(...)] for each exc_in state for q-space plot
    !     We print sum_{beta,mu} Pole_{beta,mu}(q)
    !
    use vec_operate,   ONLY:iku_v_norm
    use com,           ONLY:of_open_close
    !use YPPm,          ONLY:K_transform
    !
    implicit none
    !
    character(1) :: MODE
    character(schlen), allocatable :: headings(:)
    character(schlen)              :: file_name
    real(SP), allocatable          :: values(:)
    real(SP)                       :: tmp_q(3)
    !
    allocate(headings(n_exc_in+4))
    YAMBO_ALLOC(values,(n_exc_in+4))
    values=rZERO
    !
    headings(1)=" |q| (a.u.) "
    do i_alpha=EXCPH_states(1),EXCPH_states(2)
      headings(i_alpha+1)=' e'//intc(i_alpha)
      !
      if (MODE=="L") then
        EXC_life_of_q(i_alpha,:)=EXC_life_of_q(i_alpha,:)/SUM(EXC_life_of_q(i_alpha,:))
      endif
      !
    enddo
    headings(n_exc_in+2:n_exc_in+4)=&
&           (/"q_x (rlu)","q_y (rlu)","q_z (rlu)"/)
    !
    if (MODE=="L") file_name="LifeTimes_of_q"
    if (MODE=="P") file_name="Resonance_of_q"
    !
    call of_open_close(trim(file_name),'ot')
    call msg('o '//trim(file_name),"#",headings,INDENT=0,USE_TABS=.true.)
    call msg('o '//trim(file_name),"#")
    !
    do iq=1,q%nibz
    do i_star=1,q%nstar(iq)
      !
      iq_s=q%star(iq,i_star)
      iq_bz = q%k_table(iq, iq_s)
      !
      if (iq_bz>1) values(1)=iku_v_norm(q%ptbz(iq_bz,:3))
      !
      tmp_q = q%ptbz(iq_bz,:3)
      !
      !call K_transform(tmp_q,'iku')
      if (MODE=="L") values(1+1:n_exc_in+1)=EXC_life_of_q(:,iq_bz)
      if (MODE=="P") values(1+1:n_exc_in+1)=EXC_pole_of_q(:,iq_bz)
      values(n_exc_in+2:n_exc_in+4)=tmp_q
      call msg('o '//trim(file_name),' ',values,INDENT=0,USE_TABS=.true.)
      !
    enddo
    enddo
    call of_open_close(trim(file_name))
    !
  end subroutine print_LifeTimes_of_q
  !
end subroutine
