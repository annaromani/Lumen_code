!
! License-Identifier: GPL
!
! Copyright (C) 2020 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
subroutine ELPH_gsqF_to_Sigma(en,q,l_GF_from_CA,Sc_W)
 !
 use pars,          ONLY:SP,DP,cZERO,pi
 use functions,     ONLY:bose_f
 use units,         ONLY:HA2EV
 use LIVE_t,        ONLY:live_timing
 use R_lattice,     ONLY:nqbz
 use ELPH,          ONLY:gsqF_energy_steps,gsqF_energy_steps,gsqF_ca_corr,gsqF_fan,gsqF_life_f,&
&                        gsqF_life_bose,gsqF_dw,elph_branches,elph_nQ_used,&
&                        elph_use_q_grid,PH_freqs_sq,QP_DW,QP_OMS_Fan,elph_use_q_grid
 use QP_m,          ONLY:QP_solver,QP_Sc_steps,QP_Sc_steps,QP_n_states,QP_Sc,On_Mass_Shell_approx,&
&                        QP_table
 use parallel_m,    ONLY:PAR_nQ_bz,PAR_nQP,PAR_IND_QP,PAR_Q_bz_index,PAR_IND_Q_bz,PAR_COM_G_b_INDEX
 use parallel_int,  ONLY:PP_redux_wait
 use electrons,     ONLY:levels,spin_occ
 use R_lattice,     ONLY:bz_samp
 use frequency,     ONLY:w_samp,W_reset
 use OUTPUT,        ONLY:OUTPUT_driver
 use com,           ONLY:msg
 !
 implicit none
 !
 type(levels)  ::en
 type(bz_samp) ::q
 logical       ::l_GF_from_CA
 type(w_samp)  ::Sc_W(QP_n_states)
 !
 ! Work Space
 !
 integer            ::i_qp,iw,il,iq_loop,iq_mem,ib,ik,nq
 type(w_samp)       ::Sc_local_W
 complex(DP)        ::DP_FAN(gsqF_energy_steps),DP_DW,coefficient
 complex(SP)        ::Sc_local(gsqF_energy_steps)
 real(SP)           ::ph_freq
 complex(DP)        ::DP_FAN_mode(elph_branches(1):elph_branches(2))
 complex(DP)        ::DP_DW_mode(elph_branches(1):elph_branches(2)) 
 complex(DP),allocatable ::QP_FAN_mode(:,:),QP_DW_mode(:,:)
 integer            ::imode,n_modes
 character(32)      ::mode_file
 !
 n_modes=elph_branches(2)-elph_branches(1)+1
 mode_file='elph_mode_sigma'
 allocate(QP_FAN_mode(QP_n_states,n_modes))
 allocate(QP_DW_mode (QP_n_states,n_modes))
 QP_FAN_mode=cmplx(0.,0.,kind=DP)
 QP_DW_mode =cmplx(0.,0.,kind=DP)
 !
 call ELPH_simmetrize_and_distribute_gsqF(en)
 !
 nq=nqbz
 if (.not.elph_use_q_grid) nq=elph_nQ_used
 !
 if (trim(QP_solver)=='n') call live_timing('El-Ph Sc   [sum]',PAR_nQ_bz*PAR_nQP)
 if (trim(QP_solver)=='g') call live_timing('El-Ph Sc[sum+KK]',PAR_nQ_bz*PAR_nQP+PAR_nQP)
 !
 ! Setup the pre-coefficient for the lifetimes related gsqF. In case
 ! the QP corrections are required these are stored in the immaginary part of DP_FAN.
 ! If the self-energy has to be plotted, instead, the contribution from the
 ! lifetimes-related gsqF must be real (spectral function).
 !
 coefficient=cmplx(0.,1.,kind=DP)
 !
 ! Setup energy range for Sc/G plotting (I)
 !
 if ( trim(QP_solver)=='g' ) then
   coefficient=cmplx(1.,0.,kind=DP)
   call W_reset(Sc_local_W)
   Sc_local_W%n_freqs =QP_Sc_steps
 endif
 !
 do i_qp=1,QP_n_states   
   !
   if (.not.PAR_IND_QP%element_1D(i_qp)) cycle
   !
   ! Setup energy range for Sc/G plotting (II)
   !
   if ( trim(QP_solver)=='g' ) then
     call FREQUENCIES_Green_Function(i_qp,Sc_local_W,en%E,.FALSE.)
     !
     ! Overwrites the damping defined in FREQUENCIES_Green_Function
     !
     forall (iw=1:QP_Sc_steps) Sc_local_W%p(iw)= real( Sc_local_W%p(iw)) + cmplx(0._SP,0.00001_SP/HA2EV)
   endif
   !
   ib   =QP_table(i_qp,1)
   ik   =QP_table(i_qp,3)
   !
   DP_FAN=cmplx(0.,0.,kind=DP)
   DP_DW =cmplx(0.,0.,kind=DP)
   DP_FAN_mode=cmplx(0.,0.,kind=DP)
   DP_DW_mode =cmplx(0.,0.,kind=DP)
   do iq_loop=1,nq
     !
     if (.not.PAR_IND_Q_bz%element_1D(iq_loop)) cycle
     iq_mem=PAR_Q_bz_index(iq_loop)
     !
     do il=elph_branches(1),elph_branches(2)
       !
       if (.not.elph_use_q_grid) ph_freq=sqrt(PH_freqs_sq(iq_loop,il))
       if (     elph_use_q_grid) ph_freq=sqrt(PH_freqs_sq( q%sstar(iq_loop,1) ,il))
       !
       do iw=1,gsqF_energy_steps
         !
         ! QP lifetimes/Green's function
         !-------------------------------
         if (l_GF_from_CA.and.trim(QP_solver)=='g') then
           DP_FAN(iw)=DP_FAN(iw)+(gsqF_fan(i_qp,iq_mem,il,iw)*(2._DP*bose_f(ph_freq)/spin_occ+1._DP)/2._DP+&
&                                 gsqF_ca_corr(i_qp,iq_mem,il,iw))/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)
         else
           DP_FAN(iw)=DP_FAN(iw)+coefficient*(gsqF_life_bose(i_qp,iq_mem,il,iw)*(bose_f(ph_freq)/spin_occ+1._DP)/2._DP+&
&                                 gsqF_life_f(i_qp,iq_mem,il,iw))/real( PAR_COM_G_b_INDEX%n_CPU,kind=DP)
         endif
         !
         if ( trim(QP_solver)=='n' )  then
           !
           ! Integrate g^2 F function (QP correction).
           !
           ! QP energies: FAN+correction
           !-----------------------------
           ! 
           DP_FAN(iw)=DP_FAN(iw)+(gsqF_fan(i_qp,iq_mem,il,iw)*(2._DP*bose_f(ph_freq)/spin_occ+1._DP)/2._DP+&
&                                 gsqF_ca_corr(i_qp,iq_mem,il,iw))/real( PAR_COM_G_b_INDEX%n_CPU)
           !
         endif
         !
       enddo
       !
       ! Integrate g^2 F function (QP correction & spectral function).
       !
       ! QP energies : DW
       !------------------
       !
       DP_DW=DP_DW+gsqF_dw(i_qp,iq_mem,il)* &
        & (2._DP*bose_f(ph_freq)/spin_occ+1._DP)/ &
        & 2._DP/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)

       DP_DW_mode(il)=DP_DW_mode(il)+gsqF_dw(i_qp,iq_mem,il)* &
        & (2._DP*bose_f(ph_freq)/spin_occ+1._DP)/ &
        & 2._DP/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)

       DP_FAN_mode(il)=DP_FAN_mode(il)+coefficient* &
        & (gsqF_life_bose(i_qp,iq_mem,il,1)*(bose_f(ph_freq)/spin_occ+1._DP)/2._DP+ &
        &  gsqF_life_f(i_qp,iq_mem,il,1))/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)

       DP_FAN_mode(il)=DP_FAN_mode(il)+(gsqF_fan(i_qp,iq_mem,il,1)* &
        & (2._DP*bose_f(ph_freq)/spin_occ+1._DP)/2._DP+ &
        &  gsqF_ca_corr(i_qp,iq_mem,il,1))/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)      
      ! DP_DW=DP_DW+gsqF_dw(i_qp,iq_mem,il)* &
      !  & (2._DP*bose_f(ph_freq)/spin_occ+1._DP)/ &
      !  & 2._DP/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)
      ! DP_DW_mode(il)=DP_DW_mode(il)+gsqF_dw(i_qp,iq_mem,il)* &
      !  & (2._DP*bose_f(ph_freq)/spin_occ+1._DP)/ &
      !  & 2._DP/real(PAR_COM_G_b_INDEX%n_CPU,kind=DP)
      ! DP_FAN_mode(il)=DP_FAN_mode(il)+(gsqF_fan(i_qp,iq_mem,il,1)* &
      !  & (2._DP*bose_f(ph_freq)/spin_occ+1._DP)/2._DP+ &
      !  &  gsqF_ca_corr(i_qp,iq_mem,il,1))/real(PAR_COM_G_b_INDEX%n_CPU)
        ! 
     enddo
     !
     call live_timing(steps=1)
     !
   enddo
   do il=elph_branches(1),elph_branches(2)
     imode=il-elph_branches(1)+1
     QP_FAN_mode(i_qp,imode)=DP_FAN_mode(il)
     QP_DW_mode (i_qp,imode)=DP_DW_mode(il)
   enddo
   !
   if ( trim(QP_solver)=='n' )  then
     !
     QP_DW(i_qp)=QP_DW(i_qp)+DP_DW
     QP_OMS_Fan(i_qp)=QP_OMS_Fan(i_qp)+DP_FAN(1)
     !
     if (On_Mass_Shell_approx) then
       QP_Sc(i_qp,:)=QP_Sc(i_qp,:)+DP_FAN(1)+DP_DW
     else
       QP_Sc(i_qp,:)=QP_Sc(i_qp,:)+DP_FAN(:)+DP_DW
     endif
     !
   else if ( trim(QP_solver)=='g' )  then
     !
     if (.not.l_GF_from_CA) then
       !
       Sc_local=cZERO
       !
       call Kramers_Kronig(cmplx(DP_FAN(:))/pi,real(Sc_W(i_qp)%p),QP_Sc_steps,&
&                                    Sc_local,conjg(Sc_local_W%p),QP_Sc_steps,cZERO)
       !
       QP_Sc(i_qp,:)=QP_Sc(i_qp,:)+Sc_local(:)+DP_DW
       !
     else
       !
       QP_Sc(i_qp,:)=QP_Sc(i_qp,:)+DP_FAN(:)+DP_DW 
       !
     endif
     !
     call live_timing(steps=1)
     !
   endif
   !
 enddo
 !
 call PP_redux_wait(QP_FAN_mode)
 call PP_redux_wait(QP_DW_mode)
 !
 if (trim(QP_solver)=='n') then
   !
   call OUTPUT_driver(mode_file,action="reset")
   call OUTPUT_driver(mode_file,action="open")
   call msg('o '//trim(mode_file),'#',' Mode-resolved electron-phonon self-energy from ELPH_gsqF_to_Sigma',INDENT=0)
   call msg('o '//trim(mode_file),'#',' Values are summed over q; OUTPUT_driver converts Ha to meV.',INDENT=0)
   call msg('o '//trim(mode_file),'#',' Fan value corresponds to DP_FAN_mode(il), i.e. iw=1 in this routine.',INDENT=0)
   call msg('o '//trim(mode_file),'#',' ',INDENT=0)
   !
   do i_qp=1,QP_n_states
     ib=QP_table(i_qp,1)
     ik=QP_table(i_qp,3)
     do il=elph_branches(1),elph_branches(2)
       imode=il-elph_branches(1)+1
       !
       call OUTPUT_driver(mode_file,TITLES=(/"K-point","Band   ","Mode   "/),I_VALUES=(/ik,ib,il/))
       call OUTPUT_driver(mode_file,TITLES=(/"ReFan  ","ReDW   ","ReTot  "/),&
&        R_VALUES=(/real(real(QP_FAN_mode(i_qp,imode)),SP),&
&                   real(real(QP_DW_mode (i_qp,imode)),SP),&
&                   real(real(QP_FAN_mode(i_qp,imode)+QP_DW_mode(i_qp,imode)),SP)/),UNIT="meV")
       call OUTPUT_driver(mode_file,TITLES=(/"ImFan  ","ImDW   ","ImTot  "/),&
&        R_VALUES=(/real(aimag(QP_FAN_mode(i_qp,imode)),SP),&
&                   real(aimag(QP_DW_mode (i_qp,imode)),SP),&
&                   real(aimag(QP_FAN_mode(i_qp,imode)+QP_DW_mode(i_qp,imode)),SP)/),UNIT="meV")
       !
       call OUTPUT_driver(mode_file,ORDER=(/"K-poi","Band ","Mode "/))
       call OUTPUT_driver(mode_file,ORDER=(/"ReFan","ReDW ","ReTot"/))
       call OUTPUT_driver(mode_file,ORDER=(/"ImFan","ImDW ","ImTot"/))
       call OUTPUT_driver(mode_file,action="write")
       !
     enddo
   enddo
   !
   call OUTPUT_driver(mode_file,action="close")
   !
 endif
 !
 call PP_redux_wait(QP_DW)
 call PP_redux_wait(QP_OMS_Fan)
 call PP_redux_wait(QP_Sc)
 !
 call live_timing()
 !
 ! CLEAN
 !
 if (allocated(QP_FAN_mode)) deallocate(QP_FAN_mode)
 if (allocated(QP_DW_mode )) deallocate(QP_DW_mode )
 call W_reset(Sc_local_W)
 !
end subroutine
