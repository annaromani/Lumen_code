!
!        Copyright (C) 2000-2017 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MP DS
!
! headers
!
#include<y_memory.h>
!
! This file is distributed under the terms of the GNU
! General Public License. You can redistribute it and/or
! modify it under the terms of the GNU General Public
! License as published by the Free Software Foundation;
! either version 2, or (at your option) any later version.
!
! This program is distributed in the hope that it will
! be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU General Public
! License along with this program; if not, write to the Free
! Software Foundation, Inc., 59 Temple Place - Suite 330,Boston,
! MA 02111-1307, USA or visit http://www.gnu.org/copyleft/gpl.txt.
!
subroutine excitons_lifetimes(Xen,Xk,Xq,BS_R,BS_E,BS_E_degs)
 !
 use pars,          ONLY:SP,DP,schlen,pi,cZERO
 use units,         ONLY:HA2EV,AU2SEC,SPEED_OF_LIGHT
 use electrons,     ONLY:spin_occ,levels,spin
 use stderr,        ONLY:intc,real2ch
 use BS_solvers,    ONLY:BSS_eh_table,BS_mat,BSS_n_eig
 use BS,            ONLY:BS_bands,BS_H_dim
 use com,           ONLY:msg,of_open_close
 use R_lattice,     ONLY:q0_def_norm,nkbz,bz_samp,q_norm
 use D_lattice,     ONLY:DL_vol,non_periodic_directions,l_1D,l_2D,l_0D,l_3D
 use X_m,           ONLY:eps_2_alpha
 use DIPOLES,       ONLY:DIP_iR
 use YPPm,          ONLY:EXC_temp,EXC_mass,EXC_Phi_Grid,EXC_shift,epsi0,&
&                        DIPs,EXCITONS_user_indexes,EXCITONS_n_user_states
 use vec_operate,   ONLY:v_norm
 use IO_m,          ONLY:io_DIP
 use parallel_int,  ONLY:PARALLEL_global_indexes
 use parallel_m,    ONLY:PAR_DIPk_ibz_index
 !
 ! implicit none 
 !
 use y_memory_alloc
 !
 implicit none
 !
 type(levels)       :: Xen
 type(bz_samp)      :: Xk,Xq
 !
 complex(SP)        :: BS_E(BSS_n_eig), BS_R(BSS_n_eig) 
 integer            :: BS_E_degs(BSS_n_eig)
 !
 ! Work Space
 !
 character(schlen)  :: titles(7), ch_dummy
 !
 integer            :: i_lambda,i_grid,io_err
 integer            :: neh,i_l,i_l_grp,n_lambda_deg,i_v,i_c,i_kbz,i_kibz,i_spin 
 integer            :: Phi_grid_size
 !
 complex(SP)        :: Ir_EXC(3),rotated_iR(3),P_EXC(3)
 !
 real(SP)           :: P_EXC_square,Px_EXC_square,Py_EXC_square,Pz_EXC_square
 real(SP)           :: gamma0(BSS_n_eig)
 real(SP)           :: gammaT_aver, gamma0_grp,f_cv,alpha,tauT_aver,tau0_grp
 !
 !
 real(SP)           :: exc_energy
 real(SP)           :: tau0
 real(SP)           :: num,den,tauT_eff,gammaT_eff
 real(DP)           :: arg
 !
 !HYC add 
 integer :: ik_mem
 complex(SP)        :: iR_EXC_TRY(3)
 real(SP)           :: phi(EXC_Phi_Grid)
 !
 logical            :: is_cut(3)
 !
 call DIPOLE_dimensions(Xen,DIPs,BS_bands,(/0.,0.,0./))
 !
 call PARALLEL_global_indexes(Xen,Xk,Xq," ",RESET=.TRUE.)
 call PARALLEL_global_indexes(Xen,Xk,Xq,"DIPOLES",Dip=DIPs)
 !
 call DIPOLE_IO(Xk,Xen,DIPs,'read ',io_err,'DIP')
 !
 ! check dimensionality of the system
 is_cut(1) = .false.
 is_cut(2) = .false.
 is_cut(3) = .false.
 !
 if (index(non_periodic_directions,'x')+index(non_periodic_directions,'X')/=0) is_cut(1)=.true.
 if (index(non_periodic_directions,'y')+index(non_periodic_directions,'Y')/=0) is_cut(2)=.true.
 if (index(non_periodic_directions,'z')+index(non_periodic_directions,'Z')/=0) is_cut(3)=.true.
 call check_periodic_directions(is_cut)
 !
 call section('=','Computing Excitons radiative lifetimes')
 !
 !
 if ((l_1d) .and. (.not. is_cut(3))) then
  write (6,*) 'STOP - 1D implemented only for periodic direction along z'
  stop
 endif
 if ((l_2d) .and. (.not. is_cut(3))) then
  write (6,*) 'STOP  - 2D implemented only for non-periodic direction along z'
  stop
 endif
 if ((l_3d) .and. (epsi0(1).ne.epsi0(2)))  then
  write (6,*) 'STOP - 3D implemented only for isotropic or anisotropic axiz along z'
  stop
 endif


  titles(1)='E [eV]'
  titles(2)='Tau0(s)'
  titles(3)='tauT_aver(s)'
  titles(4)='tauT_eff(s)'
 !
 do i_lambda=1,EXCITONS_n_user_states
   !
   gammaT_aver =0._SP
   gamma0_grp=0._SP
   tauT_aver =0._SP
   tau0_grp=0._SP
   num = 0._SP
   den = 0._SP
   alpha = 0._SP
   !
   i_l_grp=EXCITONS_user_indexes(i_lambda)
   !
   if(i_lambda/=BS_E_degs(i_l_grp)) cycle
   !
   n_lambda_deg=count(BS_E_degs==BS_E_degs(i_l_grp))
   !
   if (n_lambda_deg>1) call msg('s',':: State '//trim(intc(i_l_grp))//' Merged with states '//&
&                                     trim(intc(BS_E_degs(i_l_grp)))//' -> '//&
&                                     trim(intc(BS_E_degs(i_l_grp)+n_lambda_deg-1)))
   !
   !
   ch_dummy='exc_lifetimes_at_'//trim(intc(i_l_grp))
   call of_open_close(ch_dummy,'ot')
   !
   ! First compute the excitonic dipoles <mu>  starting from indipendent-particle
   ! ones
   !==========================================
   !
   do i_l=BS_E_degs(i_l_grp),BS_E_degs(i_l_grp)+n_lambda_deg-1
     iR_EXC = cZERO
     iR_EXC_TRY= cZERO
     do neh = 1,BS_H_dim !(size BSE matrix, cvk)
       !
       i_kbz   = BSS_eh_table(neh,1)
       i_kibz  = Xk%sstar(i_kbz,1)
       i_v     = BSS_eh_table(neh,2)
       i_c     = BSS_eh_table(neh,3)
       i_spin  = spin(BSS_eh_table(neh,:))
       !
       f_cv    = (Xen%f(i_v,i_kibz,i_spin)- Xen%f(i_c,i_kibz,i_spin))/spin_occ
       !
       call DIPOLE_rotate(i_c,i_v,i_kbz,i_spin,"DIP_iR",Xk,rotated_iR)
       !
       iR_EXC(:)=iR_EXC(:) - rotated_iR(:)*BS_mat(neh,i_l)*sqrt(f_cv)
       !
     enddo
     !
     P_EXC(:) =BS_E(i_l)*iR_EXC(:)
     Px_EXC_square = real(P_EXC(1)*conjg(P_EXC(1)))
     Py_EXC_square = real(P_EXC(2)*conjg(P_EXC(2)))
     Pz_EXC_square = real(P_EXC(3)*conjg(P_EXC(3))) 
!          P_EXC_square =real(BS_E(i_l))**2*(real(iR_EXC(1)*conjg(iR_EXC(1)) + &
!&              iR_EXC(2)*conjg(iR_EXC(2)) + &
!               iR_EXC(3)*conjg(iR_EXC(3)) ))
     ! 
     ! compute the lifetime at q=0  gamma(0) and the prefactors for finite T
     !==========================================
     if (l_0D) then
      P_EXC_square=Px_EXC_square + Py_EXC_square + Pz_EXC_square
      gamma0(i_l) =  4._SP*P_EXC_square*((real(BS_E(i_l))-EXC_shift)/(3._SP*SPEED_OF_LIGHT**3))
     endif
     if (l_1D) then
      P_EXC_square = Pz_EXC_square
      gamma0(i_l) = eps_2_alpha*P_EXC_square*4.*pi/(DL_vol*real(nkbz)*SPEED_OF_LIGHT**2)
      alpha= 4._SP*(real(BS_E(i_l))-EXC_shift)/(3._SP*SPEED_OF_LIGHT*sqrt(2._SP*EXC_mass(3)*EXC_temp*PI))
      write(9,*) exc_temp,exc_energy,P_EXC_square,alpha
     endif
     if (l_2D) then
     P_EXC_square = Px_EXC_square + Py_EXC_square 
      gamma0(i_l) = 4.*pi*P_EXC_square*eps_2_alpha/((DL_VOL*real(BS_E(i_l))-EXC_shift)*SPEED_OF_LIGHT*real(nkbz))
      alpha=4.*(real(BS_E(i_l))-EXC_shift)**2/(3._SP*2._SP*sqrt(EXC_mass(1)*EXC_mass(2))*EXC_temp*SPEED_OF_LIGHT**2) 
     endif 
     if (l_3D) then
      P_EXC_square = Px_EXC_square + Py_EXC_square + Pz_EXC_square
       gamma0(i_l) =  4.*pi*(sqrt(pi*epsi0(1))/DL_vol)/((real(BS_E(i_l))-EXC_shift)**2*real(nkbz))*&
       ((2._SP/3._SP*epsi0(3)/epsi0(1)+2._SP)*(Px_EXC_square+Py_EXC_square)+8./3.*Pz_EXC_square)
      alpha=((real(BS_E(i_l))-EXC_shift)**2/(2._SP*(EXC_mass(1)**2*EXC_mass(3))**(1./3.)&
              *EXC_temp*SPEED_OF_LIGHT**2))**1.5
     endif
     gamma0_grp=gamma0_grp+gamma0(i_l)
     if (l_1D) gammaT_aver=gammaT_aver+gamma0(i_l)*alpha
     if (l_0D) gammaT_aver=gamma0(i_l) 
     if (l_2D) gammaT_aver=gammaT_aver+gamma0(i_l)*alpha
     if (l_3D) gammaT_aver=gammaT_aver+gamma0(i_l)*alpha
     !
     exc_energy = (real(BS_E(i_l)-EXC_shift)*HA2EV)
     tau0_grp= AU2SEC/gamma0_grp
     tauT_aver = AU2SEC/gammaT_aver
     !
     arg =real(BS_E(i_l)-BS_E(1))/EXC_temp
     num = num + dexp(-arg)*gammaT_aver
     den = den + dexp(-arg)
     gammaT_eff = num/den
     tauT_eff = AU2SEC/gammaT_eff
   enddo
   111 format (i4, 6d15.7)
   !
   call msg('o exc_life',"# ","exc. E(eV) = "//trim(real2ch(exc_energy)),INDENT=0)
   call msg('o exc_life',"# ","gamma zero = "//trim(real2ch(gamma0_grp)),INDENT=0)
   call msg('o exc_life',"# ","gammaT aver  = "//trim(real2ch(gammaT_aver)),INDENT=0)
   !
   call msg('o exc_life','#',titles(:4),INDENT=0,USE_TABS=.true.) 
   call msg('o exc_life',' ',(/exc_energy,tau0_grp,tauT_aver,tauT_eff/),INDENT=0,USE_TABS=.true.)
  call of_open_close(ch_dummy)
   !
 enddo
 !
end subroutine excitons_lifetimes
