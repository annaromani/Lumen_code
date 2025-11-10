!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): MG IA
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
subroutine FL_polarization(E,k,FL_V_bands,FL_Chi_pol)
 !
 use pars,           ONLY:DP,SP,cZERO,pi,cONE,cI
 use units,          ONLY:HA2FSm1,HA2EV
 use electrons,      ONLY:levels,n_sp_pol,spin_occ,n_spin
 use fields,         ONLY:Efield,Divide_by_Field
 use D_lattice,      ONLY:a,DL_vol
 use R_lattice,      ONLY:bz_samp,k_map,bz_map,nXkbz
 use vec_operate,    ONLY:get_id_perpendicular
 use electric,       ONLY:l_no_res_symms
 use QP_m,           ONLY:QP_nk
 use nl_optics,      ONLY:NL_bands,NL_initial_P,NL_damping
 use fl_optics,      ONLY:max_fl_mode,n_tot_fl_modes,FL_space_dim,FL_S_plus,FL_S_minus,&
&                         FL_states,FL_Ho_plus_Sigma,FL_time_S_det,FL_Chi_order,&
&                         max_time_step,n_tot_time_steps
 use parallel_m,     ONLY:PAR_IND_Xk_bz,PAR_COM_Xk_ibz_INDEX
 use parallel_int,   ONLY:PP_redux_wait
#if defined _TIMING
 use timing_m,       ONLY:timing
#endif
 !
 implicit none
 !
 type(levels),  intent(in)  :: E
 type(bz_samp), intent(in)  :: k
 complex(SP),   intent(in)  :: FL_V_bands(NL_bands(2),n_tot_fl_modes,maxval(E%nbf),QP_nk,n_sp_pol)
 complex(DP),   intent(inout) :: FL_Chi_pol(3,FL_Chi_order+1)  ! In cartesian coordinates
 !
 ! Work Space
 !
 complex(DP)  :: FL_time_P_Berry_RED(3,n_tot_time_steps)  ! In reduced coordinates
 complex(DP)  :: FL_time_P_Berry_CART(3,n_tot_time_steps) ! In cartesian coordinates
 complex(DP)  :: FL_Pi_pol(3,n_tot_time_steps)
 integer      :: ik_bz,i_sp,i_t,i_order,i1_fl
 integer      :: id_in,id2,id3,Nperpend,i1,i2,i3,idx(3),ikbz,id,is
 complex(DP)  :: zeta
 real(DP)     :: imag_zeta
 real(SP)     :: time,w0_sm1
 character(4) :: filename
 !
 ! Update Floquet Overlaps and determinant
 ! ===================================================================
 !
#if defined _TIMING
 call timing('FL Berry Pol NEQ',OPR='start')
#endif
 !
 FL_S_minus=cZERO
 FL_S_plus =cZERO
 FL_time_S_det =cZERO
 !
 do i_sp=1,n_sp_pol
   !
   do ik_bz=1,nXkbz
     !
     if(.not.PAR_IND_Xk_bz%element_1D(ik_bz)) cycle
     !
     call FL_overlaps(E,k,i_sp,ik_bz,FL_V_bands(:,:,:E%nbf(i_sp),:,i_sp), &
                                  &  FL_S_plus(:E%nbf(i_sp),:E%nbf(i_sp),:,:,ik_bz,i_sp),&
                                  &  FL_S_minus(:E%nbf(i_sp),:E%nbf(i_sp),:,:,ik_bz,i_sp),&
                                  &  FL_time_S_det(:,:,ik_bz,i_sp),.true.)
     !
   enddo
   !
   call PP_redux_wait(FL_S_plus(:,:,:,:,:,i_sp) , COMM=PAR_COM_Xk_ibz_INDEX%COMM)
   call PP_redux_wait(FL_S_minus(:,:,:,:,:,i_sp), COMM=PAR_COM_Xk_ibz_INDEX%COMM)
   call PP_redux_wait(FL_time_S_det(:,:,:,i_sp) , COMM=PAR_COM_Xk_ibz_INDEX%COMM)
   !
 enddo
 !
 FL_time_P_Berry_CART = cZERO
 FL_time_P_Berry_RED  = cZERO
 !
 do i_t=1,n_tot_time_steps
   !
   do i_sp=1,n_sp_pol
     !
     call Berry_polarization_EQ(FL_time_P_Berry_RED(:,i_t),FL_time_S_det(i_t,:,:,i_sp_pol),Xen,Xk,k_map)
     !
   enddo
   !
 enddo
 !
 !
 do id=1,3
   !
   do i_t=1,n_tot_time_steps
     !
     FL_time_P_Berry_CART(id,i_t)=sum(FL_time_P_Berry_RED(:,i_t)&
             &                   *a(:,id))/DL_vol
     !
   enddo
   !
 enddo
 ! 
 ! Apply residual symmetries
 !
 if(.not.l_no_res_symms) then
   do i_t=1,n_tot_time_steps      
     call Apply_symm_vec(FL_time_P_Berry_CART(:,i_t)) 
     !
     FL_time_P_Berry_CART(:,i_t) = FL_time_P_Berry_CART(:,i_t) - NL_initial_P(:3)
   enddo
 endif
 !
 FL_Chi_pol = cZERO
 FL_Pi_pol = cZERO
 !
 do id=1,3
   !
   call FL_FT_polarization(FL_time_P_Berry_CART(id,:),FL_Pi_pol(id,:))
   !
 enddo
 !
 do i_order=0,FL_Chi_order
   !
   i1_fl = i_order + (max_time_step+1)
   !
   if (i_order == 1) then
     !
     FL_Chi_pol(1,i_order+1) = 4._SP*pi*dot_product(Efield(1)%versor(:),FL_Pi_pol(:,i1_fl)) &
                                           &      * Divide_by_Field(Efield(1),i_order)
     !
   else
     !
     FL_Chi_pol(:,i_order+1) = FL_Pi_pol(:,i1_fl) * Divide_by_Field(Efield(1),i_order)
     !
   endif
   !
 enddo
 !
#if defined _TIMING
 call timing('FL Berry Pol NEQ',OPR='stop')
#endif
 !
end subroutine
