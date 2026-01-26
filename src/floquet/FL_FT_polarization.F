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
subroutine FL_FT_polarization(FL_time_pol,FL_floquet_pol)
 !
 use pars,           ONLY:DP,SP,cZERO,pi,cONE,cI
 use units,          ONLY:HA2FSm1
 use fields,         ONLY:Efield,EtoT
 use wrapper,        ONLY:M_by_V
 use interfaces,      ONLY:LINEAR_ALGEBRA_driver
 use linear_algebra,  ONLY:INV
 use fl_optics,      ONLY:max_fl_mode,n_tot_fl_modes,FL_space_dim,FL_S_plus,FL_S_minus,&
&                         FL_states,FL_Ho_plus_Sigma,FL_time_S_det,n_tot_time_steps,max_time_step
 !
 implicit none
 !
 complex(DP), intent(inout) :: FL_time_pol(n_tot_time_steps)
 complex(DP), intent(inout) :: FL_floquet_pol(n_tot_time_steps)
 !
 ! Work Space
 !
 complex(DP) :: FL_FT_Mat(n_tot_time_steps,n_tot_time_steps)
 complex(DP) :: tmp_P_t
 complex(SP) :: dummy_det
 integer     :: i_t,i1_fl,i1_shifted_fl
 real(SP)    :: time,w0_sm1
 !
 ! Build Matrix for Fourier Transform
 ! ===================================================================
 !
 FL_FT_Mat = cZERO
 FL_floquet_pol = cZERO
 !
 do i_t=1,n_tot_time_steps
   !
   do i1_fl=1,n_tot_time_steps
     !
     i1_shifted_fl = i1_fl -(max_time_step+1)
     !
     FL_FT_Mat(i_t,i1_fl) =  exp(-1._SP*cI*real(i1_shifted_fl,SP)*2._SP*pi*(i_t-1._SP)/n_tot_time_steps)
     !
   enddo
   !
 enddo
 !
 ! Invert Matrix for Fourier Transform
 ! ===================================================================
 !
 call LINEAR_ALGEBRA_driver(INV,M=FL_FT_Mat)
 !
 ! Calculate FL_floquet_pol
 ! ===================================================================
 !
  do i1_fl=(max_time_step+1),n_tot_time_steps
    do i_t=1,n_tot_time_steps
    !
    FL_floquet_pol(i1_fl)=FL_floquet_pol(i1_fl)+&
            &             FL_FT_Mat(i1_fl,i_t)*FL_time_pol(i_t)
    !
    enddo
  enddo
 !
end subroutine FL_FT_polarization
