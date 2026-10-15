!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): DS
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
subroutine RT_TRabs_rotate_exciton(DIP_kind,DIP_dir,DIP_f,Xen,Xk)
 !
 ! In Intra/Inter excitons absorption the system is initially in an excitonic state
 ! The state is chosen setting in input the excitonic state and the polarization
 ! If the excitonic state belongs to a degenerate space, the state with 
 ! excitonic dipole along the input variable needs to be selected.
 !
 ! The present subroutine does this operation. BS_mat modified
 !
 use pars,            ONLY:SP,cZERO,cONE
 use com,             ONLY:msg
 use R_lattice,       ONLY:bz_samp
 use stderr,          ONLY:intc
 use electrons,       ONLY:levels,spin,spin_occ
 use X_m,             ONLY:global_gauge
 use YPP_real_time,   ONLY:TRabs_Eeh_treshold,TRabs_Exc_Pol
 use BS,              ONLY:BS_K_dim,BS_bands
 use BS_solvers,      ONLY:BSS_eh_table,BS_mat,BSS_n_eig
 use YPPm,            ONLY:BS_E,EXCITONS_n_user_states,EXCITONS_user_indexes
 use vec_operate,     ONLY:v_norm,v_norm_c,cross_product
 use matrix_operate,  ONLY:m3inv_c,m3det_c
 use parallel_m,      ONLY:PAR_IND_DIPk_ibz,PAR_IND_VAL_BANDS_DIP,PAR_IND_CON_BANDS_DIP
 use parallel_int,    ONLY:PP_redux_wait
 !
 use y_memory_alloc
 !
 implicit none
 !
 type(levels),  intent(in) :: Xen
 type(bz_samp), intent(in) :: Xk
 complex(SP),   intent(in) :: DIP_f
 real(SP),   intent(inout) :: DIP_dir(3,3)
 !
 ! Work space
 !
 integer           :: i_cv,i_v,i_c,i_kbz,i_kibz,i_spin,i_l,i_lp,i_dir,j_dir,i_exc,i_exc_p
 character(6)      :: DIP_kind
 real(SP)          :: E_eh,DIP_proj(3,2)
 complex(SP)       :: DIP_ver(3,3),DIP_exc(3,3),DIP_inv(3,3),U(3,3),E_ver(3,3)
 complex(SP)       :: DIP_rot(3),f_eh(2),DIP_dot,DIP_tmp(3)
 complex(SP), allocatable :: A_vec_tmp(:,:)
 !
 if (EXCITONS_n_user_states>3) then
   call msg("s","No rotation performed. Too many excitons selected in input ")
   return 
 endif
 !
 ! I need to select the exciton with the dipole in the direction asked from input
 !
 ! Step 1: compute the dipoles of the degenerate excitons
 !
 do i_exc=1,EXCITONS_n_user_states
   !
   i_l=EXCITONS_user_indexes(i_exc)
   !
   DIP_tmp=0._SP
   !
   do i_cv = 1,BS_K_dim(1)
     !
     i_kbz   = BSS_eh_table(i_cv,1)
     i_v     = BSS_eh_table(i_cv,2)
     i_c     = BSS_eh_table(i_cv,3)
     i_spin  = spin(BSS_eh_table(i_cv,:))
     !
     i_kibz  = Xk%sstar(i_kbz,1)
     !
     if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
     if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
     if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
     !
     E_eh=Xen%E(i_c,i_kibz,i_spin)-Xen%E(i_v,i_kibz,i_spin)
     if (abs(E_eh) < TRabs_Eeh_treshold.and. trim(global_gauge)=="length") cycle
     !
     f_eh=(Xen%f(i_v,i_kibz,i_spin)-Xen%f(i_c,i_kibz,i_spin))/spin_occ
     !
     call DIPOLE_rotate(i_c,i_v,i_kbz,i_spin,trim(DIP_kind),Xk,DIP_rot)
     !
     DIP_tmp= DIP_tmp+DIP_f*DIP_rot(:)*BS_mat(i_cv,i_l)*sqrt(f_eh(1))
     !
   enddo
   !
   DIP_exc(:,i_exc)=DIP_tmp
   !
 enddo
 !
 call PP_redux_wait(DIP_exc)
 !
 do i_exc=1,EXCITONS_n_user_states
   DIP_ver(:,i_exc)=DIP_exc(:,i_exc)/v_norm_c(DIP_exc(:,i_exc))
 enddo
 !
 if (EXCITONS_n_user_states==1) then
   call msg("s","No rotation performed since exciton is not degenerate to anything ")
   return 
 endif
 !
 if (EXCITONS_n_user_states==2) then
   DIP_ver(:,3) = cross_product(DIP_ver(:,1),DIP_ver(:,2))
   DIP_ver(:,3) = DIP_ver(:,3)/v_norm_c(DIP_ver(:,3))
   DIP_exc(:,3) = DIP_ver(:,3)*(v_norm_c(DIP_exc(:,1))+v_norm_c(DIP_exc(:,2)))/2._SP
 endif
 !
 call msg("s","checking initial dipoles relative modulus ")
 call msg("s","Dipole modules ",(/v_norm_c(DIP_exc(:,1)),v_norm_c(DIP_exc(:,2)),v_norm_c(DIP_exc(:,3))/))
 !
 call msg("s","checking initial versors ortogonality ")
 do i_exc=1,EXCITONS_n_user_states
   i_l=EXCITONS_user_indexes(i_exc)
   do i_exc_p=i_exc+1,EXCITONS_n_user_states
     i_lp=EXCITONS_user_indexes(i_exc_p)
     DIP_proj(1,1)= real(dot_product(DIP_ver(:,i_exc),DIP_ver(:,i_exc_p)))
     DIP_proj(1,2)=aimag(dot_product(DIP_ver(:,i_exc),DIP_ver(:,i_exc_p)))
     call msg("s","States "//trim(intc(i_l))//" - "//trim(intc(i_lp)),(/DIP_proj(1,1),DIP_proj(1,2)/))
   enddo
 enddo
 call msg("s"," ")
 !
 ! Step 2: Define versors in the degenerate space
 !
 TRabs_Exc_Pol=TRabs_Exc_Pol/v_norm(TRabs_Exc_Pol)
 !
 E_ver(:,1) = cONE*TRabs_Exc_Pol
 E_ver(:,2)=(/cZERO,cONE,cZERO/)
 E_ver(:,3)=(/cZERO,cZERO,cONE/)
 !
 if (EXCITONS_n_user_states==3) then
   !
   if( (1._SP-abs(dot_product( E_ver(:,1),E_ver(:,2) )))<1.E-5_SP) E_ver(:,2)=(/cONE,cZERO,cZERO/)
   E_ver(:,2) = E_ver(:,2)-dot_product( E_ver(:,1),E_ver(:,2) )*E_ver(:,1)
   E_ver(:,2) = E_ver(:,2)/v_norm_c(E_ver(:,2))
   !
 else if (EXCITONS_n_user_states==2) then
   !
   E_ver(:,2) = DIP_ver(:,1)-dot_product(E_ver(:,1),DIP_ver(:,1))*E_ver(:,1)+ &
   &            DIP_ver(:,2)-dot_product(E_ver(:,1),DIP_ver(:,2))*E_ver(:,1)
   E_ver(:,2) = E_ver(:,2)/v_norm_c(E_ver(:,2))
   !
 endif
 !
 if( (1._SP-abs(dot_product( E_ver(:,1),E_ver(:,3) )))<1.E-5) E_ver(:,3)=(/cONE,cZERO,cZERO/)
 E_ver(:,3) = E_ver(:,3)-dot_product( E_ver(:,1),E_ver(:,3) )*E_ver(:,1) &
 &                      -dot_product( E_ver(:,2),E_ver(:,3) )*E_ver(:,2)
 E_ver(:,3) = E_ver(:,3)/v_norm_c(E_ver(:,3))
 !!
 !! Make the versors real by removing the phase
 !!
 do j_dir=1,3
   do i_dir=1,3
     if(abs(aimag(E_ver(i_dir,j_dir)))<1.E-5) cycle
     !write(*,*) j_dir,E_ver(:,j_dir)
     !write(*,*) j_dir,E_ver(:,j_dir)*conjg(E_ver(i_dir,j_dir))/abs(E_ver(i_dir,j_dir))
     E_ver(:,j_dir) = E_ver(:,j_dir)*conjg(E_ver(i_dir,j_dir))/abs(E_ver(i_dir,j_dir))
     exit
   enddo
 enddo
 !
 DIP_dir=real(transpose(E_ver),SP)
 !
 call msg("s","Projection versors ")
 do i_dir=1,3
   call msg("s","Versors "//trim(intc(i_dir))//" Re",real(E_ver(:,i_dir)))
   call msg("s","          Im",                   aimag(E_ver(:,i_dir)))
 enddo
 call msg("s"," ")
 !
 call msg("s","Excitonic dipole projections before rotation ")
 do i_exc=1,EXCITONS_n_user_states
   i_l=EXCITONS_user_indexes(i_exc)
   do i_dir=1,3
     DIP_proj(i_dir,1)= real(dot_product(DIP_ver(:,i_exc),E_ver(:,i_dir)))
     DIP_proj(i_dir,2)=aimag(dot_product(DIP_ver(:,i_exc),E_ver(:,i_dir)))
   enddo
   call msg("s","Exciton "//trim(intc(i_l))//" Re",DIP_proj(:,1))
   call msg("s","          Im",                    DIP_proj(:,2))
 enddo
 call msg("s","Exciton dipoles determinant ",abs(m3det_c(DIP_ver)))
 call msg("s"," ")
 !
 !DEBUG <
 !write(*,*) " "
 !write(*,*) "E_dir"
 !do i_dir=1,3
 !   write(*,*) E_ver(i_dir,:)
 !enddo
 !write(*,*) "m3det= ",abs(m3det_c(E_ver))
 !!
 !write(*,*) " "
 !write(*,*) "DIP_ver"
 !do i_dir=1,3
 !   write(*,*) DIP_ver(:,i_dir),dot_product(DIP_ver(:,i_dir),DIP_ver(:,1))
 !enddo
 !write(*,*) "m3det= ",abs(m3det_c(DIP_ver))
 !DEBUG <
 !
 ! Step 3: find the proper rotations in the degenerate space
 !
 DIP_inv=transpose(DIP_ver)
 call m3inv_c(DIP_inv)
 U=matmul(transpose(E_ver),DIP_inv)
 !write(*,*) "inv check ",matmul(transpose(DIP_ver),DIP_inv)
 !write(*,*) "  U check ",matmul(U,transpose(DIP_ver))
 !U=U/(abs(m3det_c(U)))**(1/3._SP)
 !
 !DEBUG <
 !write(*,*) " "
 !write(*,*) "U"
 !do i_dir=1,3
 !   write(*,*) U(i_dir,:)
 !enddo
 !write(*,*) "m3det= ",abs(m3det_c(U))
 !!
 !call m3inv_c(DIP_ver)
 !DIP_ver(:,3)=cZERO
 !E_ver=matmul(DIP_ver,U)
 !write(*,*) " "
 !write(*,*) "DIP_rot"
 !do i_dir=1,3
 !   write(*,*) E_ver(i_dir,:)
 !enddo
 !write(*,*) "m3det= ",abs(m3det_c(E_ver))
 !DEBUG >
 !
 ! Step 4: rotate the excitons
 !
 allocate(A_vec_tmp(BS_K_dim(1),EXCITONS_n_user_states))
 !
 A_vec_tmp=cZERO
 do i_exc=1,EXCITONS_n_user_states
   do i_exc_p=1,EXCITONS_n_user_states
     i_l=EXCITONS_user_indexes(i_exc_p)
     A_vec_tmp(:,i_exc)= A_vec_tmp(:,i_exc)+U(i_exc,i_exc_p)*BS_mat(:,i_l)
   enddo
 enddo
 do i_exc=1,EXCITONS_n_user_states
   i_l=EXCITONS_user_indexes(i_exc)
   BS_mat(:,i_l)=A_vec_tmp(:,i_exc)
 enddo
 !
 deallocate(A_vec_tmp)
 !
 ! Step 5: Check new excitonic dipoles
 !
 DIP_exc(:,:EXCITONS_n_user_states)=cZERO
 !
 do i_exc=1,EXCITONS_n_user_states
   !
   i_l=EXCITONS_user_indexes(i_exc)
   !
   do i_cv = 1,BS_K_dim(1)
     !
     i_kbz   = BSS_eh_table(i_cv,1)
     i_v     = BSS_eh_table(i_cv,2)
     i_c     = BSS_eh_table(i_cv,3)
     i_spin  = spin(BSS_eh_table(i_cv,:))
     !
     i_kibz  = Xk%sstar(i_kbz,1)
     !
     if (.not. PAR_IND_CON_BANDS_DIP%element_1D(i_c)) cycle
     if (.not. PAR_IND_VAL_BANDS_DIP%element_1D(i_v)) cycle
     if (.not. PAR_IND_DIPk_ibz%element_1D(i_kibz)) cycle
     !
     E_eh=Xen%E(i_c,i_kibz,i_spin)-Xen%E(i_v,i_kibz,i_spin)
     if (abs(E_eh) < TRabs_Eeh_treshold.and. trim(global_gauge)=="length") cycle
     !
     f_eh=(Xen%f(i_v,i_kibz,i_spin)-Xen%f(i_c,i_kibz,i_spin))/spin_occ
     !
     call DIPOLE_rotate(i_c,i_v,i_kbz,i_spin,trim(DIP_kind),Xk,DIP_rot)
     !
     DIP_exc(:,i_exc)= DIP_exc(:,i_exc)+DIP_f*DIP_rot(:)*BS_mat(i_cv,i_l)*sqrt(f_eh(1))
     !
   enddo
   !
 enddo
 !
 call PP_redux_wait(DIP_exc)
 !
 do i_exc=1,EXCITONS_n_user_states
   DIP_ver(:,i_exc)=DIP_exc(:,i_exc)/v_norm_c(DIP_exc(:,i_exc))
 enddo
 !
 call msg("s","Excitonic dipole projections after rotation ")
 do i_exc=1,EXCITONS_n_user_states
   i_l=EXCITONS_user_indexes(i_exc)
   do i_dir=1,3
     DIP_proj(i_dir,1)= real(dot_product(DIP_ver(:,i_exc),E_ver(:,i_dir)))
     DIP_proj(i_dir,2)=aimag(dot_product(DIP_ver(:,i_exc),E_ver(:,i_dir)))
   enddo
   call msg("s","Exciton "//trim(intc(i_l))//" Re",DIP_proj(:,1))
   call msg("s","          Im",                    DIP_proj(:,2))
 enddo
 call msg("s","Exciton dipoles determinant ",abs(m3det_c(DIP_ver)))
 call msg("s"," ")
 !
 !write(*,*) " "
 !write(*,*) "DIP_ver"
 !do i_dir=1,3
 !   write(*,*) DIP_ver(i_dir,:)
 !enddo
 !write(*,*) "m3det= ",abs(m3det_c(DIP_ver))
 !
end subroutine RT_TRabs_rotate_exciton
