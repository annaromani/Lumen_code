! 
! License-Identifier: GPL
!
! Copyright (C) 2014 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
! headers
!
#include<y_memory.h>
!
subroutine PARALLEL_Transitions_grouping()
 !
 use pars,          ONLY:SP
 use BS,            ONLY:BS_nT_at_k,BS_nT_grps,BS_nT_grps_at_k,BS_T_grp_1st_el_at_k,BS_n_eh_spaces
 use BS_solvers,    ONLY:BSS_mode
 use parallel_int,  ONLY:PP_redux_wait
 use parallel_m,    ONLY:PAR_IND_eh,PAR_BS_nT_col_grps,PAR_COM_eh_INDEX,myid,&
&                        PAR_IND_T_groups,PAR_IND_Kk_ibz,PAR_BS_T_grps_index,&
&                        PAR_COM_Xk_ibz_INDEX,PAR_IND_T_Haydock,PP_indexes_reset
 use R_lattice,     ONLY:nXkibz
 use vec_operate,   ONLY:sort
 use y_memory_alloc
 !
 implicit none
 !
 ! Work Space
 !
 logical              :: l_check
 integer              :: i_k,i_k_r,i_t_at_k,n_GRP,i_Tgrp,i_Tgrp_at_k,it_ref,i_counter,i_res_ares
 integer, allocatable :: GRP_data(:)
 integer, allocatable :: GRP_index(:)
 !
 ! E/h pairs (k resolved)
 !
 ! In this subroutine I distribute the eh transitions within each k. The COMM for this indexing is PAR_COM_eh_INDEX.
 ! Once the group of transitions are build-up the goal is to map them in the total world of transitions groups.
 !
 n_GRP=PAR_COM_eh_INDEX%n_CPU*BS_n_eh_spaces
 !
 YAMBO_ALLOC(BS_T_grp_1st_el_at_k,(n_GRP,nXkibz))
 YAMBO_ALLOC(BS_nT_grps_at_k,(nXkibz,BS_n_eh_spaces))
 allocate(GRP_data(n_GRP),GRP_index(n_GRP))
 !
 ! Here PAR_BS_nT_grps is the number of groups without taking into account the k-parallelization.
 ! This is done to define properly the BS_T_grp_1st_el_at_k array
 !
 BS_T_grp_1st_el_at_k=0
 PAR_BS_nT_col_grps =0
 !
 i_counter=0
 it_ref=PAR_COM_eh_INDEX%CPU_id+1
 !
 do i_res_ares=1,BS_n_eh_spaces
 if(i_res_ares==1) i_counter=it_ref
 if(i_res_ares==2) i_counter=it_ref+PAR_COM_eh_INDEX%n_CPU
 do i_k=1,nXkibz
   i_k_r=i_k+(i_res_ares-1)*nXkibz
   !
   if (.not.PAR_IND_Kk_ibz%element_1D(i_k) ) cycle
   !
   do i_t_at_k=1,BS_nT_at_k(i_k_r)
     !
     if (.not.PAR_IND_eh(i_k_r)%element_1D(i_t_at_k)) cycle
     !
     if ( i_t_at_k==1 ) l_check= .false.
     if ( i_t_at_k> 1 ) l_check= PAR_IND_eh(i_k_r)%element_1D(i_t_at_k-1)
     !
     if ( l_check ) cycle
     !
     PAR_BS_nT_col_grps = PAR_BS_nT_col_grps +1
     !
     BS_T_grp_1st_el_at_k(i_counter,i_k)=i_t_at_k+sum(BS_nT_at_k(:i_k_r-1))
     !
   enddo
 enddo
 enddo
 !
 call PP_redux_wait(BS_T_grp_1st_el_at_k,COMM=PAR_COM_Xk_ibz_INDEX%COMM)
 call PP_redux_wait(BS_T_grp_1st_el_at_k,COMM=PAR_COM_eh_INDEX%COMM)
 !
 ! ALL k are diveded in BS_nT_grps groups of transitions and each CPU takes PAR_BS_nT_grps of them
 !
 BS_nT_grps=0
 BS_nT_grps_at_k=0
 !
 do i_k=1,nXkibz
   !
   GRP_data(:)=BS_T_grp_1st_el_at_k(:,i_k)
   !
   call sort(GRP_data,indx=GRP_index)
   !
   BS_T_grp_1st_el_at_k(:,i_k)=0
   !
   do i_res_ares=1,BS_n_eh_spaces
   i_counter=PAR_COM_eh_INDEX%n_CPU*(i_res_ares-1)
   do i_Tgrp=1+i_counter,n_GRP/BS_n_eh_spaces+i_counter
     if (GRP_data(i_Tgrp)>0._SP) then
       BS_nT_grps=BS_nT_grps+1
       BS_nT_grps_at_k(i_k,i_res_ares)=BS_nT_grps_at_k(i_k,i_res_ares)+1
       BS_T_grp_1st_el_at_k(sum(BS_nT_grps_at_k(i_k,:)),i_k)=GRP_data(i_Tgrp)
     endif
   enddo
   enddo
   !
 enddo
 !
 ! Groups index
 !
 YAMBO_ALLOC(PAR_IND_T_groups%element_1D,(BS_nT_grps))
 YAMBO_ALLOC(PAR_IND_T_groups%n_of_elements,(PAR_COM_eh_INDEX%n_CPU))
 YAMBO_ALLOC(PAR_BS_T_grps_index,(BS_nT_grps))
 !
 if ((index(BSS_mode,'h')/=0) .or. (index(BSS_mode,'s')/=0)) then
   call PP_indexes_reset(PAR_IND_T_Haydock)
   YAMBO_ALLOC(PAR_IND_T_Haydock%element_1D,(BS_nT_grps))
   PAR_IND_T_Haydock%element_1D=.FALSE.
 endif
 !
 PAR_IND_T_groups%element_1D=.FALSE.
 PAR_IND_T_groups%n_of_elements=0
 PAR_BS_nT_col_grps = 0
 PAR_BS_T_grps_index = 0
 !
 do i_res_ares=1,BS_n_eh_spaces
 do i_k=1,nXkibz
   i_k_r=i_k+(i_res_ares-1)*nXkibz
   !
   if (.not.PAR_IND_Kk_ibz%element_1D(i_k)) cycle
   !
   do i_t_at_k=1,BS_nT_at_k(i_k_r)
     !
     if (.not.PAR_IND_eh(i_k_r)%element_1D(i_t_at_k)) cycle
     !
     if ( i_t_at_k==1 ) l_check= .false.
     if ( i_t_at_k> 1 ) l_check= PAR_IND_eh(i_k_r)%element_1D(i_t_at_k-1)
     !
     if ( l_check ) cycle
     !
     PAR_BS_nT_col_grps = PAR_BS_nT_col_grps + 1
     !
     it_ref=i_t_at_k+sum(BS_nT_at_k(:i_k_r-1))
     !
     i_counter=(i_res_ares-1)*BS_nT_grps_at_k(i_k,1)
     do i_Tgrp_at_k=1+i_counter,BS_nT_grps_at_k(i_k,i_res_ares)+i_counter
       if (BS_T_grp_1st_el_at_k(i_Tgrp_at_k,i_k)==it_ref) then
         i_Tgrp=i_Tgrp_at_k+sum(BS_nT_grps_at_k(:i_k-1,:))
         PAR_IND_T_groups%element_1D(i_Tgrp)=.true.
         PAR_IND_T_groups%n_of_elements(PAR_COM_eh_INDEX%CPU_id+1)=&
&                                         PAR_IND_T_groups%n_of_elements(PAR_COM_eh_INDEX%CPU_id+1)+1
         !
         ! Warning: this index is modified by the sorting previously done.
         !          Should be probably defined via the GRP_index or similar
         !          However it is not presently used in the code.
         PAR_BS_T_grps_index(i_Tgrp)=PAR_BS_nT_col_grps
       endif
     enddo
     !
   enddo
 enddo
 enddo
 !
 ! CLEAN UP
 !
 deallocate(GRP_data,GRP_index)
 !
end subroutine
