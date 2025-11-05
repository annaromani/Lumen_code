!
! License-Identifier: GPL
!
! Copyright (C) 2016 The Yambo Team
!
! Authors (see AUTHORS file for details): CA
!
!> @brief Apply residual symmetries if present
!!
subroutine Apply_symm_vec(vec)
 !
 ! Apply symmetries to a complex vector      
 !
 use pars,           ONLY:SP,DP,cZERO
 use D_lattice,      ONLY:nsym,dl_sop,i_time_rev
 use electric,       ONLY:l_no_res_symms
 !
 implicit none           
 !
 complex(SP),   intent(inout)  :: vec(3)
 !
 ! Work Space
 !
 complex(SP)   :: V_tmp(3)
 integer       :: is
 !
 if(nsym/=1) then
   V_tmp=cZERO
!$omp parallel do default(shared), private(is),reduction(+:V_tmp)
   do is=1,nsym/(i_time_rev+1) ! Time-reversal is never present 
     V_tmp=V_tmp+matmul(dl_sop(:,:,is),vec)
   enddo
!$omp end parallel do
   vec=V_tmp/dble(nsym/(i_time_rev+1))
 else
   return
 endif
 !
end subroutine Apply_symm_vec
