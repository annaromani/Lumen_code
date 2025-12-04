!
! License-Identifier: GPL
!
! Copyright (C) 2006 The Yambo Team
!
! Authors (see AUTHORS file for details): DS
!
!> INIT_RT_ctl_switch
!>
!> @param what
!>
subroutine INIT_RT_ctl_switch(what)
 !
 implicit none
 !
 character(1), intent(in) :: what
 !
 ! Work Space
 !
 call initactivate(1,what//'fnRTdb')
 call initactivate(1,what//'fnRT_N')
 !
end subroutine
