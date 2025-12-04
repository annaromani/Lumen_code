!
! License-Identifier: GPL
!
! Copyright (C) 2020 The Yambo Team
!
! Authors (see AUTHORS file for details): AM
!
module C_driver
 !
 implicit none
 !
 integer      :: code_version(3)
 integer      :: code_GPL_revision
 integer      :: code_revision
 character(50):: code_bin
 character(50):: code_tool
 character(50):: code_project
 character(50):: code_hash
 character(50):: code_libraries
 character(50):: code_branch
 !
 interface  
   !
!> get_version
!>
!> @param version
!> @param subversion
!> @param patchlevel
!> @param revision
!> @param hash
!>
   subroutine get_version(version,subversion,patchlevel,revision,hash)
     integer     , intent(out):: version
     integer     , intent(out):: subversion
     integer     , intent(out):: patchlevel
     integer     , intent(out):: revision
     character(*), intent(out):: hash
   end subroutine
   !
!> get_libraries
!>
!> @param libraries
!>
   subroutine get_libraries(libraries)
     character(*), intent(out):: libraries
   end subroutine
   !
!> get_running_project
!>
!> @param code_project
!>
   subroutine get_running_project(code_project)
     character(*), intent(out):: code_project
   end subroutine
   !
!> get_running_tool
!>
!> @param code_tool
!>
   subroutine get_running_tool(code_tool)
     character(*), intent(out):: code_tool
   end subroutine
   !
!> get_runlevel
!>
!> @param component
!> @param component_id
!> @param runlevel_id
!>
   subroutine get_runlevel(component,component_id,runlevel_id)
     character(*), intent(out):: component
     integer,      intent(in) :: component_id
     integer,      intent(in) :: runlevel_id
   end subroutine
   !
 end interface
 !
end module C_driver
