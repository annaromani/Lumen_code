!
! License-Identifier: GPL
!
! Copyright (C) 2025 The Yambo Team
!
! Authors (see AUTHORS file for details): CA
!
! headers
!
#include<y_memory.h>
!
integer function io_full_SYMMs(En_old,k_old,nkibz,ID)
 !
 use pars,          ONLY:schlen,SP,IP,IP_YIO,LP_YIO
 use D_lattice,     ONLY:kpoints_map,old_nsym,old_dl_sop,old_dl_tra,&
 &                       old_S_contains_TR,old_i_time_rev
 use R_lattice,     ONLY:bz_samp
 use electrons,     ONLY:levels,n_sp_pol,n_spinor
 use IO_m,          ONLY:io_sec,read_is_on,write_is_on,io_extension
 use IO_int,        ONLY:io_connect,io_disconnect,io_bulk,def_variable_elemental,io_bulk, &
&                        io_variable_elemental,io_variable_bulk,def_variable_bulk,io_header
 !
#if defined _TIMING
 use timing_m,      ONLY:timing
#endif
 use y_memory_alloc
 !
 implicit none
 !
 type(levels), intent(inout) :: En_old
 type(bz_samp),intent(inout) :: k_old
 integer, intent(in) :: ID
 integer, intent(in) :: nkibz
 !
 ! Work Space
 !
 integer :: is
 !
#if defined _TIMING
 call timing('io_full_SYMMs',OPR='start')
#endif
 !
 io_extension(ID)="full_SYMMs"
 !
 io_full_SYMMs=io_connect(desc=trim(io_extension(ID)),type=1,ID=ID)
 if (io_full_SYMMs/=0) goto 1
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   io_full_SYMMs=io_header(ID,D_LATT=.true.)
   if (io_full_SYMMs/=0) goto 1
   !
   call def_variable_elemental(ID,VAR="old_nsym",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,VAR="Old number of symmetries",I0=old_nsym)
   !
   call def_variable_elemental(ID,VAR="old_trev",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,VAR="Old time-reversal symmetry",I0=old_i_time_rev)
   !
   call def_variable_elemental(ID,VAR="old_nkibz",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,VAR="Old number of kpt in the ibz",I0=k_old%nibz)
   !
   call def_variable_elemental(ID,VAR="old_nbands",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,VAR="Old number of bands",I0=En_old%nb)
   !
 endif
 !
 if(read_is_on(ID)) then
   En_old%nk =k_old%nibz
   YAMBO_ALLOC(k_old%pt,(k_old%nibz,3))
   YAMBO_ALLOC(En_old%E,(En_old%nb,En_old%nk,n_sp_pol))
   !
   YAMBO_ALLOC(old_dl_sop,(3,3,old_nsym*n_spinor))
   YAMBO_ALLOC(old_dl_tra,(3,old_nsym*n_spinor))
   YAMBO_ALLOC(old_S_contains_TR,(old_nsym*n_spinor))
   YAMBO_ALLOC(kpoints_map,(2,nkibz))
 endif
 !
 ! K-point grid
 !
 call def_variable_bulk(ID, VAR="K-POINTS", VAR_ID=1,VAR_SZ=shape(k_old%pt),VAR_KIND=SP )
 call io_variable_bulk(ID,1,R2=k_old%pt)
 !
 ! Eigenvalues
 !
 call def_variable_bulk(ID, VAR="EIGENVALUES", VAR_ID=2, VAR_SZ=shape(En_old%E), VAR_KIND=SP )
 call io_variable_bulk(ID,2, R3=En_old%E )
 !
 call def_variable_bulk(ID,VAR="old_DL_dop",VAR_ID=3,VAR_SZ=(/3,3,old_nsym/),VAR_KIND=SP)
 call io_variable_bulk(ID,3,R3=old_dl_sop(:,:,:old_nsym))
 !
 call def_variable_bulk(ID,VAR="old_DL_tra",VAR_ID=4,VAR_SZ=(/3,old_nsym/),VAR_KIND=SP)
 call io_variable_bulk(ID,4,R2=old_dl_tra(:,:old_nsym))
 !
 if (read_is_on(ID)) then
   do is=1,old_nsym
     old_S_contains_TR(is)=is>old_nsym/(1+old_i_time_rev)
   enddo
   !
   ! Double point group
   !
   if (n_spinor==2) then
     old_dl_sop(:,:,old_nsym+1:)=old_dl_sop(:,:,:old_nsym)
     old_dl_tra(:,old_nsym+1:)=old_dl_tra(:,:old_nsym)
     old_S_contains_TR(old_nsym+1:)=old_S_contains_TR(:old_nsym)
   endif
 endif
 !
 !call def_variable_bulk(ID,VAR="kpoints_map",VAR_ID=5,VAR_SZ=shape(kpoints_map),VAR_KIND=IP_YIO)
 !call io_variable_bulk(ID,5,I2=kpoints_map)
 !
1  call io_disconnect(ID=ID)
 !
#if defined _TIMING
 call timing('io_full_SYMMs',OPR='stop')
#endif
 !
end function
