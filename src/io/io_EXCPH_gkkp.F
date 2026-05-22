!
!        Copyright (C) 2000-2019 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AM
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
integer function io_EXCPH_gkkp(ID)
 !! **io_EXCPH_gkkp - Read/write ndb.excph_gkkp database**
 !!
 !! Handles I/O for the excitonic electron-phonon coupling database,
 !! which stores the matrix elements of the lattice perturbation dV
 !! in the basis of BSE exciton states:
 !!
 !!   EXCPH_gkkp(nu, S', S)    = <S'| dV_nu(q) |S>    complex(SP)
 !!   EXCPH_gkkp_sq(nu, S', S) = |<S'| dV_nu(q) |S>|^2  real(SP)
 !!
 !! where nu is the phonon mode, S is the initial exciton state (1..EXCPH_states(2)),
 !! and S' is the scattered exciton state (1..EXCPH_sum(2)).
 !!
 !! Database layout (fragmented by q-point iq):
 !!   Section 1 (header): nqbz, EXCPH_states, EXCPH_sum, ph_modes, L_kind_in, L_kind_out
 !!   Section iq+1 (fragment iq):
 !!     EXCPH_Q{iq}                   - q-point coordinates (3)          [write only]
 !!     EXCITON_PH_GKKP_Q{iq}         - EXCPH_gkkp  (2, ph_modes, S', S) complex
 !!     EXCITON_PH_GKKP_SQUARED_Q{iq} - EXCPH_gkkp_sq (ph_modes, S', S)  real
 !!
 !! The header also checks consistency of L_kind_in/L_kind_out to ensure the
 !! BSE L matrices used during writing match those used when reading.
 !!
 use pars,           ONLY:SP,schlen
 use IO_m,           ONLY:io_sec,read_is_on,frag_ELPH,io_status,io_mode,DUMP
 use IO_int,         ONLY:io_connect,io_disconnect,io_elemental,io_bulk,io_header,&
&                         def_variable_bulk,io_variable_bulk,io_fragment
 use EXCPH,          ONLY:EXCPH_Gkkp,EXCPH_Gkkp_sq,EXCPH_sum,EXCPH_states,&
&                         BS_EXCPH_E,ph_EXCPH_E,EXCPH_q,L_kind_in,L_kind_out,Sat_alloc, &
&                         EXCPH_gkkp_alloc
 use IO_int,         ONLY:io_variable_elemental,def_variable_elemental
 use ELPH,           ONLY:ph_modes
 use stderr,         ONLY:intc
 use R_lattice,      ONLY:nqbz
 use timing_m,       ONLY:timing
 !
 use y_memory_alloc
 !
 implicit none
 !
 integer      ::ID
 !
 !Work Space
 !
 integer                :: iq,ierr,VAR_size(4),VAR_size_sq(3),ID_frag
 integer                :: states_usr(2),sum_usr(2),states_db(2),sum_db(2)
 character(schlen)      :: db_name,VAR_name,dim_names(4),dim_names_sq(3)
 !
 call timing('io_EXCPH_gkkp',OPR='start')
 !
 ID_frag=0
 io_EXCPH_gkkp=0
 !
 ! DB name
 !
 db_name='excph_gkkp'
 !
 io_EXCPH_gkkp=io_connect(desc=trim(db_name),type=2,ID=ID,FORCE_READ_MODE=(.not.any((/io_sec(ID,:)==1/))).and.frag_ELPH)
 if (io_EXCPH_gkkp/=0) goto 2
 !
 if (any((/io_sec(ID,:)==1/))) then
   !
   io_EXCPH_gkkp=io_header(ID,R_LATT=.true.,KPTS=.true.,FRAG=frag_ELPH)
   !
   ! Save user-specified ranges before the elemental read may overwrite them (DUMP mode)
   !
   if (read_is_on(ID)) then
     states_usr = EXCPH_states
     sum_usr    = EXCPH_sum
   endif
   !
   call io_elemental(ID,VAR="PARS",VAR_SZ=7,MENU=0)
   !
   call io_elemental(ID,&
&       VAR=" Number of qpoints (BZ)   :",I0=nqbz,CHECK=.true.,OP=(/"=="/))
   !
   call io_elemental(ID,&
&       VAR=" Exciton states           :",I1=EXCPH_states,CHECK=.true.,OP=(/">=","<="/))
   !
   call io_elemental(ID,&
&       VAR=" Exciton sum (scattered)  :",I1=EXCPH_sum,CHECK=.true.,OP=(/">=","<="/))
   !
   call io_elemental(ID,&
&       VAR=" Phonon modes            :",I0=ph_modes,CHECK=.true.,OP=(/"=="/))
   !
   call io_elemental(ID,VAR="",VAR_SZ=0,MENU=1)
   !
   io_EXCPH_gkkp=io_status(ID)
   if (io_EXCPH_gkkp/=0) goto 1
   !
   ! In DUMP mode the elemental read overwrites EXCPH_states/EXCPH_sum with the DB values.
   ! If the user explicitly requested a sub-range (any bound != default 1), restore their
   ! values provided the request lies within the DB range.
   !
   if (read_is_on(ID) .and. io_mode(ID)==DUMP) then
     states_db = EXCPH_states
     sum_db    = EXCPH_sum
     if ((states_usr(1) > 1 .or. states_usr(2) > 1) .and. &
&        states_usr(1) >= states_db(1) .and. states_usr(2) <= states_db(2)) &
&       EXCPH_states = states_usr
     if ((sum_usr(1) > 1 .or. sum_usr(2) > 1) .and. &
&        sum_usr(1) >= sum_db(1) .and. sum_usr(2) <= sum_db(2)) &
&       EXCPH_sum = sum_usr
   endif
   !
   call io_bulk(ID,VAR="EXCITON_STATES",VAR_SZ=(/2/))
   call io_bulk(ID,I1=EXCPH_states)
   !
   call io_bulk(ID,VAR="EXCITON_SUM",VAR_SZ=(/2/))
   call io_bulk(ID,I1=EXCPH_sum)
   !
   call io_bulk(ID,VAR="PHONON_MODES",VAR_SZ=(/1/))
   call io_bulk(ID,I0=ph_modes)
   !
   call def_variable_elemental(ID,"L_kind_in",1,0,0)
   call io_variable_elemental(ID,VAR="[K] L_kind in",CH0=L_kind_in,CHECK=.true.,OP=(/"=="/))
   !
   call def_variable_elemental(ID,"L_kind_out",1,0,0)
   call io_variable_elemental(ID,VAR="[K] L_kind out",CH0=L_kind_out,CHECK=.true.,OP=(/"=="/))
   !
 endif
 !
 io_EXCPH_gkkp=io_status(ID)
 !
 iq=maxval(io_sec(ID,:))-1
 if (iq<=0) goto 2
 !
 ! Fragmentation
 !
 call io_fragment(ID,ID_frag,i_fragment=iq,ierr=ierr)
 if (ierr<0.and.read_is_on(ID)) goto 1
 !
 ! When the DB is fragmented I allow a partial reading checking
 ! if the fragment exists or not.
 ! If the fragment does not exist (ierr<0) I return an error code
 !
 if (ierr<0.and.read_is_on(ID)) then
   io_EXCPH_gkkp=-1
   goto 1
 endif
 !
 ! Excitonic Gkkp (and SE) fragments
 !=========================
 !
 if (read_is_on(ID).and..not.allocated(EXCPH_gkkp)) call EXCPH_gkkp_alloc()
 !
 if (.not. read_is_on(ID)) then
   VAR_name="EXCPH_Q"//trim(intc(iq))
   call def_variable_bulk(ID_frag,trim(VAR_name),1,(/3/),SP)
   call io_variable_bulk(ID_frag,1,R1=EXCPH_q(iq,:))
 endif
 !
 VAR_name="EXCITON_PH_GKKP_Q"//trim(intc(iq))
 VAR_size=(/2,ph_modes,EXCPH_sum(2)-EXCPH_sum(1)+1,EXCPH_states(2)-EXCPH_states(1)+1/)
 dim_names = [character(schlen) :: 'complex','PH_modes','N_exc_sum','N_exc_states']
 call def_variable_bulk(ID_frag,trim(VAR_name),1,VAR_size,SP,dim_names)
 call io_variable_bulk(ID_frag,1,C3=EXCPH_gkkp, &
&   IPOS=(/1,1,EXCPH_sum(1),EXCPH_states(1)/))
 !
 VAR_name="EXCITON_PH_GKKP_SQUARED_Q"//trim(intc(iq))
 VAR_size_sq=(/ph_modes,EXCPH_sum(2)-EXCPH_sum(1)+1,EXCPH_states(2)-EXCPH_states(1)+1/)
 dim_names_sq = [character(schlen) :: 'PH_modes','N_exc_sum','N_exc_states']
 call def_variable_bulk(ID_frag,trim(VAR_name),1,VAR_size_sq,SP,dim_names_sq)
 call io_variable_bulk(ID_frag,1,R3=EXCPH_gkkp_sq, &
&   IPOS=(/1,EXCPH_sum(1),EXCPH_states(1)/))
 !
 io_EXCPH_gkkp=io_status(ID)
 !
1 call io_fragment_disconnect(ID,ID_frag)
 !
2 call io_disconnect(ID=ID)
 !
 call timing('io_EXCPH_gkkp',OPR='stop')
 !
end function
