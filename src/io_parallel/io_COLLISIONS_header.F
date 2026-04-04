!
! License-Identifier: GPL
!
! Copyright (C) 2021 The Yambo Team
!
! Authors (see AUTHORS file for details): AM CA DS
!
! headers
!
#include<y_memory.h>
!
integer function io_COLLISIONS_header(ID,COLL)
 !
 use pars,          ONLY:schlen,IP,IP_YIO,LP_YIO,SP,CR_YIO
 use R_lattice,     ONLY:RIM_ng,RIM_n_rand_pts
 use collision_ext, ONLY:COLLISIONS_cutoff,COLL_bands,COLLISIONS_group,&
&                        COLLISIONS_CV_only,COLL_bands_frozen,COLL_check_cutoff
 use IO_m,          ONLY:variable_exist,write_is_on,io_status,IO_OUTDATED_DB,read_is_on
 use IO_int,        ONLY:io_connect,io_header,io_variable_elemental,ver_is_gt_or_eq, &
&                        io_variable_bulk,def_variable_bulk,io_disconnect,  &
&                        def_variable_elemental,&
&                        io_variable_elemental
 use global_XC,     ONLY:QP_SE_GW_NEQ,QP_SE_GW_QED,&
&                        H_SE_FOCK,H_SE_COH,H_SE_SEX,H_SE_HARTREE
 use plasma,        ONLY:N_plasma_poles_global,PLASMA_redux_percent
 use QP_m,          ONLY:QP_ng_Sc,QP_ng_Sx,QP_ng_SH,SC_band_mixing,QP_QED_ng
 use X_m,           ONLY:X_desc
 use drivers,       ONLY:l_elel_scatt,l_elphoton_scatt
 use timing_m,      ONLY:timing
 use y_memory_alloc
 !
 implicit none
 !
 integer, intent(inout)            ::ID
 type(COLLISIONS_group)            ::COLL
 !
 ! Work Space
 !
 character(schlen)     :: db_name
 integer               :: V_kind_disk
 !
 integer, external     :: io_descriptors
 !
 call timing('io_COLLISIONS',OPR='start')
 !
 call COLLISIONS_dbname(COLL%kind,db_name)
 !
 io_COLLISIONS_header=io_connect(desc=trim(db_name)//"_header",type=2,ID=ID)
 if(io_COLLISIONS_header/=0) goto 1
 !
 io_COLLISIONS_header=io_header(ID,QPTS=.true.,R_LATT=.true.,WF=.true.,IMPOSE_SN=.true.)
 !
 if (.not.ver_is_gt_or_eq(ID,revision=5453)) then
   io_COLLISIONS_header=IO_OUTDATED_DB
   goto 1
 endif
 !
 call def_variable_elemental(ID,"COLLISIONS_name",1,CR_YIO,0)
 call io_variable_elemental(ID,CH0=COLL%name,VAR=' COLLISIONS name                 ')
 !
 call def_variable_elemental(ID,"N_Collisions",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
 call io_variable_elemental(ID,I0=COLL%N,VAR='Total number of COLLISIONS      ',CHECK=.true.,OP=(/"=="/))
 !
 call def_variable_elemental(ID,"Coll_Kind",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
 call io_variable_elemental(ID,I0=COLL%kind,      VAR='COLLISIONS kind                 ',CHECK=.true.,OP=(/"=="/))
 !
 call def_variable_elemental(ID,"Coll_Bands",VAR_SZ=2,VAR_KIND=IP_YIO,MENU=0)
 call io_variable_elemental(ID,I1=COLL_bands,VAR='COLLISIONS Bands                    ',CHECK=.true.,OP=(/"==","=="/))
 !
 if (COLL%kind==QP_SE_GW_NEQ) then
   !
   call def_variable_elemental(ID,"l_elel_scatt",1,LP_YIO,0)
   call io_variable_elemental(ID,L0=l_elel_scatt,VAR='Electron-Electron scattering    ',WARN=.true.,OP=(/"=="/))
   !
   if(l_elel_scatt) then
    !
    call def_variable_elemental(ID,"N_plasma_poles",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0) 
    call io_variable_elemental(ID,I0=N_plasma_poles_global,VAR='Max N plasma poles (q-resolved) ',CHECK=.true.,OP=(/"=="/))
    !
    call def_variable_elemental(ID,"PLASMA_redux_percent",1,SP,0)
    call io_variable_elemental(ID,R0=PLASMA_redux_percent,VAR='Plasma redux percent            ',CHECK=.true.,OP=(/"=="/))
    !
   endif
 endif
 !
 if (COLL%kind==QP_SE_GW_QED) then
   !
   call def_variable_elemental(ID,"l_elphoton_scatt",1,LP_YIO,0)
   call io_variable_elemental(ID,L0=l_elphoton_scatt,VAR='Electron-Photon   scattering    ',WARN=.true.,OP=(/"=="/))
   ! 
   call def_variable_elemental(ID,"QP_QED_ng",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,I0=QP_QED_ng,VAR='Electron-Photon   RL vectors    ',CHECK=.true.,OP=(/"=="/))
   !
 endif
 !
 call def_variable_elemental(ID,"SC_band_mixing",1,SP,0)
 call io_variable_elemental(ID,R0=SC_band_mixing,VAR='Bands mixing                    ',CHECK=.true.,OP=(/"=="/))
 !
 call def_variable_elemental(ID,"COLL_kind",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0) 
 call io_variable_elemental(ID,I0=COLL%kind,DB_I0=V_kind_disk,CHECK=.true.,OP=(/"=="/))
 !
 if (ver_is_gt_or_eq(ID,revision=21128)) then
   call def_variable_elemental(ID,"COLL_CV_only",1,LP_YIO,0)    
   call io_variable_elemental(ID,L0=COLLISIONS_CV_only,VAR='CV only scattering    ',WARN=.true.,OP=(/"=="/))
 endif
 !
 if(COLL%kind==H_SE_HARTREE.or.COLL%kind==H_SE_FOCK.or.COLL%kind==H_SE_SEX.or.COLL%kind>=H_SE_HARTREE) then
   !
   if (.not.ver_is_gt_or_eq(ID,revision=24190)) then
     call def_variable_elemental(ID,"COLLISIONS_cutoff",1,SP,0)
     call io_variable_elemental(ID,R0=COLLISIONS_cutoff,VAR='Collisions Cutoff                 ',WARN=.true.,OP=(/"=="/))
   endif
   !
   if (V_kind_disk>=H_SE_HARTREE) then
     call def_variable_elemental(ID,"QP_ng_SH",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
     call io_variable_elemental(ID,I0=QP_ng_SH,&
&         VAR='Hartree  RL vectors             ',WARN=.not.COLL_check_cutoff,CHECK=COLL_check_cutoff,OP=(/"=="/))
   endif
   if (V_kind_disk==             H_SE_FOCK.or.V_kind_disk==             H_SE_SEX.or.&
&        V_kind_disk==H_SE_HARTREE+H_SE_FOCK.or.V_kind_disk==H_SE_HARTREE+H_SE_SEX) then
     call def_variable_elemental(ID,"QP_ng_Sx",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
     call io_variable_elemental(ID,I0=QP_ng_Sx,&
&         VAR='Exchange RL vectors             ',WARN=.not.COLL_check_cutoff,CHECK=COLL_check_cutoff,OP=(/"=="/))
   endif
   if (V_kind_disk==             H_SE_SEX.or.V_kind_disk==H_SE_HARTREE+H_SE_SEX) then
     call def_variable_elemental(ID,"QP_ng_Sc",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
     call io_variable_elemental(ID,I0=QP_ng_Sc,&
&         VAR='Correlation RL vectors          ',WARN=.not.COLL_check_cutoff,CHECK=COLL_check_cutoff,OP=(/"=="/))
   endif
 endif
 !
 if (V_kind_disk==QP_SE_GW_NEQ   .or.V_kind_disk==H_SE_COH  .or.V_kind_disk==H_SE_SEX) then
   !
   ! Global Response Descriptions
   !==============================
   !
   io_COLLISIONS_header=io_descriptors(X_desc,"X",ID,.FALSE.)
   !
   call def_variable_elemental(ID,"RIM_ng",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,I0=RIM_ng,VAR=        ' RIM RL components      ',WARN=.true.)
   call def_variable_elemental(ID,"RIM_n_rand_pts",VAR_SZ=1,VAR_KIND=IP_YIO,MENU=0)
   call io_variable_elemental(ID,I0=RIM_n_rand_pts,VAR=' RIM random points      ',WARN=.true.)
   !
 endif
 !
 if (variable_exist(ID,'COLL_bands_frozen').or.write_is_on(ID)) then
   call def_variable_bulk(ID,"COLL_bands_frozen",1,shape(COLL_bands_frozen),LP_YIO)
   call  io_variable_bulk(ID,1,L1=COLL_bands_frozen)
 else
   if (allocated(COLL_bands_frozen)) COLL_bands_frozen=0
 endif
 !
 ! States
 !
 if (read_is_on(ID).and..not.allocated(COLL%state)) then
   YAMBO_ALLOC(COLL%state,(COLL%N,4))
 endif
 call def_variable_bulk(ID,'COLLISIONS_STATE',2,(/COLL%N,4/),IP_YIO)
 call  io_variable_bulk(ID,2,I2=COLL%state)
 !
 io_COLLISIONS_header=io_status(ID)
 !
1 call io_disconnect(ID=ID)
 !
 call timing('io_COLLISIONS',OPR='stop')
 !
end function io_COLLISIONS_header
