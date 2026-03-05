!
!        Copyright (C) 2000-2021 the YAMBO team
!              http://www.yambo-code.org
!
! Authors (see AUTHORS file for details): AC
!
! headers
!
#include<y_memory.h>
!
function Smear_func(W,pole)
 use units,         ONLY:HA2EV,pi
 use pars,          ONLY:SP,DP,cI
 use EXCPH,         ONLY:LORENTZIAN,GAUSSIAN,Smearing_type
 !
 implicit none
 !
 complex(SP), intent(in) :: W
 real(SP),    intent(in) :: pole
 !
 ! Working Space
 !     
 complex(SP) :: Smear_func
 real(SP)    :: real_Wp,imag_Wp
 complex(DP) :: daw
 !
 if(Smearing_type==LORENTZIAN) then
   !
   Smear_func=1._SP/(W-pole)
   !
 elseif(Smearing_type==GAUSSIAN) then
   !
   real_Wp=real(W-pole)
   imag_Wp=aimag(W-pole)
   !
   ! Full Gaussian real and imaginart compatible with KK-relation (need Dawson function)
!   Smear_func=2._SP*daw(dble(real_Wp))/(aimag_Wp*sqrt(2._SP))+cI*sqrt(pi)*(exp(-real_Wp**2/(2._SP*aimag_Wp**2)))
   Smear_func=-cI/sqrt(imag_Wp**2)*(exp(-real_Wp**2/(2._SP*imag_Wp**2)))
   !     
 endif
 !  
end function Smear_func
