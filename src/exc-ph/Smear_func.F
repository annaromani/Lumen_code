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
 use pars,          ONLY:SP,cI
 use EXCPH,         ONLY:LORENTZIAN,GAUSSIAN,Smearing_type
 !
 complex(SP), intent(in) :: W
 real(SP),    intent(in) :: pole
 !
 ! Working Space
 !     
 complex(SP) :: Smear_func
 real(SP)    :: real_Wp,imag_Wp
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
   ! Only imaginary part is implemented at present
   !
   Smear_func=cI*sqrt(pi)*(exp(-real_Wp**2/(2.0*aimag_Wp**2))
   !
   ! Full Gaussian real and imaginart compatible with KK-relation (need Dawson function)
   ! Smear_func=2.0*ddaw(real_Wp)/(aimag_Wp*sqrt(2._SP))+cI*sqrt(pi)*(exp(-real_Wp**2/(2.0*aimag_Wp**2))
   !     
 endif
 !  
end function Smear_func
