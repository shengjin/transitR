program transitr

implicit none

real*8   Rearth,Mearth,tau,gacc,Tint,Teq,kv2kth,Temp,Pressure
real*8   Hscale,Rplanet,Mplanet,Rdiff
integer*8            iline
integer*8            line
DOUBLE PRECISION            ::mH2,mHe,Y,uma,MoleculeW
DOUBLE PRECISION,PARAMETER    ::M_inearth=5.97424D27,R_inearth=6.371D8

DOUBLE PRECISION,PARAMETER    ::G=6.67D-8

!DOUBLE PRECISION,PARAMETER    ::km=1D5,RJ=7.15D4*km,Lsun=3.827D33,G=6.67D-8,Rearth=6.371D8
!DOUBLE PRECISION,PARAMETER    ::LJ=3.35D24
!LJ=8.67D-10*Lsun,gives LJ=3.318D24.  Guillot & Gautier give 3.35D24


open(unit=10,file='LineRearthPTMearthTauTintTeq.10.dat')
open (unit=20,File='LineRearthMearth_Rdiff.dat',ACCESS='append')

do iline=1,10
read(10,*) line,Rearth,Pressure,Temp,Mearth,tau,Tint,Teq

Rplanet=Rearth*R_inearth
Mplanet=Mearth*M_inearth
gacc=G*Mplanet/Rplanet**2D0

CALL kv2kthTeq(kv2kth,Teq)

uma=1.660531e-24 
mH2 = 2.*1.003 * uma
mHe = 4.008 * uma
Y=0.24
MoleculeW=mHe*Y+mH2*(1-Y)

call ScaleHight(Temp,gacc,MoleculeW,Hscale)

write(*,*)'Rplanet,Hscale,kv2kth,tau',Rplanet/R_inearth,Hscale/R_inearth,kv2kth,tau

call RTSEC(Hscale,Rplanet,kv2kth,tau,Rdiff)
write(*,*)'Rdelta',Rdiff/R_inearth

write(20,*) line,Rearth,Mearth,Rdiff/R_inearth

end do 

close(10)
close(20)


end program transitr



subroutine kv2kthTeq(kv2kth,Teq)
implicit none
Double Precision,INTENT(IN)::Teq
double Precision,INTENT(OUT)::kv2kth
Double Precision ::kv2kthtlow,kv2kthtup,Tloc
INTEGER:: iT,iTup,iTlow
LOGICAL,PARAMETER ::DEBUG=.FALSE.

Double Precision,DIMENSION(12),PARAMETER :: Tvect=(/&
260,388,584,861,1267,1460,1577,1730,1870,2015,2255,2777/)

Double Precision,DIMENSION(12),PARAMETER :: kv2kthvect=(/&
0.005,0.008,0.027,0.07,0.18,0.19,0.18,0.185,0.2,0.22,0.31,0.55/)

Tloc=Teq

!1)T>Tmax: 
IF(Tloc>=Tvect(12))THEN
   kv2kth=(Tloc/Tvect(12))**2.0*0.55
   IF(DEBUG)WRITE(*,*)'T<Tmin: scaling as T**2'
   RETURN
END IF

!2)T<Tmin:
IF(Tloc<=Tvect(1))THEN
   kv2kth=(Tloc/Tvect(1))*0.008
   IF(DEBUG)WRITE(*,*)'T<Tmin: scaling as T'
   RETURN
END IF

iT=1
do while (Tvect(iT).lt.Tloc)
   iT=iT + 1
enddo

iTup=iT
iTlow=iT-1
IF(DEBUG)WRITE(*,*)'iTlow,Tloc,iTup',iTlow,Tloc,iTup
IF(DEBUG)WRITE(*,*)'Tvect(iTlow),Tvect(iTup)',Tvect(iTlow),Tvect(iTup)
IF(DEBUG)WRITE(*,*)'kv2kthvect(iTlow),kv2kthvect(iTup)',kv2kthvect(iTlow),kv2kthvect(iTup)

!interpolate linearly in the kappa matrix

call LinIntPol(1,Tloc,kv2kthvect(iTlow),Tvect(iTlow),kv2kthvect(iTup),Tvect(iTup),kv2kth) 

IF(DEBUG)WRITE(*,*)'kv2kth,Teq', kv2kth,Tloc

RETURN
end subroutine kv2kthTeq




SUBROUTINE LinIntPol(n,x,y1,x1,y2,x2,y)
!makes linear interpolation for f(x) (Dim n) at the point x, if (x1,f(x1)=y1) and (x2,f(x2)=y2)
!are given (dimension y1, y2 also n of course) using Lagrange's formula
IMPLICIT NONE
INTEGER                         ::n
DOUBLE PRECISION                ::y1,y2,y
!DOUBLE PRECISION,DIMENSION(n)   ::y1,y2,y
DOUBLE PRECISION                ::x,x1,x2


y=((x-x2)/(x1-x2))*y1+((x-x1)/(x2-x1))*y2

END SUBROUTINE LinIntPol



SUBROUTINE ScaleHight(Teff,gacc,MoleculeW,Hscale)
IMPLICIT NONE
Double Precision,INTENT(IN)::Teff,gacc
double Precision,INTENT(OUT)::Hscale
DOUBLE PRECISION            ::MoleculeW,kB

kB=1.38062D-16

Hscale=kB*Teff/(gacc*MoleculeW)

END SUBROUTINE ScaleHight





subroutine RTSEC(Hscale,Rplanet,kv2kth,tau,Rdiff)
  IMPLICIT NONE
  INTEGER :: ISTEP
  REAL*8 :: x1, x2, xacc, DX, Func1, Func2, temp
  REAL*8 :: Hscale,Rplanet,kv2kth,tau,Rdiff
  xacc = 1.0E-06
  !x1  = 0.1*Rplanet
  x1  = -1.0E2*Rplanet
  x2  = 1.0E2*Rplanet
  DX = x1 - x2
  ISTEP = 0
  DO WHILE (ABS(DX).GT.xacc)
    Rdiff = (x1+x2)/2.0
!    write(*,*)'x1,x2,Rdiff',x1,x2,Rdiff
    call Func(x1,Hscale,Rplanet,kv2kth,tau,temp)
    Func1=temp
!    write(*,*)'Func1',Func1
    call Func(Rdiff,Hscale,Rplanet,kv2kth,tau,temp)
    Func2=temp
!    write(*,*)'Func2',Func2
!    read(*,*)
    IF ((Func1*Func2).LT.0) THEN
      x2  = Rdiff
      DX = x2-x1
    ELSE
      x1  = Rdiff
      DX = x2-x1
    END IF
    ISTEP = ISTEP+1
  END DO
end subroutine RTSEC



subroutine Func(x,Hscale,Rplanet,kv2kth,tau,temp)
IMPLICIT NONE
  REAL*8 :: temp, x
  REAL*8 :: pi
  REAL*8 :: Hscale,Rplanet,kv2kth,tau,Rdiff
  pi=3.14159265358
  !temp=Hscale*pi*kv2kth*tau*Rplanet/Hscale-x/2
  temp=Hscale*log(kv2kth*tau*(2*pi*(Rplanet+x)/Hscale)**0.5D0)-x
!  write(*,*)'x,Rplanet,Hscale,kv2kth,tau,temp',x,Rplanet,Hscale,kv2kth,tau,temp
end subroutine Func
