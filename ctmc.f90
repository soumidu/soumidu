Program orbit 
implicit none


integer,parameter:: n = 80000, p1 = 21600, p2 = 3000
real*8,parameter::x0=0.0d0,y0=0.0d0
real*8,parameter::t_max=400.d0
real*8,parameter::vz0=0.0d0,dt = 0.01d0,eip = 0.9d0,pi=acos(-1.0d0)
real*8,parameter:: dx = 0.005d0
integer::i,j,k,m
real*8::x,y,z,t,dt2,vx,vy,vz,E,E0
integer :: ix,iy,l,npoints
real*8 :: norm,part,threshold,ts,vs,dist,f,fp
real*8 ::l20,l2,lx0,ly0,lx,ly,E_fin,v0,z0,epeak,t0,Tp,vx0,dt0,dvx0,r
real*8 :: sigma, kappa, e1, F0, w, vx_fin, vy_fin, Fx, Fy, dy, vy0, Int
real*8 :: Fmax, re, beta
real*8 :: anm, v, Ax, Ay, const, vx_asymp, vy_asymp  
open(unit = 13,file="angle_new_4000_co.dat")


e1 = 1.0d0
Int = 4.5e14
F0 = dsqrt(Int/3.51E16)
Fmax = F0/dsqrt(1.d0+e1**2.d0)
w = 0.058d0
dy = 1.d0/dble(p2)

norm = 0.d0

!generation of initial conditions
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
 
  do ix = -p1, p1
     ts = ix*dx
     call sub2(ts,e1,F0,w,Fx)
     call sub3(ts,e1,F0,w,Fy)
     f = dsqrt( Fx**2.d0 + Fy**2.d0 )
     kappa = dsqrt(2.d0*eip)
     sigma = dsqrt(f/kappa)
     
    do iy = 0, p2
     
     vs = 3.d0*sigma - 6.d0*sigma*iy*dy
 
 
      
      
      call sub1(f,vs,dist)
      norm = norm + dx*dy*dist
   end do
   end do

  npoints = 25000

  threshold = 0.5d0/npoints
  part = 0.d0
  l = 0
 vx_fin = 0.0d0
 vy_fin = 0.0d0
  
  do ix = -p1,p1
     ts = ix*dx
     
  call sub2(ts,e1,F0,w,Fx)
    call sub3(ts,e1,F0,w,Fy)
     f = dsqrt( Fx**2.d0 + Fy**2.d0 )
     
    kappa = dsqrt(2.d0*eip)
    sigma = dsqrt(f/kappa)
  
    do iy = 0, p2

    vs = 3.d0*sigma - 6.d0*sigma*iy*dy


    call sub1(f,vs,dist)
     part = part + dx*dy*dist/norm
      if (part > threshold) then
        l= l+1  
t0 = ts

t = t0

     call sub2(t,e1,F0,w,Fx)
     call sub3(t,e1,F0,w,Fy)
     f = dsqrt( Fx**2.d0 + Fy**2.d0 )

  beta = 1.d0 - dsqrt(2.d0*eip)/2.d0
re = (eip + dsqrt( eip**2.d0 - 4.d0*beta*Fmax ))/( 2.d0*Fmax )

x = -(Fx/f)*re
y = -(Fy/f)*re
vx = -vs*Fy/f
vy = vs*Fx/f

z = 0.0d0

dt2=0.5d0*dt

!classical propagation of electron dynamics
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------

do while (t .le. t_max)


t = t + dt

call sub2(t,e1,F0,w,Fx)
call sub3(t,e1,F0,w,Fy)

x=x+dt2*vx
vx=vx+dt*( f_cx(x,y,z) - Fx )
x=x+dt2*vx

y=y+dt2*vy
vy=vy+dt*( f_cy(x,y,z) - Fy )
y=y+dt2*vy

enddo

if ( mod(l,1000) .eq. 0 ) write(*,*) l, vy_fin/dble(l), vx_fin/dble(l)

anm = x*vy - y*vx
r = Sqrt(x**2.d0 + y**2.d0)
v = Sqrt(vx**2.d0 + vy**2.d0-2.d0/r);

if (v .ge. 0.0d0) then
const = v/(v**2.d0*anm**2.d0+1.d0)
Ax = vy*anm - x/r
Ay = -vx*anm - y/r 
vx_asymp = -const*( Ax + Ay*anm*v )
vy_asymp = -const*( Ay - Ax*anm*v )
endif

vx_fin = vx_fin + vx_asymp
vy_fin = vy_fin + vy_asymp

threshold = threshold + 1.d0/npoints

end if
enddo 

enddo

write(13,*) vy_fin/l, vx_fin/l, abs(atan(vy_fin/vx_fin)), abs(atan(vy_fin/vx_fin))*180.d0/pi

contains

!description of Coulomb interaction
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
REAL FUNCTION f_cx(x,y,z)
implicit none
real*8,intent(in):: x,y,z

f_cx=-x/sqrt((x**2+y**2+z**2)**3)

end function f_cx

REAL FUNCTION f_cy(x,y,z)
implicit none
real*8,intent(in):: x,y,z

f_cy=-y/sqrt((x**2+y**2+z**2)**3)

end function f_cy   

REAL FUNCTION fz(x,y,z)
implicit none
real*8,intent(in):: x,y,z
real*8,parameter:: a=1.d0

fz=-z/sqrt((x**2+y**2+z**2)**3)

end function fz


end program orbit

!initial ionization rate
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
subroutine subc(x,y,dist)
    implicit none
    double precision, intent(in) :: x,y
    double precision :: dist
    dist = exp(-x**2-y**2+x*y)
   
end subroutine subc

 subroutine sub1(f,y,dist)
    implicit none
    real*8, intent(in) :: f,y
    real*8,intent(out) :: dist
real*8,parameter:: eip =  0.9d0, pi=acos(-1.0d0)
real*8:: dist1, dist2, kappa 
 
 kappa = dsqrt(2.d0*eip)
 
dist1 = exp(-2.d0*kappa**3.d0/(3.d0*f)) 
  
 dist2 = exp(-kappa*y**2.d0/f)

dist = dist1*dist2

 end subroutine sub1

!Definition of laser
!-----------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------
subroutine sub2(ts,e,F0,w,Fx)

implicit none

real*8,intent(in):: ts,e,F0,w
real*8,intent(out):: Fx
real*8::F_tilde, pw, w1
real*8,parameter:: pi = acos(-1.0d0)


F_tilde = F0/dsqrt(1.d0+e**2.d0)
w1 = w/4.d0
pw = pi/(2.0d0*w1)


              if ( abs(ts) .lt. pw ) then

              Fx = -F_tilde*cos(w1*ts)**3.d0*sin(5.d0*w*ts/4.d0)

              elseif (Abs(ts) .gt. pw) then

              Fx = 0.0d0

              endif
                    


end subroutine sub2


subroutine sub3(ts,e,F0,w,Fy)

implicit none

real*8,intent(in):: ts,e,F0,w
real*8,intent(out):: Fy
real*8::F_tilde, pw, w1
real*8,parameter:: pi = acos(-1.0d0)


F_tilde = -F0/dsqrt(1.d0+e**2.d0)
w1 = w/4.d0
pw = pi/(2.0d0*w1)


              if ( Abs(ts) .lt. pw ) then
              Fy = F_tilde*cos(w1*ts)**3.d0*cos(5.d0*w*ts/4.d0)
              elseif (Abs(ts) .gt. pw) then
              Fy = 0.0d0

              endif
                    


end subroutine sub3


