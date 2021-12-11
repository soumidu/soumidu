Program orbit 
implicit none

!parameters
!----------------------------------------------------------------------------------------------------------
integer,parameter:: n = 80000, p1 = 27082, p2 = 3000, lm = 0
real*8,parameter:: cycles = 5.d0
real*8,parameter::x0=0.0d0,y0=0.0d0
real*8::t_max = 300.d0
real*8,parameter::vz0=0.0d0,dt = 0.1d0,eip = 0.9d0,pi=acos(-1.0d0)
real*8,parameter:: dx = 0.01d0,cl = 2.67d0,din = 0.5d0
integer::i,j,k,m
real*8::x,y,z,t,dt2,vx,vy,vz,E,E0,W0
integer :: ix,iy,l,npoints
real*8 :: norm,part,threshold,ts,vs,dist,f,rho
real*8 ::l20,l2,lx0,ly0,lx,ly,E_fin,v0,z0,epeak,t0,Tp,vx0,dt0,dvx0,r
real*8 :: sigma, kappa, e1, F0, w, vx_fin, vy_fin, Fx, Fy, dy, vy0, Int
real*8 :: Fmax, re, beta, x_int, traj
real*8 :: anm, v, Ax, Ay, const, vx_asymp, vy_asymp

open(unit = 28,file="traj_spm.txt")


e1 = 0.0d0
x_int = 3.d0
Int = x_int*1e14
F0 = dsqrt(Int/3.51E16)
Fmax = F0/dsqrt(1.d0+e1**2.d0)
w = 0.058d0
dy = 1.d0/dble(p2)

!initialization of dynamics
!----------------------------------------------------------------------------------------------------------
norm = 0.d0
 
t0 = -100.d0
t = t0

  call  sub2(ts,e1,F0,w,cycles,Fx)
  call    sub3(ts,e1,F0,w,cycles,Fy)
      f = dsqrt(Fx**2.d0 + Fy**2.d0)
  beta = 1.d0 - dsqrt(2.d0*eip)/2.d0
re = (eip + dsqrt( eip**2.d0 - 4.d0*beta*Fmax ))/( 2.d0*Fmax )

x = -(Fx/f)*re
y = -(Fy/f)*re
vx = -vs*Fy/f
vy = vs*Fx/f

z = 0.0d0

dt2=0.5d0*dt

!classical propagation in the presence of electro-magnetic fields
!-----------------------------------------------------------------------------------------------------------
do while (t .le. t_max)


t = t + dt

call sub2(t,e1,F0,w,cycles,Fx)
call sub3(t,e1,F0,w,cycles,Fy)

x=x+dt2*vx
vx=vx+dt*( f_cx(x,y,z) - Fx )
x=x+dt2*vx

y=y+dt2*vy
vy=vy+dt*( f_cy(x,y,z) - Fy )
y=y+dt2*vy


!analytical free-particle solution
!----------------------------------------------------------------------------------------------------------
traj =  (F0/w)*(sin(w*t0)*(t - t0) + (cos(w*t) - cos(w*t0))/w)
!priniting only free particle solusion
!----------------------------------------------------------------------------------------------------------
write(28,*) traj,t

enddo



contains

!Coulomb interaction
!----------------------------------------------------------------------------------------------------------
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


!definition of laser
!----------------------------------------------------------------------------------------------------------
subroutine sub2(t,e,F0,w,cycles,Fx)

implicit none


real*8,intent(in):: t,e,F0,w,cycles
real*8,intent(out):: Fx
real*8::F_tilde, pw, w1
real*8,parameter:: pi = acos(-1.0d0)




F_tilde = F0/dsqrt(1.d0+e**2.d0)
w1 = w/(2.d0*cycles)
pw = pi/(2.0d0*w1)


             if ( abs(t) .le. pw ) then

              Fx = -F_tilde*cos(w1*t)**3.d0*e*sin(5.d0*w*t/4.d0)

             elseif (Abs(t) .ge. pw) then

             Fx = 0.0d0

             endif
                    


end subroutine sub2



subroutine sub3(t,e,F0,w,cycles,Fy)

implicit none


real*8,intent(in):: t,e,F0,w, cycles
real*8,intent(out):: Fy
real*8::F_tilde, pw, w1
real*8,parameter:: pi = acos(-1.0d0)




F_tilde = F0/dsqrt(1.d0+e**2.d0)
w1 = w/(2.d0*cycles)
pw = pi/(2.0d0*w1)


             if ( Abs(t) .le. pw ) then
              Fy = -F_tilde*cos(w1*t)**3.d0*cos(5.d0*w*t/4.d0)
             elseif (Abs(t) .ge. pw) then
             Fy = 0.0d0

             endif
                    


end subroutine sub3


