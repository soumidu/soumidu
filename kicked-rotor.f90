Program orbit 
implicit none

integer:: n, nmax, nrand, i
real*8:: q, p
real*8, parameter:: pi = 3.1415926535d0, K = 2.35d0


open(unit = 25,file="sm.dat")


nmax = 500
nrand = 100

call random_number(q)
call random_number(p)

do i = 1, nrand 

q = q*2.d0*pi
p = p*2.d0*pi


do n = 1 , nmax


p = p + K*sin(q)

if (p .gt. 0.0d0) then
p = mod(p, 2.d0*pi)
elseif (p .lt. 0.0d0) then
p = 2.d0*pi - abs(mod(p, 2.d0*pi))
endif


q = q + p

if (q .gt. 0.0d0) then 
q = mod(q, 2.d0*pi)
elseif (q .lt. 0.0d0) then
q = 2.d0*pi - abs(mod(q, 2.d0*pi))
endif

if (p .gt. pi) then 
p = p - 2.d0*pi 
else
p = p
endif
write(25,*) q, p

enddo 
write(25,*) 
enddo

end Program orbit
