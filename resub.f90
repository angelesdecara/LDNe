!
  ! program to generate r2_multi.txt in parts
  !
  !
  !
integer values(8)
integer N
integer ny
real*8 matriz(300,0:430), matriz_ant(300,0:430)
real*8 Nc, resul, Q, weir, feldman, Q1, weir1, feldman1, c
integer i,j,err_lectura, lunif
real*8 x1
real*8 eje_x(0:430), eje_y(210)

x1 = 0.432829272d+10
j = 0
do i = 4, 50, 2
  j = j + 1
  eje_y(j) = dble(i)
enddo 
do i = 60, 1000, 10
  j = j + 1
  eje_y(j) = dble(i)
enddo 
do i = 1100, 10000, 100
  j = j + 1
  eje_y(j) = dble(i)
enddo 
ny = j
!print *, ny
do i = 1, 430
  eje_x(i) = dble(i) / 1000
enddo
eje_x(0) = 0.5d0

do i=0,430
  print *, 'ejex ', i, eje_x(i)
enddo
do i=1,210
  print *, 'ejey ', i, eje_y(i)
enddo

!stop

matriz = 0
!do i = 1, ny
!  do j = 0, 430
do i = 1, 209
  if (eje_y(i) == 7600) then
!  do j = 0, 50
!  do j = 51, 100
!  do j = 101, 150
!  do j = 151, 200
!  do j = 201, 250
!  do j = 251, 300
!  do j = 301, 350
!  do j = 351, 400
  do j = 426, 430
    N = eje_y(i) 
    c = eje_x(j) 
    Nc = N * c
    call equil(N,Nc,resul,Q,weir,feldman,x1)
    matriz(i,j) = resul
    print '(a,2i10,2f12.5,f13.7)', 'prueba procesando ', i, j, eje_y(i), eje_x(j), resul
  enddo
  open (88,file = 'r2_multi.txt')
  do ii = 1, ny
    do jj = 0, 430
      read(88,'(3f20.10)') kk, ll, matriz_ant(ii,jj)
    enddo
  enddo
  close(88)
  !matriz_ant = 0.0d0
  do ii = 1, ny
    do jj = 0, 430
      if (matriz(ii,jj) > 0.000000001) matriz_ant(ii,jj) = matriz(ii,jj) 
    enddo
  enddo  
  open (88,file = 'r2_multi.txt')
  do ii = 1, ny
    do jj = 0, 430
      write(88,'(3f20.10)') eje_y(ii), eje_x(jj), matriz_ant(ii,jj)
    enddo
  enddo
  close(88)
  end if
enddo


stop
end

subroutine equil(N,Nc,resul,Q,weir,feldman,x1)
implicit none
integer*8 k, i, j, ii, jj, it, niter
integer N, N2
integer*8 icount
integer, allocatable, dimension(:) :: icaso2
integer, allocatable, dimension(:,:) :: icaso
real*8, allocatable, dimension(:) :: frecs
real*8, allocatable, dimension(:) :: Vlin
integer*8, allocatable, dimension(:) :: yv
integer*8, allocatable, dimension(:,:) :: xv
real*8, allocatable, dimension(:) :: V
real*8, allocatable, dimension(:) :: facto
real*8 r8_gamma_log, calc_facto, sum, c, r2tot, xtot, Q, Qt, Nc, resul, t, val
real*8 p1sum, p2sum, dt, x1, unif, r2, d, p1, p2, weir, feldman
integer x, x2  !x es ahora el caso en el que se esta en cada generacion
integer nchains, nburn, lunif, bueno

real*8 ff(4)
integer nnp(4)
N2 = 2*N
!x1 = 0.671927838d+10
!open(88, file='pepi',access='sequential', form='unformatted', status='replace')
icount = 0
k = 4
c = Nc / N

niter = 2050000
!niter =  100000
nburn =   50000
nchains = 100

allocate (frecs(k), icaso(k,nchains), icaso2(k))

!comienzo en equilibrio
do i = 1, nchains
  icaso(1,i) = N2/4
  icaso(2,i) = N2/4
  icaso(3,i) = N2/4
  icaso(4,i) = N2 - icaso(1,i) - icaso(2,i) - icaso(3,i)
enddo
frecs = 0.25d0
!print*,icaso,N2
!==============================================================

sum = 0.0d0
do ii = 1, niter
  !if (mod(ii,100000) == 0) print *, niter,ii 
  do jj = 1, nchains
    !calcular frecuencias de los gametos tras recombinacion
    do i=1,k
      frecs(i) = dble(icaso(i,jj))/(N2)
      !print *,i,icaso(i),N2, frecs(i)
    enddo
    !calcular r2 de la muestra 
    call calc_d(frecs,p1,p2,d)
    r2 = (d*d) / (p1*(1.0d0-p1)*p2*(1.0d0-p2))
    if (ii > nburn) then
      sum = sum + r2
      resul = sum / (nchains*(ii-nburn))
      !if(mod(ii,100000)==0)print '(i10,f12.4,4i10)',ii,r2,icaso
    end if
    !calcular frecuencias de los gametos tras recombinacion
    !do i=1,k
    !  frecs(i) = dble(icaso(i))/(N2)
    !  !print *,i,icaso(i),N2, frecs(i)
    !enddo
    !print '(a,4f12.4)','frecs : ',frecs
    !call calc_d(frecs,p1,p2,d)
    frecs(1) = frecs(1) - d * c 
    frecs(2) = frecs(2) + d * c 
    frecs(3) = frecs(3) + d * c 
    frecs(4) = frecs(4) - d * c  
    ! muestrear
77  call multidev(frecs,icaso2,N2,x1)
    bueno = 1
    if(icaso2(1)+icaso2(2) ==  0) bueno = 0
    if(icaso2(1)+icaso2(3) ==  0) bueno = 0
    if(icaso2(1)+icaso2(2) == N2) bueno = 0
    if(icaso2(1)+icaso2(3) == N2) bueno = 0
    if (bueno == 1) then
      icaso(:,jj) = icaso2
    else
      it = lunif(x1, nchains)
      icaso(:,jj) = icaso(:,it)
    end if
    !if(mod(ii,1)==0)print '(i10,4i10)',ii,icaso2
  enddo
enddo

Q = (1.0d0-c/2)/(1.0d0-c)**2
Q = 1.0d0/(1.0d0+4.0d0*Nc*Q)
weir = c**2+(1.0d0-c)**2
weir = weir / (2.0d0*Nc*(2.0d0-c))
feldman = 1.0d0/(1.0d0+4.0d0*Nc-2.0d0*c-2.0d0*Nc*c+c*c)

deallocate (icaso, icaso2)
deallocate (frecs)
!close(88)

return
end  

function lunif(x1,m)
real*8 x1,unif
integer lunif, m
x1 = mod (16807.0d0 * x1, 2147483647.0d0)
unif = x1/2147483647.0d0
lunif = unif * m + 1
return
end function lunif

function unif(x1)
real*8 x1,unif
x1 = mod (16807.0d0 * x1, 2147483647.0d0)
unif = x1/2147483647.0d0
return
end function unif

function calc_facto(icaso,k)
integer icaso(k)
real*8 calc_facto, r8_gamma_log
calc_facto=r8_gamma_log(dble(icaso(1)+icaso(2)+icaso(3)+icaso(4)+1))
do i = 1, 4
  calc_facto=calc_facto-r8_gamma_log(dble(icaso(i)+1))
enddo
return
end

subroutine next(k,icaso)
integer icaso(k),i,k,j
do i=k,1,-1
  if(icaso(i).gt.0)then
    j = icaso(i)
    icaso(i-1) = icaso(i-1) + 1
    icaso(i) = 0
    icaso(k) = j - 1
    return
  end if
enddo
return
end

!=================================================================
function multi(frecs, x1)
!extrae al azar un haplotipo de una multinomial de 4 frecs haplotipicas
real*8 frecs(4), x1, u, unif
integer multi
u = unif(x1)
if (u .lt. frecs(1)) then
  multi = 1
  return
end if 
if (u .lt. frecs(1)+frecs(2)) then
  multi = 2
  return
end if 
if (u .lt. frecs(1)+frecs(2)+frecs(3)) then
  multi = 3
  return
end if 
multi = 4
return
end
!=================================================================
subroutine calc_frecs(frecs,p1,p2,d)
! calcula frecuencias haplotipicas en funcion de p1, p2 y d
real*8 frecs(4), p1, p2, d
! frecs (1) es la frecuencia de AB
! frecs (2) es la frecuencia de Ab
! frecs (3) es la frecuencia de aB
! frecs (4) es la frecuencia de ab
frecs (1) = p1 * p2 + d
frecs (2) = p1 * (1.0d0 - p2) - d
frecs (3) = (1.0d0 - p1) * p2 - d
frecs (4) = (1.0d0 - p1) * (1.0d0 - p2) + d
if (abs(1.0d0 - frecs(1) - frecs(2)- frecs(3) - frecs(4)).gt.0.000001)then
   print *, 'error en calc_frecs'
   stop
end if
return
end
!=================================================================
subroutine calc_d(frecs,p1,p2,d)
! calcula frecuencias alelicas y d en funcion de frecuencias haplotipicas
real*8 frecs(4), p1, p2, d
! frecs (1) es la frecuencia de AB
! frecs (2) es la frecuencia de Ab
! frecs (3) es la frecuencia de aB
! frecs (4) es la frecuencia de ab
p1 = frecs(1) + frecs(2)
p2 = frecs(1) + frecs(3)
d = frecs(1) - p1 * p2
if (d .lt. 0.0d0 .or. d .gt. 1.0d0)then
   !print *, 'error en calc_d ', d
   !stop
end if
return
end
!==================================================================================
subroutine multidev(frecs,n_par,n,x1)
! muestreo de multinomial
real*8 frecs(4), x1, unif, bnldev
integer n_par(4), n, m
real*8 ft(4)
m = n
do i = 1, 4
  ft(i) = frecs(i)
enddo
n_par(1) = bnldev(ft(1),m,x1)
do i = 2, 3
  m = m - n_par(i-1)
  do j = i, 4
    ft(j) = ft(j) / (1.0d0 - ft(i-1))
  enddo
  n_par(i) = bnldev(ft(i),m,x1)
enddo
n_par(4) = n - n_par(1) - n_par(2) - n_par(3)
return
end
!==================================================================================
function bnldev(pp,n,x1)
INTEGER n
real*8 x1, unif
REAL*8 bnldev,pp,PI
INTEGER j,nold
REAL*8 am,em,en,g,oldg,p,pc,pclog,plog,pold,sq,t,y,gammln, alngam
SAVE nold,pold,pc,plog,pclog,en,oldg
DATA nold /-1/, pold /-1./
PI = acos(-1.0d0)
if(pp .le. 0.5d0)then
    p = pp
  else
    p = 1.0d0-pp
endif
am = n * p
if (n .lt. 25)then
  bnldev = 0.0d0
  do j = 1, n
    if (unif(x1) .lt. p) bnldev = bnldev + 1.0d0
  enddo
else if (am .lt. 1.0d0) then
  g = exp(-am)
  t=1.0d0
  do j = 0, n
    t = t * unif(x1)
    if (t.lt.g) go to 1
  enddo
  j = n
1 bnldev = j
else
  if (n.ne.nold) then
    en = n
    oldg = alngam(en+1.0d0)
    nold = n
  endif
  if (p .ne. pold) then
  pc = 1.0d0 - p
  plog = log(p)
  pclog = log(pc)
  pold = p
endif
sq = sqrt(2.0d0 * am * pc)
2 y = tan(PI * unif(x1))
em = sq * y + am
if (em .lt. 0.0d0 .or. em .ge. en+1.0d0) goto 2
em = int(em)
t = 1.2d0 * sq * (1.0d0 + y**2) * exp(oldg - alngam(em + 1.0d0) &
  - alngam (en - em + 1.0d0) + em * plog + (en - em) * pclog)
if (unif(x1) .gt. t) goto 2
  bnldev = em
endif
if (p .ne. pp) bnldev = n - bnldev
return
END

!========================================================
function gammln(xx)
real*8 gammln,xx
integer j
real*8 ser,stp,tmp,x,y,cof(6)
SAVE cof,stp
DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, 24.01409824083091d0, &
             -1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,  &
              2.5066282746310005d0/
x = xx
y = x
tmp = x + 5.5d0
tmp = (x + 0.5d0) * log(tmp) - tmp
ser = 1.000000000190015d0
do j=1,6
  y = y + 1.d0
  ser=ser+cof(j)/y
enddo
gammln = tmp + log(stp * ser/x)
return
END




FUNCTION alngam(xvalue) RESULT(fn_val)

!     ALGORITHM AS245  APPL. STATIST. (1989) VOL. 38, NO. 2

!     Calculation of the logarithm of the gamma function

! ELF90-compatible version by Alan Miller
! Latest revision - 29 November 1997

! N.B. Argument IFAULT has been removed

IMPLICIT NONE
INTEGER, PARAMETER    :: dp = SELECTED_REAL_KIND(15, 100)
REAL (dp), INTENT(IN) :: xvalue
REAL (dp)             :: fn_val

! Local variables
REAL (dp) :: x, x1, x2, y

!     Coefficients of rational functions

REAL (dp), PARAMETER :: r1(9) = (/ -2.66685511495_dp, -24.4387534237_dp,  &
                                   -21.9698958928_dp,  11.1667541262_dp,  &
                                    3.13060547623_dp,  0.607771387771_dp, &
                                    11.9400905721_dp,  31.4690115749_dp,  &
                                    15.2346874070_dp /)
REAL (dp), PARAMETER :: r2(9) = (/ -78.3359299449_dp, -142.046296688_dp,  &
                                    137.519416416_dp,  78.6994924154_dp,  &
                                    4.16438922228_dp,  47.0668766060_dp,  &
                                    313.399215894_dp,  263.505074721_dp,  &
                                    43.3400022514_dp /)
REAL (dp), PARAMETER :: r3(9) = (/ -2.12159572323E5_dp,  2.30661510616E5_dp,  &
                                    2.74647644705E4_dp, -4.02621119975E4_dp,  &
                                   -2.29660729780E3_dp, -1.16328495004E5_dp,  &
                                   -1.46025937511E5_dp, -2.42357409629E4_dp,  &
                                   -5.70691009324E2_dp /)
REAL (dp), PARAMETER :: r4(5) = (/ 0.279195317918525_dp, 0.4917317610505968_dp, &
                                   0.0692910599291889_dp, 3.350343815022304_dp, &
                                   6.012459259764103_dp /)

!     Fixed constants

REAL (dp), PARAMETER :: alr2pi = 0.918938533204673_dp, four = 4._dp,  &
                        half = 0.5_dp, one = 1._dp, onep5 = 1.5_dp,   &
                        twelve = 12._dp, zero = 0._dp

!     Machine-dependant constants.
!     A table of values is given at the top of page 399 of the paper.
!     These values are for the IEEE double-precision format for which
!     B = 2, t = 53 and U = 1023 in the notation of the paper.

REAL (dp), PARAMETER :: xlge = 5.10E6_dp, xlgst = HUGE(1.0_dp)

x = xvalue
fn_val = zero

!     Test for valid function argument

IF (x >= xlgst) THEN
  WRITE(*, *) 'AS 245: Argument x too large'
  RETURN
END IF
IF (x <= zero) THEN
  WRITE(*, *) 'AS 245: Argument x <= 0'
  RETURN
END IF

!     Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined

IF (x < onep5) THEN
  IF (x < half) THEN
    fn_val = -LOG(x)
    y = x + one
    
!     Test whether X < machine epsilon
    
    IF (y == one) RETURN
  ELSE
    fn_val = zero
    y = x
    x = (x - half) - half
  END IF
  fn_val = fn_val + x * ((((r1(5)*y + r1(4))*y + r1(3))*y + r1(2))*y + r1(1)) / &
                    ((((y + r1(9))*y + r1(8))*y+ r1(7))*y + r1(6))
  RETURN
END IF

!     Calculation for 1.5 <= X < 4.0

IF (x < four) THEN
  y = (x - one) - one
  fn_val = y * ((((r2(5)*x + r2(4))*x + r2(3))*x + r2(2))*x + r2(1)) /  &
               ((((x + r2(9))*x + r2(8))*x + r2(7))*x+ r2(6))
  RETURN
END IF

!     Calculation for 4.0 <= X < 12.0

IF (x < twelve) THEN
  fn_val = ((((r3(5)*x + r3(4))*x + r3(3))*x + r3(2))*x + r3(1)) /  &
           ((((x + r3(9))*x + r3(8))*x + r3(7))*x + r3(6))
  RETURN
END IF

!     Calculation for X >= 12.0

y = LOG(x)
fn_val = x * (y - one) - half * y + alr2pi
IF (x > xlge) RETURN
x1 = one / x
x2 = x1 * x1
fn_val = fn_val + x1 * ((r4(3)*x2 + r4(2))*x2 + r4(1)) /  &
         ((x2 + r4(5))*x2 + r4(4))
RETURN
END FUNCTION alngam

include 'subru.f90'

