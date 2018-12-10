! program to generate r^2 as in Fig 1 of Sved 1971
! nind individuals, 1000 loci over 100 cM or 10000 over 1000cM
!
! BP at nloci/2-nloci/2+1 for 2 chromosomes
! at nloci/10 for 10
!
! save r2 every 100 generations, as it's very time consuming
program r2Sved
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: i,j,k,l, gametetemp, n, igen, irep
  integer, parameter :: nind=5000, ngen=5000, nloci=10000, nrep=1, ngam=2
  real(dp), parameter :: length=10
  integer :: genotipo(nind,2,nloci), gen1(2*nind),gen2(2*nind)
  integer :: genotipotemp(nind,2,nloci), locusstart, geneal
  integer :: crossover, poisson, contador(0:nloci), cruz(nloci), ncross
  integer::ldcontador(0:nloci),kl
  real (dp):: averr2(0:ngen,0:nloci), r2dummy(0:nloci), r2, dprime, pk,pl
  real(dp):: ldr2dummy(0:nloci),ldr2,noDr2, noDr2dummy(0:nloci) !! noD r2 is ~correlation
  !! i.e. denominatori of r2 correlation is pq_A pq_B instead of product of variances
  real(dp)::varr2(0:ngen,0:nloci),ldvarr2(0:ngen,0:nloci),ldaverr2(0:ngen,0:nloci)
  real(dp)::noDvarr2(0:ngen,0:nloci),noDaverr2(0:ngen,0:nloci)
  ! ld are r2 ld based, while without ld are correlatoon based measures
  real:: ranmar
  character(len=10)::aind,alength,afix
  character(len=20)::aformat,agen,aformat2
  call rmarin(6589,25468)

  write(alength,'(i3)')int(length)
  write(aind,'(i6)')nind
  write(agen,'(i6)')ngen+1
  write(aformat2,'(i6)')(ngen+1)*3
  afix='N'//trim(adjustl(aind))//'L'//trim(adjustl(alength))
  aformat='(i5,'//trim(adjustl(agen))//'f12.7)'
  aformat2='(i5,'//trim(adjustl(aformat2))//'f12.7)'



  !//'.dat'
!
! start from random...
!
  contador =0; averr2=0.d0; varr2=0.d0; noDr2=0.d0; noDaverr2=0.d0;noDvarr2=0.d0
  ldcontador =0; ldaverr2=0.d0; ldvarr2=0.d0; noDr2dummy=0.d0
  do irep=1,nrep

     do i=1,nind
        do k=1,ngam
           do l=1,nloci
              genotipo(i,k,l)=ranmar()*2+1
           enddo
        enddo
     enddo
!
! at t=0
!

     r2dummy=0; contador = 0; ldr2dummy=0
     do l=1,nloci-1

        do k=l+1,nloci
           r2=0.d0; ldr2=0.d0; noDr2=0.d0
!
! first, test its not fixed
!
           pk=count(genotipo(:,:,k)==1)/dble(2*nind)
           pl=count(genotipo(:,:,l)==1)/dble(2*nind)
           if ((pk.gt.1.d-3).and.(pk.lt.(1-1.d-3)).and.((pl.gt.1.d-3).and.(pl.lt.(1-1.d-3))))then
              gen1(1:nind)=genotipo(:,1,k)
              gen1(1+nind:2*nind)=genotipo(:,2,k)
!
              gen2(1:nind)=genotipo(:,1,l)
              gen2(1+nind:2*nind)=genotipo(:,2,l)

              kl=k-l
              if(k/1000.ne.l/1000)kl=0
              contador(kl) = contador(kl)+1
              call mkdprime (dprime, r2,ldr2,noDr2, gen1, gen2,2*nind)

              r2dummy(kl)=r2dummy(kl)+r2
              ldr2dummy(kl)=ldr2dummy(kl)+ldr2
              noDr2dummy(kl)=noDr2dummy(kl)+noDr2
        
           endif
        enddo
     enddo

     do l=0,nloci/10
!        contador(l) = contador(l)+1
        averr2(0,l)=averr2(0,l)+r2dummy(l)/dble(contador(l))
        ldaverr2(0,l)=ldaverr2(0,l)+ldr2dummy(l)/dble(contador(l))
        ! don't think I need sum+=
        noDaverr2(0,l)=noDr2dummy(l)/dble(contador(l))
        !           noDaverr2(igen,l)=noDr2dummy(l)/dble(contador(l))
     enddo

!                                                                                                              
! neutral evolution for ngen generations                                                                       
!                                                                                                              
     do igen=1,ngen
        !if(mod(igen,100).eq.0)
        write(6,*)'gen ',igen
        do i=1,nind
        ! select parental gametes                                                                              
           do k=1,2
              geneal=ranmar()*nind+1
!              tt=0.d0                                                                                         

              ncross=poisson(length)!crossover(length*dble(ranmar()))                                          
              do while(ncross.gt.nloci)
                 ncross=poisson(length)!!crossover(dble(ranmar()))                                             
              enddo
              cruz=0
              do l=1,ncross
                 cruz(l)=int(ranmar()*nloci)+1
              enddo
              do l=1,9 ! 10 chr by forcing breakpoints at nloci*l/10 every generation
                 if(ranmar().lt..5)then
                    ncross=ncross+1
                    cruz(ncross)=nloci*l/10.d0
                 endif
              enddo
              call sortincrease(cruz,ncross)
! Copy with overcrossing at positions given by cruz(ncross)                                                    
              locusstart=1
              gametetemp=int(ranmar()*ngam)+1
!                                                                                                              
              genotipotemp(i,k,:)=genotipo(geneal,k,:)
              do n=1,ncross
!                      write(6,*)i,'crossing',n,k,gametetemp                                                   
                 genotipotemp(i,k,locusstart:cruz(n))=&
                      genotipo(geneal,gametetemp,locusstart:cruz(n))
                 locusstart=cruz(n)+1
                 gametetemp=3-gametetemp
              enddo
! out of the loop at the last overcrossing, let's copy the remainder of the chromosome                         
              genotipotemp(i,k,locusstart:nloci)=&
                   genotipo(geneal,gametetemp,locusstart:nloci)
           enddo
        enddo ! all guys done                                                                                  
        genotipo=genotipotemp

!        if(mod(igen,100).eq.0)then ! to save time                                                             
        r2dummy=0; contador = 0; ldr2dummy=0; noDr2dummy=0.d0
        do l=1,nloci-1
           do k=l+1,nloci
              r2=0.d0;ldr2=0.d0
!                                                                                                              
! first, test its not fixed                                                                                    
!                                                                                                              
              if ((count(genotipo(:,:,k)==1).ne.2*nind) &
                  .and.(count(genotipo(:,:,k)==1).ne.0) &
                  .and.(count(genotipo(:,:,l)==1).ne.2*nind) &
                  .and.(count(genotipo(:,:,l)==1).ne.0))then

                 gen1(1:nind)=genotipo(:,1,k)
                 gen1(1+nind:2*nind)=genotipo(:,2,k)
!                                                                                                              
                 gen2(1:nind)=genotipo(:,1,l)
                 gen2(1+nind:2*nind)=genotipo(:,2,l)
!             
                 kl=k-l
                 if(k/1000.ne.l/1000)kl=0
!                 if((k.gt.nloci/2).and.(l.le.nloci/2))kl=0
                 contador(kl) = contador(kl)+1
! save results every 100 gens
                 if(mod(igen,100).eq.0)then
                    call mkdprime (dprime, r2,ldr2,noDr2, gen1, gen2,2*nind)
                 endif
                r2dummy(kl)=r2dummy(kl)+r2
                 ldr2dummy(kl)=ldr2dummy(kl)+ldr2
                 noDr2dummy(kl)=noDr2dummy(kl)+noDr2
              endif
           enddo
        enddo


        do l=0,nloci/10
!                                                                                                              
           averr2(igen,l)=r2dummy(l)/dble(contador(l))
           varr2(igen,l)=(r2dummy(l)*r2dummy(l))/dble(contador(l))
!
           ldaverr2(igen,l)=ldr2dummy(l)/dble(contador(l))
           ldvarr2(igen,l)=(ldr2dummy(l)*ldr2dummy(l))/dble(contador(l))
           !                                                                                                   
           noDaverr2(igen,l)=noDr2dummy(l)/dble(contador(l))
           noDvarr2(igen,l)=(noDr2dummy(l)*noDr2dummy(l))/dble(contador(l))

        enddo

!     endif                                                                                                    
!                 print *,irep,igen,l,k,dprime,r2                                                              


     enddo


     open(65,file=trim(adjustl(afix))//'corr2chr.dat')
     open(75,file=trim(adjustl(afix))//'ld2chr.dat')
     open(95,file=trim(adjustl(afix))//'noDcorr2chr.dat')
     open(85,file=trim(adjustl(afix))//'var2chr.dat') !! file to write variances
     do l=0,nloci/10
        !format
        write(65,aformat)l,(averr2(igen,l)/dble(irep),igen=0,ngen)
        write(75,aformat)l,(ldaverr2(igen,l)/dble(irep),igen=0,ngen)
        write(95,aformat)l,(noDaverr2(igen,l)/dble(irep),igen=0,ngen)
     enddo
     write(65,*)'irep',irep
     close(65)
     close(75)
     close(85)
     close(95)
  enddo
!
! averages over reps and how to write them ...
!
  open(65,file=trim(afix)//'corr2chr.dat')
  open(75,file=trim(afix)//'ld2chr.dat')
  open(85,file=trim(afix)//'var2chr.dat') !! file to write variances
  open(95,file=trim(adjustl(afix))//'noDcorr2chr.dat')
  do l=0,nloci/10
     !     write(65,'(i4,201f12.7)')l,(averr2(igen,l)/dble(nrep),igen=0,ngen)
     write(65,aformat)l,(averr2(igen,l)/dble(nrep),igen=0,ngen)
     write(75,aformat)l,(ldaverr2(igen,l)/dble(nrep),igen=0,ngen)
     write(95,aformat)l,(noDaverr2(igen,l)/dble(nrep),igen=0,ngen)
     varr2(igen,l)=(varr2(igen,l)-averr2(igen,l)/dble(nrep))/dble(nrep)
     ldvarr2(igen,l)=(ldvarr2(igen,l)-ldaverr2(igen,l)/dble(nrep))/dble(nrep)
     write(85,aformat2)l,(ldaverr2(igen,l),varr2(igen,l),noDvarr2(igen,l),igen=0,ngen)
  enddo
  close(65)

end program r2Sved
!------ subroutine to order array in increasing order --
subroutine sortincrease(vectortosort,ncross)
  implicit none
  integer :: ncross, i,j, temp
  integer:: vectortosort(ncross), vectorlength

  vectorlength=ncross
  do i=1,vectorlength-1
     do j=i,vectorlength
        if (vectortosort(j).lt.vectortosort(i))then
           temp=vectortosort(i)
           vectortosort(i)=vectortosort(j)
           vectortosort(j)=temp
        endif
     enddo
  enddo
end subroutine sortincrease
!------------------
!--------------
function poisson(media)
  implicit none
  integer :: k, poisson, sem1,sem2
  real :: ranmar
  double precision :: L, p, media
!  call rmarin(sem1,sem2)

  L=exp(-dble(media))
  k=0
  p=1
  do while(p>L)
     k=k+1
     p=p*dble(ranmar())
  enddo


  poisson=k-1
end function poisson
!====================
function crossover(u) !samples number of crossovers from a poisson distribution
!====================
! Poisson distribution: Prob k occurrences, when the expected value is l
! is f(k,l)= ( l^k exp(-l)) / k!
! mean is l, standard dev is sqrt(l)
!
  implicit none
  INTEGER :: crossover,nc
  double precision :: u,x,y
  x=EXP(-1.0d0) 
!variance of a poisson is equal to its mean, mean is 1 recombination per Morgan
  y=x
  nc=0
  DO
     IF(x>u) exit
     nc=nc+1
     y=y/dble(nc)
     x=x+y
  END do
  crossover=nc
END function crossover
!----------------------------------------------
SUBROUTINE mkdprime (dprime, r2,ldr2,noDr2, gen_1, gen_2,n_ind)
!-------------------
! Computes D prime and R2 linkage disequilibrium measures
  ! between loci gen_1 and gen_2
  implicit none
  real*8, parameter :: tol=0.0001d0
  integer :: i, j, n_ind, gen_1, gen_2, nmedios
  real*8:: pi, pj, dij, dmax, dprime, r2,noDr2
  dimension gen_1(n_ind),gen_2(n_ind)
  real*8 g1(n_ind/2), g2(n_ind/2),ldr2
  real*8 sa1,sa2,sb1,sb2,sab,var1,var2,cov,pa,pb
  r2=0.d0;ldr2=0.d0; noDr2=0.d0
  
  g1 = 0
  g2 = 0
  pa=0.d0;pb=0.d0
  pa=dble(count(gen_1(:).eq.1))/dble(n_ind)
  pb=dble(count(gen_2(:).eq.1))/dble(n_ind)
!  write(6,*)pa,pb,maxval(gen_1(:)),maxval(gen_2(:))
  nmedios = n_ind/2
  do i = 1, nmedios
    ! codificar genotipos de los n/2 individuos
    if (gen_1(i)*gen_1(i+nmedios) == 1) then
      g1(i) = 1
    elseif (gen_1(i)*gen_1(i+nmedios) == 2) then
      g1(i) = 2
    elseif ((gen_1(i)*gen_1(i+nmedios) == 4)) then
      g1(i) = 3
    end if
    if (gen_2(i)*gen_2(i+nmedios) == 1) then
      g2(i) = 1
    elseif (gen_2(i)*gen_2(i+nmedios) == 2) then
      g2(i) = 2
    elseif ((gen_2(i)*gen_2(i+nmedios) == 4)) then
      g2(i) = 3
    end if
  enddo

  ldr2=0.d0
  r2 = 0.d0
  dprime = 0.d0

  do i=minval(gen_1(:)), maxval(gen_1(:))
     pi = count(gen_1(:)==i)
     pi = pi/n_ind
     do j=minval(gen_2(:)), maxval(gen_2(:))
        pj = count(gen_2(:)==j)
        pj = pj / n_ind
        !print *,pi,pj
        if (pi>tol .and. pj>tol .and. (1-pi)>tol .and. (1-pj)>tol) then
           dij = count((gen_1(:)==i) .and. (gen_2(:)==j))
           dij = dij/n_ind - pi*pj  
           if (dij<0.0) then
              dmax = min (pi*pj,(1-pi)*(1-pj) )
           else
              dmax = min (pi*(1-pj), pj*(1-pi) )
           endif
           if (dmax>0.) dprime = dprime + pi*pj*abs(dij)/dmax
           ldr2 = ldr2 + dij**2 / ((1.d0-pi)*(1.d0-pj))
       !    write(6,'(3f10.7,2i3)')dij,pi,pj,i,j
        !   pause
        endif
     enddo
  enddo


  !calcular r2 de weir (machaca el otro)
  sa1=0
  sa2=0
  sb1=0
  sb2=0
  sab=0
  do i = 1, nmedios
    sa1=sa1+g1(i)
    sa2=sa2+g1(i)**2
    sb1=sb1+g2(i)
    sb2=sb2+g2(i)**2
    sab=sab+g1(i)*g2(i)
  enddo
  var1 = (sa2-sa1*sa1/nmedios)/nmedios
  var2 = (sb2-sb1*sb1/nmedios)/nmedios
  cov  = (sab-sa1*sb1/nmedios)/nmedios
  r2 = (cov*cov)/(var1*var2)
  noDr2= (cov*cov)/(4*pa*pb*(1.d0-pa)*(1.d0-pb))

 ! write(6,*)cov**2,pa,pb,ldr2,r2
 ! write(6,*)var1,var2
 ! write(6,*)pa*(1-pa),pb*(1-pb)
 ! pause


!--------------
 ! print *,'========================='
end SUBROUTINE mkdprime

!----------------------------------------------
! ranmar y rmarin
!---------------------------------------------
! I received this from STUART@ADS.COM,
! "Toward a Universal Random Number Generator" by George Marsaglia and Arif
! Zaman. Florida State University Report: FSU-SCRI-87-50 (1987).
! It was later modified by F. James and  published in "A Review of
! Pseudo-random Number Generators"
! Stuart says this is the BEST KNOWN RANDOM NUMBER GENERATOR available.
! It passes all tests for a random number generator, and has a period
! of 2^144, and will give bit-identical results on all machines with
! at least 24-bit manitissas in the floating point representation.
! The algorithm is a combination of Fibonacci sequence (with lags of 97
! and 33, and operation "subtraction plus one, modulo one") and an
! "arithmetic sequence" (using subtraction.
! There are three routines contained herein. I have separated them and the
! notes with double lines of asterisks (just delete all between them and the
! lines of asterisks).
! The first program, called TstRAN is a test for correct operation.
! The second part, RANMAR, is the function itself, which must be included
! in any program where you wish to call it from.
! The third part is the initialization routine, which asks for the number seeds.
! Note the way the function is called, for anyone not familiar. The value
! returned by the RANMAR function is fully scalable: I use it between 0 & 100
! by saying
!           x=100*RANMAR()
! Now, RMARIN and RANMAR share initialization variables in a common block,
! so they can be separated between subroutines, etc (again, for anyone
! unfamiliar, I call RMARIN as a subroutine from my main, and then include
! RANMAR as a function in some of my subroutines where I need to generate
! random numbers).
! **************************************************************************
! **************************************************************************
!      PROGRAM tstran
!      INTEGER ij, kl, i
!C These are the seeds needed to produce the test case results
!      ij = 1802
!      kl = 9373
!C Do the initialization
!      CALL rmarin(ij,kl)
!C Generate 20000 random numbers
!      DO 10 i = 1, 20000
!        x = ranmar()
!   10 CONTINUE
!      PRINT *,'GOT THIS FAR!'
!c If the random number generator is working properly, the next six random
!C numbers should be:
!C           6533892.0  14220222.0  7275067.0
!C           6172232.0  8354498.0   10633180.0
!      WRITE(*,*)'Correct:     6533892.0  14220222.0  7275067.0 '
!      WRITE(*,*)'Correct:     6172232.0  8354498.0   10633180.0'
!      WRITE(*,20) (4096.0*4096.0*ranmar(), i=1,6)
!   20 FORMAT (3f12.1)
!      END
! **************************************************************************
! **************************************************************************
      SUBROUTINE rmarin(ij,kl)
! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
      REAL u(97), c, cd, cm
      INTEGER i97, j97
      LOGICAL test
      COMMON /raset1/ u, c, cd, cm, i97, j97, test
      DATA test /.false./
      IF( ij .LT. 0  .OR.  ij .GT. 31328  .OR.&
         kl .LT. 0  .OR.  kl .GT. 30081 ) THEN
        PRINT '(A)', ' The first random number seed must have a value &
    & between 0 and 31328'
        PRINT '(A)',' The second seed must have a value between 0 and &
     &30081'
        STOP
      ENDIF
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      DO 2 ii = 1, 97
        s = 0.0
        t = 0.5
        DO 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          IF (mod(l*m, 64) .ge. 32) then
            s = s + t
          ENDIF
          t = 0.5 * t
    3   CONTINUE
        u(ii) = s
    2 CONTINUE
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      test = .true.
      RETURN
      END
! **************************************************************************
! **************************************************************************
      real FUNCTION ranmar()
      REAL u(97), c, cd, cm
      INTEGER i97, j97
      LOGICAL test
      COMMON /raset1/ u, c, cd, cm, i97, j97, test

      IF( .NOT. test ) THEN
        PRINT '(A)',' Call the init routine (RMARIN) before calling RANMAR'
        STOP
      ENDIF
      uni = u(i97) - u(j97)
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      IF(i97 .EQ. 0) i97 = 97
      j97 = j97 - 1
      IF(j97 .EQ. 0) j97 = 97
      c = c - cd
      IF( c .LT. 0.0 ) c = c + cm
      uni = uni - c
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      ranmar = uni
      RETURN
    END FUNCTION ranmar
! **************************************************************************
