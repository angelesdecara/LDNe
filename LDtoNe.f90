!
! program to see SF, WH or Sved depending on r2 measured from ld
! or r2 from correlations
!
program r2eqs
  implicit none
  integer,parameter:: dp=KIND(1.d0)
  integer:: i,j,k,l,idummy,igen,io,jant,jpost,m,cmin,cmax
  integer,parameter::nloci=1000,nind=50,ngen=200 ! within chromosome
  real(dp):: r2ld(0:ngen,0:nloci+1),r2corr(0:ngen,0:nloci),Q,weir,feldman,c
  real(dp)::Nc,length, N_Q, N_weir, N_feldman,Ntpm
  real(dp)::tabla(0:430,209),ejec(0:430),ejen(209),r2dummy(209),p
  character(len=10)::aind,alength,afix
  character(len=20)::aformat,agen,aformat2

  length=10.d0
  write(aind,'(i6)')nind
  write(alength,'(i3)')int(length)
  afix='N'//trim(adjustl(aind))//'L'//trim(adjustl(alength))
!  aformat='(i5,'//trim(adjustl(agen))//'f12.7)'
!  aformat2='(i5,'//trim(adjustl(aformat2))//'f12.7)'

  open(21,file='r2_multi.txt')
  do i=1,209
     do j=0,430
        read(21,*)ejen(i),ejec(j),tabla(j,i)
        !write(6,*)ejen(i),ejec(j),tabla(j,i)
     enddo
  enddo
  open(65,file=trim(afix)//'corr2chr.dat')
  open(75,file=trim(afix)//'ld2chr.dat')

  io=0; j=0
!  do while(io.eq.0)!i=0,nloci
  do j=0,nloci
     read(65,*)i,(r2corr(l,j),l=0,ngen)
     read(75,*)i,(r2ld(l,j),l=0,ngen)
!     write(6,*)i,j,r2corr(j,20)
     write(6,*)i,j,r2ld(j,200)
 !    pause
     !j=j+1
  enddo
  close(65)
  close(75)
  write(6,*)'aqui'
  open(11,file=trim(adjustl(afix))//'N_Qwf.dat')
  do i=0,nloci
     N_Q=0.d0; N_weir=0.d0; N_feldman=0.d0
!
! Haldane's function to pass i within chromosome
! to rec fraction r=(1/2)(1-exp(-2d)), d = i/nloci for nloci in 1M
!
     !igen=100
     igen=100
     c=0.5d0
     if(i.gt.0) c=0.5d0*(1.d0-dexp(-2.d0*dble(i)/dble(nloci)))
     Nc=nind*c
     N_Q = (1.0d0/r2ld(igen,i)-1.d0)
     N_Q = N_Q/(4.0d0*c*((1.0d0-c/2.0d0)/(1.0d0-c)**2))
     N_weir = c**2+(1.0d0-c)**2
     N_weir = N_weir / (2.0d0*r2corr(igen,i)*c*(2.0d0-c))
     N_feldman = (1.0d0/r2ld(igen,i)-(1.0d0-c)**2)/(4.0d0*c-2.0d0*c**2)
     !
     !      write(6,*)'r2ld ',i,r2ld(igen,i)
     m=c*1000; cmin=m;cmax=m+1
     if(c.lt..5d0)then
        do j=1,209
           if(tabla(cmin,j)>0.d0) then
              p=(ejec(cmax)-ejec(cmin))/(c-ejec(cmin))
              r2dummy(j)=(tabla(cmax,j)-tabla(cmin,j)+p*tabla(cmin,j))/p
           endif
        enddo
     else
        do j=1,209
           r2dummy(j)=tabla(0,j)
           !write(6,*)'aqui',j,r2dummy(j),r2ld(igen,i)
        enddo
     endif
!     write(6,*)'aqui',c,p
     do j=1,209
        if(r2dummy(j)>1.d-5.and.r2dummy(j)<r2ld(igen,i))then
           jpost=j
           exit
        endif
     enddo
     do j=209,1,-1
        if(r2dummy(j)>1.d-5.and.r2dummy(j)>r2ld(igen,i))then
           jant=j
           exit
        endif
     enddo
     p=(r2dummy(jpost)-r2dummy(jant))/(r2ld(igen,i)-r2dummy(jant))
     Ntpm=(ejen(jpost)-ejen(jant)+p*ejen(jant))/p
     
     write(11,'(i5,7f15.7)')i,c,r2ld(igen,i),r2corr(igen,i),&
          N_Q,N_weir,N_feldman,Ntpm
     write(*,'(i5,7f15.7)')i,c,r2ld(igen,i),r2corr(igen,i),&
          N_Q,N_weir,N_feldman,Ntpm
  enddo



end program r2eqs
