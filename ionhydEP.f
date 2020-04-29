      program lyapunov  
      implicit none
      integer n,niter,maxit,numit,i,j,int_time,nx,ny,m1,m2
      parameter (n=6, niter=10000,maxit=1000000,nx=1,ny=1)
      integer ion(nx,ny),nk,fli,nok,nbad
      double precision F,alpha,B,E, t0,tf,ttout(maxit),
     :yout(maxit,n),
     :dn(maxit),lce(maxit),x1(nx),x2(ny),x3,x6,
     :x(n),yy(n),w,nrm,oldy(maxit),oldx(maxit),df(maxit),
     :time(maxit),tin(maxit),oldxx(maxit),oldyy(maxit),
	:oldl(maxit),oldll(maxit),llce(maxit),oldpx(maxit),
	:oldpy(maxit),oldppx(maxit),oldppy(maxit),dfnew(maxit),h,
	:k(maxit),q(maxit),eps,h1,ystart(6),
	:sub1,sub2,sub3,sub4,sub5,sub6,sub7,tt,pw(maxit),ppw(maxit),
	:p2(maxit),r2(maxit),pr(maxit),avlz2(maxit),av_n(maxit),
	:a2(maxit),Tx(maxit),tion
      external impampl,norm,derivs,rkqs,odeint,adm2DEP,eqnjac2DCP
c
      t0=0.
      tf=800.
	h1=tf/float(niter)
c      h1=tf/10000.
      call impampl(F,B,E,alpha)
      x(3)=0.
      x(6)=-E
	eps=1.e-8
c------Initial conditions---------------------------------------
      do 97, i=1,nx
      x2(i)=-0.01
   97 continue
      do 98, i=1,ny
      x1(i)=3.
   98 continue
      open(unit=8,file='indlyap.dat')
      open(unit=3,file='timelyap.dat') 
      open(unit=2,file='coord1hyd.dat')
      open(unit=4,file='coord2hyd.dat')
	open(unit=5,file='momen1hyd.dat')
	open(unit=6,file='momen2hyd.dat')
	open(unit=10,file='Q_jac.dat')
	open(unit=12,file='anmomenthyd.dat')


c	open(unit=11,file='FLI1.dat')
c     do 31, m1=1,291
c     do 32, m2=1,ny
c     read(9,112) ion(m1,m2)
c  32 continue
c  31 continue
      do 94, m1=1,nx
c      print*,m1
      do 95, m2=1,ny
c      print*,m2
c     read(9,112) ion(m1,m2)
c    if (ion(m1,m2).eq.0) then
      sub1=x1(m1)**2+x2(m2)**2
	x(1)=dsqrt(dsqrt(sub1)+x1(m1)*dcos(x(3))+x2(m2)*dsin(x(3)))
      x(2)=sign(dsqrt(dsqrt(sub1)-x1(m1)*dcos(x(3))-x2(m2)*dsin(x(3))),
	:x2(m2))   
	sub1=x(1)**2+x(2)**2
      x(5)=0.5*(1+0.5*B)*sub1*x(1)
      if (abs(x(1)**2).le.(0.1e-02)) then
      fli=0
	print*,'abs(x(1)**2)<0.1e-02' 
      write(8,112) fli
      goto 95
      endif
      if (((1+B)*sub1**2-8*x(6)+16/sub1
     :-4*F*(x(1)**2-x(2)**2)-(alpha-1)*8*F*(0.5*(x(1)**2-x(2)**2)
     :*dsin(x(3))+x(1)*x(2)*dcos(x(3)))*dsin(x(3)).lt.0.)) then
      fli=0
      print*,'Choose some other initial conditions'
      write(8,112) fli
      goto 95
      endif
      x(5)=0.5*x(1)*((1+0.5*B)*sub1-dsqrt((1+B)*sub1**2-
     :8*x(6)+16/sub1-4*F*(x(1)**2-x(2)**2)-
     :(alpha-1)*8*F*(0.5*(x(1)**2-x(2)**2)*dsin(x(3))+x(1)*x(2)
     :*dcos(x(3)))*dsin(x(3))))
      x(4)=-x(5)*x(2)/x(1)
c      do 16, i=n+1,2*n
c      x(i)=1./sqrt(6.)
c   16 continue
     	do 15, i=1,n
       ystart(i)=x(i)
   15 continue  
      nk=niter
	  df(1)=1.
c      call adm2DEP(n,t0,tf,ttout,x,yout,nk)
c	call RK(n,t0,tf,ttout,x,yout,niter,nk,maxit);
	call odeint(ystart,n,t0,tf,eps,h1,nk,ttout,yout,derivs,
	:rkqs)
	print*,nk 
c      if (nk.lt.niter) then
c	  print*,'lyapunov indicator > 0.1e+15'
c     write(8,112) nk
c      goto 95
c      endif
c------ Original coordinates and momenta expressed in terms of parabolic 
      do 56, i=1,nk	
	sub1=yout(i,1)**2-yout(i,2)**2
	sub2=yout(i,1)*yout(i,2)
	sub3=yout(i,1)**2+yout(i,2)**2
	sub4=yout(i,1)*yout(i,4)-yout(i,2)*yout(i,5)
	sub5=yout(i,2)*yout(i,4)+yout(i,1)*yout(i,5)
	if (alpha.eq.1.) then
	tin(i)=ttout(i)
	else
	tin(i)=yout(i,3)
	endif
c	sub6=dcos(tin(i))
c	sub7=dsin(tin(i))
	sub6=1.
	sub7=0.
      oldx(i)=0.5*sub1*sub6-sub2*sub7
      oldy(i)=0.5*sub1*sub7+sub2*sub6
	if (dsqrt(oldx(i)**2+oldy(i)**2).ge.150.) then
	tion=tin(i)
	else
	tion=tin(nk)
	endif 
	oldpx(i)=(sub4*sub6-sub5*sub7)/sub3
	oldpy(i)=(-sub4*sub7+sub5*sub6)/sub3
   56 continue
c-------------------------------------------------------------------     
	print*,tion
      call interp(nk,niter,tin,time,oldx,oldxx)
	call interp(nk,niter,tin,time,oldy,oldyy)
	call interp(nk,niter,tin,time,oldpy,oldppy)
      call interp(nk,niter,tin,time,oldpx,oldppx)
 	call interp(nk,niter,tin,time,df,dfnew)
	do 76, i=1,niter
	oldl(i)=oldxx(i)*oldppy(i)-oldyy(i)*oldppx(i)	
	p2(i)=oldppx(i)**2+oldppy(i)**2
	r2(i)=oldxx(i)**2+oldyy(i)**2
	pr(i)=oldppx(i)*oldxx(i)+oldppy(i)*oldyy(i)
c	a2(i)=1/(2/dsqrt(r2(i))-p2(i))
c	:+pr(i)**2-p2(i)*r2(i)
	avlz2(i)=oldl(i)*oldl(i)
	av_n(i)=((2/dsqrt(r2(i))-p2(i)))
	Tx(i)=av_n(i)**3
	k(i)=0.5*p2(i)-(1+0.5*B)*(oldxx(i)*oldppy(i)
	:-oldyy(i)*oldppx(i))+F*oldxx(i)-1/dsqrt(r2(i))
     :+F*(alpha-1)*(oldxx(i)*dsin(time(i))
     :+oldyy(i)*dcos(time(i)))*dsin(time(i))+B**2/float(8)*(r2(i))
	q(i)=(k(i)-F**2/(2*(1+0.5*B)**2))/(oldxx(i)*oldppy(i)-
	:oldyy(i)*oldppx(i)-F/((1+0.5*B)**2)*(oldppy(i)+(1+0.5*B)*oldxx(i))
     :+F**2/((1+0.5*B)**3)) 
   76	continue
c      do 18, j=1,n
c      yy(j)=yout(i,j+n)
c   18 continue 
c      call norm(yy,n,nrm)
c      lce(i)=dlog(nrm)
cc	write(8,111) lce(i)
c   17 continue
c	subroutine simpson(av_n,niter,sum_n)
c	subroutine simpson(avlz2,niter,sum_Lz)

cc      fli(m1,m2)=lce(1)
c      do 19, i=1,niter-1
c      if (lce(i+1).gt.lce(i)) then  
cc      fli(m1,m2)=lce(i+1)
c      endif
c   19 continue
c      write(8,112) nk
c	print*, fli(m1,m2)
c     else
c     fli(m1,m2)=0.
c     write(8,111) fli(m1,m2)
c     endif     
   95 continue
   94 continue
	do 22, i=1,niter
	write(3,111) time(i)
      write(4,111) oldyy(i)
	write(2,111) oldxx(i)
	write(5,111) oldppx(i)
	write(6,111) oldppy(i)
	write(12,111) oldl(i)
	write(10,111) q(i)
c	dfnew(i)=dlog(dfnew(i))
c	write(8,111) dfnew(i)
   22 continue

      close(unit=3)
	close(unit=2)
      close(unit=5)
	close(unit=6) 
      close(unit=12)
      close(unit=4)
	close(unit=10)
   2  close(unit=8) 
  111 format(1x,e13.7/)
  112 format(1x,i8/) 
      end

