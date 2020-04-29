      program ionhydEP 
	use wavelet
      implicit none	 
      integer n,niter,maxit,numit,i,j,int_time,nx,ny,m1,m2,l,k1,k2,
	:nsteps,l2,l1,k3,lk
      parameter (n=6, niter=100000,maxit=1000000,nsteps=10,nx=1,
	:ny=1)
      integer ion(nx,ny),nk,fli,nok,nbad
      double precision F,alpha,B,E, t0,tf,
	:tout(maxit),yout(maxit,n),yyout(maxit,n),ttout(maxit),
     :dn(maxit),lce(maxit),x1(nx),x2(nx),x3,x6,
     :x(n),yy(n),w,nrm,oldy(maxit),oldx(maxit),df,
     :ttime(maxit),tin(maxit*nsteps),oldxx(maxit),oldyy(maxit),
	:oldl(maxit),oldll(maxit),llce(maxit),oldpx(maxit),
	:oldpy(maxit),oldppx(maxit),oldppy(maxit),dfnew(maxit),h,
	:k(maxit),q(maxit),eps,h1,ystart(n),hmin,
	:sub1,sub2,sub3,sub4,sub5,sub6,sub7,tt,pw(maxit),ppw(maxit),
	:p2(maxit),r2(maxit),pr(maxit),avlz2(maxit),av_n(maxit),
	:a2(maxit),Tx(maxit),tion,px(nx,ny),py(nx,ny),tea(nx,ny),Ai(nx,ny),
	:J3(nx,ny),Lz(nx,ny),d,el_time,qj(nx),pj,jy,jx,newx(N2),
     :newtime(N2),newx1(niter), ep(nx,ny),
	:t1(niter),t2(niter),newy(N2),newy1(niter),tn(nsteps)
	real fx1(N2),fx2(N2),fx3(N2),freqx1(nscale),
	:freqx2(nscale),freqx3(nscale),nscalex1(nscale,N2),
     :nscalex2(nscale,N2),time1(N2)
	complex cyout1(N2),cyout2(N2) 
      external impampl,norm,derivs,rkqs,odeint,adm2DEP1,eqnjac2DCP,
	:interp
c
      t0=0.
      tf=100.
	do i=1,nsteps-1
	tn(i)=t0+(i-1)*(tf-t0)/nsteps
	enddo
	tn(nsteps)=tf
	h1=tf/float(niter)
c      h1=tf/10000.
      call impampl(F,B,E,alpha)
      x(3)=0.
	eps=1.e-8
	x(6)=-E+eps
c------Initial conditions---------------------------------------
c	do 98, i=1,nx
c      x2(i)=-10.+0.09*i
c	x1(i)=-10.+0.09*i
c   98	continue
c      J3(i)=0.1
c   97 continue
c      do 98, i=1,nx
c	qj(i)=-0.1
      pj=-1.
      x1(1)=5.2
	x2(1)=-0.1
	open(unit=1,file='mainfreq.dat')
      open(unit=8,file='indlyap.dat')
      open(unit=3,file='timelyap.dat') 
      open(unit=2,file='coord1hyd.dat')
      open(unit=4,file='coord2hyd.dat')
	open(unit=5,file='momen1hyd.dat')
	open(unit=6,file='momen2hyd.dat')
	open(unit=7,file='freq1hyd.dat')
	open(unit=9,file='scalogram.dat')
      open(unit=10,file='Q_jac.dat')
	open(unit=12,file='anmomenthyd.dat')
	open(unit=11,file='L_z.dat')
	open(unit=13,file='J_3.dat')
	open(unit=14,file='e_p.dat')
	

c	open(unit=11,file='FLI1.dat')
c     do 31, m1=1,291
c     do 32, m2=1,ny
c     read(9,112) ion(m1,m2)
c  32 continue
c  31 continue
      do 94, m1=1,nx
c      print*,m1
c	Lz(m1)=0.5
      do 95, m2=1,ny

c      print*,m2
c     read(9,112) ion(m1,m2)
c    if (ion(m1,m2).eq.0) then
c	if (Lz(m2)**2.le.J3(m1)**2) then
c      tea(m1,m2)=J3(m1)**4*(dsqrt(1-Lz(m2)**2/J3(m1)**2)+1)**2
c      Ai(m1,m2)=E+0.5/J3(m1)**2+(1+0.5*B)*Lz(m2)
c      if (tea(m1,m2)-1/F**2*(Ai(m1,m2)-B**2/8.*tea(m1,m2))**2.ge.0) then
c      x1(m1,m2)=-1/F*(B**2/8*tea(m1,m2)-Ai(m1,m2))
c	write(*,111) x1(m1,m2)
c      x2(m1,m2)=dsqrt(tea(m1,m2)-x1(m1,m2)**2)
c	write(*,111) x2(m1,m2)
c	py(m1,m2)=Lz(m2)*x1(m1,m2)/tea(m1,m2)
c	write(*,111) py(m1,m2)
c	px(m1,m2)=-x2(m1,m2)/x1(m1,m2)*py(m1,m2) 
c	que(m2)=2*(py(m1,m2)*Lz(m2)-1/dsqrt(x1(m1,m2)**2+x2(m1,m2)**2))
c	:/dsqrt(J3(m1)-Lz(m2))
c	d=0.5*(px(m1,m2)**2+py(m1,m2)**2)-1/dsqrt(tea(m1,m2))-
c	:(1+0.5*B)*Lz(m2)+F*x1(m1,m2)+
c     :B**2/8.*tea(m1,m2)
c      write(*,111) d 
c	else
c	goto 95
c	fli=0
c	write(8,111) fli
c 	endif
c	else
c	goto 95
c	fli=0
c	write(8,111) fli
c	endif
c      if (2*J3(m2).ge.qj(m1)**2+pj) then
c 	Lz=qj(m1)**2+pj**2-J3(m2)
c	jx=qj(m1)*0.5*dsqrt(2*J3(m2)-qj(m1)-pj)
c	jy=-pj*0.5*dsqrt(2*J3(m2)-qj(m1)-pj)
c	px(m1,m2)=2*jy/Lz
c	x1(m1,m2)=Lz**2/(2*jx+1)
c	x2(m1,m2)=0
c	x(6)=0.5/J3(m2)**2+(1+0.5*B)*Lz-B**2/8.*
c	:(x1(m1,m2)+x2(m1,m2))-F*x1(m1,m2)
c	print*,x(6) 
c	py(m1,m2)=Lz/x1(m1,m2)      
c	else
c	goto 95
c	endif
	sub1=x1(m1)**2+x2(m2)**2
	x(1)=dsqrt(dsqrt(sub1)+x1(m1)*dcos(x(3))+x2(m2)*dsin(x(3)))
      x(2)=dsqrt(dsqrt(sub1)-x1(m1)*dcos(x(3))
	:-x2(m2)*dsin(x(3)))   
c	sub1=x(1)**2+x(2)**2
c      x(5)=0.5*(1+0.5*B)*sub1*x(1)
      if (abs(x(1)**2).le.(0.1e-07)) then
c     fli=0
	J3(m1,m2)=70
	ep(m1,m2)=70
	Lz(m1,m2)=70
c	print*,'abs(x(1)**2)<0.1e-04' 
c      write(8,112) fli 
	write(11,111) J3(m1,m2)
	write(13,111) Lz(m1,m2)
	write(14,111) ep(m1,m2) 
      goto 95
      endif
      if (((1+B)*(x(1)**2+x(2)**2)**2-8*x(6)+16/(x(1)**2+x(2)**2)
     :-4*F*(x(1)**2-x(2)**2)-(alpha-1)*8*F*(0.5*(x(1)**2-x(2)**2)
     :*dsin(x(3))+x(1)*x(2)*dcos(x(3)))*dsin(x(3)).lt.0.)) then
c      fli=0	
	J3(m1,m2)=70
 	ep(m1,m2)=70
	Lz(m1,m2)=70
c     print*,'Choose some other initial conditions'
c      write(8,112) fli
	write(11,111) J3(m1,m2)
	write(13,111) Lz(m1,m2)
	write(14,111) ep(m1,m2) 
	print*,ep(m1,m2),x(2)
      goto 95
      endif
      x(5)=0.5*x(1)*((1+0.5*B)*(x(1)**2+x(2)**2)-dsqrt((1+B)*
	:(x(1)**2+x(2)**2)**2-
     :8*x(6)+16/(x(1)**2+x(2)**2)-4*F*(x(1)**2-x(2)**2)-
     :(alpha-1)*8*F*(0.5*(x(1)**2-x(2)**2)*dsin(x(3))+x(1)*x(2)
     :*dcos(x(3)))*dsin(x(3)))) 
	 x(4)=-x(2)/x(1)*x(5)
c	write(*,111) x(4)
c	write(*,111) x(5)
c	write(*,111) (x(1)*x(4)-x(2)*x(5))/(x(1)**2+x(2)**2)
c	write(*,111) (x(1)*x(5)+x(2)*x(4))/(x(1)**2+x(2)**2)
	if	((2/sqrt(x1(m1)**2+x2(m2)**2)-
	:((x(1)*x(4)-x(2)*x(5))**2+(x(2)*x(4)+x(1)*x(5))**2)/
     :((x(1)**2+x(2)**2)**2)).lt.0.) then
	J3(m1,m2)=70
 	ep(m1,m2)=70
	Lz(m1,m2)=70
	write(11,111) J3(m1,m2)
 	write(13,111) Lz(m1,m2)
	write(14,111) ep(m1,m2)
	goto 95
	endif 
	J3(m1,m2)=sqrt(1/(2/sqrt(x1(m1)**2+x2(m2)**2)-		
	:((x(1)*x(4)-x(2)*x(5))**2+(x(2)*x(4)+x(1)*x(5))**2)/
     :((x(1)**2+x(2)**2)**2)))
c	 print*,J3(m1)
	Lz(m1,m2)=0.5*(x(1)*x(5)-x(2)*x(4))
c	print*,Lz(m1)
	ep(m1,m2)=sqrt(1-Lz(m1,m2)**2/J3(m1,m2)**2)
	write(13,111) J3(m1,m2)
	write(11,111)  Lz(m1,m2)
	write(14,111)  ep(m1,m2)
cc      do 16, i=n+1,n
cc      x(i)=1./sqrt(6.)
cc   16 continue
cc      x(4)=px(m1,m2)*x(1)+py(m1,m2)*x(2)
cc	x(5)=py(m1,m2)*x(1)-px(m1,m2)*x(2)
cc	write(*,111) x(4)
cc	write(*,111) x(5)
     	do 15, i=1,n							 
       ystart(i)=x(i)
   15 continue  
cc      nk=niter
cc	  df=1.
cc	el_time=TIMEF()
cc      call adm2DEP(n,t0,tf,ttout,x,yout,nk,df)
cc	el_time=TIMEF()
cc	write(*,111) el_time 
cc	call RK(n,t0,tf,ttout,x,yout,niter,nk,maxit);
cc	print*,nk
	l2=0   
      do l=1,nsteps-1
	h=(tn(2)-tn(1))/float(niter)
	hmin=h*niter/float(maxit)
	nk=0
	call odeint(ystart,n,tn(l),tn(l+1),eps,h,hmin,nk,ttout,yyout,
	*derivs,rkqs)  
c	print*,nk 
cc	write(1,112) dlg(df)
cc	time(l+1)=tn(l+1)
cc	print*,l,tn(l+1)/dt
	do l1=1,n
ccc	y(l1)=yyout(l1)
cc      ystart(l1)=yyout(nk,l1)
	do k3=1,nk
	yout(k3,l1)=yyout(k3,l1)
ccc	yout(l+1,l1)=yyout(l1)
      enddo
	enddo
	do k3=1,nk
	tout(k3)=ttout(k3)
	enddo 
	l2=l2+nk
c	l2=nk
cc  53	write(1,112) dlog(df)
	enddo  	

cc	call odeint(ystart,n,t0,tf,eps,h1,nk,ttout,yout,
cc	:derivs,rkqs) 
cc	print*, nk
cc      if (nk.lt.niter) then
cc	print*,nk
cc	  print*,'lyapunov indicator > 0.1e+15'
cc      write(8,111) df
cc      goto 95
cc      endif
cc	goto 95
c------ Original coordinates and momenta expressed in terms of parabolic 
      do 56, i=1,l2	
cc	write(*,111) yout(i,3) 
	sub1=yout(i,1)**2-yout(i,2)**2
	sub2=yout(i,1)*yout(i,2)
	sub3=yout(i,1)**2+yout(i,2)**2
	sub4=yout(i,1)*yout(i,4)-yout(i,2)*yout(i,5)
	sub5=yout(i,2)*yout(i,4)+yout(i,1)*yout(i,5)
	if (alpha.eq.1.) then
	tin(i)=tout(i)
cc	print*, tin(i)
	else
	tin(i)=yout(i,3)
cc	print*, yout(i,3)
	endif
cc	sub6=dcos(tin(i))
cc	sub7=dsin(tin(i))
	sub6=1.
	sub7=0.
      oldx(i)=0.5*sub1*sub6-sub2*sub7
      oldy(i)=0.5*sub1*sub7+sub2*sub6
	if (dsqrt(oldx(i)**2+oldy(i)**2).ge.150.) then
	tion=tin(i)
	else
	tion=tin(l2)
	endif 
	oldpx(i)=(sub4*sub6-sub5*sub7)/sub3
	oldpy(i)=(-sub4*sub7+sub5*sub6)/sub3
  56  continue
c-------------------------------------------------------------------     
      call interp(l2,niter,tin,ttime,oldx,oldxx)
	call interp(l2,niter,tin,ttime,oldy,oldyy)
	call interp(l2,niter,tin,ttime,oldpy,oldppy)
      call interp(l2,niter,tin,ttime,oldpx,oldppx)
	do 76, i=1,niter
	oldl(i)=oldxx(i)*oldppy(i)-oldyy(i)*oldppx(i)	
	p2(i)=oldppx(i)**2+oldppy(i)**2
	r2(i)=oldxx(i)**2+oldyy(i)**2
	pr(i)=oldppx(i)*oldxx(i)+oldppy(i)*oldyy(i)
cc	a2(i)=1/(2/dsqrt(r2(i))-p2(i))
cc	:+pr(i)**2-p2(i)*r2(i)
	avlz2(i)=oldl(i)*oldl(i)
cc	av_n(i)=(1/(2/dsqrt(r2(i))-p2(i))) 
c
c	Tx(i)=av_n(i)**3
c	k(i)=0.5*p2(i)-(1+0.5*B)*oldl(i)+F*oldxx(i)-1/dsqrt(r2(i))
c    :+F*(alpha-1)*(oldxx(i)*dsin(ttime(i))
c     :+oldyy(i)*dcos(ttime(i)))*dsin(ttime(i))+B**2/float(8)*(r2(i))
c	av_n(i)=1/(0.5*(-k(i)-(1+0.5*B)*oldl(i)+F*oldxx(i)
c    :+F*(alpha-1)*(oldxx(i)*dsin(ttime(i))
c     :+oldyy(i)*dcos(ttime(i)))*dsin(ttime(i))+B**2/float(8)*(r2(i))))
cc	print*, av_n(i) 
c	q(i)=(k(i)-F**2/(2*(1+0.5*B)**2))/(oldxx(i)*oldppy(i)-
c	:oldyy(i)*oldppx(i)-F/((1+0.5*B)**2)*(oldppy(i)+(1+0.5*B)*oldxx(i))
c    :+F**2/((1+0.5*B)**3))
c	ep(i)=(1-Lz(i)**2/J3(i)**2)
c	J3(i)=sqrt(1/(2/sqrt(r2(i))-p2(i))) 
   76	continue
c      print*, av_n(1) 
c      do 17, i=1,nk
c      do 18, j=1,n
c      yy(j)=yout(i,j+n)
c   18 continue 
c      call norm(yy,n,nrm)
c      lce(i)=dlog(nrm)
c	write(8,111) lce(i)
c   17 continue
c	subroutine simpson(av_n,niter,sum_n)
c	subroutine simpson(avlz2,niter,sum_Lz) 
c
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
c 
c      do k2=1,niter
c	newx1(k2)=oldxx(k2)
c	newy1(k2)=oldyy(k2)
c	t1(k2)=ttime(k2)
c	enddo 
c	do k2=1,niter
c	write(3,111) t1(k2)
c	enddo
c--------Linear Interpolation of time sequence of action coordinates----     
c      call interp(niter,N2,t1,newtime,newx1,newx)
cc	do i4=1,niter/10
cc	write(3,111) t1(i4)
cc	enddo
c	call interp(niter,N2,t1,newtime,newy1,newy)
cc	call interp(nk,N2,ttout,time,I_z,newI_z)
cc	call interp(nk,N2,ttout,time,q_j,newq_j)
c------Time-frequency wavelet analysis----------------------------------
c	do 89, lk=1,N2
c	time1(lk)=newtime(lk)
cc	write(2,111) newx(l)
c      cyout1(lk)=cmplx(newx(lk),0)
cc      cyout1(l)=dsin(0.7*time(l))+dsin(0.1*time(l)) 
c      cyout2(lk)=cmplx(newy(lk),0.)
cc	cyout3(l)=cmplx(newI_z(l),0.)
cc      cyout2(l)=dcos(0.7*time(l))+dcos(0.01*time(l))	 
c      write(2,113) newx(lk),newy(lk)
cc	write(3,111) time1(l)
c   89 continue 
c      print*,'zdes'
c	call t_f(N2,time1,cyout1,fx1,freqx1,nscalex1)
c	print*,'zdes'
c	call t_f(N2,time1,cyout2,fx2,freqx2,nscalex2)
cc	call t_f(N2,time1,cyout3,fx3,freqx,nscalex)
cc	print*,'zdes'
cc      do 449, i=1,N2
cc	write(9,111) fx2
cc  449 continue
cc      if (kk.lt.7) then
cc	do 31, i=1,nscale/2
cc	print*,i
c	write(9,113) freqx1(i),freqx2(i)
c	do 29, j=1,N2/64
c	write(7,111) nscalex1(i*2,j*64)
c  29	continue
c   31	continue
c	do i=1,N2/64
c	write(1,113) fx1(i),fx2(i)
c	enddo
cc      do i=1,N2
cc	write(7,111) newx(i)
cc	enddo
	do 22, i=1,niter
c	write(3,111) time1(i)
cc	write(*,111) ttime(i)
      write(4,111) oldyy(i)
	write(2,111) oldxx(i)

cc 	write(*,111) oldyy(i)
c	write(5,111) r2(i)
c	write(6,111) p2(i)
c	write(12,111) oldl(i)
cc	write(13,111) av_n(i)
cc	write(11,111)  oldl(i)
cc	write(14,111)  ep(i)

c	write(10,111) k(i)
cc	dfnew(i)=dlog(dfnew(i))
cc	write(8,111) dfnew(i)
   22 continue
c	enddo
   95 continue
   94 continue
	close(unit=1)
      close(unit=3)
	close(unit=2)
      close(unit=5)
	close(unit=6) 
	close(unit=7)
	close(unit=9)
      close(unit=12)
      close(unit=4)
	close(unit=10)
   2  close(unit=8) 
	close(unit=11)
	close(unit=13)
	close(unit=14)
  111 format(1x,e13.7/)
  112 format(1x,i8/)
  113 format(1x,e13.7,1x,e13.7/) 
      end

