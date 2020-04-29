	subroutine  RK(n,t0,tf,time,y0,yout,niter,numit,maxit)
c	 niter desired number of iteration
c	 numit actual number of iteration 
c	 y0  initial value
c      [t0 tf] interval of time
c	 yout,time output values
c	 h step size
C	 hmin minimum step size
c	 n number of equations		 
c	 tol required tolerance
c
      integer numit,n,maxit,niter,status
      double precision y0(n),t0,tf,tol,hmin,yout(maxit,n),
     :time(maxit),
     :told,yold(n),tfin1,tfin2,tfin3,yfin1(n),yfin2(n), 
     :yfin3(n),h0,h,h1
c	 
c	 t0=0.
c	 tf=10.
      hmin=(tf-t0)/float(maxit)
c	 niter=1000
c	 y0(1)=0.
c	 y0(2)=1.
      time(1)=t0
      do 14, i=1,n
   14 yout(1,i)=y0(i)
      h0=(tf-t0)/float(niter)
      tol=0.00001
      told=t0
      do 13, i=1,n
   13	yold(i)=y0(i)
      h=h0
      numit=1
      h=h0
      h1=h
      do while ((tf-told).gt.hmin)
      status=1
      call r_k_t(h,n,told,yold,tfin1,yfin1)
c	if ((abs(yfin1(1))+abs(yfin1(2))).gt.1.e+10) then
c	goto 5
c	endif
      call r_k_t(h*0.5,n,told,yold,tfin2,yfin2)
c	if ((abs(yfin1(1))+abs(yfin1(2))).gt.1.e+10) then
c	goto 5
c	endif
      call r_k_t(h*0.5,n,tfin2,yfin2,tfin3,yfin3)	 
c      do 15, i=1,n
      if (abs(yfin3(i)-yfin1(i))*16/float(15).gt.tol) then
      status=0
      endif
c   15	continue
      if (status.eq.0) then
      h=h*0.5
      h1=h
      endif 
      if (status.eq.1) then
      numit=numit+1
      time(numit)=told+h
      do 46, i=1,n
   46	yout(numit,i)=yfin1(i)
c      if ((abs(yfin1(1))+abs(yfin1(2))).gt.1.e+10) then
c	goto 5
c	endif
      told=told+h
      do 76 i=1,n
   76	yold(i)=yfin1(i)
      h=h0
      h1=h
      if ((tf-told).lt.h0 .and. (tf-told).gt.0.) then
      h=(tf-told)
      endif
      endif
      enddo 
      if (status.eq.0) then
      print*, "The required tolerance cann't be achieved"
      else
      print*,"ok"
c	 open(unit=10,file='timeOCS3.dat') 
c       do 67, i=0,numit
c67	 write(10,116) time(i)
c	 close(unit=10)
c	 open(unit=11,file='trajOCS3.dat')
c	 do 68, i=0,numit
c68	 write(11,116) yout(i,1)
c	 close(unit=11)
c	 open(unit=12,file='trajrkutt2.dat')
c	 do 69, i=0,numit
c69	 write(12,116) yout(i,2)
c 	 close(unit=12)
      endif
	numit=numit-1
  116	format(1x,e13.7/)
  114	format(1x,e13.7/)
    5 return
      end subroutine
c
      subroutine r_k_t(hh,n,t,yinit,tfinal,y1)
      integer n
      double precision k1(1:n),k2(1:n),k3(1:n),k4(1:n),ydh1(1:n),
     :ydh2(1:n),ydh3(1:n),ydh4(1:n),gg(1:n),yinit(1:n),
     :y1(1:n),hh,t,tfinal
	external derivs 
c
c      call eqn_H2D(n,t,yinit,ydh1)
      call derivs(t,yinit,ydh1)
      do 34, i=1,n
      k1(i)=hh*ydh1(i)
      gg(i)=yinit(i)+k1(i)*0.5
   34	continue
c      call eqn_H2D(n,t+hh*0.5,gg,ydh2)
      call derivs(t+hh*0.5,gg,ydh2)
      do 35, i=1,n
      k2(i)=hh*ydh2(i)
      gg(i)=yinit(i)+k2(i)*0.5
   35	continue
c      call eqn_H2D(n,t+hh*0.5,gg,ydh3)
      call derivs(t+hh*0.5,gg,ydh3)	 
      do 63, i=1,n
      k3(i)=hh*ydh3(i)
      gg(i)=yinit(i)+k3(i)
   63	continue
c      call eqn_H2D(n,t+hh,gg,ydh4)
      call derivs(t+hh,gg,ydh4)
      do 73, i=1,n
      k4(i)=hh*ydh4(i)
      y1(i)=yinit(i)+1/6.*(k1(i)+2.*k2(i)+2.*k3(i)+k4(i))
   73 continue	  
      tfinal=t+hh
      end subroutine
c
      subroutine eqn_H2D(neq,t,x,yprime)
	integer neq
	double precision B,alpha,F,E,x(neq),yprime(neq),t,
	:qn_n,qn_lz,sigma_xf
	external impampl
	call impampl(F,B,E,alpha)
	yprime(1)=x(4)+0.5*(1+0.5*B)*(x(1)**2+x(2)**2)*x(2)
      yprime(2)=x(5)-0.5*(1+0.5*B)*x(1)*(x(1)**2+x(2)**2)
      yprime(3)=x(1)**2+x(2)**2
      yprime(4)=-(2*x(6)*x(1)-(1+0.5*B)*x(1)*
     :(x(1)*x(5)-x(2)*x(4))-0.5*(1+0.5*B)*
     :(x(1)**2+x(2)**2)*x(5)+2*F*x(1)**3+
     :6./32.*B**2*x(1)*(x(1)**2+x(2)**2)**2
     :+2*(alpha-1)*F*x(1)*(1/2*
     :(x(1)**2-x(2)**2)*dsin(x(3))+x(2)*x(1)*dcos(x(3)))*dsin(x(3))+
     :(alpha-1)*F*(x(1)**2+x(2)**2)*
     :(x(1)*dsin(x(3))+x(2)*dcos(x(3)))*dsin(x(3)))
      yprime(5)=-(2*x(2)*x(6)-(1+0.5*B)*x(2)*(x(1)*x(5)-x(2)*x(4))
     :+0.5*(1+0.5*B)*
     :(x(1)**2+x(2)**2)*x(4)-2.*F*x(2)**3+
     :6./32.*B**2*x(2)*(x(1)**2+x(2)**2)**2
     :+2*(alpha-1)*F*x(2)*(1/2*(x(1)**2-x(2)**2)
     :*dsin(x(3))+x(1)*x(2)*dcos(x(3)))*dsin(x(3))
     :+(alpha-1)*F*(x(1)**2+x(2)**2)*
     :(-x(2)*dsin(x(3))+x(1)*dcos(x(3)))*dsin(x(3)))
      yprime(6)=-((alpha-1)*F*(x(1)**2+x(2)**2)*
     :((x(1)**2-x(2)**2)*dsin(x(3))*dcos(x(3))+
     :x(1)*x(2)*(dcos(x(3))**2-dsin(x(3))**2))) 
      
      end subroutine eqn_H2D	
