      subroutine adm2DEP(n,t0,tf,tm,y0,yout,nit,dff)
c     Constant step Adams method
c     niter number of iteration
c     y0  initial value
c     [t0 tf] interval of time
c     yout,time output values
c     h step size
c     n number of equations		
c
      integer nn,nit,n,l
      double precision y0(1:2*n),t0,tf,yout(nit,2*n),tm(1:nit),
     :h,ti(1:nit),yy(1:nit,1:2*n),
     :y2(1:2*n),y3(1:2*n),y4(1:2*n),y5(1:2*n),y7(1:2*n),
     :y6(1:2*n),yk(1:2*n),ydh1(1:2*n),
     :ydh2(1:2*n),ydh3(1:2*n),ydh4(1:2*n),ydh5(1:2*n),
     :ydh6(1:2*n),
     :k1(2*n),k2(2*n),k3(2*n),k4(2*n),y1(2*n),dh(4,2*n),
     :ydh7(1:2*n),dff 
	external derivs,eqnjac2DCP
 
      nn=2*n
      tm(1)=t0
      h=(tf-t0)/float(nit)
      do 110, i=1,7
      ti(i)=t0+(i-1)*h
  110 continue
c      dff(1)=1.
      do 44, i=1,nn 
      yy(1,i)=y0(i)
      y1(i)=yy(1,i)
   44 continue
c     call funk(n,ti(1),y0,ydh1)
c-------------Method Runge-Kutta--------------
      do 12, l=2,4
      call derivs(ti(l-1),y1,k1)
c	call eqnjac2DCP(n,ti(l-1),y1,k1,dff)
 

      do 13, i=1,nn
      dh(l-1,i)=k1(i)
   13 continue
      do 45, i=1,nn
c     yy(2,i)=yy(1,i)+h*ydh1(i)
      yy(l,i)=y1(i)+h*0.5*k1(i)
      y2(i)=yy(l,i)
   45 continue
      call derivs(ti(l-1)+h*0.5,y2,k2)
c      call eqnjac2DCP(n,ti(l-1)+h*0.5,y2,k2,dff)

	do 145, i=1,nn
      yy(l,i)=yy(l-1,i)+0.5*h*k2(i)
      y2(i)=yy(l,i)
  145 continue
      call derivs(ti(l-1)+h*0.5,y2,k3)
c	call eqnjac2DCP(n,ti(l-1)+h*0.5,y2,k3,dff)

      do 146, i=1,nn
      yy(l,i)=yy(l-1,i)+k3(i)*h
      y2(i)=yy(l,i)
  146 continue
      call derivs(ti(l),y2,k4)
c	call eqnjac2DCP(n,ti(l),y2,k4,dff)

      do 147, i=1,nn
      yy(l,i)=yy(l-1,i)+1/6.*h*(k1(i)+2*k2(i)+2*k3(i)+k4(i)) 
      y1(i)=yy(l,i)
  147 continue
   12 continue
      call derivs(ti(4),y1,k1)
c	call eqnjac2DCP(n,ti(4),y1,k1,dff)

      do 77, i=1,nn
      dh(4,i)=k1(i)
   77 continue
      do 14, i=1,nn
      ydh1(i)=dh(1,i)
      ydh2(i)=dh(2,i)
      ydh3(i)=dh(3,i)
      ydh4(i)=dh(4,i)
   14 continue
c-----------Adams-Bashforth Two-Step Explicit Method------
c     call funk(n,ti(2),y2,ydh2)
c     call funk(n,ti(1),y0,ydh1)
c     do 46, i=1,n
c     yy(3,i)=y2(i)+h*(3.*ydh2(i)/2.-ydh1(i)/2.)
c     y3(i)=yy(3,i)
c  46 continue
cc-----------Adams-Bashforth Two-Step Explicit Method-----
c     call funk(n,ti(3),y3,ydh3)
c     do 11, i=1,n
c     yy(3,i)=
cc-----------Adams-Bashforth Three-Step Explicit Method-----
c     do 47, i=1,n
c     yy(4,i)=y3(i)+h*(23.*ydh3(i)-16.*ydh2(i)+5.*ydh1(i))/12.
c     y4(i)=yy(4,i)
c  47 continue
c     call funk(n,ti(4),y4,ydh4)
cc----------A-B Four-Step Explicit Method--------------------
      do 48, i=1,nn
      yy(5,i)=y1(i)+h*(55.*ydh4(i)-59.*ydh3(i)+37.*ydh2(i)-9.
     :*ydh1(i))/24.
      y5(i)=yy(5,i)
   48 continue
      call derivs(ti(5),y5,ydh5)
c	call eqnjac2DCP(n,ti(5),y5,ydh5,dff)

c---------Adams-Moulton Three Step Implicit Method-----------
      do 49, i=1,nn
      yy(5,i)=y1(i)+h*(9.*ydh5(i)+19.*ydh4(i)-5.*ydh3(i)
     :+ydh2(i))/24.
      y5(i)=yy(5,i)
   49 continue
      call derivs(ti(5),y5,ydh5)
c	call eqnjac2DCP(n,ti(5),y5,ydh5,dff)

c----------A-B Five-Step Explicit Method---------------------
      do 33, i=1,nn
c     yy(6,i)=y5(i)+h*(55.*ydh5(i)-59.*ydh4(i)+37.*ydh3(i)-9.
c    :*ydh2(i))/24.
      yy(6,i)=y5(i)+h/720.*(1901.*ydh5(i)-2774.*ydh4(i)+2616.*ydh3(i)
     :-1274.*ydh2(i)+251.*ydh1(i))
c     yy(6,i)=y5(i)+h*(941.*ydh6(i)/2400.+431.*ydh5(i)/450.+
c    :53.*ydh4(i)/150.-1213.*ydh3(i)/1200.+173.*ydh2(i)/1440.
c    :+3.*ydh1(i)/16.)
c
      y6(i)=yy(6,i) 
   33 continue
      call derivs(ti(6),y6,ydh6)
c	call eqnjac2DCP(n,ti(6),y6,ydh6,dff)

c----------A-M Four Step Implicit Method---------------------
      do 50, i=1,nn
      yy(6,i)=y5(i)+h*(251.*ydh6(i)+646.*ydh5(i)-264.*ydh4(i)
     :+106.*ydh3(i)-19.*ydh2(i))/720.
      y6(i)=yy(6,i)
   50 continue
      call derivs(ti(6),y6,ydh6)
c	call eqnjac2DCP(n,ti(6),y6,ydh6,dff)

c---------A-B Six Step Explicit Method-----------------------
      do 51, i=1,nn
      yy(7,i)=y6(i)+h*(4277.*ydh6(i)-7923.*ydh5(i)+9982.*ydh4(i)
     :-7298.*ydh3(i)+2877.*ydh2(i)-475.*ydh1(i))/1440.
      y7(i)=yy(7,i)
   51 continue
      call derivs(ti(7),y7,ydh7)
c	call eqnjac2DCP(n,ti(6),y7,ydh7,dff)

c--------A-M Five Step Implicit Method-------------------------
      do 52, i=1,nn
      yy(7,i)=y6(i)+h*(475.*ydh7(i)+1427.*ydh6(i)-798.*ydh5(i)
     :+482.*ydh4(i)-173.*ydh3(i)+27.*ydh2(i))/1440.
      y7(i)=yy(7,i)
   52 continue
      do 170, k=8,nit
      ti(k)=ti(k-1)+h
      do 171, i=1,nn
      yy(k,i)=yy(k-1,i)+h*(4277.*ydh7(i)-7923.*ydh6(i)+9982.*ydh5(i)-
     :7298.*ydh4(i)+2877.*ydh3(i)-475.*ydh2(i))/1440.
      yk(i)=yy(k,i)
  171 continue
      call derivs(ti(k),yk,ydh7)
c	call eqnjac2DCP(n,ti(k),yk,ydh7,dff)

      do 172, i=1,nn
      yy(k,i)=yy(k-1,i)+h*(475.*ydh7(i)+1427.*ydh6(i)-798.*ydh5(i)
     :+482.*ydh4(i)-173.*ydh3(i)+27.*ydh2(i))/1440.
c     yy(k,i)=yy(k-1,i)+h*(941.*ydh6(i)/2400.+431.*ydh5(i)/450.+
c    :53.*ydh4(i)/150.-1213.*ydh3(i)/1200.+173.*ydh2(i)/1440.
c    :+3.*ydh1(i)/16.)
      yk(i)=yy(k,i)
      ydh1(i)=ydh2(i)
      ydh2(i)=ydh3(i)
      ydh3(i)=ydh4(i)
      ydh4(i)=ydh5(i)
      ydh5(i)=ydh6(i)
      ydh6(i)=ydh7(i)
  172 continue
      call derivs(ti(k),yk,ydh7)
c	call eqnjac2DCP(n,ti(k),yk,ydh7,dff)
c      if ((abs(yk(1))+abs(yk(2))).ge.(0.1e+15)) then
	if (dff.ge.(0.1e+15)) then
      nit=k
      goto 173
      endif
      if (ti(k).ge.tf) goto 173
  170 continue
c	
  173 do 36, i=1,nit
      tm(i)=ti(i)
   36 continue
      do 38, k=1,nit
      do 37, j=1,nn
      yout(k,j)=yy(k,j)
   37 continue 
   38 continue
c      open(unit=10,file='adamstime.dat')
c      do 67, i=1,nit
c      write(10,117) tm(i)
c   67 continue
c      close(unit=10)
c      open(unit=11,file='coor1adams1.dat')
c      do 66, i=1,nit
c   66 write(11,117)	yout(i,1)
c      close(unit=11)
c      open(unit=12,file='coor2adams1.dat')
c     do 68, i=1,nit
c  68 write(12,117) yout(i,2)
c     close(unit=12)
  117 format(1x,e13.7/)
      return
      end
